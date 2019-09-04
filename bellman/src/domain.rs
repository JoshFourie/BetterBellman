//! This module contains an `EvaluationDomain` abstraction for
//! performing various kinds of polynomial arithmetic on top of
//! the scalar field.
//!
//! In pairing-based SNARKs like Groth16, we need to calculate
//! a quotient polynomial over a target polynomial with roots
//! at distinct points associated with each constraint of the
//! constraint system. In order to be efficient, we choose these
//! roots to be the powers of a 2^n root of unity in the field.
//! This allows us to perform polynomial operations in O(n)
//! by performing an O(n log n) FFT over such a domain.

use ff::{Field, PrimeField, ScalarEngine};
use group::CurveProjective;

use crate::error::{SynthesisError, Result};
use crate::multicore::Worker;

pub struct EvaluationDomain<'a,E,G> 
where
    E: ScalarEngine,
    G: Group<E>
{
    coeffs: Vec<G>,
    exp: u32,
    omega: E::Fr,
    omegainv: E::Fr,
    geninv: E::Fr,
    minv: E::Fr,
    worker: &'a Worker
}

impl<'a,E,G> EvaluationDomain<'a,E,G> 
where
    E: ScalarEngine,
    G: Group<E>
{

    pub fn new(mut coeffs: Vec<G>, worker: &'a Worker) -> Result<Self> {
        let (m,exp): (usize,u32) = Self::size_of(&coeffs)?;
        let omega: E::Fr = Self::square_primitive_root_of_unity_to_degree(exp);

        let omegainv: _ = omega.inverse()?;
        let geninv: _ = E::Fr::multiplicative_generator().inverse()?;

        let casted_m: _ = format!("{}",m);
        let minv: _ = E::Fr::from_str(&casted_m)?.inverse()?;

        Self::maybe_pad_with_zeroes(&mut coeffs, m);

        let domain: _ = EvaluationDomain {
            coeffs,
            exp,
            omega,
            omegainv,
            geninv,
            minv,
            worker
        };
        Ok(domain)
    }

    fn maybe_pad_with_zeroes(coeffs: &mut Vec<G>, size: usize) {
        coeffs.resize(size, G::group_zero());
    }

    // Compute omega, the 2^exp primitive root of unity
    fn square_primitive_root_of_unity_to_degree(degree: u32) -> E::Fr {
        let mut omega: _ = E::Fr::root_of_unity();
        for _ in degree..E::Fr::S {
            omega.square();
        }
        omega
    }

    fn size_of(coeffs: &Vec<G>) -> Result<(usize,u32)> {
        let mut m: usize = 1;
        let mut exp: u32 = 0;
        while m < coeffs.len() {
            m *= 2;
            exp += 1;

            let upper_bound: _ = E::Fr::S;
            if exp >= upper_bound {
                return Err(SynthesisError::PolynomialDegreeTooLarge);
            }
        }
        Ok((m,exp))
    }

    pub fn as_mut(&mut self) -> &mut [G] {
        &mut self.coeffs
    }

    pub fn into_coeffs(self) -> Vec<G> {
        self.coeffs
    }

    pub fn fft(&mut self) {
        best_fft(&mut self.coeffs, self.worker, &self.omega, self.exp);
    }

    pub fn ifft(&mut self) {
        best_fft(&mut self.coeffs, self.worker, &self.omegainv, self.exp);

        self.worker
            .scope(self.coeffs.len(), |scope, chunk| {
                let minv = self.minv;

                for v in self.coeffs.chunks_mut(chunk) {
                    scope.spawn(move || {
                        for v in v {
                            v.group_mul_assign(&minv);
                        }
                    });
                }
            });
    }

    pub fn distribute_powers(&mut self, g: E::Fr) {
        self.worker
            .scope(self.coeffs.len(), |scope, chunk| {
                for (i, v) in self.coeffs
                    .chunks_mut(chunk)
                    .enumerate() 
                {
                    scope.spawn(move || {
                        let mut u = g.pow(&[(i * chunk) as u64]);
                        for v in v.iter_mut() {
                            v.group_mul_assign(&u);
                            u.mul_assign(&g);
                        }
                    });
                }
            });
    }

    pub fn coset_fft(&mut self) {
        self.distribute_powers(E::Fr::multiplicative_generator());
        self.fft();
    }

    pub fn icoset_fft(&mut self) {
        let geninv = self.geninv;

        self.ifft();
        self.distribute_powers(geninv);
    }

    pub fn raise_tau_to_size(&self, tau: E::Fr) -> E::Fr {
        let size: u64 = self.coeffs.len() as u64;
        let mut tmp: E::Fr = tau.pow(&[size]);
        tmp.sub_assign(&E::Fr::one());
        tmp
    }

    pub fn divide_over_coset(&mut self) -> Result<()> {
        let tau: _ = E::Fr::multiplicative_generator();
        let i = self.raise_tau_to_size(tau).inverse()?;

        self.worker
            .scope(self.coeffs.len(), |scope, chunk| {
                for v in self.coeffs.chunks_mut(chunk) {
                    scope.spawn(move || {
                        for v in v {
                            v.group_mul_assign(&i);
                        }
                    });
                }
            });
        Ok(())
    }
}

impl<'a,'b,E,G> std::ops::SubAssign<&'a Self> for EvaluationDomain<'b,E,G> 
where
    E: ScalarEngine,
    G: Group<E>
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        assert_eq!(self.coeffs.len(), rhs.coeffs.len());

        self.worker
            .scope(self.coeffs.len(), |scope, chunk| {
                for (a, b) in self
                    .coeffs
                    .chunks_mut(chunk)
                    .zip(rhs.coeffs.chunks(chunk))
                {
                    scope.spawn(move || {
                        for (a, b) in a.iter_mut().zip(b.iter()) {
                            a.group_sub_assign(&b);
                        }
                    });
                }
            });
    }
}

impl<'a,'b,'c,E,G> std::ops::MulAssign<&'a EvaluationDomain<'b,E,Scalar<E>>> for EvaluationDomain<'c,E,G>
where
    E: ScalarEngine,
    G: Group<E>
{
    fn mul_assign(&mut self, rhs: &'a EvaluationDomain<'b,E,Scalar<E>>) {
        assert_eq!(self.coeffs.len(), rhs.coeffs.len());

        self.worker
            .scope(self.coeffs.len(), |scope, chunk| {
                for (a, b) in self
                    .coeffs
                    .chunks_mut(chunk)
                    .zip(rhs.coeffs.chunks(chunk))
                {
                    scope.spawn(move || {
                        for (a, b) in a.iter_mut().zip(b.iter()) {
                            a.group_mul_assign(&b.0);
                        }
                    });
                }
            });
    }       
}

impl<'a,E,G> AsRef<[G]> for EvaluationDomain<'a,E,G> 
where
    E: ScalarEngine,
    G: Group<E>
{
    fn as_ref(&self) -> &[G] {
        &self.coeffs
    }
}

pub trait Group<E: ScalarEngine>: Sized 
    + Copy 
    + Clone 
    + Send 
    + Sync 
{
    fn group_zero() -> Self;
    fn group_mul_assign(&mut self, by: &E::Fr);
    fn group_add_assign(&mut self, other: &Self);
    fn group_sub_assign(&mut self, other: &Self);
}

pub struct Point<G: CurveProjective>(pub G);

impl<G: CurveProjective> PartialEq for Point<G> {
    fn eq(&self, other: &Point<G>) -> bool {
        self.0 == other.0
    }
}

impl<G: CurveProjective> Copy for Point<G> {}

impl<G: CurveProjective> Clone for Point<G> {
    fn clone(&self) -> Point<G> {
        *self
    }
}

impl<G: CurveProjective> Group<G::Engine> for Point<G> {
    fn group_zero() -> Self {
        Point(G::zero())
    }
    fn group_mul_assign(&mut self, by: &G::Scalar) {
        self.0.mul_assign(by.into_repr());
    }
    fn group_add_assign(&mut self, other: &Self) {
        self.0.add_assign(&other.0);
    }
    fn group_sub_assign(&mut self, other: &Self) {
        self.0.sub_assign(&other.0);
    }
}

pub struct Scalar<E: ScalarEngine>(pub E::Fr);

impl<E: ScalarEngine> PartialEq for Scalar<E> {
    fn eq(&self, other: &Scalar<E>) -> bool {
        self.0 == other.0
    }
}

impl<E: ScalarEngine> Copy for Scalar<E> {}

impl<E: ScalarEngine> Clone for Scalar<E> {
    fn clone(&self) -> Scalar<E> {
        *self
    }
}

impl<E: ScalarEngine> Group<E> for Scalar<E> {
    fn group_zero() -> Self {
        Scalar(E::Fr::zero())
    }
    fn group_mul_assign(&mut self, by: &E::Fr) {
        self.0.mul_assign(by);
    }
    fn group_add_assign(&mut self, other: &Self) {
        self.0.add_assign(&other.0);
    }
    fn group_sub_assign(&mut self, other: &Self) {
        self.0.sub_assign(&other.0);
    }
}

fn best_fft<E: ScalarEngine, T: Group<E>>(a: &mut [T], worker: &Worker, omega: &E::Fr, log_n: u32) {
    let log_cpus = worker.log_num_cpus();

    if log_n <= log_cpus {
        serial_fft(a, omega, log_n);
    } else {
        parallel_fft(a, worker, omega, log_n, log_cpus);
    }
}

fn serial_fft<E: ScalarEngine, T: Group<E>>(a: &mut [T], omega: &E::Fr, log_n: u32) {
    fn bitreverse(mut n: u32, l: u32) -> u32 {
        let mut r = 0;
        for _ in 0..l {
            r = (r << 1) | (n & 1);
            n >>= 1;
        }
        r
    }

    let n = a.len() as u32;
    assert_eq!(n, 1 << log_n);

    for k in 0..n {
        let rk = bitreverse(k, log_n);
        if k < rk {
            a.swap(rk as usize, k as usize);
        }
    }

    let mut m = 1;
    for _ in 0..log_n {
        let w_m = omega.pow(&[(n / (2 * m)) as u64]);

        let mut k = 0;
        while k < n {
            let mut w = E::Fr::one();
            for j in 0..m {
                let mut t = a[(k + j + m) as usize];
                t.group_mul_assign(&w);
                let mut tmp = a[(k + j) as usize];
                tmp.group_sub_assign(&t);
                a[(k + j + m) as usize] = tmp;
                a[(k + j) as usize].group_add_assign(&t);
                w.mul_assign(&w_m);
            }

            k += 2 * m;
        }

        m *= 2;
    }
}

fn parallel_fft<E: ScalarEngine, T: Group<E>>(
    a: &mut [T],
    worker: &Worker,
    omega: &E::Fr,
    log_n: u32,
    log_cpus: u32,
) {
    assert!(log_n >= log_cpus);

    let num_cpus = 1 << log_cpus;
    let log_new_n = log_n - log_cpus;
    let mut tmp = vec![vec![T::group_zero(); 1 << log_new_n]; num_cpus];
    let new_omega = omega.pow(&[num_cpus as u64]);

    worker.scope(0, |scope, _| {
        let a = &*a;

        for (j, tmp) in tmp.iter_mut().enumerate() {
            scope.spawn(move || {
                // Shuffle into a sub-FFT
                let omega_j = omega.pow(&[j as u64]);
                let omega_step = omega.pow(&[(j as u64) << log_new_n]);

                let mut elt = E::Fr::one();
                for i in 0..(1 << log_new_n) {
                    for s in 0..num_cpus {
                        let idx = (i + (s << log_new_n)) % (1 << log_n);
                        let mut t = a[idx];
                        t.group_mul_assign(&elt);
                        tmp[i].group_add_assign(&t);
                        elt.mul_assign(&omega_step);
                    }
                    elt.mul_assign(&omega_j);
                }

                // Perform sub-FFT
                serial_fft(tmp, &new_omega, log_new_n);
            });
        }
    });

    // TODO: does this hurt or help?
    worker.scope(a.len(), |scope, chunk| {
        let tmp = &tmp;

        for (idx, a) in a.chunks_mut(chunk).enumerate() {
            scope.spawn(move || {
                let mut idx = idx * chunk;
                let mask = (1 << log_cpus) - 1;
                for a in a {
                    *a = tmp[idx & mask][idx >> log_cpus];
                    idx += 1;
                }
            });
        }
    });
}

// Test multiplying various (low degree) polynomials together and
// comparing with naive evaluations.
#[cfg(feature = "pairing")]
#[test]
fn polynomial_arith() {
    use pairing::bls12_381::Bls12;
    use rand_core::RngCore;

    fn test_mul<E: ScalarEngine, R: RngCore>(rng: &mut R) {
        let worker = Worker::new();

        for coeffs_a in 0..70 {
            for coeffs_b in 0..70 {
                let mut a: Vec<_> = (0..coeffs_a)
                    .map(|_| Scalar::<E>(E::Fr::random(rng)))
                    .collect();
                let mut b: Vec<_> = (0..coeffs_b)
                    .map(|_| Scalar::<E>(E::Fr::random(rng)))
                    .collect();

                // naive evaluation
                let mut naive = vec![Scalar(E::Fr::zero()); coeffs_a + coeffs_b];
                for (i1, a) in a.iter().enumerate() {
                    for (i2, b) in b.iter().enumerate() {
                        let mut prod = *a;
                        prod.group_mul_assign(&b.0);
                        naive[i1 + i2].group_add_assign(&prod);
                    }
                }

                a.resize(coeffs_a + coeffs_b, Scalar(E::Fr::zero()));
                b.resize(coeffs_a + coeffs_b, Scalar(E::Fr::zero()));

                let mut a = EvaluationDomain::new(a, &worker).unwrap();
                let mut b = EvaluationDomain::new(b, &worker).unwrap();

                a.fft();
                b.fft();
                a *= &b;
                a.ifft();

                for (naive, fft) in naive.iter().zip(a.coeffs.iter()) {
                    assert!(naive == fft);
                }
            }
        }
    }

    let rng = &mut rand::thread_rng();

    test_mul::<Bls12, _>(rng);
}

#[cfg(feature = "pairing")]
#[test]
fn fft_composition() {
    use pairing::bls12_381::Bls12;
    use rand_core::RngCore;

    fn test_comp<E: ScalarEngine, R: RngCore>(rng: &mut R) {
        let worker = Worker::new();

        for coeffs in 0..10 {
            let coeffs = 1 << coeffs;

            let mut v = vec![];
            for _ in 0..coeffs {
                v.push(Scalar::<E>(E::Fr::random(rng)));
            }

            let mut domain = EvaluationDomain::new(v.clone(), &worker).unwrap();

            domain.ifft();
            domain.fft();
            assert!(v == domain.coeffs);

            domain.fft();
            domain.ifft();
            assert!(v == domain.coeffs);

            domain.icoset_fft();
            domain.coset_fft();
            assert!(v == domain.coeffs);

            domain.coset_fft();
            domain.icoset_fft();
            assert!(v == domain.coeffs);
        }
    }

    let rng = &mut rand::thread_rng();

    test_comp::<Bls12, _>(rng);
}

#[cfg(feature = "pairing")]
#[test]
fn parallel_fft_consistency() {
    use pairing::bls12_381::Bls12;
    use rand_core::RngCore;
    use std::cmp::min;

    fn test_consistency<E: ScalarEngine, R: RngCore>(rng: &mut R) {
        let worker = Worker::new();

        for _ in 0..5 {
            for log_d in 0..10 {
                let d = 1 << log_d;

                let v1 = (0..d)
                    .map(|_| Scalar::<E>(E::Fr::random(rng)))
                    .collect::<Vec<_>>();
                let mut v1 = EvaluationDomain::new(v1, &worker).unwrap();
                let mut v2 = EvaluationDomain::new(v1.coeffs.clone(), &worker).unwrap();

                for log_cpus in log_d..min(log_d + 1, 3) {
                    parallel_fft(&mut v1.coeffs, &worker, &v1.omega, log_d, log_cpus);
                    serial_fft(&mut v2.coeffs, &v2.omega, log_d);

                    assert!(v1.coeffs == v2.coeffs);
                }
            }
        }
    }

    let rng = &mut rand::thread_rng();

    test_consistency::<Bls12, _>(rng);
}
