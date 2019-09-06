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

use crate::{error, multicore, fft, arith};
use error::{SynthesisError, Result};
use multicore::Worker;
use arith::{Scalar, Group};

use std::ops;

pub struct EvaluationDomain<'a,E,G> 
where
    E: ScalarEngine
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
    for <'b> G: Group<'b,E>
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
        fft::run_optimal_fft(&mut self.coeffs, self.worker, &self.omega, self.exp);
    }

    pub fn ifft(&mut self) {
        fft::run_optimal_fft(&mut self.coeffs, self.worker, &self.omegainv, self.exp);

        self.worker
            .scope(self.coeffs.len(), |scope, chunk| {
                let minv = self.minv;

                for v in self.coeffs.chunks_mut(chunk) {
                    scope.spawn(move || {
                        for v in v {
                            *v *= &minv;
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
                            *v *= &u;
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
                for v in self.coeffs
                    .chunks_mut(chunk) 
                {
                    scope.spawn(move || {
                        for v in v {
                            *v *= &i;
                        }
                    });
                }
            });
        Ok(())
    }
}

impl<'a,'b,E,G> ops::SubAssign<&'a Self> for EvaluationDomain<'b,E,G> 
where
    E: ScalarEngine,
    G: Group<'a,E>
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
                        for (a, b) in a.iter_mut()
                            .zip(b.iter()) 
                        {
                            *a -= b;
                        }
                    });
                }
            });
    }
}

impl<'a,'b,'c,E,G> ops::MulAssign<&'a EvaluationDomain<'b,E,Scalar<E>>> for EvaluationDomain<'c,E,G>
where
    E: ScalarEngine,
    G: Group<'a,E>
{
    fn mul_assign(&mut self, rhs: &'a EvaluationDomain<'b,E,Scalar<E>>) {
        assert_eq!(self.coeffs.len(), rhs.coeffs.len());

        self.worker
            .scope(self.coeffs.len(), |scope, chunk| {
                for (lhs, rhs) in self
                    .coeffs
                    .chunks_mut(chunk)
                    .zip(rhs.coeffs.chunks(chunk))
                {
                    scope.spawn(move || {
                        for (l,r) in lhs.iter_mut()
                            .zip(rhs.iter()) 
                        {
                            *l *= &r.0;
                        }
                    });
                }
            });
    }       
}

impl<'a,E,G> AsRef<[G]> for EvaluationDomain<'a,E,G> 
where
    E: ScalarEngine,
    G: Group<'a,E>
{
    fn as_ref(&self) -> &[G] {
        &self.coeffs
    }
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
                        prod *= &b.0;
                        naive[i1 + i2] += &prod;
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
                let mut v2 = EvaluationDomain::new(v1.into_coeffs().clone(), &worker).unwrap();

                for log_cpus in log_d..min(log_d + 1, 3) {
                    fft::parallel_fft(&mut v1.into_coeffs(), &worker, &v1.omega, log_d, log_cpus);
                    fft::serial_fft(&mut v2.into_coeffs(), &v2.omega, log_d);

                    assert!(v1.into_coeffs() == v2.into_coeffs());
                }
            }
        }
    }

    let rng = &mut rand::thread_rng();

    test_consistency::<Bls12, _>(rng);
}

