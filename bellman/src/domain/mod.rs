pub mod arith;
pub mod linear;
pub mod fft;
pub mod multiexp;

pub use arith::*;
pub use linear::*;
pub use fft::*;
pub use multiexp::*;

use ff::{Field, PrimeField, ScalarEngine};

use crate::{error, multicore, multi_thread};
use error::{SynthesisError, Result};
use arith::{Scalar, Group};

use std::ops;

/// A `Domain` abstraction for
/// performing various kinds of polynomial arithmetic on top of
/// the scalar field.
///
/// In pairing-based SNARKs like Groth16, we need to calculate
/// a quotient polynomial over a target polynomial with roots
/// at distinct points associated with each constraint of the
/// constraint system. In order to be efficient, we choose these
/// roots to be the powers of a 2^n root of unity in the field.
/// This allows us to perform polynomial operations in O(n)
/// by performing an O(n log n) FFT over such a domain.
pub struct Domain<E,G> 
where
    E: ScalarEngine
{
    coeffs: Vec<G>,
    exp: u32,
    omega: E::Fr,
    omegainv: E::Fr,
    geninv: E::Fr,
    minv: E::Fr
}

impl<E,G> Domain<E,G> 
where
    E: ScalarEngine,
    for <'a> G: Group<'a,E>
{
    pub fn new(mut coeffs: Vec<G>) -> Result<Self> {
        let (m,exp): (usize,u32) = Self::size_of(&coeffs)?;
        let omega: E::Fr = Self::square_primitive_root_of_unity_to_degree(exp);

        let omegainv: _ = omega.inverse()?;
        let geninv: _ = E::Fr::multiplicative_generator().inverse()?;

        let casted_m: _ = format!("{}",m);
        let minv: _ = E::Fr::from_str(&casted_m)?.inverse()?;

        Self::maybe_pad_with_zeroes(&mut coeffs, m);

        let domain: _ = Domain {
            coeffs,
            exp,
            omega,
            omegainv,
            geninv,
            minv
        };
        Ok(domain)
    }

    fn maybe_pad_with_zeroes(coeffs: &mut Vec<G>, size: usize) {
        coeffs.resize(size, G::zero());
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

    pub fn as_coeffs(&self) -> &[G] {
        &self.coeffs
    }

    pub fn fft(&mut self) {
        fft::run_optimal_fft(&mut self.coeffs, &self.omega, self.exp);
    }

    pub fn ifft(&mut self) {
        fft::run_optimal_fft(&mut self.coeffs, &self.omegainv, self.exp);
        let coeff_len: usize = self.coeffs.len();
        let mul_inv: E::Fr = self.minv;
        multi_thread!(coeff_len, iter(self.coeffs) => {
            for v in iter => {
                *v *= &mul_inv;
            }
        });
    }

    pub fn distribute_powers(&mut self, g: E::Fr) {
        multi_thread!(self.coeffs.len(), enumerate(self.coeffs) => {
            for (i, v) in coeffs => {
                *v *= &g.pow(&[i as u64]);
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

    pub fn raise_tau_to_size(&self, tau: &mut E::Fr) {
        let size: u64 = self.coeffs.len() as u64;
        let mut tmp: E::Fr = tau.pow(&[size]);
        tmp.sub_assign(&E::Fr::one());
        *tau = tmp;
    }

    pub fn divide_over_coset(&mut self) -> Result<()> {
        let mut tau: _ = E::Fr::multiplicative_generator();
        self.raise_tau_to_size(&mut tau);
        let tau_inv: E::Fr = tau.inverse()?;
        
        multi_thread!(self.coeffs.len(), iter(self.coeffs) => {
            for val in coeffs => {
                *val *= &tau_inv
            }
        });

        Ok(())
    }
}

impl<'a,E,G> ops::SubAssign<&'a Self> for Domain<E,G> 
where
    E: ScalarEngine,
    G: Group<'a,E>
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        assert_eq!(self.coeffs.len(), rhs.coeffs.len());

        multi_thread!(self.coeffs.len(), iter(self.coeffs, rhs.coeffs) => {
            for (l,r) in lhs_coeffs, rhs_coeffs => {
                *l -= &r
            }
        });
    }
}

impl<'a,E,G> ops::MulAssign<&'a Domain<E,Scalar<E>>> for Domain<E,G>
where
    E: ScalarEngine,
    G: Group<'a,E>
{
    fn mul_assign(&mut self, rhs: &'a Domain<E,Scalar<E>>) {
        assert_eq!(self.coeffs.len(), rhs.coeffs.len());

        multi_thread!(self.coeffs.len(), iter(self.coeffs, rhs.coeffs) => {
            for (l,r) in lhs_coeffs, rhs_coeffs => {
                *l *= &r.0
            }
        });
    }       
}

impl<E,G> AsRef<[G]> for Domain<E,G> 
where
    E: ScalarEngine
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

                let mut a = Domain::new(a).unwrap();
                let mut b = Domain::new(b).unwrap();

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
        for _ in 0..5 {
            for log_d in 0..10 {
                let d = 1 << log_d;

                let v1 = (0..d)
                    .map(|_| Scalar::<E>(E::Fr::random(rng)))
                    .collect::<Vec<_>>();
                let mut v1 = Domain::new(v1).unwrap();
                let mut v2 = Domain::new(v1.coeffs.clone()).unwrap();

                for log_cpus in log_d..min(log_d + 1, 3) {
                    fft::parallel_fft(&mut v1.coeffs, &v1.omega, log_d, log_cpus);
                    fft::serial_fft(&mut v2.coeffs, &v2.omega, log_d);

                    assert!(v1.coeffs == v2.coeffs);
                }
            }
        }
    }

    let rng = &mut rand::thread_rng();

    test_consistency::<Bls12, _>(rng);
}
