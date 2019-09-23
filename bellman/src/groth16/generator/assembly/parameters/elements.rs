use pairing::Engine;
use ff::{PrimeField, Field};
use group::{CurveProjective, Wnaf};

use crate::{multi_thread, domain, arith, error};
use arith::Scalar;
use domain::Domain;
use error::{SynthesisError, Result};

pub struct ParameterGroups<E>
where
    E: Engine
{
    pub g1: E::G1,
    pub g2: E::G2
}

impl<E> ParameterGroups<E> 
where
    E: Engine
{
    pub fn new(g1: E::G1, g2: E::G2) -> Self {
        Self { g1, g2 }
    }
}

pub struct Elements<E>
where
    E: Engine
{
    pub alpha: E::Fr,
    pub beta: E::Fr,
    pub gamma: E::Fr,
    pub delta: E::Fr,
    tau: E::Fr
}

impl<E> Elements<E>
where
    E: Engine
{
    pub fn new(alpha: E::Fr, beta: E::Fr, gamma: E::Fr, delta: E::Fr, tau: E::Fr) -> Self {
        Self { alpha, beta, gamma, delta, tau }
    }

    pub fn map_powers_of_tau(&self, domain: &mut [Scalar<E>]) {
        multi_thread!(domain.len(), enumerate(domain) => {
            for (i, power) in powers => {
                let exp: E::Fr = self.tau.pow(&[i as u64]);
                *power = Scalar(exp)
            }
        });
    }

    // Set coeff := t(x)/delta
    pub fn set_tau_over_delta(&mut self, domain: &Domain<E, Scalar<E>>, inverse: &InverseElements<E>) {
        domain.raise_tau_to_size(&mut self.tau);
        self.tau.mul_assign(&inverse.delta);   
    }

    // Set values of the H query to g1^{(tau^i * t(tau)) / delta}
    pub fn map_exponent_of_tau(&self, h: &mut Vec<E::G1>, domain: &[Scalar<E>],  based_g1: &Wnaf<usize, &[E::G1], &mut Vec<i64>>) {
        multi_thread!(h.len(), iter(h, domain) => {
            for (value, power) in h_iter, domain_iter => {
                let mut g1_wnaf: _ = based_g1.shared();
                let exponent: _ = self.exponentiate_tau(power).into_repr();
                *value = g1_wnaf.scalar(exponent);
            }
        });
        E::G1::batch_normalization(h);
    }

    fn exponentiate_tau(&self, power: &Scalar<E>) -> E::Fr {
        let Scalar(mut exp): Scalar<E> = *power;
        exp.mul_assign(&self.tau);   
        exp
    }
}

pub struct InverseElements<E>
where
    E: Engine
{
    pub delta: E::Fr,
    pub gamma: E::Fr
}

impl<E> InverseElements<E> 
where
    E: Engine
{
    pub fn new(delta: &E::Fr, gamma: &E::Fr) -> Result<Self> {
        let gamma_inv: E::Fr = gamma.inverse().ok_or(SynthesisError::UnexpectedIdentity)?;
        let delta_inv: E::Fr = delta.inverse().ok_or(SynthesisError::UnexpectedIdentity)?; 

        Ok(Self { 
            delta: delta_inv, 
            gamma: gamma_inv 
        })
    }
}