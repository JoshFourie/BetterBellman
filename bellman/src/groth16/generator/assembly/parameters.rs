use pairing::Engine;
use group::{Wnaf, CurveProjective, CurveAffine};
use ff::{Field, PrimeField};

use crate::{error, arith, domain, groth16, multi_thread};
use crate::Circuit;
use error::{SynthesisError, Result};
use arith::Scalar;
use domain::Domain;
use groth16::VerifyingKey;

use super::{eval, key_pair, windows};
use eval::WireEvaluation;
use key_pair::{KeyPairAssembly, KeyPairWires};
use windows::BasedWindowTables;

pub struct ParameterAssembly<E,C> 
where
    E: Engine
{
    circuit: Option<C>,
    pub g1: E::G1,
    pub g2: E::G2,
    alpha: E::Fr,
    beta: E::Fr,
    gamma: E::Fr,
    delta: E::Fr,
    tau: E::Fr,
    delta_inv: E::Fr,
    gamma_inv: E::Fr
}

impl<E,C> ParameterAssembly<E,C> 
where
    E: Engine,
    C: Circuit<E>
{
    pub fn new(
        circuit: C, 
        g1: E::G1, 
        g2: E::G2, 
        alpha: E::Fr, 
        beta: E::Fr, 
        gamma: E::Fr,
        delta: E::Fr, 
        tau: E::Fr
    ) -> Result<Self> {
        
        let gamma_inv: E::Fr = gamma.inverse().ok_or(SynthesisError::UnexpectedIdentity)?;
        let delta_inv: E::Fr = delta.inverse().ok_or(SynthesisError::UnexpectedIdentity)?; 

        Ok(ParameterAssembly {
            circuit: Some(circuit),
            g1,
            g2,
            alpha,
            beta,
            gamma,
            delta,
            tau,
            delta_inv,
            gamma_inv
        })
    }

    pub fn build_key_pair_assembly(&mut self) -> Result<KeyPairAssembly<E>> {
        let mut key_assembly: _ = KeyPairAssembly::default();
        key_assembly.allocate_input_one()?;
        key_assembly.synthesize_circuit(self.circuit.take()?)?;
        key_assembly.enforce_full_density()?;
        Ok(key_assembly)
    }

    pub fn compute_h(&mut self, domain: &mut Domain<E, Scalar<E>>, based_g1: &Wnaf<usize, &[E::G1], &mut Vec<i64>>) -> Result<Vec<E::G1Affine>> {   
        into_powers_of_tau(domain, &self.tau);
        let mut h: Vec<E::G1> = vec![E::G1::zero(); domain.as_ref().len() - 1];
        // Set coeff := t(x)/delta
        let tau_delta: _ = self.set_tau_over_delta(domain)?;
        // Set values of the H query to g1^{(tau^i * t(tau)) / delta}
        set_tau_and_exponentiate_g1_into_h(&mut h, domain.as_ref(), based_g1, &tau_delta);

        let into_affine: Vec<_> = h.into_iter()
            .map(|e| e.into_affine())
            .collect();

        Ok(into_affine)
    }

    // Set coeff := t(x)/delta
    fn set_tau_over_delta(&mut self, domain: &Domain<E, Scalar<E>>) -> Result<&E::Fr> {
        domain.raise_tau_to_size(&mut self.tau);
        self.tau.mul_assign(&self.delta_inv);   
        Ok(&self.tau)
    }


    pub fn evaluate(&self, result: &mut WireEvaluation<E>, kp: KeyPairAssembly<E>, win: &BasedWindowTables<'_,E>, coeffs: &[Scalar<E>]) -> Result<()> {
        self.input_eval(result, kp.inputs, kp.num.inputs, win, coeffs)?;
        self.aux_eval(result, kp.aux, kp.num.aux, win, coeffs)?;
        Ok(())
    }

    fn input_eval(&self, result: &mut WireEvaluation<E>, inputs: KeyPairWires<E>, input_len: usize, win: &BasedWindowTables<'_,E>, coeffs: &[Scalar<E>]) -> Result<()> {
        eval::eval(
            win,
            coeffs,
            inputs,
            result.as_mut_inputs(input_len)?,
            &self.gamma_inv,
            &self.alpha,
            &self.beta
        )?;
        Ok(())
    }

    fn aux_eval(&self, result: &mut WireEvaluation<E>, aux: KeyPairWires<E>, aux_len: usize, win: &BasedWindowTables<'_,E>, coeffs: &[Scalar<E>]) -> Result<()> {
        eval::eval(
            win,
            coeffs,
            aux,
            result.as_mut_auxilliaries(aux_len),
            &self.delta_inv,
            &self.alpha,
            &self.beta
        )?;
        Ok(())
    }

    pub fn build_verifying_key(self, writer: &mut WireEvaluation<E>) -> VerifyingKey<E> {
        let g1: E::G1Affine = self.g1.into_affine();
        let g2: E::G2Affine = self.g2.into_affine();

        let alpha_g1: E::G1Affine = g1.mul(self.alpha).into_affine();
        let beta_g1: E::G1Affine = g1.mul(self.beta).into_affine();
        let beta_g2: E::G2Affine = g2.mul(self.beta).into_affine();

        let gamma_g2: E::G2Affine = g2.mul(self.gamma).into_affine();
        let delta_g1: E::G1Affine = g1.mul(self.delta).into_affine();
        let delta_g2: E::G2Affine = g2.mul(self.delta).into_affine();

        let ic: Vec<E::G1Affine> = writer.ic
            .take()
            .expect("ic value should not have been taken before building verification key")
            .into_iter()
            .map(|e| e.into_affine())
            .collect();

        VerifyingKey {
            alpha_g1,
            beta_g1,
            beta_g2,
            gamma_g2,
            delta_g1,
            delta_g2,
            ic
        }
    }
}

pub fn into_powers_of_tau<E>(domain: &mut Domain<E,Scalar<E>>, tau: &E::Fr)
where
    E: Engine
{
    let powers: &mut [Scalar<E>] = domain.as_mut();
    multi_thread!(powers.len(), enumerate(powers) => {
        for (i, power) in powers => {
            *power = Scalar(tau.pow(&[i as u64]))
        }
    });
}

pub fn into_lagrange_coefficients<E>(mut domain: Domain<E, Scalar<E>>) -> Vec<Scalar<E>>
where
    E: Engine
{
    domain.ifft();
    domain.into_coeffs()
}

fn set_tau_and_exponentiate_g1_into_h<E: Engine>(
    h: &mut Vec<E::G1>,
    domain: &[Scalar<E>], 
    based_g1: &Wnaf<usize, &[E::G1], &mut Vec<i64>>, 
    tau_delta: &E::Fr, 
) {
    multi_thread!(h.len(), iter(h, domain.as_ref()) => {
        for (value, power) in h_iter, domain_iter => {
            let mut g1_wnaf: _ = based_g1.shared();
            let exponent: <E::Fr as PrimeField>::Repr = raise_tau_delta_to_exponent(power, &tau_delta).into_repr();
            *value = g1_wnaf.scalar(exponent);
        }
    });
    E::G1::batch_normalization(h);

}

fn raise_tau_delta_to_exponent<E: Engine>(power: &Scalar<E>, tau_delta: &E::Fr) -> E::Fr {
    let mut exp = power.0;
    exp.mul_assign(tau_delta);   
    exp
}
