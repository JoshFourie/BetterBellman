use pairing::Engine;
use group::{Wnaf, CurveProjective, CurveAffine};
use ff::{Field, PrimeField};

use crate::{multicore, error, arith, domain, groth16};
use crate::Circuit;
use multicore::Worker;
use error::{SynthesisError, Result};
use arith::Scalar;
use domain::Domain;
use groth16::VerifyingKey;

use super::{eval, key_pair, windows};
use eval::EvaluationWriter;
use key_pair::KeyPairAssembly;
use windows::BasedWindowTables;

pub struct ParameterAssembly<'a,E,C> 
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
    gamma_inv: E::Fr,
    worker: &'a Worker
}

impl<'a,E,C> ParameterAssembly<'a,E,C> 
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
        tau: E::Fr, 
        worker: &'a Worker
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
            gamma_inv,
            worker
        })
    }

    pub fn build_key_pair_assembly(&mut self) -> Result<KeyPairAssembly<E>> {
        let mut key_assembly: _ = KeyPairAssembly::default();
        key_assembly.allocate_input_one()?;
        key_assembly.synthesize_circuit(self.circuit.take()?)?;
        key_assembly.enforce_full_density()?;
        Ok(key_assembly)
    }

    pub fn compute_h(&mut self, domain: &mut Domain<E, Scalar<E>>, based_g1: &Wnaf<usize, &[E::G1], &mut Vec<i64>>, worker: &Worker) -> Result<Vec<E::G1Affine>> {   
        into_powers_of_tau(domain, &self.tau, worker);
        let mut h: Vec<E::G1> = vec![E::G1::zero(); domain.as_ref().len() - 1];
        // Set coeff := t(x)/delta
        let tau_delta: _ = self.set_tau_over_delta(domain)?;
        // Set values of the H query to g1^{(tau^i * t(tau)) / delta}
        set_tau_and_exponentiate_g1_into_h(&mut h, domain.as_ref(), based_g1, &tau_delta, worker);

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


    pub fn evaluate(
        &self, 
        writer: &mut EvaluationWriter<E>, 
        key_pair: &KeyPairAssembly<E>,
        based: BasedWindowTables<'_,E>,
        lagrange_coeffs: &[Scalar<E>]
    ) {
        self.evaluate_inputs(writer, key_pair, &based, lagrange_coeffs);
        self.evaluate_auxilliaries(writer, key_pair, based, lagrange_coeffs);
    }

    fn evaluate_inputs(
        &self, 
        writer: &mut EvaluationWriter<E>, 
        key_pair: &KeyPairAssembly<E>,
        based: &BasedWindowTables<'_,E>,
        lagrange_coeffs: &[Scalar<E>]
    ) {
        let ic: &mut _ = writer.ic
            .as_mut()
            .expect("ic value should not have been taken before evaluation");

        eval::eval(
            &based.g1,
            &based.g2,
            &lagrange_coeffs,
            &key_pair.at_inputs,
            &key_pair.bt_inputs,
            &key_pair.ct_inputs,
            &mut writer.a[0..key_pair.num_inputs],
            &mut writer.b_g1[0..key_pair.num_inputs],
            &mut writer.b_g2[0..key_pair.num_inputs],
            ic,
            &self.gamma_inv,
            &self.alpha,
            &self.beta,
            self.worker,
        );
    }

    fn evaluate_auxilliaries(
        &self, 
        writer: &mut EvaluationWriter<E>, 
        key_pair: &KeyPairAssembly<E>,
        based: BasedWindowTables<'_,E>,
        lagrange_coeffs: &[Scalar<E>]
    ) {
        eval::eval(
            &based.g1,
            &based.g2,
            &lagrange_coeffs,
            &key_pair.at_aux,
            &key_pair.bt_aux,
            &key_pair.ct_aux,
            &mut writer.a[key_pair.num_inputs..],
            &mut writer.b_g1[key_pair.num_inputs..],
            &mut writer.b_g2[key_pair.num_inputs..],
            &mut writer.l,
            &self.delta_inv,
            &self.alpha,
            &self.beta,
            self.worker,
        );
    }

    pub fn build_verifying_key(self, writer: &mut EvaluationWriter<E>) -> VerifyingKey<E> {
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

pub fn into_powers_of_tau<E>(domain: &mut Domain<E,Scalar<E>>, tau: &E::Fr, worker: &Worker)
where
    E: Engine
{
    let powers: &mut [Scalar<E>] = domain.as_mut();

    worker.scope(powers.len(), |scope,chunk| {
        for (i, chunk_of_powers) in powers.chunks_mut(chunk)
            .enumerate()
        {
            scope.spawn(move || {
                let idx: _ = (i * chunk) as u64;
                let mut new = tau.pow(&[idx]);

                for power in chunk_of_powers {
                    *power = Scalar(new);
                    new.mul_assign(&tau)
                }
            });
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
    worker: &Worker
) {
    worker.scope(h.len(), |scope, chunk_size| {
        for (chunk_of_h, chunk_of_powers) in h
            .chunks_mut(chunk_size)
            .zip(domain.as_ref().chunks(chunk_size))
        {
            let mut g1_wnaf = based_g1.shared();

            scope.spawn(move || {
                for (value, power) in chunk_of_h.iter_mut()
                    .zip(chunk_of_powers.iter()) 
                {
                    let exponent: <E::Fr as PrimeField>::Repr = raise_tau_delta_to_exponent(power, &tau_delta).into_repr();

                    // Exponentiate
                    *value = g1_wnaf.scalar(exponent);
                }

                E::G1::batch_normalization(chunk_of_h);
            });
        }
    });
}

fn raise_tau_delta_to_exponent<E: Engine>(power: &Scalar<E>, tau_delta: &E::Fr) -> E::Fr {
    let mut exp = power.0;
    exp.mul_assign(tau_delta);   
    exp
}
