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
use eval::{Evaluation, Writer};
use key_pair::{KeyPairWires, KeyPairAssembly};
use windows::BasedWindows;

pub struct ParameterAssembly<E,C> 
where
    E: Engine
{
    circuit: Option<C>,
    pub groups: ParameterGroups<E>,
    elements: Elements<E>,
    inverse: InverseElements<E>
}

impl<E,C> ParameterAssembly<E,C> 
where
    E: Engine,
    C: Circuit<E>
{
    pub fn new(circuit: C, g1: E::G1, g2: E::G2, alpha: E::Fr, beta: E::Fr, gamma: E::Fr, delta: E::Fr, tau: E::Fr) -> Result<Self> {
        let groups: _ = ParameterGroups::new(g1, g2);
        let elements: _ = Elements::new(alpha, beta, gamma, delta, tau);
        let inverse: _ = InverseElements::new(&delta, &gamma)?;

        Ok(ParameterAssembly {
            circuit: Some(circuit),
            groups,
            elements,
            inverse
        })
    }

    pub fn key_assembly(&mut self) -> Result<KeyPairAssembly<E>> {
        let mut key_assembly: _ = KeyPairAssembly::default();

        key_assembly.allocate_input_one()?;
        key_assembly.synthesize_circuit(self.circuit.take()?)?;
        key_assembly.enforce_full_density()?;

        Ok(key_assembly)
    }

    pub fn h(&mut self, domain: &mut Domain<E, Scalar<E>>, based_g1: &Wnaf<usize, &[E::G1], &mut Vec<i64>>) -> Result<Vec<E::G1Affine>> {  
        let mut h: Vec<E::G1> = vec![E::G1::zero(); domain.as_ref().len() - 1];

        self.elements.map_powers_of_tau(domain.as_mut());        
        self.elements.set_tau_over_delta(&domain, &self.inverse);
        self.elements.map_exponent_of_tau(&mut h, domain.as_ref(), based_g1);

        let into_affine: Vec<_> = h.into_iter()
            .map(|e| e.into_affine())
            .collect();

        Ok(into_affine)
    }

    pub fn evaluate(&self, result: &mut Evaluation<E>, kp: KeyPairAssembly<E>, win: &BasedWindows<'_,E>, coeffs: &[Scalar<E>]) -> Result<()> {
        let input_result_writer: _ = result.as_inputs(kp.num.inputs)?;
        self.input_eval(input_result_writer, kp.inputs, win, coeffs)?;

        let aux_result_writer: _ = result.as_aux(kp.num.aux);
        self.aux_eval(aux_result_writer, kp.aux, win, coeffs)?;

        Ok(())
    }

    fn input_eval(&self, input_results: Writer<E>, input_wires: KeyPairWires<E>, win: &BasedWindows<'_,E>, coeffs: &[Scalar<E>]) -> Result<()> {
        if input_results.sanity_check(&input_wires) {
            input_results.eval(win, coeffs, input_wires, &self.inverse.gamma, &self.elements);
            Ok(())  
        } else {
            Err(SynthesisError::MalformedWireSize) 
        }
    }

    fn aux_eval(&self, aux_results: Writer<E>, aux_wires: KeyPairWires<E>, win: &BasedWindows<'_,E>, coeffs: &[Scalar<E>]) -> Result<()> {
        if aux_results.sanity_check(&aux_wires) {
            aux_results.eval(win, coeffs, aux_wires, &self.inverse.delta, &self.elements);
            Ok(())
        } else {
            Err(SynthesisError::MalformedWireSize)
        }
    }

    pub fn into_verifying_key(self, writer: &Evaluation<E>) -> VerifyingKey<E> {
        let g1: E::G1Affine = self.groups.g1.into_affine();
        let g2: E::G2Affine = self.groups.g2.into_affine();
        let elements: _ = self.elements;
        
        let alpha_g1: E::G1Affine = g1.mul(elements.alpha).into_affine();
        let beta_g1: E::G1Affine = g1.mul(elements.beta).into_affine();
        let beta_g2: E::G2Affine = g2.mul(elements.beta).into_affine();

        let gamma_g2: E::G2Affine = g2.mul(elements.gamma).into_affine();
        let delta_g1: E::G1Affine = g1.mul(elements.delta).into_affine();
        let delta_g2: E::G2Affine = g2.mul(elements.delta).into_affine();

        let ic: Vec<E::G1Affine> = writer.ic
            .iter()
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
    fn new(g1: E::G1, g2: E::G2) -> Self {
        Self { g1, g2 }
    }
}

pub struct Elements<E>
where
    E: Engine
{
    pub alpha: E::Fr,
    pub beta: E::Fr,
    gamma: E::Fr,
    delta: E::Fr,
    tau: E::Fr
}

impl<E> Elements<E>
where
    E: Engine
{
    fn new(alpha: E::Fr, beta: E::Fr, gamma: E::Fr, delta: E::Fr, tau: E::Fr) -> Self {
        Self { alpha, beta, gamma, delta, tau }
    }

    fn map_powers_of_tau(&self, domain: &mut [Scalar<E>]) {
        multi_thread!(domain.len(), enumerate(domain) => {
            for (i, power) in powers => {
                let exp: E::Fr = self.tau.pow(&[i as u64]);
                *power = Scalar(exp)
            }
        });
    }

    // Set coeff := t(x)/delta
    fn set_tau_over_delta(&mut self, domain: &Domain<E, Scalar<E>>, inverse: &InverseElements<E>) {
        domain.raise_tau_to_size(&mut self.tau);
        self.tau.mul_assign(&inverse.delta);   
    }

    // Set values of the H query to g1^{(tau^i * t(tau)) / delta}
    fn map_exponent_of_tau(&self, h: &mut Vec<E::G1>, domain: &[Scalar<E>],  based_g1: &Wnaf<usize, &[E::G1], &mut Vec<i64>>) {
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
    delta: E::Fr,
    gamma: E::Fr
}

impl<E> InverseElements<E> 
where
    E: Engine
{
    fn new(delta: &E::Fr, gamma: &E::Fr) -> Result<Self> {
        let gamma_inv: E::Fr = gamma.inverse().ok_or(SynthesisError::UnexpectedIdentity)?;
        let delta_inv: E::Fr = delta.inverse().ok_or(SynthesisError::UnexpectedIdentity)?; 

        Ok(Self { 
            delta: delta_inv, 
            gamma: gamma_inv 
        })
    }
}

pub fn into_lagrange_coefficients<E>(mut domain: Domain<E, Scalar<E>>) -> Vec<Scalar<E>>
where
    E: Engine
{
    domain.ifft();
    domain.into_coeffs()
}
