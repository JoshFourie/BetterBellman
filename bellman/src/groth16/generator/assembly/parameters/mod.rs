use pairing::Engine;
use group::{Wnaf, CurveProjective, CurveAffine};

use crate::{error, domain, groth16};
use crate::Circuit;
use error::{SynthesisError, Result};
use domain::{Domain, Scalar};
use groth16::VerifyingKey;

use super::{eval, key_pair, windows};
use eval::{Evaluation, Writer};
use key_pair::{KeyPairWires, KeyPairAssembly};
use windows::BasedWindows;

mod elements;
pub use elements::{Elements, InverseElements, ParameterGroups};

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
        let input_size: usize = kp.num.inputs;
        let input_result_writer: _ = result.as_inputs(input_size);
        self.input_eval(input_result_writer, kp.inputs, win, coeffs)?;

        let aux_result_writer: _ = result.as_aux(input_size);
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

pub fn into_lagrange_coefficients<E>(mut domain: Domain<E, Scalar<E>>) -> Vec<Scalar<E>>
where
    E: Engine
{
    domain.ifft();
    domain.into_coeffs()
}
