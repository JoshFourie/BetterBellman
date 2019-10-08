use rand_core::RngCore;

use std::sync::Arc;

use ff::Field;
use group::CurveProjective;
use pairing::Engine;

use super::{Parameters, VerifyingKey};
use crate::{Circuit, SynthesisError};
use crate::domain::Domain;
use crate::error::Result;

mod assembly;
use assembly::Assembly;

/// Generates a random common reference string for
/// a circuit.
pub fn generate_random_parameters<E,C,R>(circuit: C, rng: &mut R) -> Result<Parameters<E>>
where
    E: Engine,
    C: Circuit<E>,
    R: RngCore,
{
    let g1 = E::G1::random(rng);
    let g2 = E::G2::random(rng);
    let alpha = E::Fr::random(rng);
    let beta = E::Fr::random(rng);
    let gamma = E::Fr::random(rng);
    let delta = E::Fr::random(rng);
    let tau = E::Fr::random(rng);

    generate_parameters(circuit, g1, g2, alpha, beta, gamma, delta, tau)
}

/// Create parameters for a circuit, given some toxic waste.
pub fn generate_parameters<E,C>(
    circuit: C,
    g1: E::G1,
    g2: E::G2,
    alpha: E::Fr,
    beta: E::Fr,
    gamma: E::Fr,
    delta: E::Fr,
    tau: E::Fr,
) -> Result<Parameters<E>>
where
    E: Engine,
    C: Circuit<E>,
{
    let mut assembly: _ = Assembly::new(circuit, g1, g2, alpha, beta, gamma, delta, tau)?;
    let mut evaluation_domain: Domain<_,_> = assembly.evaluation_domain()?; 

    let mut windows: _ = assembly::Windows::default();
    let based: _ = windows.as_based(&assembly, &evaluation_domain)?;

    let h: Vec<E::G1Affine> = assembly.h(&mut evaluation_domain, &based.g1)?;

    let lagrange_coeffs = assembly::into_lagrange_coefficients(evaluation_domain);

    assembly.evaluate(&based, &lagrange_coeffs)?;
    
    if assembly.result_is_unconstrained()? {
        return Err(SynthesisError::UnconstrainedVariable)
    }

    let vk: VerifyingKey<E> = assembly.verifying_key()?;
    
    let (l, a, b_g1, b_g2): _ = assembly.results().filter_into_affine();

    Ok(Parameters {
        vk,
        h: Arc::new(h),
        l: Arc::new(l),
        a: Arc::new(a),
        b_g1: Arc::new(b_g1),
        b_g2: Arc::new(b_g2)
    })
}
