use rand_core::RngCore;

use std::sync::Arc;

use ff::Field;
use group::CurveProjective;
use pairing::Engine;

use super::{Parameters, VerifyingKey};
use crate::{Circuit, SynthesisError};
use crate::domain::Domain;
use crate::multicore::Worker;
use crate::error::Result;

mod assembly;
use assembly::*;

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
    let worker = Worker::new();
    let mut assembly: ParameterAssembly<E,C> = ParameterAssembly::new(circuit, g1, g2, alpha, beta, gamma, delta, tau, &worker)?;
    let key_pair: KeyPairAssembly<E> = assembly.build_key_pair_assembly()?;
    let mut evaluation_domain: Domain<_,_> = key_pair.blind_evaluation_base(&worker)?; 

    let mut windows: _ = WindowTables::default();
    let based: _ = windows.into_based(&key_pair, &assembly, &evaluation_domain);

    let h: Vec<E::G1Affine> = assembly.compute_h(&mut evaluation_domain, &based.g1, &worker)?;
    let lagrange_coeffs = into_lagrange_coefficients(evaluation_domain);

    let mut writer: _ = EvaluationWriter::new(&key_pair);
    assembly.evaluate(&mut writer, &key_pair, based, &lagrange_coeffs);
    
    if writer.is_unconstrained() {
        return Err(SynthesisError::UnconstrainedVariable)
    }

    let vk: VerifyingKey<E> = assembly.build_verifying_key(&mut writer);
    let (l, a, b_g1, b_g2): _ = writer.filter_non_zero_and_map_to_affine();

    Ok(Parameters {
        vk,
        h: Arc::new(h),
        l: Arc::new(l),
        a: Arc::new(a),
        b_g1: Arc::new(b_g1),
        b_g2: Arc::new(b_g2)
    })
}
