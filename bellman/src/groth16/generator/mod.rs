use rand_core::RngCore;

use std::sync::Arc;

use ff::{Field, PrimeField};
use group::{CurveAffine, CurveProjective, Wnaf};
use pairing::Engine;

use super::{Parameters, VerifyingKey};
use crate::{Circuit, SynthesisError};
use crate::domain::EvaluationDomain;
use crate::multicore::Worker;
use crate::arith::Scalar;
use crate::error::Result;

mod assembly;
use assembly::*;

mod eval;
use eval::*;

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
pub fn generate_parameters<E, C>(
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
    let mut assembly: _ = KeypairAssembly::default();

    assembly.allocate_one_input()?;

    // Synthesize the circuit.
    circuit.synthesize(&mut assembly)?;

    assembly.enforce_full_density()?;

    let worker = Worker::new();

    // Create bases for blind evaluation of polynomials at tau
    let powers_of_tau = vec![Scalar::<E>(E::Fr::zero()); assembly.num_constraints];
    let mut evaluation_domain = EvaluationDomain::new(powers_of_tau, &worker)?;

    // Compute G1 window table
    let mut g1_wnaf = Wnaf::new();
    let g1_wnaf = g1_wnaf.base(g1, {
        // H query
        (evaluation_domain.as_ref().len() - 1)
        // IC/L queries
        + assembly.num_inputs + assembly.num_aux
        // A query
        + assembly.num_inputs + assembly.num_aux
        // B query
        + assembly.num_inputs + assembly.num_aux
    });

    // Compute G2 window table
    let mut g2_wnaf = Wnaf::new();
    let g2_wnaf = g2_wnaf.base(g2, {
        // B query
        assembly.num_inputs + assembly.num_aux
    });

    let gamma_inverse = gamma.inverse().ok_or(SynthesisError::UnexpectedIdentity)?;
    let delta_inverse = delta.inverse().ok_or(SynthesisError::UnexpectedIdentity)?;

    let mut h = vec![E::G1::zero(); evaluation_domain.as_ref().len() - 1];
    {
        // Compute powers of tau
        {
            let powers_of_tau = evaluation_domain.as_mut();
            worker.scope(powers_of_tau.len(), |scope, chunk| {
                for (i, powers_of_tau) in powers_of_tau.chunks_mut(chunk).enumerate() {
                    scope.spawn(move || {
                        let mut current_tau_power = tau.pow(&[(i * chunk) as u64]);

                        for p in powers_of_tau {
                            p.0 = current_tau_power;
                            current_tau_power.mul_assign(&tau);
                        }
                    });
                }
            });
        }

        // coeff = t(x) / delta
        let mut coeff = evaluation_domain.raise_tau_to_size(tau);
        coeff.mul_assign(&delta_inverse);

        // Compute the H query with multiple threads
        worker.scope(h.len(), |scope, chunk| {
            for (h, p) in h
                .chunks_mut(chunk)
                .zip(evaluation_domain.as_ref().chunks(chunk))
            {
                let mut g1_wnaf = g1_wnaf.shared();

                scope.spawn(move || {
                    // Set values of the H query to g1^{(tau^i * t(tau)) / delta}
                    for (h, p) in h.iter_mut().zip(p.iter()) {
                        // Compute final exponent
                        let mut exp = p.0;
                        exp.mul_assign(&coeff);

                        // Exponentiate
                        *h = g1_wnaf.scalar(exp.into_repr());
                    }

                    // Batch normalize
                    E::G1::batch_normalization(h);
                });
            }
        });
    }

    // Use inverse FFT to convert powers of tau to Lagrange coefficients
    evaluation_domain.ifft();
    let powers_of_tau = evaluation_domain.into_coeffs();

    let mut a = vec![E::G1::zero(); assembly.num_inputs + assembly.num_aux];
    let mut b_g1 = vec![E::G1::zero(); assembly.num_inputs + assembly.num_aux];
    let mut b_g2 = vec![E::G2::zero(); assembly.num_inputs + assembly.num_aux];
    let mut ic = vec![E::G1::zero(); assembly.num_inputs];
    let mut l = vec![E::G1::zero(); assembly.num_aux];

    // Evaluate for inputs.
    eval(
        &g1_wnaf,
        &g2_wnaf,
        &powers_of_tau,
        &assembly.at_inputs,
        &assembly.bt_inputs,
        &assembly.ct_inputs,
        &mut a[0..assembly.num_inputs],
        &mut b_g1[0..assembly.num_inputs],
        &mut b_g2[0..assembly.num_inputs],
        &mut ic,
        &gamma_inverse,
        &alpha,
        &beta,
        &worker,
    );

    // Evaluate for auxiliary variables.
    eval(
        &g1_wnaf,
        &g2_wnaf,
        &powers_of_tau,
        &assembly.at_aux,
        &assembly.bt_aux,
        &assembly.ct_aux,
        &mut a[assembly.num_inputs..],
        &mut b_g1[assembly.num_inputs..],
        &mut b_g2[assembly.num_inputs..],
        &mut l,
        &delta_inverse,
        &alpha,
        &beta,
        &worker,
    );

    // Don't allow any elements be unconstrained, so that
    // the L query is always fully dense.
    for e in l.iter() {
        if e.is_zero() {
            return Err(SynthesisError::UnconstrainedVariable);
        }
    }

    let g1 = g1.into_affine();
    let g2 = g2.into_affine();

    let vk = VerifyingKey::<E> {
        alpha_g1: g1.mul(alpha).into_affine(),
        beta_g1: g1.mul(beta).into_affine(),
        beta_g2: g2.mul(beta).into_affine(),
        gamma_g2: g2.mul(gamma).into_affine(),
        delta_g1: g1.mul(delta).into_affine(),
        delta_g2: g2.mul(delta).into_affine(),
        ic: ic.into_iter().map(|e| e.into_affine()).collect(),
    };

    Ok(Parameters {
        vk: vk,
        h: Arc::new(h.into_iter().map(|e| e.into_affine()).collect()),
        l: Arc::new(l.into_iter().map(|e| e.into_affine()).collect()),

        // Filter points at infinity away from A/B queries
        a: Arc::new(
            a.into_iter()
                .filter(|e| !e.is_zero())
                .map(|e| e.into_affine())
                .collect(),
        ),
        b_g1: Arc::new(
            b_g1.into_iter()
                .filter(|e| !e.is_zero())
                .map(|e| e.into_affine())
                .collect(),
        ),
        b_g2: Arc::new(
            b_g2.into_iter()
                .filter(|e| !e.is_zero())
                .map(|e| e.into_affine())
                .collect(),
        ),
    })
}
