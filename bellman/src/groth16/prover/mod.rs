use std::sync::Arc;

use futures::Future;

use ff::{Field, PrimeField};
use group::{CurveAffine, CurveProjective};
use pairing::Engine;

use super::{ParameterSource, Proof, Result};

use crate::{Circuit, ConstraintSystem, Index, SynthesisError, Variable};
use crate::multiexp::{multiexp, FullDensity};
use crate::multicore::Worker;

pub mod system;
pub use system::*;

pub fn create_proof<E, C, P: ParameterSource<E>>(
    circuit: C,
    mut params: P,
    r: E::Fr,
    s: E::Fr,
) -> Result<Proof<E>>
where
    E: Engine,
    C: Circuit<E>,
{
    let mut prover: _ = ProvingSystem::default();

    prover.alloc_input(
        || "", 
        || Ok(E::Fr::one())
    )?;

    circuit.synthesize(&mut prover)?;

    for i in 0..prover.assignment.input.len() {
        prover.enforce(
            || "", 
            |lc| lc + Variable(Index::Input(i)), 
            |lc| lc, 
            |lc| lc
        );
    }

    let worker = Worker::new();

    let vk = params.get_vk(prover.assignment.input.len())?;

    let h = {
        let a: _ = prover.primefield(&worker)?;
        multiexp(&worker, params.get_h(a.len())?, FullDensity, a)
    };

    // TODO: parallelize if it's even helpful
    let input_assignment = Arc::new(
        prover.assignment
            .input
            .into_iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>(),
    );
    let aux_assignment = Arc::new(
        prover.assignment
            .aux
            .into_iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>(),
    );

    let l = multiexp(
        &worker,
        params.get_l(aux_assignment.len())?,
        FullDensity,
        aux_assignment.clone(),
    );

    let a_aux_density_total = prover.density
        .a_aux
        .get_total_density();

    let (a_inputs_source, a_aux_source) = params.get_a(
        input_assignment.len(), 
        a_aux_density_total
    )?;

    let a_inputs = multiexp(
        &worker,
        a_inputs_source,
        FullDensity,
        input_assignment.clone(),
    );
    let a_aux = multiexp(
        &worker,
        a_aux_source,
        Arc::new(prover.density.a_aux),
        aux_assignment.clone(),
    );

    let b_input = Arc::new(prover.density.b_input);
    let b_input_total = b_input.get_total_density();
    let b_aux_density = Arc::new(prover.density.b_aux);
    let b_aux_density_total = b_aux_density.get_total_density();

    let (b_g1_inputs_source, b_g1_aux_source) =
        params.get_b_g1(b_input_total, b_aux_density_total)?;

    let b_g1_inputs = multiexp(
        &worker,
        b_g1_inputs_source,
        b_input.clone(),
        input_assignment.clone(),
    );
    let b_g1_aux = multiexp(
        &worker,
        b_g1_aux_source,
        b_aux_density.clone(),
        aux_assignment.clone(),
    );

    let (b_g2_inputs_source, b_g2_aux_source) =
        params.get_b_g2(b_input_total, b_aux_density_total)?;

    let b_g2_inputs = multiexp(
        &worker,
        b_g2_inputs_source,
        b_input,
        input_assignment,
    );
    let b_g2_aux = multiexp(&worker, b_g2_aux_source, b_aux_density, aux_assignment);

    if vk.delta_g1.is_zero() || vk.delta_g2.is_zero() {
        // If this element is zero, someone is trying to perform a
        // subversion-CRS attack.
        return Err(SynthesisError::UnexpectedIdentity);
    }

    let mut g_a = vk.delta_g1.mul(r);
    g_a.add_assign_mixed(&vk.alpha_g1);
    let mut g_b = vk.delta_g2.mul(s);
    g_b.add_assign_mixed(&vk.beta_g2);
    let mut g_c;
    {
        let mut rs = r;
        rs.mul_assign(&s);

        g_c = vk.delta_g1.mul(rs);
        g_c.add_assign(&vk.alpha_g1.mul(s));
        g_c.add_assign(&vk.beta_g1.mul(r));
    }
    let mut a_answer = a_inputs.wait()?;
    a_answer.add_assign(&a_aux.wait()?);
    g_a.add_assign(&a_answer);
    a_answer.mul_assign(s);
    g_c.add_assign(&a_answer);

    let mut b1_answer = b_g1_inputs.wait()?;
    b1_answer.add_assign(&b_g1_aux.wait()?);
    let mut b2_answer = b_g2_inputs.wait()?;
    b2_answer.add_assign(&b_g2_aux.wait()?);

    g_b.add_assign(&b2_answer);
    b1_answer.mul_assign(r);
    g_c.add_assign(&b1_answer);
    g_c.add_assign(&h.wait()?);
    g_c.add_assign(&l.wait()?);

    Ok(Proof {
        a: g_a.into_affine(),
        b: g_b.into_affine(),
        c: g_c.into_affine(),
    })
}
