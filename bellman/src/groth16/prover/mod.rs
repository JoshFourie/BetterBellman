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
        let a: _ = prover.fft_shortcut(&worker)?;
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

    let b_input_density = Arc::new(prover.density.b_input);
    let b_aux_density = Arc::new(prover.density.b_aux);

    let b_input_total = b_input_density.get_total_density();
    let b_aux_density_total = b_aux_density.get_total_density();

    let (b1_input_src, b1_aux_src): _ = ProvingSystem::src_gb1(
        b_input_density.clone(),
        b_aux_density.clone(), 
        &mut params
    )?;

    let (b2_input_src, b2_aux_src): _ = ProvingSystem::src_gb2(
        b_input_total,
        b_aux_density_total,
        &mut params
    )?;

    if vk.delta_g1.is_zero() || vk.delta_g2.is_zero() {
        // If this element is zero, someone is trying to perform a
        // subversion-CRS attack.
        return Err(SynthesisError::UnexpectedIdentity);
    }
    
    let answer: _ = AnswerContext::try_new::<P>(
        &worker, 
        a_inputs_source,
        b1_input_src,
        b2_input_src,
        b_input_density,
        input_assignment
    )?;

    let aux: _ = AuxiliaryContext::try_new::<P>(
        &worker,
        a_aux_source,
        b1_aux_src,
        b2_aux_src,
        prover.density.a_aux,
        b_aux_density,
        aux_assignment
    )?;

    let builder: _ = GroupBuilder {
        vk: &vk,
        r,
        s,
        answer,
        aux,
        h: h.wait()?,
        l: l.wait()?
    };

    let (ga,gb,gc): _ = builder.try_build()?;

    Ok(Proof {
        a: ga.into_affine(),
        b: gb.into_affine(),
        c: gc.into_affine(),
    })
}
