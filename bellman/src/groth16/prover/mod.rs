use futures::Future;
use rand_core::RngCore;

use ff::{Field};
use group::{CurveProjective};
use pairing::Engine;

use super::{ParameterSource, Proof, Result};

use crate::{Circuit, ConstraintSystem, Index, SynthesisError, Variable};

mod system;
use system::*;

pub fn create_random_proof<E,C,R,P>(circuit: C, params: P, rng: &mut R) -> Result<Proof<E>>
where
    E: Engine,
    C: Circuit<E>,
    P: ParameterSource<E>,
    R: RngCore,
{
    let r = E::Fr::random(rng);
    let s = E::Fr::random(rng);

    create_proof::<E, C, P>(circuit, params, r, s)
}

pub fn create_proof<E, C, P>(circuit: C, mut params: P, r: E::Fr, s: E::Fr) -> Result<Proof<E>>
where
    E: Engine,
    C: Circuit<E>,
    P: ParameterSource<E>
{
    let mut prover: _ = ProvingSystem::default();
    prover.alloc_input(
        || "", 
        || Ok(E::Fr::one())
    )?;
    circuit.synthesize(&mut prover)?;
    
    let (ga,gb,gc): _ = prover.prepare(&mut params, r, s)?.try_build()?;

    Ok(Proof {
        a: ga.into_affine(),
        b: gb.into_affine(),
        c: gc.into_affine(),
    })
}
