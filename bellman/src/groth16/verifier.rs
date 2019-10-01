use ff::PrimeField;
use group::{CurveAffine, CurveProjective};
use pairing::{Engine, PairingCurveAffine};

use super::{PreparedVerifyingKey, Proof, VerifyingKey, Result};

use crate::SynthesisError;

pub fn prepare_verifying_key<E>(vk: &VerifyingKey<E>) -> PreparedVerifyingKey<E> 
where
    E: Engine
{
    let mut gamma = vk.gamma_g2;
    gamma.negate();
    let mut delta = vk.delta_g2;
    delta.negate();

    PreparedVerifyingKey {
        alpha_g1_beta_g2: E::pairing(vk.alpha_g1, vk.beta_g2),
        neg_gamma_g2: gamma.prepare(),
        neg_delta_g2: delta.prepare(),
        ic: vk.ic.clone(),
    }
}

pub fn verify_proof<E>(pvk: &PreparedVerifyingKey<E>, proof: &Proof<E>, public_inputs: &[E::Fr]) -> Result<bool> 
where
    E: Engine
{
    if (public_inputs.len() + 1) != pvk.ic.len() {
        return Err(SynthesisError::MalformedVerifyingKey);
    }

    let acc: _ = public_inputs.iter()
        .zip(pvk.ic.iter().skip(1))
        .fold(pvk.ic[0].into_projective(), |mut acc, (i,b)| {
            acc.add_assign( &b.mul(i.into_repr()) );
            acc
        });

    // let mut acc = pvk.ic[0].into_projective();

    // for (i, b) in public_inputs.iter().zip(pvk.ic.iter().skip(1)) {
    //     acc.add_assign(&b.mul(i.into_repr()));
    // }

    // The original verification equation is:
    // A * B = alpha * beta + inputs * gamma + C * delta
    // ... however, we rearrange it so that it is:
    // A * B - inputs * gamma - C * delta = alpha * beta
    // or equivalently:
    // A * B + inputs * (-gamma) + C * (-delta) = alpha * beta
    // which allows us to do a single final exponentiation.
    let extension_field: _ = E::miller_loop(&[
        (&proof.a.prepare(), &proof.b.prepare()),
        (&acc.into_affine().prepare(), &pvk.neg_gamma_g2),
        (&proof.c.prepare(), &pvk.neg_delta_g2),
    ]);
    let exponentiation: _ = E::final_exponentiation(&extension_field)?;

    Ok(exponentiation == pvk.alpha_g1_beta_g2)
}
