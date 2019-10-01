use group::{CurveAffine, CurveProjective};
use ff::{PrimeField, Field, ScalarEngine};
use std::sync::Arc;
use futures::Future;

use crate::error::{Result, SynthesisError};

mod density;
mod inner;
mod settings;
mod source;

pub use density::*;
use source::SourceIter;
use settings::MultiExpSettings;

type Exponents<G: CurveAffine> = Vec<<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr>;

/// An object that builds a source of bases.
pub trait SourceBuilder<G: CurveAffine>: Send 
    + 'static
    + Sync 
    + Clone 
{
    fn new<'a>(&'a self) -> SourceIter<'_,G>; 
}

impl<G> SourceBuilder<G> for (Arc<Vec<G>>, usize) 
where
    G: CurveAffine
{
    fn new(&self) -> SourceIter<'_,G> {
        SourceIter::new(&self.0, self.1)
    }
}

/// Perform multi-exponentiation. The caller is responsible for ensuring the
/// query size is the same as the number of exponents.
pub fn multiexp<Q,D,G,S>(bases: S, density_map: D, exponents: Arc<Exponents<G>>) -> Box<dyn Future<Item=G::Projective, Error=SynthesisError>>
where
    for<'a> &'a Q: QueryDensity,
    D: Send + Sync + 'static + Clone + AsRef<Q>,
    G: CurveAffine,
    S: SourceBuilder<G>,
{
    let settings: _ = MultiExpSettings::try_new::<G,Q>(&exponents, density_map.as_ref())
        .expect("could not build settings for multi-exponentiation");     
    inner::multiexp_inner(bases, density_map, exponents, settings)
}

#[cfg(feature = "pairing")]
#[test]
fn test_with_bls12() {
    fn naive_multiexp<G: CurveAffine>(
        bases: Arc<Vec<G>>,
        exponents: Arc<Vec<<G::Scalar as PrimeField>::Repr>>,
    ) -> G::Projective {
        assert_eq!(bases.len(), exponents.len());

        let mut acc = G::Projective::zero();

        for (base, exp) in bases.iter().zip(exponents.iter()) {
            acc.add_assign(&base.mul(*exp));
        }

        acc
    }

    use pairing::{bls12_381::Bls12, Engine};
    use rand;

    const SAMPLES: usize = 1 << 14;

    let rng = &mut rand::thread_rng();
    let v = Arc::new(
        (0..SAMPLES)
            .map(|_| <Bls12 as ScalarEngine>::Fr::random(rng).into_repr())
            .collect::<Vec<_>>(),
    );
    let g = Arc::new(
        (0..SAMPLES)
            .map(|_| <Bls12 as Engine>::G1::random(rng).into_affine())
            .collect::<Vec<_>>(),
    );

    let naive = naive_multiexp(g.clone(), v.clone());

    let fast = multiexp((g, 0), FullDensity, v).wait().unwrap();

    assert_eq!(naive, fast);
}
