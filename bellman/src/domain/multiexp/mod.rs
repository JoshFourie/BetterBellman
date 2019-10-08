use group::CurveAffine;
use ff::{PrimeField, ScalarEngine};
use std::sync::Arc;
use futures::Future;

use crate::error::SynthesisError;

mod density;
mod inner;
mod region;
mod source;

pub use density::*;
use source::SourceIter;
use region::RegionCounter;

type Exponents<G> = Vec<<<<G as CurveAffine>::Engine as ScalarEngine>::Fr as PrimeField>::Repr>;

/// Perform multi-exponentiation. The thread will panic if the
/// query size is the not the same as the number of exponents.
pub fn multiexp<Q,D,G,S>(bases: S, density_map: D, exponents: Arc<Exponents<G>>) -> Box<dyn Future<Item=G::Projective, Error=SynthesisError>>
where
    for<'a> &'a Q: QueryDensity,
    D: Send + Sync + 'static + Clone + AsRef<Q>,
    G: CurveAffine,
    S: SourceBuilder<G>,
{
    let region: _ = RegionCounter::try_new::<G,Q>(&exponents, density_map.as_ref())
        .expect("could not build region for multi-exponentiation");     
    inner::multiexp_inner(bases, density_map, exponents, region)
}

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

#[cfg(feature = "pairing")]
#[test]
fn test_with_bls12() {
    use ff::{Field, PrimeField};
    use group::CurveProjective;

    fn naive_multiexp<G: CurveAffine>(bases: Arc<Vec<G>>, exponents: Arc<Exponents<G>>) -> G::Projective {
        assert_eq!(bases.len(), exponents.len());
        let mut acc = G::Projective::zero();
        for (base, exp) in bases.iter()
            .zip(exponents.iter()) 
        {
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
