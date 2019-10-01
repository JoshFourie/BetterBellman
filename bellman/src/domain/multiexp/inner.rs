use crate::multicore::MULTI_THREAD;
use ff::{Field, PrimeField, PrimeFieldRepr, ScalarEngine};
use futures::Future;
use group::{CurveAffine, CurveProjective};
use std::io;
use std::sync::Arc;

use crate::error;
use error::{SynthesisError, Result};
use super::{SourceBuilder, QueryDensity, Exponents, MultiExpSettings, SourceIter};

pub fn multiexp_inner<Q,D,G,S>(bases: S, density_map: D, exponents: Arc<Exponents<G>>, mut settings: MultiExpSettings) -> Box<dyn Future<Item=G::Projective, Error=SynthesisError>>
where
    for<'a> &'a Q: QueryDensity,
    D: Send + Sync + 'static + Clone + AsRef<Q>,
    G: CurveAffine,
    S: SourceBuilder<G> + Send,
{
    // Perform this region of the multiexp
    let this = {
        let exponents = exponents.clone();
        let density_map = density_map.clone();
        let bases: _ = bases.clone();

        MULTI_THREAD.compute(move || {
            let mut bases: SourceIter<_> = bases.new();
            bases.settings(settings);

            let mut buckets = vec![<G as CurveAffine>::Projective::zero(); (1 << settings.get_cpu()) - 1];
            let density_iter: _ = density_map.as_ref();
            let mut forward_total: G::Projective = exponents.iter()
                .zip(density_iter)
                .fold(Ok(G::Projective::zero()), |accumulator, (exp, density)| {
                    if density {
                        accumulator.and_then(|acc| bases.try_sort(acc, &mut buckets, exp))
                    } else { accumulator }
                })?;

            add_assign_by_parts::<G>(&mut forward_total, buckets);
            Ok(forward_total)
        })
    };
    
    settings.next_region();
    if settings.still_more_regions::<G>() {
        // There's another region more significant. Calculate and join it with
        // this region recursively.
        Box::new(
            this.join(multiexp_inner(
                bases,
                density_map,
                exponents,
                settings
            )).map(move |(this, mut higher)| {
                for _ in 0..settings.get_cpu() {
                    higher.double();
                }

                higher.add_assign(&this);

                higher
            }),
        )
    } else {
        Box::new(this)
    }
}

// Summation by parts
// e.g. 3a + 2b + 1c = a +
//                    (a) + b +
//                    ((a) + b) + c
fn add_assign_by_parts<G>(lhs: &mut G::Projective, buckets: Vec<G::Projective>) 
where
    G: CurveAffine
{
    let mut sigma: _ = G::Projective::zero();
    for exponent in buckets.iter()
        .rev()
    {
        sigma.add_assign(exponent);
        lhs.add_assign(&sigma)
    }
}
