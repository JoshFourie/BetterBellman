use crate::multicore::MULTI_THREAD;
use futures::Future;
use group::{CurveAffine, CurveProjective};
use std::sync::Arc;

use crate::error::SynthesisError;
use super::{SourceBuilder, QueryDensity, Exponents, RegionCounter, SourceIter};

pub fn multiexp_inner<Q,D,G,S>(bases: S, density_map: D, exponents: Arc<Exponents<G>>, mut rc: RegionCounter) -> Box<dyn Future<Item=G::Projective, Error=SynthesisError>>
where
    for<'a> &'a Q: QueryDensity,
    D: Send + Sync + 'static + Clone + AsRef<Q>,
    G: CurveAffine,
    S: SourceBuilder<G> + Send,
{
    let this_region: _ = {
        let exponents = exponents.clone();
        let density_map = density_map.clone();
        let bases: _ = bases.clone();

        MULTI_THREAD.compute(move || {
            let mut bases: SourceIter<_> = bases.new();
            bases.configure(rc);

            let mut buckets = vec![<G as CurveAffine>::Projective::zero(); (1 << rc.get_cpu()) - 1];
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
    
    rc.next_region();
    if rc.still_more_regions::<G>() {
        Box::new(
            this_region.join(
                multiexp_inner(bases, density_map, exponents, rc)
            ).map(move |(this, mut higher)| {
                for _ in 0..rc.get_cpu() {
                    higher.double();
                }
                higher.add_assign(&this);
                higher
            }),
        )
    } else { Box::new(this_region) }
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
