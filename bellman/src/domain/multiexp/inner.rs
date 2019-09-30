use crate::multicore::MULTI_THREAD;
use ff::{Field, PrimeField, PrimeFieldRepr, ScalarEngine};
use futures::Future;
use group::{CurveAffine, CurveProjective};
use std::io;
use std::sync::Arc;

use crate::error::{SynthesisError, Result};
use super::{SourceBuilder, QueryDensity, Exponents};

pub struct SourceIter<'a,G> {
    elements: &'a Arc<Vec<G>>,
    _count: usize,
}

impl<'a,G> SourceIter<'a,G> {
    pub fn new(elements: &'a Arc<Vec<G>>, _count: usize) -> Self {
        Self { 
            elements, 
            _count
        }
    }

    fn skip_forward(&mut self, amt: usize) {
        self._count += amt;
    }
}

impl<'a,G> Iterator for SourceIter<'a,G> {
    type Item = &'a G;

    fn next(&mut self) -> Option<Self::Item> {
        if self._count < self.elements.len() {
            let item: _ = &self.elements[self._count];
            self._count += 1;
            Some(item)
        } else { None }
    }
}

pub fn multiexp_inner<Q,D,G,S>(
    bases: S,
    density_map: D,
    exponents: Arc<Exponents<G>>,
    mut skip: u32,
    c: u32,
    handle_trivial: bool,
) -> Box<dyn Future<Item = G::Projective, Error = SynthesisError>>
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

            // Create space for the buckets
            let mut buckets = vec![<G as CurveAffine>::Projective::zero(); (1 << c) - 1];

            let ref zero: <G::Scalar as ff::PrimeField>::Repr = <G::Engine as ScalarEngine>::Fr::zero().into_repr();
            let ref one: <G::Scalar as ff::PrimeField>::Repr = <G::Engine as ScalarEngine>::Fr::one().into_repr();

            let mut bases: SourceIter<_> = bases.new();
            let density_iter: _ = density_map.as_ref();
            let mut forward_total: _ = exponents.iter()
                .zip(density_iter)
                .fold(Ok(G::Projective::zero()), |wrapped_acc: Result<G::Projective>, (exp, density)| {
                    wrapped_acc.and_then(|mut acc| {
                        if density {
                            if exp == zero || (exp == one && !handle_trivial) {
                                bases.skip_forward(1)
                            } else if exp == one && handle_trivial {
                                let base: Option<&G> = bases.next();
                                try_add_assign_mixed(&mut acc, base)?
                            } else {
                                let exponent: _ = {
                                    let mut buf: _ = exp.clone();
                                    buf.shr(skip);
                                    let exp_inner: _ = buf.as_ref();
                                    exp_inner[0] % (1 << c)
                                };
                                if exponent != 0 {
                                    let bucket: _ = &mut buckets[(exponent - 1) as usize];
                                    let base: Option<&G> = bases.next();
                                    try_add_assign_mixed(bucket, base)?
                                } else {
                                    bases.skip_forward(1)
                                }
                            }
                        }
                        Ok(acc)
                    })
                })?;

            // Summation by parts
            // e.g. 3a + 2b + 1c = a +
            //                    (a) + b +
            //                    ((a) + b) + c
            let mut running_sum = G::Projective::zero();
            for exp in buckets.iter()
                .rev() 
            {
                running_sum.add_assign(exp);
                forward_total.add_assign(&running_sum);
            }

            Ok(forward_total)
        })
    };

    skip += c;

    if skip >= <G::Engine as ScalarEngine>::Fr::NUM_BITS {
        // There isn't another region.
        Box::new(this)
    } else {
        // There's another region more significant. Calculate and join it with
        // this region recursively.
        Box::new(
            this.join(multiexp_inner(
                bases,
                density_map,
                exponents,
                skip,
                c,
                false,
            )).map(move |(this, mut higher)| {
                for _ in 0..c {
                    higher.double();
                }

                higher.add_assign(&this);

                higher
            }),
        )
    }
}

fn try_add_assign_mixed<G>(lhs: &mut G::Projective, rhs: Option<&G>) -> Result<()> 
where
    G: CurveAffine
{
    match rhs {
        Some(item) => {
            if item.is_zero() {
                Err(SynthesisError::UnexpectedIdentity)
            } else {
                lhs.add_assign_mixed(item);
                Ok(())
            }
        },
        None => Err(io::Error::new(io::ErrorKind::UnexpectedEof, "expected more bases from source").into())
    }
}
