use super::multicore::MULTI_THREAD;
use bit_vec::{self, BitVec};
use ff::{Field, PrimeField, PrimeFieldRepr, ScalarEngine};
use futures::Future;
use group::{CurveAffine, CurveProjective};
use std::io;
use std::iter;
use std::sync::Arc;

use crate::error::{SynthesisError, Result};

/// An object that builds a source of bases.
pub trait SourceBuilder<G: CurveAffine>: Send 
    + 'static
    + Sync 
    + Clone 
{
    fn new<'a>(&'a self) -> SourceIter<'_,G>; 
}

/// A source of bases, like an iterator.
pub trait Source<G> 
where
    G: CurveAffine
{
    /// Parses the element from the source. Fails if the point is at infinity.
    fn add_assign_mixed(&mut self, to: &mut G::Projective) -> Result<()>;

    /// Skips `amt` elements from the source, avoiding deserialization.
    fn skip(&mut self, amt: usize) -> Result<()>;
}

impl<G> SourceBuilder<G> for (Arc<Vec<G>>, usize) 
where
    G: CurveAffine
{
    fn new(&self) -> SourceIter<'_,G> {
        SourceIter::new(&self.0, self.1)
        // self.clone()
    }
}

impl<G> Source<G> for (Arc<Vec<G>>, usize) 
where
    G: CurveAffine
{
    fn add_assign_mixed(&mut self, to: &mut G::Projective) -> Result<()> {
        if self.0.len() > self.1 {
            to.add_assign_mixed(&self.0[self.1]);
            self.1 += 1;
            Ok(())
        } else if self.0[self.1].is_zero() {
            Err(SynthesisError::UnexpectedIdentity)
        } else {
            Err(io::Error::new(
                io::ErrorKind::UnexpectedEof, 
                "expected more bases from source"
            ).into())
        }
    }

    fn skip(&mut self, amt: usize) -> Result<()> {
        if self.0.len() > self.1 {
            self.1 += amt;
            Ok(())
        } else {
            Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "expected more bases from source",
            ).into())
        } 
    }
}

pub trait QueryDensity: IntoIterator<Item = bool> {
    fn get_query_size(self) -> Option<usize>;
}

#[derive(Clone)]
pub struct FullDensity;

impl AsRef<FullDensity> for FullDensity {
    fn as_ref(&self) -> &FullDensity {
        self
    }
}

impl<'a> QueryDensity for &'a FullDensity {
    fn get_query_size(self) -> Option<usize> {
        None
    }
}

impl<'a> IntoIterator for &'a FullDensity {
    type Item = bool;
    type IntoIter = iter::Repeat<bool>;

    fn into_iter(self) -> Self::IntoIter {
        iter::repeat(true)
    }
}

pub struct DensityTracker {
    bv: BitVec,
    total_density: usize,
}

impl<'a> QueryDensity for &'a DensityTracker {
    fn get_query_size(self) -> Option<usize> {
        Some(self.bv.len())
    }
}

impl<'a> IntoIterator for &'a DensityTracker {
    type Item = bool;
    type IntoIter = bit_vec::Iter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.bv.iter()
    }
}

impl DensityTracker {
    pub fn new() -> DensityTracker {
        DensityTracker {
            bv: BitVec::new(),
            total_density: 0,
        }
    }

    pub fn add_element(&mut self) {
        self.bv.push(false);
    }

    pub fn inc(&mut self, idx: usize) {
        if !self.bv.get(idx).unwrap() {
            self.bv.set(idx, true);
            self.total_density += 1;
        }
    }

    pub fn get_total_density(&self) -> usize {
        self.total_density
    }
}

type Exponents<G: CurveAffine> = Vec<<<G::Engine as ScalarEngine>::Fr as PrimeField>::Repr>;

fn multiexp_inner<Q,D,G,S>(
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
            // Accumulate the result
            let mut acc = G::Projective::zero();

            // Create space for the buckets
            let mut buckets = vec![<G as CurveAffine>::Projective::zero(); (1 << c) - 1];

            let zero = <G::Engine as ScalarEngine>::Fr::zero().into_repr();
            let one = <G::Engine as ScalarEngine>::Fr::one().into_repr();

            let bases: _ = bases.new();

            // Sort the bases into buckets
            for ((&exp, density), base) in exponents.iter()
            // for (&exp, density) in exponents.iter()
                .zip(density_map.as_ref()) 
                .zip(bases)
            {
                if density {
                    if exp == zero {
                        // <SourceIter<G> as Source<G>>::skip(&mut bases, 1);
                    } else if exp == one {
                        if handle_trivial {
                            acc.add_assign_mixed(base)
                            // bases.add_assign_mixed(&mut acc)?;
                        } else {
                            // <SourceIter<G> as Source<G>>::skip(&mut bases, 1);
                        }
                    } else {
                        let mut exp = exp;
                        exp.shr(skip);
                        let exp = exp.as_ref()[0] % (1 << c);

                        if exp != 0 {
                            let bucket: _ = &mut buckets[(exp - 1) as usize];
                            // bases.add_assign_mixed(bucket)?;
                            bucket.add_assign_mixed(base);
                        } else {

                        }
                    }
                }
            }

            // Summation by parts
            // e.g. 3a + 2b + 1c = a +
            //                    (a) + b +
            //                    ((a) + b) + c
            let mut running_sum = G::Projective::zero();
            for exp in buckets.into_iter().rev() {
                running_sum.add_assign(&exp);
                acc.add_assign(&running_sum);
            }

            Ok(acc)
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
            ))
            .map(move |(this, mut higher)| {
                for _ in 0..c {
                    higher.double();
                }

                higher.add_assign(&this);

                higher
            }),
        )
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
    let c = if exponents.len() < 32 {
        3u32
    } else {
        (f64::from(exponents.len() as u32)).ln().ceil() as u32
    };

    if let Some(query_size) = density_map.as_ref().get_query_size() {
        // If the density map has a known query size, it should not be
        // inconsistent with the number of exponents.

        assert!(query_size == exponents.len());
    }

    multiexp_inner(bases, density_map, exponents, 0, c, true)
}

pub struct SourceIter<'a,G> {
    elements: &'a Arc<Vec<G>>,
    _count: usize,
}

impl<'a,G> SourceIter<'a,G> {
    fn new(elements: &'a Arc<Vec<G>>, _count: usize) -> Self {
        Self { 
            elements, 
            _count
        }
    }
}

impl<'a,G> Source<G> for SourceIter<'a,G> 
where
    G: CurveAffine
{
    fn add_assign_mixed(&mut self, to: &mut G::Projective) -> Result<()> {
        match self.next() {
            Some(item) => {
                if !item.is_zero() {
                    to.add_assign_mixed(&item);
                    Ok(())
                } else { Err(SynthesisError::UnexpectedIdentity) }
            },
            None => {
                Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof, 
                    "expected more bases from source"
                ).into())
            }
        }
    }

    fn skip(&mut self, amt: usize) -> Result<()> {
        self._count += amt;
        Ok(())
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

/// A source is 'like an iterator' is stupid and I need to check its behaviour
#[test]
fn temp_compare_iterator() {
    use pairing::{bls12_381::Bls12, Engine};
    use rand;

    const SAMPLES: usize = 10;
    let rng = &mut rand::thread_rng();
    
    let mut tuple_source: (Arc<Vec<_>>, usize) = {
        let g = Arc::new(
            (0..SAMPLES)
                .map(|_| <Bls12 as Engine>::G1::random(rng).into_affine())
                .collect::<Vec<_>>(),
        );
        (g, 0)
    };

    let cloned_source: _ = tuple_source.clone();
    let iter_source: _ = cloned_source.new();

    for base in iter_source {
        assert_eq!(base, &tuple_source.0[tuple_source.1]);
        tuple_source.1 += 1;
    }
}
