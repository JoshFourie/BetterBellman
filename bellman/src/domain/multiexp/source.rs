use crate::multicore::MULTI_THREAD;
use ff::{Field, PrimeField, PrimeFieldRepr, ScalarEngine};
use futures::Future;
use group::{CurveAffine, CurveProjective};
use std::io;
use std::sync::Arc;

use crate::error;
use error::{SynthesisError, Result};
use super::{SourceBuilder, QueryDensity, Exponents, MultiExpSettings};

pub struct SourceIter<'a,G> {
    elements: &'a Arc<Vec<G>>,
    settings: MultiExpSettings,
    _count: usize,
}

impl<'a,G> SourceIter<'a,G> {
    pub fn new(elements: &'a Arc<Vec<G>>, _count: usize) -> Self {
        Self { 
            elements,
            settings: MultiExpSettings::default(),
            _count
        }
    }

    pub fn skip_forward(&mut self, amt: usize) {
        self._count += amt;
    }

    pub fn settings(&mut self, settings: MultiExpSettings) {
        self.settings = settings
    }
}

impl<'a,G> SourceIter<'a,G> 
where
    G: CurveAffine
{
    pub fn try_sort(&mut self, mut acc: G::Projective, buckets: &mut Vec<G::Projective>, exp: &<G::Scalar as ff::PrimeField>::Repr) -> Result<G::Projective> {
        let settings: _ = &mut self.settings;
        let ref zero: _ = Self::repr_zero();
        let ref one: _ = Self::repr_one();

        if exp == zero || (exp == one && !settings.handle_trivial()) {
            self.skip_forward(1)
        } else if exp == one && settings.handle_trivial() {
            try_add_assign_mixed(&mut acc, self)?
        } else {
            self.try_into_bucket(buckets, exp)?
        }
        Ok(acc)
    }

    fn try_into_bucket(&mut self, buckets: &mut Vec<G::Projective>, exp: &<G::Scalar as ff::PrimeField>::Repr) -> Result<()> {
        let adjustment_source: _ = exp.clone();
        let adjusted_exponent: _ = self.settings.adjust_for_multithreading::<G>(adjustment_source);

        if adjusted_exponent != 0 {
            let bucket: _ = &mut buckets[(adjusted_exponent - 1) as usize];
            try_add_assign_mixed(bucket, self)?
        } else {
            self.skip_forward(1)
        };
        Ok(())
    }

    fn repr_zero() -> <G::Scalar as ff::PrimeField>::Repr {
        <G::Engine as ScalarEngine>::Fr::zero().into_repr()
    }

    fn repr_one() -> <G::Scalar as ff::PrimeField>::Repr {
        <G::Engine as ScalarEngine>::Fr::one().into_repr()
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

fn try_add_assign_mixed<G>(lhs: &mut G::Projective, bases: &mut SourceIter<'_,G>) -> Result<()> 
where
    G: CurveAffine
{
    bases.next()
        .ok_or(io::Error::new(io::ErrorKind::UnexpectedEof, "expected more bases from source").into())
        .and_then(|base| {
            if !base.is_zero() {
                lhs.add_assign_mixed(base);
                Ok(())
            } else { Err(SynthesisError::UnexpectedIdentity) }
        })
}
