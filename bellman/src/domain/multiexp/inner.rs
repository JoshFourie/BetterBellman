use crate::multicore::MULTI_THREAD;
use ff::{Field, PrimeField, PrimeFieldRepr, ScalarEngine};
use futures::Future;
use group::{CurveAffine, CurveProjective};
use std::io;
use std::sync::Arc;

use crate::error;
use error::{SynthesisError, Result};
use super::{SourceBuilder, QueryDensity, Exponents};

pub fn multiexp_inner<Q,D,G,S>(
    bases: S,
    density_map: D,
    exponents: Arc<Exponents<G>>,
    skip: u32,
    c: u32,
    handle_trivial: bool,
) -> Box<dyn Future<Item = G::Projective, Error = SynthesisError>>
where
    for<'a> &'a Q: QueryDensity,
    D: Send + Sync + 'static + Clone + AsRef<Q>,
    G: CurveAffine,
    S: SourceBuilder<G> + Send,
{
    let mut settings: _ = MultiExpSettings::default();
    settings.set_cpu(c);
    settings.set_skip_by(skip);
    settings.set_handle_trivial(handle_trivial);    

    // Perform this region of the multiexp
    let this = {
        let exponents = exponents.clone();
        let density_map = density_map.clone();
        let bases: _ = bases.clone();

        MULTI_THREAD.compute(move || {
            let mut bases: SourceIter<_> = bases.new();
            bases.settings(settings);

            // Create space for the buckets
            let mut buckets = vec![<G as CurveAffine>::Projective::zero(); (1 << c) - 1];

            let density_iter: _ = density_map.as_ref();
            let mut forward_total: G::Projective = exponents.iter()
                .zip(density_iter)
                .fold(Ok(G::Projective::zero()), |accumulator, (exp, density)| {
                    if density {
                        accumulator.and_then(|acc| bases.try_sort(acc, &mut buckets, exp))
                    } else { accumulator }
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
    
    settings.next_region();
    if settings.still_more_regions::<G>() {
        // There's another region more significant. Calculate and join it with
        // this region recursively.
        Box::new(
            this.join(multiexp_inner(
                bases,
                density_map,
                exponents,
                settings.get_skip(),
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
    } else {
        Box::new(this)
    }
}

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

    fn skip_forward(&mut self, amt: usize) {
        self._count += amt;
    }

    fn settings(&mut self, settings: MultiExpSettings) {
        self.settings = settings
    }
}

impl<'a,G> SourceIter<'a,G> 
where
    G: CurveAffine
{
    fn try_sort(&mut self, mut acc: G::Projective, buckets: &mut Vec<G::Projective>, exp: &<G::Scalar as ff::PrimeField>::Repr) -> Result<G::Projective> {
        let settings: MultiExpSettings = self.settings;
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
        let adjusted_exp: _ = adjust_for_multithreading(exp, self.settings);
        if adjusted_exp != 0 {
            let bucket: _ = &mut buckets[(adjusted_exp - 1) as usize];
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

#[derive(Copy, Clone)]
struct MultiExpSettings {
    skip: u32,
    cpu_count: u32,
    handle_trivial: bool
}

impl MultiExpSettings {
    fn set_handle_trivial(&mut self, setting: bool) {
        self.handle_trivial = setting
    }

    fn set_cpu(&mut self, cpu: u32) {
        self.cpu_count = cpu
    }

    fn set_skip_by(&mut self, skip: u32) {
        self.skip = skip
    }

    fn get_skip(&self) -> u32 {
        self.skip
    }

    fn get_cpu(&self) -> u32 {
        self.cpu_count
    }

    fn handle_trivial(&self) -> bool {
        self.handle_trivial
    }

    fn next_region(&mut self) {
        let new_skip: _ = self.get_skip() + self.get_cpu();
        self.set_skip_by(new_skip)
    }

    fn still_more_regions<G>(&self) -> bool 
    where
        G: CurveAffine
    {
        self.skip < <G::Engine as ScalarEngine>::Fr::NUM_BITS
    }
}

impl Default for MultiExpSettings {
    fn default() -> Self {
        MultiExpSettings {
            skip: 0,
            cpu_count: 1,
            handle_trivial: true
        }
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

fn adjust_for_multithreading<G>(exp: &G, settings: MultiExpSettings) -> u64
where
    G: PrimeFieldRepr
{
    let mut buf: _ = exp.clone();
    let skip: _ = settings.get_skip();
    buf.shr(skip);

    let adjusted: _ = buf.as_ref();
    let cpu_count: _ = settings.get_cpu();
    let cpu_adjusted: _ = adjusted[0] % (1 << cpu_count);
    cpu_adjusted
}
