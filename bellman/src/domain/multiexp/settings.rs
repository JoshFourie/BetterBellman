use crate::multicore::MULTI_THREAD;
use ff::{Field, PrimeField, PrimeFieldRepr, ScalarEngine};
use futures::Future;
use group::{CurveAffine, CurveProjective};
use std::io;
use std::sync::Arc;

use crate::error;
use error::{SynthesisError, Result};
use super::{SourceBuilder, QueryDensity, Exponents};

#[derive(Copy, Clone)]
pub struct MultiExpSettings {
    skip: u32,
    cpu_count: u32,
    handle_trivial: bool
}

impl MultiExpSettings {
    pub fn try_new<G,Q>(exponents: &Arc<Exponents<G>>, density_query: &Q) -> Result<Self> 
    where
        G: CurveAffine,
        for <'a> &'a Q: QueryDensity
    {
        if let Some(query_size) = density_query.get_query_size() {
            if query_size != exponents.len() {
                Err(io::Error::new(io::ErrorKind::InvalidInput, "exected length of exponents to match query size"))?
            }
        }

        let cpu_count = if exponents.len() < 32 {
            3_u32
        } else {
            (f64::from(exponents.len() as u32)).ln().ceil() as u32
        };

        Ok(Self {
            cpu_count,
            handle_trivial: true,
            skip: 0
        })
    }

    pub fn adjust_for_multithreading<T>(&self, mut exp: <T::Scalar as ff::PrimeField>::Repr) -> u64 
    where
        T: CurveAffine
    {
        let skip: _ = self.get_skip();
        exp.shr(skip);

        let adjusted: _ = exp.as_ref();
        let cpu_count: _ = self.get_cpu();
        let cpu_adjusted: _ = adjusted[0] % (1 << cpu_count);
        cpu_adjusted
    }

    pub fn still_more_regions<U>(&self) -> bool 
    where
        U: CurveAffine
    {
        self.skip < <U::Engine as ScalarEngine>::Fr::NUM_BITS
    }

    pub fn set_handle_trivial(&mut self, setting: bool) {
        self.handle_trivial = setting
    }

    pub fn set_cpu(&mut self, cpu: u32) {
        self.cpu_count = cpu
    }

    pub fn set_skip_by(&mut self, skip: u32) {
        self.skip = skip
    }

    pub fn get_skip(&self) -> u32 {
        self.skip
    }

    pub fn get_cpu(&self) -> u32 {
        self.cpu_count
    }

    pub fn handle_trivial(&self) -> bool {
        self.handle_trivial
    }

    pub fn next_region(&mut self) {
        let new_skip: _ = self.get_skip() + self.get_cpu();
        self.set_handle_trivial(false);
        self.set_skip_by(new_skip)
    }
}

impl Default for MultiExpSettings {
    fn default() -> Self {
        Self {
            handle_trivial: true,
            skip: 0,
            cpu_count: 1
        }
    }
}
