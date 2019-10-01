use ff::{PrimeField, PrimeFieldRepr, ScalarEngine};
use group::CurveAffine;
use std::io;
use std::sync::Arc;

use crate::error::Result;
use super::{QueryDensity, Exponents};

/// A structure to maintain awareness of the current region
/// so that the multithreading will consider separate parts 
/// of the source of bases.
#[derive(Copy, Clone)]
pub struct RegionCounter {
    count: u32,
    cpu: u32,
    handle_trivial: bool
}

impl RegionCounter {
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

        let cpu = if exponents.len() < 32 {
            3_u32
        } else {
            let casted_size_of_exp: _ = f64::from(exponents.len() as u32);
            let log_n: _ = casted_size_of_exp.ln();
            let casted_cpu: u32 = log_n.ceil() as u32;
            casted_cpu            
        };

        Ok(RegionCounter {
            cpu,
            handle_trivial: true,
            count: 0
        })
    }

    pub fn adjust_exponent_by_region<T>(&self, mut exp: <T::Scalar as ff::PrimeField>::Repr) -> u64 
    where
        T: CurveAffine
    {
        let skip: _ = self.get_count();
        exp.shr(skip);

        let adjusted: _ = exp.as_ref();
        let cpu: _ = self.get_cpu();
        let cpu_adjusted: _ = adjusted[0] % (1 << cpu);
        cpu_adjusted
    }

    pub fn still_more_regions<U>(&self) -> bool 
    where
        U: CurveAffine
    {
        self.count < <U::Engine as ScalarEngine>::Fr::NUM_BITS
    }

    pub fn get_count(&self) -> u32 {
        self.count
    }

    pub fn get_cpu(&self) -> u32 {
        self.cpu
    }

    pub fn handle_trivial(&self) -> bool {
        self.handle_trivial
    }

    pub fn next_region(&mut self) {
        self.handle_trivial = false;
        self.count += self.cpu;
    }
}

impl Default for RegionCounter {
    fn default() -> Self {
        RegionCounter {
            handle_trivial: true,
            count: 0,
            cpu: 1
        }
    }
}
