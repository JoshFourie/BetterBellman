use bit_vec::{self, BitVec};
use std::iter;

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