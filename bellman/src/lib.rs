#![feature(try_trait)]

#[cfg(feature = "multicore")]
extern crate crossbeam;

#[cfg(feature = "multicore")]
extern crate num_cpus;

#[cfg(test)]
#[macro_use]
extern crate hex_literal;

#[cfg(test)]
extern crate rand;

#[cfg(feature = "groth16")] 
pub mod groth16;

pub mod domain;
pub mod gadgets;
pub mod error;
pub mod namespace;
pub mod constraint;
pub mod multicore;

pub use error::{Result, SynthesisError};
pub use namespace::Namespace;
pub use constraint::ConstraintSystem;
pub use domain::*;

/// Computations are expressed in terms of arithmetic circuits, in particular
/// rank-1 quadratic constraint systems. The `Circuit` trait represents a
/// circuit that can be synthesized. The `synthesize` method is called during
/// CRS generation and during proving.
pub trait Circuit<E> 
where
    E: ff::ScalarEngine
{
    /// Synthesize the circuit into a rank-1 quadratic constraint system
    fn synthesize<CS>(self, cs: &mut CS) -> Result<()>
    where
        CS: ConstraintSystem<E>;
}
