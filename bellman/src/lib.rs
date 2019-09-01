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

pub mod domain;
pub mod gadgets;
pub mod error;
pub mod linear;
pub mod namespace;
#[cfg(feature = "groth16")] pub mod groth16;
pub mod multicore;
mod multiexp;

use ff::{ScalarEngine};

use std::marker::PhantomData;

pub use error::{Result, SynthesisError};
pub use linear::{Variable, LinearCombination, Index};
pub use namespace::Namespace;

/// Computations are expressed in terms of arithmetic circuits, in particular
/// rank-1 quadratic constraint systems. The `Circuit` trait represents a
/// circuit that can be synthesized. The `synthesize` method is called during
/// CRS generation and during proving.
pub trait Circuit<E> 
where
    E: ScalarEngine
{
    /// Synthesize the circuit into a rank-1 quadratic constraint system
    fn synthesize<CS>(self, cs: &mut CS) -> Result<()>
    where
        CS: ConstraintSystem<E>;
}

/// Represents a constraint system which can have new variables
/// allocated and constrains between them formed.
pub trait ConstraintSystem<E>: Sized 
where
    E: ScalarEngine
{
    /// Represents the type of the "root" of this constraint system
    /// so that nested namespaces can minimize indirection.
    type Root: ConstraintSystem<E>;

    /// Return the "one" input variable
    fn one() -> Variable {
        Variable::new_unchecked(Index::Input(0))
    }

    /// Allocate a private variable in the constraint system. The provided function is used to
    /// determine the assignment of the variable. The given `annotation` function is invoked
    /// in testing contexts in order to derive a unique name for this variable in the current
    /// namespace.
    fn alloc<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>;

    /// Allocate a public variable in the constraint system. The provided function is used to
    /// determine the assignment of the variable.
    fn alloc_input<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>;

    /// Enforce that `A` * `B` = `C`. The `annotation` function is invoked in testing contexts
    /// in order to derive a unique name for the constraint in the current namespace.
    fn enforce<A, AR, LA, LB, LC>(&mut self, annotation: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LB: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LC: FnOnce(LinearCombination<E>) -> LinearCombination<E>;

    /// Create a new (sub)namespace and enter into it. Not intended
    /// for downstream use; use `namespace` instead.
    fn push_namespace<NR, N>(&mut self, name_fn: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR;

    /// Exit out of the existing namespace. Not intended for
    /// downstream use; use `namespace` instead.
    fn pop_namespace(&mut self);

    /// Gets the "root" constraint system, bypassing the namespacing.
    /// Not intended for downstream use; use `namespace` instead.
    fn get_root(&mut self) -> &mut Self::Root;

    /// Begin a namespace for this constraint system.
    fn namespace<'a, NR, N>(&'a mut self, name_fn: N) -> Namespace<'a, E, Self::Root>
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        self.get_root().push_namespace(name_fn);

        Namespace::new(self.get_root())
    }
}


/// Convenience implementation of ConstraintSystem<E> for mutable references to
/// constraint systems.
impl<'cs, E: ScalarEngine, CS: ConstraintSystem<E>> ConstraintSystem<E> for &'cs mut CS {
    type Root = CS::Root;

    fn one() -> Variable {
        CS::one()
    }

    fn alloc<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        (**self).alloc(annotation, f)
    }

    fn alloc_input<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        (**self).alloc_input(annotation, f)
    }

    fn enforce<A, AR, LA, LB, LC>(&mut self, annotation: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LB: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LC: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
    {
        (**self).enforce(annotation, a, b, c)
    }

    fn push_namespace<NR, N>(&mut self, name_fn: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        (**self).push_namespace(name_fn)
    }

    fn pop_namespace(&mut self) {
        (**self).pop_namespace()
    }

    fn get_root(&mut self) -> &mut Self::Root {
        (**self).get_root()
    }
}
