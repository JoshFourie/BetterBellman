use ff::ScalarEngine;
use std::marker;

use super::{ConstraintSystem, domain, error};
use error::Result;
use domain::{Coefficient, LinearCombination};

/// This is a "namespaced" constraint system which borrows a constraint system (pushing
/// a namespace context) and, when dropped, pops out of the namespace context.

pub struct Namespace<'a,E,CS> 
where
    E: ScalarEngine,
    CS: ConstraintSystem<E>
{
    cs: &'a mut CS,
    _marker: marker::PhantomData<E>
}

impl<'a,E,CS> Namespace<'a,E,CS> 
where
    E: ScalarEngine,
    CS: ConstraintSystem<E>
{
    pub fn new(cs: &'a mut CS) -> Self {
        Self {
            cs,
            _marker: marker::PhantomData
        }
    }
}

impl<'a,E,CS> ConstraintSystem<E> for Namespace<'a,E,CS> 
where
    E: ScalarEngine,
    CS: ConstraintSystem<E>
{
    type Root = CS::Root;

    fn one() -> Coefficient {
        CS::one()
    }

    fn alloc<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Coefficient>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        self.cs.alloc(annotation, f)
    }

    fn alloc_input<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Coefficient>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        self.cs.alloc_input(annotation, f)
    }

    fn enforce<A, AR, LA, LB, LC>(&mut self, annotation: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LB: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LC: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
    {
        self.cs.enforce(annotation, a, b, c)
    }

    // Downstream users who use `namespace` will never interact with these
    // functions and they will never be invoked because the namespace is
    // never a root constraint system.
    fn push_namespace<NR, N>(&mut self, _: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        panic!("only the root's push_namespace should be called");
    }

    fn pop_namespace(&mut self) {
        panic!("only the root's pop_namespace should be called");
    }

    fn get_root(&mut self) -> &mut Self::Root {
        self.cs.get_root()
    }
}

impl<'a,E,CS> Drop for Namespace<'a,E,CS> 
where
    E: ScalarEngine,
    CS: ConstraintSystem<E>
{
    fn drop(&mut self) {
        self.get_root().pop_namespace()
    }
}

