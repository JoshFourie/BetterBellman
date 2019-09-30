use std::sync::Arc;
use super::{Future, SynthesisError};

use ff::{Field, PrimeField, ScalarEngine};
use pairing::Engine;

use super::{ParameterSource, Result};

use crate::{ConstraintSystem, Index, LinearCombination, Coefficient};
use crate::{arith, multiexp};
use arith::Scalar;
use multiexp::DensityTracker;

mod builder;
mod source;
mod fourier;

type AssignmentField<E> = Arc<Vec<<<E as ScalarEngine>::Fr as PrimeField>::Repr>>;

pub struct ProvingSystem<E: Engine> {
    density: QueryDensity,
    eval: PolynomialEvaluation<E>,
    pub assignment: ProvingAssignment<E>
}

impl<E: Engine> ProvingSystem<E> {
    pub fn prepare<T>(mut self, params: &mut T, r: E::Fr, s: E::Fr) -> Result<builder::Builder<E>>
    where
        T: ParameterSource<E>
    {
        for i in 0..self.assignment.input.len() {
            self.enforce(
                || "", 
                |lc| lc + Coefficient::new_unchecked(Index::Input(i)), 
                |lc| lc, 
                |lc| lc
            );
        }

        builder::Builder::try_new(self, params, r, s)
    }
}

impl<E> ConstraintSystem<E> for ProvingSystem<E> 
where
    E: Engine
{
    type Root = Self;

    fn alloc<F, A, AR>(&mut self, _: A, f: F) -> Result<Coefficient>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        self.assignment.aux.push(f()?);
        self.density.a_aux.add_element();
        self.density.b_aux.add_element();

        Ok(Coefficient::new_unchecked(Index::Aux(self.assignment.aux.len() - 1)))
    }

    fn alloc_input<F, A, AR>(&mut self, _: A, f: F) -> Result<Coefficient>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        self.assignment.input.push(f()?);
        self.density.b_input.add_element();

        let index: _ = Index::Input(self.assignment.input.len() - 1);
        Ok(Coefficient::new_unchecked(index))
    }

    fn enforce<A, AR, LA, LB, LC>(&mut self, _: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LB: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LC: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
    {
        let a = a(LinearCombination::zero());
        let b = b(LinearCombination::zero());
        let c = c(LinearCombination::zero());

        let eval_a: E::Fr = evalaluate_linear_combination(
            &a,
            |index| match index {
                Index::Input(i) => self.assignment.input[i],
                Index::Aux(i) => {
                    self.density.a_aux.inc(i);
                    self.assignment.aux[i]
                }
            }
        );
        self.eval
            .a
            .as_mut()
            .expect("tried calling an operation on Linear Combination A, but found None")
            .push(Scalar(eval_a));

        let eval_b: E::Fr = evalaluate_linear_combination(
            &b,
            |index| match index {
                Index::Input(i) => {
                    self.density.b_input.inc(i);
                    self.assignment.input[i]
                },
                Index::Aux(i) => {
                    self.density.b_aux.inc(i);
                    self.assignment.aux[i]
                }
            }
        );
        self.eval
            .b
            .as_mut()
            .expect("tried calling an operation on Linear Combination B, but found None")
            .push(Scalar(eval_b));

        let eval_c: E::Fr = evalaluate_linear_combination(
            &c,
            |index| match index {
                Index::Input(i) => self.assignment.input[i],
                Index::Aux(i) => self.assignment.aux[i]
            }
        );
        self.eval
            .c
            .as_mut()
            .expect("tried calling an operation on Linear Combination C, but found None")
            .push(Scalar(eval_c));
    }

    fn push_namespace<NR, N>(&mut self, _: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        // Do nothing; we don't care about namespaces in this context.
    }

    fn pop_namespace(&mut self) {
        // Do nothing; we don't care about namespaces in this context.
    }

    fn get_root(&mut self) -> &mut Self::Root {
        self
    }
}

fn evalaluate_linear_combination<E,F>(linear: &LinearCombination<E>, mut func: F) -> E::Fr 
where
    E: Engine,
    F: FnMut(Index) -> E::Fr    
{
    linear.0
        .iter()
        .fold(E::Fr::zero(), |mut acc, (idx, coeff)| {
            let mut buf: _ = func(idx.get_unchecked());
            if coeff != &E::Fr::one() {
                buf.mul_assign(&coeff)
            }
            acc.add_assign(&buf);
            acc 
        })
}

impl<E: Engine> Default for ProvingSystem<E> {
    fn default() -> Self {
        ProvingSystem {
            density: QueryDensity::default(),           
            eval: PolynomialEvaluation::default(),
            assignment: ProvingAssignment::default()
        }
    }
}

// todo: resolve naming collision with object in multiexp.
struct QueryDensity {
    a_aux: DensityTracker,
    b_input: DensityTracker,
    b_aux: DensityTracker,
}

impl Default for QueryDensity {
    fn default() -> Self {
        QueryDensity {
            a_aux: DensityTracker::new(),
            b_input: DensityTracker::new(),
            b_aux: DensityTracker::new()
        }
    }
}

pub struct PolynomialEvaluation<E: Engine> {
    a: Option<Vec<Scalar<E>>>,
    b: Option<Vec<Scalar<E>>>,
    c: Option<Vec<Scalar<E>>>,
}

impl<E: Engine> Default for PolynomialEvaluation<E> {
    fn default() -> Self {
        PolynomialEvaluation {
            a: Some(Vec::new()),
            b: Some(Vec::new()),
            c: Some(Vec::new())
        }
    }
}

pub struct ProvingAssignment<E: Engine> {
    pub input: Vec<E::Fr>,
    pub aux: Vec<E::Fr>,
}

impl<E: Engine> Default for ProvingAssignment<E> {
    fn default() -> Self {
        ProvingAssignment {
            input: Vec::new(),
            aux: Vec::new()
        }
    }
} 
