use ff::Field;
use pairing::Engine;

use crate::{ConstraintSystem, Circuit, Index, LinearCombination, Coefficient};
use crate::{domain, error, arith};
use domain::Domain;
use error::Result;
use arith::Scalar;

pub mod wires;
pub use wires::*;

/// This is our assembly structure that we'll use to synthesize the
/// circuit into a QAP.
pub struct KeyPairAssembly<E: Engine> {
    pub num: KeyPairNum,
    pub inputs: KeyPairWires<E>,
    pub aux: KeyPairWires<E>
}

impl<E> KeyPairAssembly<E>
where
    E: Engine
{
    pub fn allocate_input_one(&mut self) -> Result<()> {
        self.alloc_input(
            || "", 
            || Ok(E::Fr::one())
        )?;
        Ok(())
    }

    pub fn enforce_full_density(&mut self) -> Result<()> {
        for i in 0..self.num.inputs {
            self.enforce(
                || "", 
                |lc| lc + Coefficient::new_unchecked(Index::Input(i)), 
                |lc| lc, 
                |lc| lc
            );
        }
        Ok(())
    }

    pub fn synthesize_circuit<C>(&mut self, circuit: C) -> Result<()>
    where
        C: Circuit<E>
    {
        circuit.synthesize(self)
    }

    pub fn blind_evaluation_base(&self) -> Result<Domain<E,Scalar<E>>> {
        let powers_of_tau = vec![Scalar(E::Fr::zero()); self.num.constraints];
        Domain::new(powers_of_tau)
    } 
}

impl<E> ConstraintSystem<E> for KeyPairAssembly<E> 
where
    E: Engine
{
    type Root = Self;

    fn alloc<F, A, AR>(&mut self, _: A, _: F) -> Result<Coefficient>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        // There is no assignment, so we don't even invoke the
        // function for obtaining one.

        let index = self.num.aux;
        self.num.aux += 1;

        self.aux.at.push(vec![]);
        self.aux.bt.push(vec![]);
        self.aux.ct.push(vec![]);

        Ok(Coefficient::new_unchecked(Index::Aux(index)))
    }

    fn alloc_input<F, A, AR>(&mut self, _: A, _: F) -> Result<Coefficient>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        // There is no assignment, so we don't even invoke the
        // function for obtaining one.

        let index = self.num.inputs;
        self.num.inputs += 1;

        self.inputs.at.push(Vec::new());
        self.inputs.bt.push(Vec::new());
        self.inputs.ct.push(Vec::new());

        Ok(Coefficient::new_unchecked(Index::Input(index)))
    }

    fn enforce<A, AR, LA, LB, LC>(&mut self, _: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LB: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LC: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
    {
        fn eval<E: Engine>(
            l: LinearCombination<E>,
            inputs: &mut [Vec<(E::Fr, usize)>],
            aux: &mut [Vec<(E::Fr, usize)>],
            this_constraint: usize,
        ) {
            for (index, coeff) in l.0 {
                match index.get_unchecked() {
                    Index::Input(id) => inputs[id].push((coeff, this_constraint)),
                    Index::Aux(id) => aux[id].push((coeff, this_constraint)),
                }
            }
        }

        eval(
            a(LinearCombination::zero()),
            &mut self.inputs.at,
            &mut self.aux.at,
            self.num.constraints,
        );
        eval(
            b(LinearCombination::zero()),
            &mut self.inputs.bt,
            &mut self.aux.bt,
            self.num.constraints,
        );
        eval(
            c(LinearCombination::zero()),
            &mut self.inputs.ct,
            &mut self.aux.ct,
            self.num.constraints,
        );

        self.num.constraints += 1;
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

impl<E> Default for KeyPairAssembly<E>
where
    E: Engine
{
    fn default() -> Self {
        KeyPairAssembly {
            num: KeyPairNum::default(),
            inputs: KeyPairWires::default(),
            aux: KeyPairWires::default()
        }
    }
}
