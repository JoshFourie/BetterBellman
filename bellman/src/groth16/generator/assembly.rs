use ff::Field;
use pairing::Engine;

use crate::{ConstraintSystem, Index, LinearCombination, Variable};
use crate::error::Result;

/// This is our assembly structure that we'll use to synthesize the
/// circuit into a QAP.
pub struct KeypairAssembly<E: Engine> {
    pub num_inputs: usize,
    pub num_aux: usize,
    pub num_constraints: usize,
    pub at_inputs: Vec<Vec<(E::Fr, usize)>>,
    pub bt_inputs: Vec<Vec<(E::Fr, usize)>>,
    pub ct_inputs: Vec<Vec<(E::Fr, usize)>>,
    pub at_aux: Vec<Vec<(E::Fr, usize)>>,
    pub bt_aux: Vec<Vec<(E::Fr, usize)>>,
    pub ct_aux: Vec<Vec<(E::Fr, usize)>>,
}

impl<E> KeypairAssembly<E>
where
    E: Engine
{
    pub fn allocate_one_input(&mut self) -> Result<()>
    where
        E: Engine
    {
        self.alloc_input(
            || "", 
            || Ok(E::Fr::one())
        )?;
        Ok(())
    }

    pub fn enforce_full_density(&mut self) -> Result<()>
    where
        E: Engine
    {
        for i in 0..self.num_inputs {
            self.enforce(
                || "", 
                |lc| lc + Variable::new_unchecked(Index::Input(i)), 
                |lc| lc, 
                |lc| lc
            );
        }
        Ok(())
    }
}

impl<E> ConstraintSystem<E> for KeypairAssembly<E> 
where
    E: Engine
{
    type Root = Self;

    fn alloc<F, A, AR>(&mut self, _: A, _: F) -> Result<Variable>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        // There is no assignment, so we don't even invoke the
        // function for obtaining one.

        let index = self.num_aux;
        self.num_aux += 1;

        self.at_aux.push(vec![]);
        self.bt_aux.push(vec![]);
        self.ct_aux.push(vec![]);

        Ok(Variable::new_unchecked(Index::Aux(index)))
    }

    fn alloc_input<F, A, AR>(&mut self, _: A, _: F) -> Result<Variable>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        // There is no assignment, so we don't even invoke the
        // function for obtaining one.

        let index = self.num_inputs;
        self.num_inputs += 1;

        self.at_inputs.push(vec![]);
        self.bt_inputs.push(vec![]);
        self.ct_inputs.push(vec![]);

        Ok(Variable::new_unchecked(Index::Input(index)))
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
            &mut self.at_inputs,
            &mut self.at_aux,
            self.num_constraints,
        );
        eval(
            b(LinearCombination::zero()),
            &mut self.bt_inputs,
            &mut self.bt_aux,
            self.num_constraints,
        );
        eval(
            c(LinearCombination::zero()),
            &mut self.ct_inputs,
            &mut self.ct_aux,
            self.num_constraints,
        );

        self.num_constraints += 1;
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

impl<E> Default for KeypairAssembly<E>
where
    E: Engine
{
    fn default() -> Self {
        KeypairAssembly {
            num_inputs: 0,
            num_aux: 0,
            num_constraints: 0,
            at_inputs: vec![],
            bt_inputs: vec![],
            ct_inputs: vec![],
            at_aux: vec![],
            bt_aux: vec![],
            ct_aux: vec![],
        }
    }
}
