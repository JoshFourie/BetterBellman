use rand_core::RngCore;

use std::sync::Arc;

use ff::{Field, PrimeField};
use pairing::Engine;

use super::{ParameterSource, Proof, Result};

use crate::{Circuit, ConstraintSystem, Index, LinearCombination, Variable};
use crate::domain::{EvaluationDomain, Scalar};
use crate::multiexp::{multiexp, DensityTracker, FullDensity};
use crate::multicore::Worker;
use crate::groth16::VerifyingKey;
use group::CurveAffine;

pub struct ProvingSystem<E: Engine> {
    pub density: QueryDensity,
    pub eval: PolynomialEvaluation<E>,
    pub assignment: ProvingAssignment<E>
}

impl<E: Engine> ProvingSystem<E> {
    fn eval<F>(linear: &LinearCombination<E>, mut func: F) -> E::Fr 
    where
        F: FnMut(Index) -> E::Fr
    {
        linear.0
            .iter()
            .fold(E::Fr::zero(), |mut acc, (idx, coeff)| {
                let mut buf: _ = func(idx.0);
                if coeff != &E::Fr::one() {
                    buf.mul_assign(&coeff)
                }
                acc.add_assign(&buf);
                acc 
            })
    }
}

impl<E> ConstraintSystem<E> for ProvingSystem<E> 
where
    E: Engine
{
    type Root = Self;

    fn alloc<F, A, AR>(&mut self, _: A, f: F) -> Result<Variable>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        self.assignment.aux.push(f()?);
        self.density.a_aux.add_element();
        self.density.b_aux.add_element();

        Ok(Variable(Index::Aux(self.assignment.aux.len() - 1)))
    }

    fn alloc_input<F, A, AR>(&mut self, _: A, f: F) -> Result<Variable>
    where
        F: FnOnce() -> Result<E::Fr>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        self.assignment.input.push(f()?);
        self.density.b_input.add_element();

        Ok(Variable(Index::Input(self.assignment.input.len() - 1)))
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

        let eval_a: E::Fr = ProvingSystem::eval(
            &a,
            |index| match index {
                Index::Input(i) => self.assignment.input[i],
                Index::Aux(i) => {
                    self.density.a_aux.inc(i);
                    self.assignment.aux[i]
                }
            }
        );
        self.eval.a
            .as_mut()
            .unwrap()
            .push(Scalar(eval_a));

        let eval_b: E::Fr = ProvingSystem::eval(
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
            .unwrap()
            .push(Scalar(eval_b));

        let eval_c: E::Fr = ProvingSystem::eval(
            &c,
            |index| match index {
                Index::Input(i) => self.assignment.input[i],
                Index::Aux(i) => self.assignment.aux[i]
            }
        );
        self.eval
            .c
            .as_mut()
            .unwrap()
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

impl<E: Engine> ProvingSystem<E> {
    pub fn random<C, R, P>(circuit: C, params: P, rng: &mut R) -> Result<Proof<E>>
    where
        C: Circuit<E>,
        P: ParameterSource<E>,
        R: RngCore,
    {
        let r = E::Fr::random(rng);
        let s = E::Fr::random(rng);

        super::create_proof::<E, C, P>(circuit, params, r, s)
    }

    pub fn primefield(&mut self, worker: &Worker) ->  Result<Arc<Vec<<<E as ff::ScalarEngine>::Fr as PrimeField>::Repr>>> {
        let field: _ = ProvingField::new(&mut self.eval, worker)?;
        let mut a: _ = field.transform().into_coefficient();

        let new_len = a.len() - 1;
        a.truncate(new_len);

        let repr: Vec<_> =  a.into_iter()
            .map(|s| s.0.into_repr())
            .collect();

        Ok(Arc::new(repr))
    }   

    pub fn into_variables<P>(self, worker: &Worker, mut params: P) -> Result<()> 
    where
        P: ParameterSource<E>
    {
        let aux: Arc<Vec<_>> = Arc::new(
            self.assignment
                .aux
                .into_iter()
                .map(|s| s.into_repr())
                .collect()
        );
        let input: Arc<_> = Arc::new(
            self.assignment
                .input
                .into_iter()
                .map(|s| s.into_repr())
                .collect::<Vec<_>>()
        );
        let l = multiexp(
            &worker,
            params.get_l(aux.len())?,
            FullDensity,
            aux.clone(),
        );

        unimplemented!()
    }

    fn build_ga(vk: VerifyingKey<E>, r: E::Fr) -> () {
        let deltar: _ = vk.delta_g1.mul(r);

    }

    
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

pub struct QueryDensity {
    pub a_aux: DensityTracker,
    pub b_input: DensityTracker,
    pub b_aux: DensityTracker,
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
    pub a: Option<Vec<Scalar<E>>>,
    pub b: Option<Vec<Scalar<E>>>,
    pub c: Option<Vec<Scalar<E>>>,
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

struct ProvingField<'a, E: Engine> {
    a: EvaluationDomain<E,Scalar<E>>,
    b: EvaluationDomain<E,Scalar<E>>,
    c: EvaluationDomain<E,Scalar<E>>,
    work: &'a Worker
}

impl<'a,E> ProvingField<'a,E> 
where
    E: Engine
{
    fn new(eval: &'a mut PolynomialEvaluation<E>, work: &'a Worker) -> Result<Self> {
        let a = EvaluationDomain::from_coeffs(eval.a.take()?)?;
        let b = EvaluationDomain::from_coeffs(eval.b.take()?)?;
        let c = EvaluationDomain::from_coeffs(eval.c.take()?)?;
        Ok(ProvingField {a, b, c, work})
    }

    fn transform(mut self) -> Self {
        self.a.ifft(self.work);
        self.a.coset_fft(self.work);
        self.b.ifft(self.work);
        self.b.coset_fft(self.work);
        self.c.ifft(self.work);
        self.c.coset_fft(self.work);
        self
    }

    fn into_coefficient(mut self) -> Vec<Scalar<E>> {
        self.a.mul_assign(self.work, &self.b);
        drop(self.b);

        self.a.sub_assign(self.work, &self.c);
        drop(self.c);

        self.a.divide_by_z_on_coset(self.work);
        self.a.icoset_fft(self.work);

        self.a.into_coeffs()
    }
}