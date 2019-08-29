use rand_core::RngCore;

use std::sync::Arc;
use super::{Future, SynthesisError};

use ff::{Field, PrimeField, ScalarEngine};
use pairing::Engine;

use super::{ParameterSource, Proof, Result};

use crate::{Circuit, ConstraintSystem, Index, LinearCombination, Variable};
use crate::domain::{EvaluationDomain, Scalar};
use crate::multiexp::{multiexp, DensityTracker, FullDensity};
use crate::multicore::Worker;
use crate::groth16::VerifyingKey;
use group::{CurveAffine, CurveProjective};

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

    pub fn fft_shortcut(&mut self, worker: &Worker) ->  Result<Arc<Vec<<<E as ff::ScalarEngine>::Fr as PrimeField>::Repr>>> {
        let field: _ = FourierField::new(&mut self.eval, worker)?;
        let mut a: _ = field.transform().into_coefficient();

        let new_len = a.len() - 1;
        a.truncate(new_len);

        let repr: Vec<_> =  a.into_iter()
            .map(|s| s.0.into_repr())
            .collect();


        Ok(Arc::new(repr))
    }   

    pub fn src_gb1<T>(
        input_density: Arc<DensityTracker>, 
        aux_density: Arc<DensityTracker>, 
        params: &mut T
    ) -> Result<(T::G1Builder, T::G1Builder)> 
    where
        T: ParameterSource<E>
    {
        let b_input_total = input_density.get_total_density();
        let b_aux_density_total = aux_density.get_total_density();

        Ok(params.get_b_g1(b_input_total, b_aux_density_total)?)
    }

    pub fn src_gb2<U>(
        input_density: usize, 
        aux_density: usize, 
        params: &mut U
    ) -> Result<(U::G2Builder, U::G2Builder)> 
    where
        U: ParameterSource<E>
    {
        Ok(params.get_b_g2(input_density, aux_density)?)
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

pub struct GroupBuilder<'a, E: Engine> {
    pub vk: &'a VerifyingKey<E>, 
    pub r: E::Fr, 
    pub s: E::Fr, 
    pub answer: AnswerContext<E>,
    pub aux: AuxiliaryContext<E>,
    pub h: E::G1,
    pub l: E::G1
}

impl<'a,E> GroupBuilder<'a,E>
where
    E: Engine
{
    pub fn try_new<P>(prover: ProvingSystem<E>, params: P) -> Result<Self> 
    where
        P: ParameterSource<E>
    {
        unimplemented!()
    }

    pub fn try_build(mut self) -> Result<(E::G1, E::G2, E::G1)> {
        let ga: _ = self.try_ga()?;
        let gb: _ = self.try_gb()?;
        let gc: _ = self.try_gc()?;
        Ok((ga, gb, gc))
    }

    fn try_ga(&mut self) -> Result<E::G1> {
        let mut ga: _ = self.vk.delta_g1.mul(self.r);
        ga.add_assign_mixed(&self.vk.alpha_g1);

        self.answer.a.add_assign(&self.aux.a);
        ga.add_assign(&self.answer.a);
        
        Ok(ga)
    }

    fn try_gb(&mut self) -> Result<E::G2> {
        let mut gb: _ = self.vk.delta_g2.mul(self.s);
        gb.add_assign_mixed(&self.vk.beta_g2);

        self.answer.b2.add_assign(&self.aux.b2);
        gb.add_assign(&self.answer.b2);

        Ok(gb)
    }   

    fn try_gc(mut self) -> Result<E::G1> {
        let delta_rs: E::G1 = {
            let mut rs: _ = self.r; 
            rs.mul_assign(&self.s);
            self.vk.delta_g1.mul(rs)
        };
        let As: _ = self.vk.alpha_g1.mul(self.s);
        let Br: _ = self.vk.beta_g1.mul(self.r);

        let mut gc: _ = delta_rs;
        gc.add_assign(&As);
        gc.add_assign(&Br);

        self.answer.a.mul_assign(self.s);
        gc.add_assign(&self.answer.a);

        self.answer.b1.add_assign(&self.aux.b1);
        self.answer.b1.mul_assign(self.r);
        gc.add_assign(&self.answer.b1);

        gc.add_assign(&self.h);
        gc.add_assign(&self.l); 

        Ok(gc)
    } 
}

pub struct AnswerContext<E: Engine> {
    a: E::G1,
    b1: E::G1,
    b2: E::G2
}

impl<E> AnswerContext<E>
where
    E: Engine
{
    pub fn try_new<P: ParameterSource<E>>(
        worker: &Worker, 
        a_src: P::G1Builder, 
        b1_src: P::G1Builder, 
        b2_src: P::G2Builder,
        b_in: Arc<DensityTracker>,
        input: Arc<Vec<<<E as ff::ScalarEngine>::Fr as PrimeField>::Repr>>
    ) -> Result<Self> {
        let a: E::G1 = multiexp(
            &worker,
            a_src,
            FullDensity,
            input.clone(),
        ).wait()?;

        let b1: E::G1 = multiexp(
            &worker,
            b1_src,
            b_in.clone(),
            input.clone(),
        ).wait()?;

        let b2: E::G2 = multiexp(
            &worker,
            b2_src,
            b_in,
            input
        ).wait()?;

        Ok(AnswerContext { a, b1, b2 })
    }
}

pub struct AuxiliaryContext<E: Engine> {
    a: E::G1,
    b1: E::G1,
    b2: E::G2,  
}

impl<E> AuxiliaryContext<E> 
where
    E: Engine
{
    pub fn try_new<P>(
        worker: &Worker,
        a_src: P::G1Builder,
        b1_src: P::G1Builder,
        b2_src: P::G2Builder,
        a_density: DensityTracker,
        b_density: Arc<DensityTracker>,
        assignment: Arc<Vec<<<E as ScalarEngine>::Fr as PrimeField>::Repr>>
    ) -> Result<Self> 
    where
        P: ParameterSource<E>
    {
        let a: _ = multiexp(
            &worker,
            a_src,
            Arc::new(a_density),
            assignment.clone(),
        ).wait()?;

        let b1: _ = multiexp(
            &worker,
            b1_src,
            b_density.clone(),
            assignment.clone(),
        ).wait()?;

        let b2 = multiexp(
            &worker, 
            b2_src, 
            b_density, 
            assignment
        ).wait()?;

        Ok(AuxiliaryContext{ a, b1, b2 })
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

struct FourierField<'a, E: Engine> {
    a: EvaluationDomain<E,Scalar<E>>,
    b: EvaluationDomain<E,Scalar<E>>,
    c: EvaluationDomain<E,Scalar<E>>,
    work: &'a Worker
}

impl<'a,E> FourierField<'a,E> 
where
    E: Engine
{
    fn new(eval: &'a mut PolynomialEvaluation<E>, work: &'a Worker) -> Result<Self> {
        let a = EvaluationDomain::from_coeffs(eval.a.take()?)?;
        let b = EvaluationDomain::from_coeffs(eval.b.take()?)?;
        let c = EvaluationDomain::from_coeffs(eval.c.take()?)?;
        Ok(FourierField {a, b, c, work})
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