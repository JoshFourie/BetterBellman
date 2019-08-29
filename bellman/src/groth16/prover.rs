use rand_core::RngCore;

use std::sync::Arc;

use futures::Future;

use ff::{Field, PrimeField};
use group::{CurveAffine, CurveProjective};
use pairing::Engine;

use super::{ParameterSource, Proof, Result};

use crate::{Circuit, ConstraintSystem, Index, LinearCombination, SynthesisError, Variable};
use crate::domain::{EvaluationDomain, Scalar};
use crate::multiexp::{multiexp, DensityTracker, FullDensity};
use crate::multicore::Worker;

pub struct ProvingSystem<E: Engine> {
    density: QueryDensity,
    eval: PolynomialEvaluation<E>,
    assignment: ProvingAssignment<E>
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

        create_proof::<E, C, P>(circuit, params, r, s)
    }

    fn as_primefield(&mut self, worker: &Worker) ->  Result<Arc<Vec<<<E as ff::ScalarEngine>::Fr as PrimeField>::Repr>>> {
        let mut a = EvaluationDomain::from_coeffs(self.eval.a.take()?)?;
        let mut b = EvaluationDomain::from_coeffs(self.eval.b.take()?)?;
        let mut c = EvaluationDomain::from_coeffs(self.eval.c.take()?)?;
        
        a.ifft(&worker);
        a.coset_fft(&worker);
        b.ifft(&worker);
        b.coset_fft(&worker);
        c.ifft(&worker);
        c.coset_fft(&worker);

        a.mul_assign(&worker, &b);
        drop(b);
        a.sub_assign(&worker, &c);
        drop(c);
        a.divide_by_z_on_coset(&worker);
        a.icoset_fft(&worker);

        let mut a = a.into_coeffs();
        let new_len = a.len() - 1;
        a.truncate(new_len);

        let repr: Vec<_> =  a.into_iter()
            .map(|s| s.0.into_repr())
            .collect();

        Ok(Arc::new(repr))
    }   
}

pub fn create_proof<E, C, P: ParameterSource<E>>(
    circuit: C,
    mut params: P,
    r: E::Fr,
    s: E::Fr,
) -> Result<Proof<E>>
where
    E: Engine,
    C: Circuit<E>,
{
    let mut prover: _ = ProvingSystem::default();

    prover.alloc_input(
        || "", 
        || Ok(E::Fr::one())
    )?;

    circuit.synthesize(&mut prover)?;

    for i in 0..prover.assignment.input.len() {
        prover.enforce(
            || "", 
            |lc| lc + Variable(Index::Input(i)), 
            |lc| lc, |lc| lc
        );
    }

    let worker = Worker::new();

    let vk = params.get_vk(prover.assignment.input.len())?;

    let h = {
        let a: _ = prover.as_primefield(&worker)?;
        multiexp(&worker, params.get_h(a.len())?, FullDensity, a)
    };

    // TODO: parallelize if it's even helpful
    let input_assignment = Arc::new(
        prover.assignment
            .input
            .into_iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>(),
    );
    let aux_assignment = Arc::new(
        prover.assignment
            .aux
            .into_iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>(),
    );

    let l = multiexp(
        &worker,
        params.get_l(aux_assignment.len())?,
        FullDensity,
        aux_assignment.clone(),
    );

    let a_aux_density_total = prover.density
        .a_aux
        .get_total_density();

    let (a_inputs_source, a_aux_source) = params.get_a(
        input_assignment.len(), 
        a_aux_density_total
    )?;

    let a_inputs = multiexp(
        &worker,
        a_inputs_source,
        FullDensity,
        input_assignment.clone(),
    );
    let a_aux = multiexp(
        &worker,
        a_aux_source,
        Arc::new(prover.density.a_aux),
        aux_assignment.clone(),
    );

    let b_input = Arc::new(prover.density.b_input);
    let b_input_total = b_input.get_total_density();
    let b_aux_density = Arc::new(prover.density.b_aux);
    let b_aux_density_total = b_aux_density.get_total_density();

    let (b_g1_inputs_source, b_g1_aux_source) =
        params.get_b_g1(b_input_total, b_aux_density_total)?;

    let b_g1_inputs = multiexp(
        &worker,
        b_g1_inputs_source,
        b_input.clone(),
        input_assignment.clone(),
    );
    let b_g1_aux = multiexp(
        &worker,
        b_g1_aux_source,
        b_aux_density.clone(),
        aux_assignment.clone(),
    );

    let (b_g2_inputs_source, b_g2_aux_source) =
        params.get_b_g2(b_input_total, b_aux_density_total)?;

    let b_g2_inputs = multiexp(
        &worker,
        b_g2_inputs_source,
        b_input,
        input_assignment,
    );
    let b_g2_aux = multiexp(&worker, b_g2_aux_source, b_aux_density, aux_assignment);

    if vk.delta_g1.is_zero() || vk.delta_g2.is_zero() {
        // If this element is zero, someone is trying to perform a
        // subversion-CRS attack.
        return Err(SynthesisError::UnexpectedIdentity);
    }

    let mut g_a = vk.delta_g1.mul(r);
    g_a.add_assign_mixed(&vk.alpha_g1);
    let mut g_b = vk.delta_g2.mul(s);
    g_b.add_assign_mixed(&vk.beta_g2);
    let mut g_c;
    {
        let mut rs = r;
        rs.mul_assign(&s);

        g_c = vk.delta_g1.mul(rs);
        g_c.add_assign(&vk.alpha_g1.mul(s));
        g_c.add_assign(&vk.beta_g1.mul(r));
    }
    let mut a_answer = a_inputs.wait()?;
    a_answer.add_assign(&a_aux.wait()?);
    g_a.add_assign(&a_answer);
    a_answer.mul_assign(s);
    g_c.add_assign(&a_answer);

    let mut b1_answer = b_g1_inputs.wait()?;
    b1_answer.add_assign(&b_g1_aux.wait()?);
    let mut b2_answer = b_g2_inputs.wait()?;
    b2_answer.add_assign(&b_g2_aux.wait()?);

    g_b.add_assign(&b2_answer);
    b1_answer.mul_assign(r);
    g_c.add_assign(&b1_answer);
    g_c.add_assign(&h.wait()?);
    g_c.add_assign(&l.wait()?);

    Ok(Proof {
        a: g_a.into_affine(),
        b: g_b.into_affine(),
        c: g_c.into_affine(),
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

struct PolynomialEvaluation<E: Engine> {
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

struct ProvingAssignment<E: Engine> {
    input: Vec<E::Fr>,
    aux: Vec<E::Fr>,
}

impl<E: Engine> Default for ProvingAssignment<E> {
    fn default() -> Self {
        ProvingAssignment {
            input: Vec::new(),
            aux: Vec::new()
        }
    }
} 
