use std::sync::Arc;
use super::{
    PolynomialEvaluation, ParameterSource, Result, 
    ProvingSystem, Future, SynthesisError, 
    AssignmentField, ProvingAssignment,
    source, context, fourier, 
};

use ff::{Field, PrimeField};
use pairing::Engine;

use crate::multiexp::{multiexp, FullDensity};
use crate::multicore::Worker;
use crate::groth16::VerifyingKey;
use group::{CurveAffine, CurveProjective};

pub struct Builder<E: Engine> {
    vk: VerifyingKey<E>, 
    r: E::Fr, 
    s: E::Fr, 
    h: E::G1,
    l: E::G1,
    answer: context::Answer<E>,
    aux: context::Auxiliary<E>,
}

impl<E> Builder<E>
where
    E: Engine
{
    pub fn try_new<P>(mut prover: ProvingSystem<E>, params: &mut P, r: E::Fr, s: E::Fr) -> Result<Self> 
    where
        P: ParameterSource<E>
    {
        let worker: _ = Worker::new();

        let vk: _ = try_vk(params)?;

        let h: _ = try_h(&mut prover.eval, &worker, params)?;

        let (input, aux): _ = into_primefield(prover.assignment);

        let l = try_l(&aux, &worker, params)?;

        let mut src: _ = source::Source::try_new(prover.density, input.len(), params)?;
        let answer: _ = src.into_answer(&worker, input)?;
        let aux: _ = src.into_auxiliary(&worker, aux)?;

        let builder: _ = Self {
            vk,
            r,
            s,
            answer,
            aux,
            h: h.wait()?,
            l: l.wait()?
        };
        Ok(builder)
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
        let a_mul_s: _ = self.vk.alpha_g1.mul(self.s);
        let b_mul_r: _ = self.vk.beta_g1.mul(self.r);

        let mut gc: _ = delta_rs;
        gc.add_assign(&a_mul_s);
        gc.add_assign(&b_mul_r);

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

fn into_primefield<E>(assignment: ProvingAssignment<E>) -> (AssignmentField<E>, AssignmentField<E>) 
where
    E: Engine
{
    let input = Arc::new(
        assignment.input
            .into_iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>(),
    );

    let aux = Arc::new(
        assignment.aux
            .into_iter()
            .map(|s| s.into_repr())
            .collect::<Vec<_>>(),
    );

    (input, aux)
}

fn try_h<E,P>(eval: &mut PolynomialEvaluation<E>, worker: &Worker, params: &mut P) -> Result<impl Future<Item=E::G1, Error=SynthesisError>>
where
    E: Engine,
    P: ParameterSource<E>
{
    let linear_coeffs: _ = fourier::evaluate_coefficients(eval, worker)?;
    let multi_exponentiated_coeffs: _ = multiexp(&worker, params.get_h()?, FullDensity, linear_coeffs);
    Ok(multi_exponentiated_coeffs)
}

fn try_l<E,P>(aux: &AssignmentField<E>, worker: &Worker, params: &mut P) -> Result<impl Future<Item=E::G1, Error=SynthesisError>> 
where
    E: Engine,
    P: ParameterSource<E>
{
    let l: _ = multiexp(
        &worker,
        params.get_l()?,
        FullDensity,
        aux.clone(),
    );
    Ok(l)
}

fn try_vk<E,P>(params: &mut P) -> Result<VerifyingKey<E>> 
where
    E: Engine,
    P: ParameterSource<E>
{
    let vk = params.get_vk()?;
    if vk.delta_g1.is_zero() || vk.delta_g2.is_zero() {
        // If this element is zero, someone is trying to perform a
        // subversion-CRS attack.
        return Err(SynthesisError::UnexpectedIdentity);
    } else { Ok(vk) }
}
