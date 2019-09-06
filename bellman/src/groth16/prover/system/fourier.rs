use std::sync::Arc;

use ff::PrimeField;
use pairing::Engine;

use crate::domain::EvaluationDomain;
use crate::arith::Scalar;
use super::{PolynomialEvaluation, Worker, AssignmentField, Result};

pub fn evaluate_coefficients<E>(eval: &mut PolynomialEvaluation<E>, worker: &Worker) -> Result<AssignmentField<E>>
where
    E: Engine
{
    let fourier_eval_domain: _ = FourierEvaluationDomain::new(eval, worker)?;
    fourier_eval_domain.coeffs_by_fft()
}

// fn into_domains(eval: &'a mut PolynomialEvaluation<E>, worker: &'a Worker) 

struct FourierEvaluationDomain<'a, E: Engine> {
    a: EvaluationDomain<'a,E,Scalar<E>>,
    b: EvaluationDomain<'a,E,Scalar<E>>,
    c: EvaluationDomain<'a,E,Scalar<E>>,
}

impl<'a,E> FourierEvaluationDomain<'a,E> 
where
    E: Engine
{
    fn new(eval: &'a mut PolynomialEvaluation<E>, worker: &'a Worker) -> Result<Self> {
        let a = EvaluationDomain::new(eval.a.take()?, worker)?;
        let b = EvaluationDomain::new(eval.b.take()?, worker)?;
        let c = EvaluationDomain::new(eval.c.take()?, worker)?;
        Ok(FourierEvaluationDomain {a, b, c})
    }

    // The efficiency shortcut for building coefficients from the groth16 paper.
    fn coeffs_by_fft(self) ->  Result<AssignmentField<E>> {
        let mut a: _ = self.transform_field().into_coefficients()?;
        let new_len = a.len() - 1;
        a.truncate(new_len);

        let repr: Vec<_> =  a.into_iter()
            .map(|s| s.0.into_repr())
            .collect();
            
        Ok(Arc::new(repr))
    }   

    fn transform_field(mut self) -> Self {
        self.a.ifft();
        self.a.coset_fft();

        self.b.ifft();
        self.b.coset_fft();

        self.c.ifft();
        self.c.coset_fft();

        self
    }

    fn into_coefficients(mut self) -> Result<Vec<Scalar<E>>> {
        self.a *= &self.b;
        drop(self.b);

        self.a -= &self.c;
        drop(self.c);

        self.a.divide_over_coset()?;
        self.a.icoset_fft();

        let coeffs: _ = self.a.into_coeffs();
        Ok(coeffs)
    }
}
