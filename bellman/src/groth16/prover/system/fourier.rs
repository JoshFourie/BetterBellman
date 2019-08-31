use std::sync::Arc;

use ff::PrimeField;
use pairing::Engine;

use crate::domain::{EvaluationDomain, Scalar};
use super::{PolynomialEvaluation, Worker, AssignmentField, Result};

pub struct FourierField<'a, E: Engine> {
    a: EvaluationDomain<E,Scalar<E>>,
    b: EvaluationDomain<E,Scalar<E>>,
    c: EvaluationDomain<E,Scalar<E>>,
    worker: &'a Worker
}

impl<'a,E> FourierField<'a,E> 
where
    E: Engine
{
    pub fn new(eval: &'a mut PolynomialEvaluation<E>, worker: &'a Worker) -> Result<Self> {
        let a = EvaluationDomain::from_coeffs(eval.a.take()?)?;
        let b = EvaluationDomain::from_coeffs(eval.b.take()?)?;
        let c = EvaluationDomain::from_coeffs(eval.c.take()?)?;
        Ok(FourierField {a, b, c, worker})
    }

    // The efficiency shortcut for building coefficients from the groth16 paper.
    pub fn fft_shortcut(self) ->  Result<AssignmentField<E>> {
        let mut a: _ = self.transform().into_coefficient();
        let new_len = a.len() - 1;
        a.truncate(new_len);

        let repr: Vec<_> =  a.into_iter()
            .map(|s| s.0.into_repr())
            .collect();
            
        Ok(Arc::new(repr))
    }   

    fn transform(mut self) -> Self {
        let worker: _ = self.worker;

        self.a.ifft(worker);
        self.a.coset_fft(worker);
        self.b.ifft(worker);
        self.b.coset_fft(worker);
        self.c.ifft(worker);
        self.c.coset_fft(worker);
        self
    }

    fn into_coefficient(mut self) -> Vec<Scalar<E>> {
        let worker: _ = self.worker;

        self.a.mul_assign(worker, &self.b);
        drop(self.b);

        self.a.sub_assign(worker, &self.c);
        drop(self.c);

        self.a.divide_by_z_on_coset(worker);
        self.a.icoset_fft(worker);
        self.a.into_coeffs()
    }
}
