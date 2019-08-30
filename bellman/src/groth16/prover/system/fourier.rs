use std::sync::Arc;

use ff::PrimeField;
use pairing::Engine;

use crate::domain::{EvaluationDomain, Scalar};
use super::{PolynomialEvaluation, Worker, AssignmentField, Result};

pub struct FourierField<'a, E: Engine> {
    a: EvaluationDomain<E,Scalar<E>>,
    b: EvaluationDomain<E,Scalar<E>>,
    c: EvaluationDomain<E,Scalar<E>>,
    work: &'a Worker
}

impl<'a,E> FourierField<'a,E> 
where
    E: Engine
{
    pub fn new(eval: &'a mut PolynomialEvaluation<E>, work: &'a Worker) -> Result<Self> {
        let a = EvaluationDomain::from_coeffs(eval.a.take()?)?;
        let b = EvaluationDomain::from_coeffs(eval.b.take()?)?;
        let c = EvaluationDomain::from_coeffs(eval.c.take()?)?;
        Ok(FourierField {a, b, c, work})
    }

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
