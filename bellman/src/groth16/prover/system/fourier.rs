use std::sync::Arc;

use ff::PrimeField;
use pairing::Engine;

use crate::domain::{Domain, Scalar};
use super::{PolynomialEvaluation, AssignmentField, Result};

pub fn evaluate_coefficients<E>(eval: &mut PolynomialEvaluation<E>) -> Result<AssignmentField<E>>
where
    E: Engine
{
    let fourier_eval_domain: _ = FourierEvaluationDomain::new(eval)?;
    fourier_eval_domain.coeffs_by_fft()
}

struct FourierEvaluationDomain<E: Engine> {
    a: Domain<E,Scalar<E>>,
    b: Domain<E,Scalar<E>>,
    c: Domain<E,Scalar<E>>,
}

impl<E> FourierEvaluationDomain<E> 
where
    E: Engine
{
    fn new(eval: &mut PolynomialEvaluation<E>) -> Result<Self> {
        let a = Domain::new(eval.a.take()?)?;
        let b = Domain::new(eval.b.take()?)?;
        let c = Domain::new(eval.c.take()?)?;
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
