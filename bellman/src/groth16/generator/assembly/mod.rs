mod parameters;
mod key_pair;
mod eval;
mod windows;

use pairing::Engine;
use group::Wnaf;

use crate::{domain, Circuit};
use domain::Scalar;

use parameters::ParameterAssembly;
use key_pair::KeyPairAssembly;
use eval::Evaluation;
use windows::BasedWindows;

use super::*;

pub use windows::Windows;
pub use parameters::into_lagrange_coefficients;

pub struct Assembly<E,C> 
where
    E: Engine,
    C: Circuit<E>
{
    param: Option<ParameterAssembly<E,C>>,
    key_pair: Option<KeyPairAssembly<E>>,
    result: Evaluation<E>
}

impl<E,C> Assembly<E,C>
where
    E: Engine,
    C: Circuit<E>
{
    pub fn new(circuit: C, g1: E::G1, g2: E::G2, alpha: E::Fr, beta: E::Fr, gamma: E::Fr, delta: E::Fr, tau: E::Fr) -> Result<Self> {
        let mut param: _ = ParameterAssembly::new(circuit, g1, g2, alpha, beta, gamma, delta, tau)?;
        let key_pair: KeyPairAssembly<E> = param.key_assembly()?;
        let result: _ = Evaluation::new(&key_pair);

        Ok(Self { 
            param: Some(param), 
            key_pair: Some(key_pair), 
            result: result 
        })
    }

    pub fn evaluation_domain(&self) -> Result<Domain<E, Scalar<E>>> {
        self.key_pair
            .as_ref()?
            .blind_evaluation_base()
    }

    pub fn h(&mut self, domain: &mut Domain<E, Scalar<E>>, based_g1: &Wnaf<usize, &[E::G1], &mut Vec<i64>>) -> Result<Vec<E::G1Affine>> {
        self.param
            .as_mut()?
            .h(domain, based_g1)
    }

    pub fn evaluate(&mut self, win: &BasedWindows<'_,E>, coeffs: &[Scalar<E>]) -> Result<()> {
        self.param
            .as_mut()?
            .evaluate(
                &mut self.result, 
                self.key_pair.take()?, 
                win, 
                coeffs
            )
    }

    pub fn verifying_key(&mut self) -> Result<VerifyingKey<E>> {
        let vk: _ = self.param
            .take()?
            .into_verifying_key(&self.result);
        Ok(vk)
    }

    pub fn result_is_unconstrained(&self) -> Result<bool> {
        Ok(self.result.is_unconstrained())
    }

    pub fn results(self) -> Evaluation<E> {
        self.result
    }
}
