use std::sync::Arc;
use super::{ParameterSource, Result, ProvingSystem, Future, source, context, fourier};

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
    answer: context::Answer<E>,
    aux: context::Auxiliary<E>,
    h: E::G1,
    l: E::G1
}

impl<E> Builder<E>
where
    E: Engine
{
    pub fn try_new<P>(
        mut prover: ProvingSystem<E>, 
        worker: Worker, 
        params: &mut P,
        vk: VerifyingKey<E>,
        r: E::Fr, 
        s: E::Fr
    ) -> Result<Self> 
    where
        P: ParameterSource<E>
    {
        let h = {
            let field: _ = fourier::FourierField::new(&mut prover.eval, &worker)?;
            let linear_coeff: _ = field.fft_shortcut()?;
            multiexp(&worker, params.get_h(linear_coeff.len())?, FullDensity, linear_coeff)
        };

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

        let mut src: _ = source::Source::try_new(
            prover.density, 
            input_assignment.len(), 
            params
        )?;

        let answer: _ = src.into_answer(&worker, input_assignment)?;
        let aux: _ = src.into_auxiliary(&worker, aux_assignment)?;

        Ok(Builder {
            vk,
            r,
            s,
            answer,
            aux,
            h: h.wait()?,
            l: l.wait()?
        })
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
