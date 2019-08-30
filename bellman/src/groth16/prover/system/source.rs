use std::sync::Arc;
use pairing::Engine;

use super::{ParameterSource, QueryDensity, ArcAssignment, Result, source, context};

use crate::multiexp::DensityTracker;
use crate::multicore::Worker;

pub struct Source<P: ParameterSource<E>, E: Engine> {
    answer: Option<source::Answer<P,E>>,
    aux: Option<source::Auxiliary<P,E>>,
}

impl<P,E> Source<P,E>
where
    P: ParameterSource<E>, 
    E: Engine
{
    pub(super) fn try_new(density: QueryDensity, a_input: usize, params: &mut P) -> Result<Self> {
        let a_aux_density: _ = density.a_aux;
        let b_input_density: _ = Arc::new(density.b_input);
        let b_aux_density: _ = Arc::new(density.b_aux);

        let a_aux_total: usize = a_aux_density.get_total_density();
        let b_input_total: usize = b_input_density.get_total_density();
        let b_aux_total: usize = b_aux_density.get_total_density();       

        let (a_input_src, a_aux_src): _ = params.get_a(a_input, a_aux_total)?;
        let (b1_input_src, b1_aux_src): _ = params.get_b_g1(b_input_total, b_aux_total)?;
        let (b2_input_src, b2_aux_src): _ = params.get_b_g2(b_input_total, b_aux_total)?;

        let answer_src: _ = Answer::new(
            a_input_src, 
            b1_input_src, 
            b2_input_src, 
            b_input_density
        );

        let aux_src: _ = Auxiliary::new(
            a_aux_src,
            b1_aux_src,
            b2_aux_src,
            a_aux_density,
            b_aux_density
        );
        
        Ok(Source {
            answer: Some(answer_src),
            aux: Some(aux_src)
        })
    }        

    pub fn into_answer(&mut self, worker: &Worker, input: ArcAssignment<E>) -> Result<context::Answer<E>> {
        context::Answer::try_new(
            worker, 
            self.answer.take()?, 
            input
        )
    }

    pub fn into_auxiliary(&mut self, worker: &Worker, aux: ArcAssignment<E>) -> Result<context::Auxiliary<E>> {
        context::Auxiliary::try_new(
            worker,
            self.aux.take()?, 
            aux
        )
    }
}

pub struct Answer<P: ParameterSource<E>, E: Engine> {
    pub a_input_src: P::G1Builder,
    pub b1_input_src: P::G1Builder,
    pub b2_input_src: P::G2Builder,
    pub b_input_density: Arc<DensityTracker>,
    _marker: std::marker::PhantomData<E>
}

impl<P,E> Answer<P,E>
where
    P: ParameterSource<E>,
    E: Engine
{
    fn new(
        a_input_src: P::G1Builder,
        b1_input_src: P::G1Builder,
        b2_input_src: P::G2Builder,
        b_input_density: Arc<DensityTracker>
    ) -> Self {
        Answer {
            a_input_src,
            b1_input_src,
            b2_input_src,
            b_input_density,
            _marker: std::marker::PhantomData
        }
    }
}

pub struct Auxiliary<P: ParameterSource<E>, E: Engine> {
    pub a_aux_src: P::G1Builder,
    pub b1_aux_src: P::G1Builder,
    pub b2_aux_src: P::G2Builder,
    pub a_aux_density: DensityTracker,
    pub b_aux_density: Arc<DensityTracker>,
    _marker: std::marker::PhantomData<E>
}

impl<P,E> Auxiliary<P,E>
where
    P: ParameterSource<E>,
    E: Engine
{   
    fn new(
        a_aux_src: P::G1Builder,
        b1_aux_src: P::G1Builder,
        b2_aux_src: P::G2Builder,
        a_aux_density: DensityTracker,
        b_aux_density: Arc<DensityTracker>
    ) -> Self {
        Auxiliary {
            a_aux_src,
            b1_aux_src,
            b2_aux_src,
            a_aux_density,
            b_aux_density,
            _marker: std::marker::PhantomData
        }
    }
}