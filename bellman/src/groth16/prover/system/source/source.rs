use std::sync::Arc;
use pairing::Engine;

use super::{ParameterSource, QueryDensity, AssignmentField, Result, source};
use crate::multiexp::DensityTracker;

pub struct Source<P: ParameterSource<E>, E: Engine> {
    answer: Option<source::AnswerSource<P,E>>,
    aux: Option<source::AuxiliarySource<P,E>>,
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

        let b_input_total: usize = b_input_density.get_total_density();
        
        let (a_input_src, a_aux_src): _ = params.a(a_input)?;
        let (b1_input_src, b1_aux_src): _ = params.b_g1(b_input_total)?;
        let (b2_input_src, b2_aux_src): _ = params.b_g2(b_input_total)?;

        let answer_src: _ = AnswerSource::new(
            a_input_src, 
            b1_input_src, 
            b2_input_src, 
            b_input_density
        );

        let aux_src: _ = AuxiliarySource::new(
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

    pub fn into_answer(&mut self, input: AssignmentField<E>) -> Result<super::Answer<E>> {
        super::Answer::try_new(self.answer.take()?, input)
    }

    pub fn into_auxiliary(&mut self, aux: AssignmentField<E>) -> Result<super::Auxiliary<E>> {
        super::Auxiliary::try_new(self.aux.take()?, aux)
    }
}

pub struct AnswerSource<P: ParameterSource<E>, E: Engine> {
    pub a_input_src: P::G1Builder,
    pub b1_input_src: P::G1Builder,
    pub b2_input_src: P::G2Builder,
    pub b_input_density: Arc<DensityTracker>,
    _marker: std::marker::PhantomData<E>
}

impl<P,E> AnswerSource<P,E>
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
        AnswerSource {
            a_input_src,
            b1_input_src,
            b2_input_src,
            b_input_density,
            _marker: std::marker::PhantomData
        }
    }
}

pub struct AuxiliarySource<P: ParameterSource<E>, E: Engine> {
    pub a_aux_src: P::G1Builder,
    pub b1_aux_src: P::G1Builder,
    pub b2_aux_src: P::G2Builder,
    pub a_aux_density: DensityTracker,
    pub b_aux_density: Arc<DensityTracker>,
    _marker: std::marker::PhantomData<E>
}

impl<P,E> AuxiliarySource<P,E>
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
        AuxiliarySource {
            a_aux_src,
            b1_aux_src,
            b2_aux_src,
            a_aux_density,
            b_aux_density,
            _marker: std::marker::PhantomData
        }
    }
}