use std::sync::Arc;
use pairing::Engine;

use super::{Future, QueryDensity, ParameterSource, AssignmentField, Result};
use crate::multiexp::{multiexp, FullDensity};

mod source;

pub struct SourceFactory<E>
where
    E: Engine
{
    answer: Answer<E>,
    auxiliary: Auxiliary<E>
}

impl<E> SourceFactory<E>
where
    E: Engine
{
    pub(super) fn try_new<P>(density: QueryDensity, input: AssignmentField<E>, aux: AssignmentField<E>, params: &mut P) -> Result<Self> 
    where
        P: ParameterSource<E>
    {
        let mut src: _ = source::Source::try_new(density, input.len(), params)?;
        Ok(SourceFactory {
            answer: src.into_answer(input)?,
            auxiliary: src.into_auxiliary(aux)?
        })
    }

    pub fn unpack(self) -> (Answer<E>, Auxiliary<E>) {
        (self.answer, self.auxiliary)
    }
}

pub struct Answer<E: Engine> {
    pub a: E::G1,
    pub b1: E::G1,
    pub b2: E::G2
}

impl<E> Answer<E>
where
    E: Engine
{
    pub fn try_new<P>(src: source::AnswerSource<P,E>, input: AssignmentField<E>) -> Result<Self> 
    where
        P: ParameterSource<E>
    {
        let a: E::G1 = multiexp(src.a_input_src, FullDensity, input.clone()).wait()?;

        let b1: E::G1 = multiexp(
            src.b1_input_src,
            src.b_input_density.clone(),
            input.clone(),
        ).wait()?;

        let b2: E::G2 = multiexp(
            src.b2_input_src,
            src.b_input_density,
            input
        ).wait()?;

        Ok(Answer { a, b1, b2 })
    }
}

pub struct Auxiliary<E: Engine> {
    pub a: E::G1,
    pub b1: E::G1,
    pub b2: E::G2,  
}

impl<E> Auxiliary<E> 
where
    E: Engine
{
    pub fn try_new<P>(src: source::AuxiliarySource<P,E>, assignment: AssignmentField<E>) -> Result<Self> 
    where
        P: ParameterSource<E>
    {
        let a: _ = multiexp(
            src.a_aux_src,
            Arc::new(src.a_aux_density),
            assignment.clone(),
        ).wait()?;

        let b1: _ = multiexp(
            src.b1_aux_src,
            src.b_aux_density.clone(),
            assignment.clone(),
        ).wait()?;

        let b2 = multiexp(
            src.b2_aux_src, 
            src.b_aux_density, 
            assignment
        ).wait()?;

        Ok(Auxiliary{ a, b1, b2 })
    }
}
