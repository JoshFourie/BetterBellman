use std::ops::{Add, Sub};
use ff::{ScalarEngine, Field};

/// Represents a variable in our constraint system.
#[derive(Copy, Clone, Debug)]
pub struct Variable(Index);

impl Variable {
    /// This constructs a variable with an arbitrary index.
    /// Circuit implementations are not recommended to use this.
    pub fn new_unchecked(idx: Index) -> Variable {
        Variable(idx)
    }

    /// This returns the index underlying the variable.
    /// Circuit implementations are not recommended to use this.
    pub fn get_unchecked(&self) -> Index {
        self.0
    }
}

/// Represents the index of either an input variable or
/// auxiliary variable.
#[derive(Copy, Clone, PartialEq, Debug)]
pub enum Index {
    Input(usize),
    Aux(usize),
}

/// This represents a linear combination of some variables, with coefficients
/// in the scalar field of a pairing-friendly elliptic curve group.
#[derive(Clone)]
pub struct LinearCombination<E: ScalarEngine>(pub Vec<(Variable, E::Fr)>);

impl<E> AsRef<[(Variable, E::Fr)]> for LinearCombination<E> 
where
    E: ScalarEngine
{
    fn as_ref(&self) -> &[(Variable, E::Fr)] {
        &self.0
    }
}

impl<E> LinearCombination<E> 
where
    E: ScalarEngine
{
    pub fn zero() -> Self {
        LinearCombination(vec![])
    }
}

impl<E> Add<(E::Fr, Variable)> for LinearCombination<E> 
where
    E: ScalarEngine
{
    type Output = Self;

    fn add(mut self, (coeff, var): (E::Fr, Variable)) -> LinearCombination<E> {
        self.0.push((var, coeff));

        self
    }
}

impl<E> Sub<(E::Fr, Variable)> for LinearCombination<E> 
where
    E: ScalarEngine
{
    type Output = Self;

    fn sub(self, (mut coeff, var): (E::Fr, Variable)) -> LinearCombination<E> {
        coeff.negate();

        self + (coeff, var)
    }
}

impl<E> Add<Variable> for LinearCombination<E> 
where
    E: ScalarEngine
{
    type Output = Self;

    fn add(self, other: Variable) -> LinearCombination<E> {
        self + (E::Fr::one(), other)
    }
}

impl<E> Sub<Variable> for LinearCombination<E> 
where
    E: ScalarEngine
{
    type Output = Self;

    fn sub(self, other: Variable) -> LinearCombination<E> {
        self - (E::Fr::one(), other)
    }
}

impl<'a, E> Add<&'a LinearCombination<E>> for LinearCombination<E> 
where
    E: ScalarEngine
{
    type Output = Self;

    fn add(mut self, other: &'a LinearCombination<E>) -> LinearCombination<E> {
        for s in &other.0 {
            self = self + (s.1, s.0);
        }

        self
    }
}

impl<'a, E> Sub<&'a LinearCombination<E>> for LinearCombination<E> 
where
    E: ScalarEngine
{
    type Output = Self;

    fn sub(mut self, other: &'a LinearCombination<E>) -> LinearCombination<E> {
        for s in &other.0 {
            self = self - (s.1, s.0);
        }

        self
    }
}

impl<'a, E> Add<(E::Fr, &'a LinearCombination<E>)> for LinearCombination<E> 
where
    E: ScalarEngine
{
    type Output = Self;

    fn add(mut self, (coeff, other): (E::Fr, &'a LinearCombination<E>)) -> LinearCombination<E> {
        for s in &other.0 {
            let mut tmp = s.1;
            tmp.mul_assign(&coeff);
            self = self + (tmp, s.0);
        }

        self
    }
}

impl<'a, E> Sub<(E::Fr, &'a LinearCombination<E>)> for LinearCombination<E> 
where
    E: ScalarEngine
{
    type Output = Self;

    fn sub(mut self, (coeff, other): (E::Fr, &'a LinearCombination<E>)) -> LinearCombination<E> {
        for s in &other.0 {
            let mut tmp = s.1;
            tmp.mul_assign(&coeff);
            self = self - (tmp, s.0);
        }

        self
    }
}
