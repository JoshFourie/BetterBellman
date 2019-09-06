use std::ops;

use ff::{ScalarEngine, Field, PrimeField};
use group::CurveProjective;

pub trait Group<'a,E>: Sized 
    + Copy 
    + Clone 
    + Send 
    + Sync 
    + ops::MulAssign<&'a E::Fr>
    + ops::SubAssign<&'a Self>
    + ops::AddAssign<&'a Self>
where
    Self: 'a,
    E: ScalarEngine    
{
    fn zero() -> Self;
}

pub struct Point<G>(pub G);

impl<'a,G> ops::MulAssign<&'a G::Scalar> for Point<G> 
where
    G: CurveProjective
{
    fn mul_assign(&mut self, rhs: &'a G::Scalar) {
        self.0.mul_assign(rhs.into_repr());
    }
} 

impl<'a,G> ops::AddAssign<&'a Self> for Point<G> 
where
    G: CurveProjective
{
    fn add_assign(&mut self, rhs: &'a Self) {
        self.0.add_assign(&rhs.0);   
    }
}

impl<'a,G> ops::SubAssign<&'a Self> for Point<G> 
where
    G: CurveProjective
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        self.0.sub_assign(&rhs.0);
    }
}

impl<G: CurveProjective> PartialEq for Point<G> {
    fn eq(&self, other: &Point<G>) -> bool {
        self.0 == other.0
    }
}

impl<G: CurveProjective> Copy for Point<G> {}

impl<G: CurveProjective> Clone for Point<G> {
    fn clone(&self) -> Point<G> {
        *self
    }
}

impl<'a,G> Group<'a,G::Engine> for Point<G> 
where
    G: CurveProjective
{
    fn zero() -> Self {
        Point(G::zero())
    }
}

pub struct Scalar<E: ScalarEngine>(pub E::Fr);

impl<E: ScalarEngine> PartialEq for Scalar<E> {
    fn eq(&self, other: &Scalar<E>) -> bool {
        self.0 == other.0
    }
}

impl<E: ScalarEngine> Copy for Scalar<E> {}

impl<E: ScalarEngine> Clone for Scalar<E> {
    fn clone(&self) -> Scalar<E> {
        *self
    }
}

impl<'a,E> Group<'a,E> for Scalar<E>
where
    E: ScalarEngine
{
    fn zero() -> Self {
        Scalar(E::Fr::zero())
    }
}

impl<'a,E> ops::MulAssign<&'a E::Fr> for Scalar<E> 
where
    E: ScalarEngine
{
    fn mul_assign(&mut self, rhs: &'a E::Fr) {
        self.0.mul_assign(&rhs);
    }
} 

impl<'a,E> ops::AddAssign<&'a Self> for Scalar<E> 
where
    E: ScalarEngine
{
    fn add_assign(&mut self, rhs: &'a Self) {
        self.0.add_assign(&rhs.0);   
    }
}

impl<'a,E> ops::SubAssign<&'a Self> for Scalar<E> 
where
    E: ScalarEngine
{
    fn sub_assign(&mut self, rhs: &'a Self) {
        self.0.sub_assign(&rhs.0);
    }
}