use group::CurveProjective;
use pairing::Engine;

use super::key_pair::KeyPairAssembly;

mod writer;
pub use writer::*;

pub struct Evaluation<E>
where
    E: Engine
{
    a: Vec<E::G1>,
    b_g1: Vec<E::G1>,
    b_g2: Vec<E::G2>,
    l: Vec<E::G1>,
    pub ic: Vec<E::G1>
}

impl<E> Evaluation<E> 
where
    E: Engine
{
    pub fn new(key_pair: &KeyPairAssembly<E>) -> Self {
        let size_of_key_pair: usize = key_pair.num.inputs + key_pair.num.aux; 

        let a = vec![E::G1::zero(); size_of_key_pair];
        let b_g1 = vec![E::G1::zero(); size_of_key_pair];
        let b_g2 = vec![E::G2::zero(); size_of_key_pair];
        let ic = vec![E::G1::zero(); key_pair.num.inputs];
        let l = vec![E::G1::zero(); key_pair.num.aux];
        
        Evaluation { a, b_g1, b_g2, ic, l }
    }

    pub fn as_aux(&mut self, aux_bound: usize) -> Writer<E> {
        Writer::new(
            &mut self.a[aux_bound..],
            &mut self.b_g1[aux_bound..],
            &mut self.b_g2[aux_bound..],
            &mut self.l
        )
    }

    pub fn as_inputs<'a>(&'a mut self, input_bound: usize) -> Writer<'a,E> {
        Writer::new(
            &mut self.a[0..input_bound],
            &mut self.b_g1[0..input_bound],
            &mut self.b_g2[0..input_bound],
            &mut self.ic
        )
    }

    pub fn is_unconstrained(&self) -> bool {
        for e in self.l.iter() {
            if e.is_zero() {
                return true
            }
        }
        false
    }
    
    pub fn filter_into_affine(self) -> (Vec<E::G1Affine>, Vec<E::G1Affine>, Vec<E::G1Affine>, Vec<E::G2Affine>) {
        macro_rules! map_to_affine {
            ( $($field:expr),+ ) => {
                (
                    $(
                        $field.into_iter()
                            .filter(|e| !e.is_zero())
                            .map(|e| e.into_affine())
                            .collect() 
                    ),+
                )
            }
        }
        map_to_affine!(self.l, self.a, self.b_g1, self.b_g2)
    }
}
