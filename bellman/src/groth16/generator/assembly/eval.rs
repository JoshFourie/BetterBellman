use ff::{Field, PrimeField};
use group::{CurveProjective, Wnaf};
use pairing::Engine;

use crate::{multicore, arith};
use multicore::Worker;
use arith::Scalar;

use super::key_pair::KeyPairAssembly;

pub fn eval<E: Engine>(
    // wNAF window tables
    g1_wnaf: &Wnaf<usize, &[E::G1], &mut Vec<i64>>,
    g2_wnaf: &Wnaf<usize, &[E::G2], &mut Vec<i64>>,

    // Lagrange coefficients for tau
    powers_of_tau: &[Scalar<E>],

    // QAP polynomials
    at: &[Vec<(E::Fr, usize)>],
    bt: &[Vec<(E::Fr, usize)>],
    ct: &[Vec<(E::Fr, usize)>],

    // Resulting evaluated QAP polynomials
    a: &mut [E::G1],
    b_g1: &mut [E::G1],
    b_g2: &mut [E::G2],
    ext: &mut [E::G1],

    // Inverse coefficient for ext elements
    inv: &E::Fr,

    // Trapdoors
    alpha: &E::Fr,
    beta: &E::Fr,

    // Worker
    worker: &Worker,
) {
    // Sanity check
    assert_eq!(a.len(), at.len());
    assert_eq!(a.len(), bt.len());
    assert_eq!(a.len(), ct.len());
    assert_eq!(a.len(), b_g1.len());
    assert_eq!(a.len(), b_g2.len());
    assert_eq!(a.len(), ext.len());

    // Evaluate polynomials in multiple threads
    worker.scope(a.len(), |scope, chunk| {
        for ((((((a, b_g1), b_g2), ext), at), bt), ct) in a
            .chunks_mut(chunk)
            .zip(b_g1.chunks_mut(chunk))
            .zip(b_g2.chunks_mut(chunk))
            .zip(ext.chunks_mut(chunk))
            .zip(at.chunks(chunk))
            .zip(bt.chunks(chunk))
            .zip(ct.chunks(chunk))
        {
            let mut g1_wnaf = g1_wnaf.shared();
            let mut g2_wnaf = g2_wnaf.shared();

            scope.spawn(move || {
                for ((((((a, b_g1), b_g2), ext), at), bt), ct) in a
                    .iter_mut()
                    .zip(b_g1.iter_mut())
                    .zip(b_g2.iter_mut())
                    .zip(ext.iter_mut())
                    .zip(at.iter())
                    .zip(bt.iter())
                    .zip(ct.iter())
                {
                    fn eval_at_tau<E: Engine>(
                        powers_of_tau: &[Scalar<E>],
                        p: &[(E::Fr, usize)],
                    ) -> E::Fr {
                        let mut acc = E::Fr::zero();

                        for &(ref coeff, index) in p {
                            let mut n = powers_of_tau[index].0;
                            n.mul_assign(coeff);
                            acc.add_assign(&n);
                        }

                        acc
                    }

                    // Evaluate QAP polynomials at tau
                    let mut at = eval_at_tau(powers_of_tau, at);
                    let mut bt = eval_at_tau(powers_of_tau, bt);
                    let ct = eval_at_tau(powers_of_tau, ct);

                    // Compute A query (in G1)
                    if !at.is_zero() {
                        *a = g1_wnaf.scalar(at.into_repr());
                    }

                    // Compute B query (in G1/G2)
                    if !bt.is_zero() {
                        let bt_repr = bt.into_repr();
                        *b_g1 = g1_wnaf.scalar(bt_repr);
                        *b_g2 = g2_wnaf.scalar(bt_repr);
                    }

                    at.mul_assign(&beta);
                    bt.mul_assign(&alpha);

                    let mut e = at;
                    e.add_assign(&bt);
                    e.add_assign(&ct);
                    e.mul_assign(inv);

                    *ext = g1_wnaf.scalar(e.into_repr());
                }

                // Batch normalize
                E::G1::batch_normalization(a);
                E::G1::batch_normalization(b_g1);
                E::G2::batch_normalization(b_g2);
                E::G1::batch_normalization(ext);
            });
        }
    });
}

pub struct EvaluationWriter<E>
where
    E: Engine
{
    pub a: Vec<E::G1>,
    pub b_g1: Vec<E::G1>,
    pub b_g2: Vec<E::G2>,
    pub ic: Option<Vec<E::G1>>,
    pub l: Vec<E::G1>
}

impl<E> EvaluationWriter<E> 
where
    E: Engine
{
    pub fn new(key_pair: &KeyPairAssembly<E>) -> Self {
        let a = vec![E::G1::zero(); key_pair.num_inputs + key_pair.num_aux];
        let b_g1 = vec![E::G1::zero(); key_pair.num_inputs + key_pair.num_aux];
        let b_g2 = vec![E::G2::zero(); key_pair.num_inputs + key_pair.num_aux];
        let ic = vec![E::G1::zero(); key_pair.num_inputs];
        let l = vec![E::G1::zero(); key_pair.num_aux];
        
        Self { 
            a, 
            b_g1, 
            b_g2, 
            ic: Some(ic), 
            l 
        }
    }

    pub fn is_unconstrained(&self) -> bool {
        for e in self.l.iter() {
            if e.is_zero() {
                return true
            }
        }
        false
    }
    
    pub fn filter_non_zero_and_map_to_affine(self) -> (Vec<E::G1Affine>, Vec<E::G1Affine>, Vec<E::G1Affine>, Vec<E::G2Affine>) {
        let l: _ = map_g1_to_affine::<_,E>(self.l);
        let a: _ = filter_and_map_g1_to_affine::<_,E>(self.a);
        let b_g1: _ = filter_and_map_g1_to_affine::<_,E>(self.b_g1);
        let b_g2: _ = filter_and_map_g2_to_affine::<_,E>(self.b_g2);
            
        (l, a, b_g1, b_g2)
    }
}

fn map_g1_to_affine<I,E>(iter: I) -> Vec<E::G1Affine>
where
    I: IntoIterator<Item = E::G1>,
    E: Engine
{
    iter.into_iter()
        .map(|e| e.into_affine())
        .collect()
}

fn filter_and_map_g1_to_affine<I,E>(iter: I) -> Vec<E::G1Affine>
where
    I: IntoIterator<Item = E::G1>,
    E: Engine
{
    iter.into_iter()
        .filter(|e| !e.is_zero())
        .map(|e| e.into_affine())
        .collect()
}

fn filter_and_map_g2_to_affine<I,E>(iter: I) -> Vec<E::G2Affine>
where
    I: IntoIterator<Item = E::G2>,
    E: Engine
{
    iter.into_iter()
        .filter(|e| !e.is_zero())
        .map(|e| e.into_affine())
        .collect()
}
