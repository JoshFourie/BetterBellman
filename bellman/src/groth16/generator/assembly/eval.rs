use ff::{Field, PrimeField};
use group::CurveProjective;
use pairing::Engine;

use crate::{arith, error, multi_thread};
use arith::Scalar;
use error::Result;

use super::{key_pair, windows};
use key_pair::{KeyPairAssembly, KeyPairWires, FlatKeyPairWires};
use windows::BasedWindowTables;

pub fn eval<E: Engine>(
    wnaf: &BasedWindowTables<'_,E>,
    lagrange_coeffs: &[Scalar<E>],
    qap_polynomials: KeyPairWires<E>,
    writer: EvaluationWriter<'_,E>,

    // Inverse coefficient for ext elements
    inv: &E::Fr,

    // Trapdoors
    alpha: &E::Fr,
    beta: &E::Fr,
) {
    // // Sanity check
    // assert_eq!(writer.a.len(), qap_polynomials.at.len());
    // assert_eq!(writer.a.len(), qap_polynomials.bt.len());
    // assert_eq!(writer.a.len(), qap_polynomials.ct.len());
    // assert_eq!(writer.a.len(), writer.b_g1.len());
    // assert_eq!(writer.a.len(), writer.b_g2.len());
    // assert_eq!(writer.a.len(), writer.ext.len());

    let coeff_len: usize = writer.a.len();
    let mut flat_writer: FlatEvaluationWriter<E> = writer.flatten();
    let flat_poly: FlatKeyPairWires<E> = qap_polynomials.flatten();

    multi_thread!(coeff_len, iter(flat_writer, flat_poly) => {
        for ((a, b_g1, b_g2, ext), (at, bt, ct)) in writer, poly => {

            let mut g1_wnaf = wnaf.g1.shared();
            let mut g2_wnaf = wnaf.g2.shared();

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
            let mut at = eval_at_tau(lagrange_coeffs, at);
            let mut bt = eval_at_tau(lagrange_coeffs, bt);
            let ct = eval_at_tau(lagrange_coeffs, ct);

            // Compute A query (in G1)
            if !at.is_zero() {
                **a = g1_wnaf.scalar(at.into_repr());
            }

            // Compute B query (in G1/G2)
            if !bt.is_zero() {
                let bt_repr = bt.into_repr();
                **b_g1 = g1_wnaf.scalar(bt_repr);
                **b_g2 = g2_wnaf.scalar(bt_repr);
            }

            at.mul_assign(&beta);
            bt.mul_assign(&alpha);

            let mut e = at;
            e.add_assign(&bt);
            e.add_assign(&ct);
            e.mul_assign(inv);

            **ext = g1_wnaf.scalar(e.into_repr());
        }

        // Batch normalize
        map_to_chunk!{
            // E::G1::batch_normalization(writer.a);
            // E::G1::batch_normalization(writer.b_g1);
            // E::G2::batch_normalization(writer.b_g2);
            // E::G1::batch_normalization(writer.ext);
        }
    });
}

pub struct WireEvaluation<E>
where
    E: Engine
{
    pub a: Vec<E::G1>,
    pub b_g1: Vec<E::G1>,
    pub b_g2: Vec<E::G2>,
    pub ic: Option<Vec<E::G1>>,
    pub l: Vec<E::G1>
}

impl<E> WireEvaluation<E> 
where
    E: Engine
{
    pub fn new(key_pair: &KeyPairAssembly<E>) -> Self {
        let a = vec![E::G1::zero(); key_pair.num.inputs + key_pair.num.aux];
        let b_g1 = vec![E::G1::zero(); key_pair.num.inputs + key_pair.num.aux];
        let b_g2 = vec![E::G2::zero(); key_pair.num.inputs + key_pair.num.aux];
        let ic = vec![E::G1::zero(); key_pair.num.inputs];
        let l = vec![E::G1::zero(); key_pair.num.aux];
        
        WireEvaluation { 
            a, 
            b_g1, 
            b_g2, 
            ic: Some(ic), 
            l 
        }
    }

    pub fn as_mut_auxilliaries(&mut self, aux_bound: usize) -> EvaluationWriter<E> {
        EvaluationWriter::new(
            &mut self.a[aux_bound..],
            &mut self.b_g1[aux_bound..],
            &mut self.b_g2[aux_bound..],
            &mut self.l
        )
    }

    pub fn as_mut_inputs<'a>(&'a mut self, input_bound: usize) -> Result<EvaluationWriter<'a,E>> {
        Ok(EvaluationWriter::new(
            &mut self.a[0..input_bound],
            &mut self.b_g1[0..input_bound],
            &mut self.b_g2[0..input_bound],
            self.ic.as_mut()?
        ))
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

pub struct EvaluationWriter<'a, E: Engine> {
    a: &'a mut [E::G1],
    b_g1: &'a mut [E::G1],
    b_g2: &'a mut [E::G2],
    ext: &'a mut [E::G1]
} 

impl<'a,E> EvaluationWriter<'a,E> 
where
    E: Engine
{
    fn new(a: &'a mut [E::G1], b_g1: &'a mut [E::G1], b_g2: &'a mut [E::G2], ext: &'a mut [E::G1]) -> Self {
        Self {
            a,
            b_g1,
            b_g2,
            ext   
        }
    }

    fn flatten(self) -> FlatEvaluationWriter<'a,E> {
        FlatEvaluationWriter::from(self)
    }
}

struct FlatEvaluationWriter<'a,E: Engine>(Vec<(&'a mut E::G1, &'a mut E::G1, &'a mut E::G2, &'a mut E::G1)>);

impl<'a,E> FlatEvaluationWriter<'a,E> 
where
    E: Engine
{
    fn chunks_mut(&mut self, chunk_size: usize) -> std::slice::ChunksMut<'_, (&'a mut E::G1, &'a mut E::G1, &'a mut E::G2, &'a mut E::G1)> {
        self.0.chunks_mut(chunk_size)
    }
}

impl<'a,E> From <EvaluationWriter<'a,E>> for FlatEvaluationWriter<'a,E> 
where
    E: Engine
{
    fn from(writer: EvaluationWriter<'a, E>) -> Self {
        let flattened: Vec<_> = writer.a.into_iter()
            .zip(writer.b_g1.into_iter())
            .zip(writer.b_g2.into_iter())
            .zip(writer.ext.into_iter())
            .map(|(((a, b_g1), b_g2), ext)| {
                (a, b_g1, b_g2, ext)
            }).collect();
        FlatEvaluationWriter(flattened)
    }
}
