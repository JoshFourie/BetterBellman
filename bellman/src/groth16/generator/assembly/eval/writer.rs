use ff::{Field, PrimeField};
use group::CurveProjective;
use pairing::Engine;

use crate::{domain, multi_thread};
use domain::Scalar;

use super::super::{key_pair, windows, parameters};
use parameters::Elements;
use key_pair::{KeyPairWires, FlatKeyPairWires};
use windows::BasedWindows;

pub struct Writer<'a, E: Engine> {
    a: &'a mut [E::G1],
    b_g1: &'a mut [E::G1],
    b_g2: &'a mut [E::G2],
    ext: &'a mut [E::G1]
} 

impl<'a,E> Writer<'a,E> 
where
    E: Engine
{
    pub fn new(a: &'a mut [E::G1], b_g1: &'a mut [E::G1], b_g2: &'a mut [E::G2], ext: &'a mut [E::G1]) -> Self {
        Writer { a, b_g1, b_g2, ext }
    }

    pub fn sanity_check(&self, qap_polynomials: &KeyPairWires<E>) -> bool 
    where
        E: Engine
    {
        self.a.len() == qap_polynomials.at.len() &&
        self.a.len() == qap_polynomials.bt.len() &&
        self.a.len() == qap_polynomials.ct.len() &&
        self.a.len() == self.b_g1.len() &&
        self.a.len() == self.b_g2.len() &&
        self.a.len() == self.ext.len()
    }

    pub fn eval(self, wnaf: &BasedWindows<'_,E>, coeffs: &[Scalar<E>], qap: KeyPairWires<E>, inverse_coeff: &E::Fr, trapdoors: &Elements<E>) {

        let coeff_len: usize = self.a.len();
        let mut flat_writer: FlatWriter<E> = self.flatten();
        let flat_poly: FlatKeyPairWires<E> = qap.flatten();

        multi_thread!(coeff_len, iter(flat_writer, flat_poly) => {
            for ((a, b_g1, b_g2, ext), (at, bt, ct)) in writer, poly => {

                let mut g1_wnaf = wnaf.g1.shared();
                let mut g2_wnaf = wnaf.g2.shared();

                // Evaluate QAP polynomials at tau
                let mut at = eval_at_tau(coeffs, at);
                let mut bt = eval_at_tau(coeffs, bt);
                let ct = eval_at_tau(coeffs, ct);

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

                at.mul_assign(&trapdoors.beta);
                bt.mul_assign(&trapdoors.alpha);

                let mut e = at;
                e.add_assign(&bt);
                e.add_assign(&ct);
                e.mul_assign(inverse_coeff);

                **ext = g1_wnaf.scalar(e.into_repr());
            }
        });

        flat_writer.batch_normalization();
    }


    fn flatten(self) -> FlatWriter<'a,E> { FlatWriter::from(self) }
}

fn eval_at_tau<E>(powers_of_tau: &[Scalar<E>], wires: &[(E::Fr, usize)]) -> E::Fr 
where
    E: Engine
{
    wires.iter()
        .fold(E::Fr::zero(), |mut acc, (coeff, idx)| {
            let Scalar(mut exp): Scalar<E> = powers_of_tau[*idx];
            exp.mul_assign(coeff);
            acc.add_assign(&exp);
            acc
        })
}

struct FlatWriter<'a,E: Engine>(Vec<(&'a mut E::G1, &'a mut E::G1, &'a mut E::G2, &'a mut E::G1)>);

impl<'a,E> FlatWriter<'a,E> 
where
    E: Engine
{
    fn chunks_mut(&mut self, chunk_size: usize) -> std::slice::ChunksMut<'_, (&'a mut E::G1, &'a mut E::G1, &'a mut E::G2, &'a mut E::G1)> {
        self.0.chunks_mut(chunk_size)
    }

    fn batch_normalization(self) {
        let mut buf_a: Vec<E::G1> = Vec::new();
        let mut buf_b_g1: Vec<E::G1> = Vec::new();
        let mut buf_b_g2: Vec<E::G2> = Vec::new();
        let mut buf_ext: Vec<E::G1> = Vec::new();

        for (a, b_g1, b_g2, ext) in self.0.into_iter() {
            buf_a.push(*a);
            buf_b_g1.push(*b_g1);
            buf_b_g2.push(*b_g2);
            buf_ext.push(*ext);
        }

        E::G1::batch_normalization(&mut buf_a);
        E::G1::batch_normalization(&mut buf_b_g1);
        E::G2::batch_normalization(&mut buf_b_g2);
        E::G1::batch_normalization(&mut buf_ext);
    }
}

impl<'a,E> From <Writer<'a,E>> for FlatWriter<'a,E> 
where
    E: Engine
{
    fn from(writer: Writer<'a, E>) -> Self {
        let flattened: Vec<_> = writer.a.into_iter()
            .zip(writer.b_g1.into_iter())
            .zip(writer.b_g2.into_iter())
            .zip(writer.ext.into_iter())
            .map(|(((a, b_g1), b_g2), ext)| (a, b_g1, b_g2, ext))
            .collect();
        FlatWriter(flattened)
    }
}