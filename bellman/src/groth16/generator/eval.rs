use ff::{Field, PrimeField};
use group::{CurveProjective, Wnaf};
use pairing::Engine;

use crate::multicore::Worker;
use crate::arith::Scalar;

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