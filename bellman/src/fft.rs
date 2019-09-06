use crate::{multicore, arith, domain};
use multicore::Worker;
use arith::{Scalar, Group};
use domain::EvaluationDomain;

use ff::{Field, ScalarEngine};

pub fn run_optimal_fft<E,T>(a: &mut [T], worker: &Worker, omega: &E::Fr, log_n: u32) 
where
    E: ScalarEngine,
    for <'a> T: Group<'a,E> 
{
    let log_cpus = worker.log_num_cpus();

    if log_n <= log_cpus {
        serial_fft(a, omega, log_n);
    } else {
        parallel_fft(a, worker, omega, log_n, log_cpus);
    }
}

pub fn serial_fft<E,T>(a: &mut [T], omega: &E::Fr, log_n: u32) 
where 
    E: ScalarEngine, 
    for <'a> T: Group<'a,E>
{
    fn bitreverse(mut n: u32, l: u32) -> u32 {
        let mut r = 0;
        for _ in 0..l {
            r = (r << 1) | (n & 1);
            n >>= 1;
        }
        r
    }

    let n = a.len() as u32;
    assert_eq!(n, 1 << log_n);

    for k in 0..n {
        let rk = bitreverse(k, log_n);
        if k < rk {
            a.swap(rk as usize, k as usize);
        }
    }

    let mut m = 1;
    for _ in 0..log_n {
        let w_m = omega.pow(&[(n / (2 * m)) as u64]);

        let mut k = 0;
        while k < n {
            let mut w = E::Fr::one();
            for j in 0..m {
                let mut t = a[(k + j + m) as usize];
                t *= &w;
                let mut tmp = a[(k + j) as usize];
                tmp -= &t;
                a[(k + j + m) as usize] = tmp;
                a[(k + j) as usize] += &t;
                w.mul_assign(&w_m);
            }

            k += 2 * m;
        }

        m *= 2;
    }
}

pub fn parallel_fft<E,T>(a: &mut [T], worker: &Worker, omega: &E::Fr, log_n: u32, log_cpus: u32) 
where
    E: ScalarEngine,
    for <'a> T: Group<'a,E>
{
    assert!(log_n >= log_cpus);

    let num_cpus = 1 << log_cpus;
    let log_new_n = log_n - log_cpus;
    let mut tmp = vec![vec![T::group_zero(); 1 << log_new_n]; num_cpus];
    let new_omega = omega.pow(&[num_cpus as u64]);

    worker.scope(0, |scope, _| {
        let a = &*a;

        for (j, tmp) in tmp.iter_mut().enumerate() {
            scope.spawn(move || {
                // Shuffle into a sub-FFT
                let omega_j = omega.pow(&[j as u64]);
                let omega_step = omega.pow(&[(j as u64) << log_new_n]);

                let mut elt = E::Fr::one();
                for i in 0..(1 << log_new_n) {
                    for s in 0..num_cpus {
                        let idx = (i + (s << log_new_n)) % (1 << log_n);
                        let mut t = a[idx];
                        t *= &elt;
                        tmp[i] += &t;
                        elt.mul_assign(&omega_step);
                    }
                    elt.mul_assign(&omega_j);
                }

                // Perform sub-FFT
                serial_fft(tmp, &new_omega, log_new_n);
            });
        }
    });

    // TODO: does this hurt or help?
    worker.scope(a.len(), |scope, chunk| {
        let tmp = &tmp;

        for (idx, a) in a.chunks_mut(chunk).enumerate() {
            scope.spawn(move || {
                let mut idx = idx * chunk;
                let mask = (1 << log_cpus) - 1;
                for a in a {
                    *a = tmp[idx & mask][idx >> log_cpus];
                    idx += 1;
                }
            });
        }
    });
}

#[cfg(feature = "pairing")]
#[test]
fn fft_composition() {
    use pairing::bls12_381::Bls12;
    use rand_core::RngCore;

    fn test_comp<E: ScalarEngine, R: RngCore>(rng: &mut R) {
        let worker = Worker::new();

        for coeffs in 0..10 {
            let coeffs = 1 << coeffs;

            let mut v = vec![];
            for _ in 0..coeffs {
                v.push(Scalar::<E>(E::Fr::random(rng)));
            }

            let mut domain = EvaluationDomain::new(v.clone(), &worker).unwrap();

            domain.ifft();
            domain.fft();
            assert!(v == domain.into_coeffs());

            domain.fft();
            domain.ifft();
            assert!(v == domain.into_coeffs());

            domain.icoset_fft();
            domain.coset_fft();
            assert!(v == domain.into_coeffs());

            domain.coset_fft();
            domain.icoset_fft();
            assert!(v == domain.into_coeffs());
        }
    }

    let rng = &mut rand::thread_rng();

    test_comp::<Bls12, _>(rng);
}
