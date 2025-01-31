use crate::{domain, multi_thread, multicore};
use multicore::MULTI_THREAD;
use domain::Group;

use ff::{Field, ScalarEngine};

pub fn run_optimal_fft<E,T>(a: &mut [T], omega: &E::Fr, log_n: u32) 
where
    E: ScalarEngine,
    for <'a> T: Group<'a,E> 
{
    let log_cpus = MULTI_THREAD.log_num_cpus();

    if log_n <= log_cpus {
        serial_fft(a, omega, log_n);
    } else {
        parallel_fft(a, omega, log_n, log_cpus);
    }
}

pub fn serial_fft<E,T>(series: &mut [T], omega: &E::Fr, log_n: u32) 
where 
    E: ScalarEngine, 
    for <'a> T: Group<'a,E>
{
    let len: u32 = series.len() as u32;
    assert_eq!(len, 1 << log_n);

    for k in 0..len {
        let rk: u32 = bitreverse(k, log_n);
        if k < rk {
            series.swap(rk as usize, k as usize);
        }
    }

    let mut m: u32 = 1;
    for _ in 0..log_n {
        let w_m = omega.pow(&[(len / (2 * m)) as u64]);

        let mut k = 0;
        while k < len {
            let mut w: _ = E::Fr::one();
            for j in 0..m {
                let mut t: T = series[(k + j + m) as usize];
                t *= &w;

                let mut tmp: T = series[(k + j) as usize];
                tmp -= &t;

                series[(k + j + m) as usize] = tmp;
                series[(k + j) as usize] += &t;

                w.mul_assign(&w_m);
            }
            k += 2 * m;
        }
        m *= 2;
    }
}

pub fn parallel_fft<E,T>(a: &mut [T], omega: &E::Fr, log_n: u32, log_cpus: u32) 
where
    E: ScalarEngine,
    for <'a> T: Group<'a,E>
{
    assert!(log_n >= log_cpus);

    let num_cpus = 1 << log_cpus;
    let log_new_n = log_n - log_cpus;
    let mut tmp = vec![vec![T::zero(); 1 << log_new_n]; num_cpus];
    let new_omega = omega.pow(&[num_cpus as u64]);

    let ref_a: &_ = a; 
    multi_thread!(tmp.len(), enumerate(tmp) => {
        for (j, tmp) in iter => {
            // Shuffle into a sub-FFT
            let omega_j: _ = omega.pow(&[j as u64]);
            let omega_step: _ = omega.pow(&[(j as u64) << log_new_n]);

            let mut elt: _ = E::Fr::one();
            for i in 0..(1 << log_new_n) {
                for s in 0..num_cpus {
                    let idx = (i + (s << log_new_n)) % (1 << log_n);
                    let mut t = ref_a[idx];
                    t *= &elt;
                    tmp[i] += &t;
                    elt.mul_assign(&omega_step);
                }
                elt.mul_assign(&omega_j);
            }

            // Perform sub-FFT
            serial_fft(tmp, &new_omega, log_new_n);
        }
    });

    let tmp: &[_] = &tmp;
    let mask: _ = (1 << log_cpus) - 1;
    multi_thread!(a.len(), enumerate(a) => {
        for (i, val) in a => {
            let idx: usize = (i as u32 & mask) as usize;
            *val = tmp[idx][i >> log_cpus]
        } 
    });
}

fn bitreverse(mut n: u32, l: u32) -> u32 {
    let mut r = 0;
    for _ in 0..l {
        r = (r << 1) | (n & 1);
        n >>= 1;
    }
    r
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::domain::{Domain, Scalar};

    #[cfg(feature = "pairing")]
    #[test]
    fn fft_composition() {
        use pairing::bls12_381::Bls12;
        use rand_core::RngCore;

        fn test_comp<E: ScalarEngine, R: RngCore>(rng: &mut R) {
            for coeffs in 0..10 {
                let coeffs = 1 << coeffs;

                let mut v = vec![];
                for _ in 0..coeffs {
                    v.push(Scalar::<E>(E::Fr::random(rng)));
                }

                let mut domain = Domain::new(v.clone()).unwrap();

                domain.ifft();
                domain.fft();
                assert!(v == domain.as_coeffs());

                domain.fft();
                domain.ifft();
                assert!(v == domain.as_coeffs());

                domain.icoset_fft();
                domain.coset_fft();
                assert!(v == domain.as_coeffs());

                domain.coset_fft();
                domain.icoset_fft();
                assert!(v == domain.as_coeffs());
            }
        }

        let rng = &mut rand::thread_rng();

        test_comp::<Bls12, _>(rng);
    }
}
