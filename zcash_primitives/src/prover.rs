use crate::jubjub::{edwards, Unknown};
use pairing::bls12_381::Bls12;

use crate::{
    redjubjub::{PublicKey, Signature},
    transaction::components::{Amount, GROTH_PROOF_SIZE},
};

/// Interface for creating zero-knowledge proofs for shielded transactions.
pub trait TxProver {
    
    fn new() -> Self;

    fn spend_proof<T>(&mut self, order: T) -> Option<SpendProof>;

    fn output_proof<U>(&mut self, order: U) -> Option<OutputProof>;

    fn binding_sig(&mut self, value_balance: Amount, sighash: &[u8; 32]) -> Option<Signature>;
}

pub struct SpendProof {
    buf: [u8; GROTH_PROOF_SIZE], 
    point: edwards::Point<Bls12, Unknown>, 
    pk: PublicKey<Bls12>
}

impl SpendProof {
    pub fn new(buf: [u8; GROTH_PROOF_SIZE], point: edwards::Point<Bls12, Unknown>, pk: PublicKey<Bls12>) -> Self {
        Self { buf, point, pk }
    }

    pub fn into_tuple(self) -> ([u8; GROTH_PROOF_SIZE], edwards::Point<Bls12, Unknown>, PublicKey<Bls12>) {
        (self.buf, self.point, self.pk)
    }
}

pub struct OutputProof {
    buf: [u8; GROTH_PROOF_SIZE], 
    point: edwards::Point<Bls12, Unknown>, 
}

impl OutputProof {
    pub fn new(buf: [u8; GROTH_PROOF_SIZE], point: edwards::Point<Bls12, Unknown>) -> Self {
        Self { buf, point }
    }

    pub fn into_tuple(self) -> ([u8; GROTH_PROOF_SIZE], edwards::Point<Bls12, Unknown>) {
        (self.buf, self.point)
    }
}


#[cfg(test)]
pub(crate) mod mock {
    use ff::Field;
    use pairing::bls12_381::{Bls12, Fr};
    use rand_os::OsRng;
    use super::*;

    use crate::{
        jubjub::{edwards, fs::Fs, FixedGenerators, Unknown},
        primitives::{Diversifier, PaymentAddress, ProofGenerationKey, ValueCommitment},
    };

    use crate::{
        merkle_tree::CommitmentTreeWitness,
        redjubjub::{PublicKey, Signature},
        sapling::Node,
        transaction::components::{Amount, GROTH_PROOF_SIZE},
        JubjubBls12,
    };

    pub(crate) struct MockTxProver;

    #[cfg(test)]
    impl TxProver for MockTxProver {

        fn new() -> Self { MockTxProver }

        fn spend_proof<T>(&mut self, spend_proof_details: T) -> Option<SpendProof> {
            let mut rng = OsRng;
            let value: u64 = unsafe { std::mem::MaybeUninit::zeroed().assume_init() };
            let proof_generation_key: ProofGenerationKey<Bls12> = unsafe { std::mem::MaybeUninit::zeroed().assume_init() };
            let ar: Fs = unsafe{ std::mem::MaybeUninit::zeroed().assume_init() };
            
            let entropy: _ = Fs::random(&mut rng);
            let cv: edwards::Point<Bls12, Unknown> = ValueCommitment::new(value, entropy)
                .cm(&JubjubBls12::new())
                .into();

            let rk = PublicKey::<Bls12>(proof_generation_key.ak.clone().into()).randomize(
                ar,
                FixedGenerators::SpendingKeyGenerator,
                &JubjubBls12::new(),
            );

            Some(SpendProof::new([0u8; GROTH_PROOF_SIZE], cv, rk))
        }

        fn output_proof<U>(&mut self, output_proof_details: U) -> Option<OutputProof> {

            let mut rng = OsRng;

            let cv = ValueCommitment::<Bls12> {
                value: unsafe { std::mem::MaybeUninit::zeroed().assume_init() },
                randomness: Fs::random(&mut rng),
            }
            .cm(&JubjubBls12::new())
            .into();

            Some(OutputProof::new([0u8; GROTH_PROOF_SIZE], cv))
        }

        fn binding_sig(&mut self, value_balance: Amount, sighash: &[u8; 32]) -> Option<Signature> {
            None
        }
    }
}
