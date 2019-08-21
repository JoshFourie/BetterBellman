//! Abstractions over the proving system and parameters.

use crate::{
    jubjub::{edwards, fs::Fs, Unknown},
    primitives::{Diversifier, PaymentAddress, ProofGenerationKey},
};
use pairing::bls12_381::{Bls12, Fr};

use crate::{
    merkle_tree::CommitmentTreeWitness,
    redjubjub::{PublicKey, Signature},
    sapling::Node,
    transaction::components::{Amount, GROTH_PROOF_SIZE},
};

/// Interface for creating zero-knowledge proofs for shielded transactions.
pub trait TxProver {
    /// Type for persisting any necessary context across multiple Sapling proofs.
    type Context;

    type Order;

    /// Instantiate a new Sapling proving context.
    fn new_sapling_proving_context(&self) -> Self::Context;

    fn order(
        proof_generation_key: ProofGenerationKey<Bls12>,
        diversifier: Diversifier,
        rcm: Fs,
        ar: Fs,
        value: u64,
        anchor: Fr,
        witness: CommitmentTreeWitness<Node>,
    ) -> Self::Order;

    /// Create the value commitment, re-randomized key, and proof for a Sapling
    /// [`SpendDescription`], while accumulating its value commitment randomness inside
    /// the context for later use.
    ///
    /// [`SpendDescription`]: crate::transaction::components::SpendDescription
    fn spend_proof(&self, ctx: &mut Self::Context, spend_order: Self::Order) -> Option<SpendProof>;

    /// Create the value commitment and proof for a Sapling [`OutputDescription`],
    /// while accumulating its value commitment randomness inside the context for later
    /// use.
    ///
    /// [`OutputDescription`]: crate::transaction::components::OutputDescription
    fn output_proof(
        &self,
        ctx: &mut Self::Context,
        esk: Fs,
        payment_address: PaymentAddress<Bls12>,
        rcm: Fs,
        value: u64,
    ) -> ([u8; GROTH_PROOF_SIZE], edwards::Point<Bls12, Unknown>);

    /// Create the `bindingSig` for a Sapling transaction. All calls to
    /// [`TxProver::spend_proof`] and [`TxProver::output_proof`] must be completed before
    /// calling this function.
    fn binding_sig(
        &self,
        ctx: &mut Self::Context,
        value_balance: Amount,
        sighash: &[u8; 32],
    ) -> Result<Signature, ()>;
}

pub struct SpendProof {
    buf: [u8; GROTH_PROOF_SIZE], 
    point: edwards::Point<Bls12, Unknown>, 
    pk: Option<PublicKey<Bls12>>
}

impl SpendProof {
    pub fn new(
        buf: [u8; GROTH_PROOF_SIZE], 
        point: edwards::Point<Bls12, Unknown>, 
        pk: Option<PublicKey<Bls12>>
    ) -> Self {
        Self { buf, point, pk }
    }
}

impl From<SpendProof> for ([u8; GROTH_PROOF_SIZE], edwards::Point<Bls12, Unknown>, Option<PublicKey<Bls12>>) {
    fn from(proof: SpendProof) -> Self {
        (proof.buf, proof.point, proof.pk)
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

        type Context = ();

        type Order = ();

        fn new_sapling_proving_context(&self) -> Self::Context {}

        fn spend_proof(&self, _ctx: &mut Self::Context, dep: Self::Order) -> Option<SpendProof> {
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

            Some(SpendProof::new([0u8; GROTH_PROOF_SIZE], cv, Some(rk)))
        }

        fn output_proof(
            &self,
            _ctx: &mut Self::Context,
            _esk: Fs,
            _payment_address: PaymentAddress<Bls12>,
            _rcm: Fs,
            value: u64,
        ) -> ([u8; GROTH_PROOF_SIZE], edwards::Point<Bls12, Unknown>) {
            let mut rng = OsRng;

            let cv = ValueCommitment::<Bls12> {
                value,
                randomness: Fs::random(&mut rng),
            }
            .cm(&JubjubBls12::new())
            .into();

            ([0u8; GROTH_PROOF_SIZE], cv)
        }

        fn binding_sig(
            &self,
            _ctx: &mut Self::Context,
            _value_balance: Amount,
            _sighash: &[u8; 32],
        ) -> Result<Signature, ()> {
            Err(())
        }

        fn order(
            proof_generation_key: ProofGenerationKey<Bls12>,
            diversifier: Diversifier,
            rcm: Fs,
            ar: Fs,
            value: u64,
            anchor: Fr,
            witness: CommitmentTreeWitness<Node>
        ) -> Self::Order { }
    }
}
