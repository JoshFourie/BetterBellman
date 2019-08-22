mod output;

use crate::{
    jubjub::fs::Fs,
    primitives::{Diversifier, Note, PaymentAddress},
};
use ff::Field;
use pairing::bls12_381::{Bls12, Fr};
use rand::{rngs::OsRng, seq::SliceRandom, CryptoRng, RngCore};
use zip32::ExtendedSpendingKey;

use crate::{
    keys::OutgoingViewingKey,
    legacy::TransparentAddress,
    merkle_tree::{CommitmentTreeWitness, IncrementalWitness},
    note_encryption::{generate_esk, Memo, SaplingNoteEncryption},
    prover::TxProver,
    redjubjub::PrivateKey,
    sapling::{spend_sig, Node},
    transaction::{
        components::{amount::DEFAULT_FEE, Amount, OutputDescription, SpendDescription, TxOut},
        signature_hash_data, Transaction, TransactionData, SIGHASH_ALL,
    },
    JUBJUB,
};

const DEFAULT_TX_EXPIRY_DELTA: u32 = 20;

/// If there are any shielded inputs, always have at least two shielded outputs, padding
/// with dummy outputs if necessary. See https://github.com/zcash/zcash/issues/3615
const MIN_SHIELDED_OUTPUTS: usize = 2;

type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, PartialEq)]
pub enum Error {
    AnchorMismatch,
    BindingSig,
    ChangeIsNegative(Amount),
    InvalidAddress,
    InvalidAmount,
    InvalidWitness,
    NoChangeAddress,
    SpendProof,
    OperationOnNull(&'static str)
}

/// Generates a [`Transaction`] from its inputs and outputs.
pub struct Builder<R: RngCore + CryptoRng> {
    rng: R,
    mtx: TransactionData,
    fee: Amount,
    anchor: Option<Fr>,
    spends: Option<Vec<SpendDescriptionInfo>>,
    outputs: Option<Vec<output::SaplingOutput>>,
    change_address: Option<(OutgoingViewingKey, PaymentAddress<Bls12>)>,
}

impl Builder<OsRng> {
    /// Creates a new `Builder` targeted for inclusion in the block with the given height,
    /// using default values for general transaction fields and the default OS random.
    ///
    /// # Default values
    ///
    /// The expiry height will be set to the given height plus the default transaction
    /// expiry delta (20 blocks).
    ///
    /// The fee will be set to the default fee (0.0001 ZEC).
    pub fn new(height: u32) -> Self {
        Builder::new_with_rng(height, OsRng)
    }
}

impl<R> Builder<R> 
where
    R: RngCore + CryptoRng
{
    /// Creates a new `Builder` targeted for inclusion in the block with the given height
    /// and randomness source, using default values for general transaction fields.
    ///
    /// # Default values
    ///
    /// The expiry height will be set to the given height plus the default transaction
    /// expiry delta (20 blocks).
    ///
    /// The fee will be set to the default fee (0.0001 ZEC).
    pub fn new_with_rng(height: u32, rng: R) -> Self {
        let mut mtx = TransactionData::new();
        mtx.expiry_height = height + DEFAULT_TX_EXPIRY_DELTA;

        Builder {
            rng,
            mtx,
            fee: DEFAULT_FEE,
            anchor: None,
            spends: Some(Vec::new()),
            outputs: Some(Vec::new()),
            change_address: None,
        }
    }

    /// Adds a Sapling note to be spent in this transaction.
    ///
    /// Returns an error if the given witness does not have the same anchor as previous
    /// witnesses, or has no path.
    pub fn add_sapling_spend(
        &mut self,
        extsk: ExtendedSpendingKey,
        diversifier: Diversifier,
        note: Note<Bls12>,
        witness: IncrementalWitness<Node>,
    ) -> Result<()> {
        // Consistency check: all anchors must equal the first one
        if let Some(anchor) = self.anchor {
            let witness_root: Fr = witness.root().into();
            if witness_root != anchor {
                return Err(Error::AnchorMismatch);
            }
        } else {
            self.anchor = Some(witness.root().into())
        }
        let witness = witness.path().ok_or(Error::InvalidWitness)?;

        let alpha = Fs::random(&mut self.rng);

        self.mtx.value_balance += Amount::from_u64(note.value).map_err(|_| Error::InvalidAmount)?;

        // F LAG
        self.spends
            .as_mut()
            .map(|spends| {
                let spend_description: _ = SpendDescriptionInfo {
                    extsk,
                    diversifier,
                    note,
                    alpha,
                    witness,
                };
                spends.push(spend_description)
            }).ok_or(Error::OperationOnNull("Expected Empty Vec to Push SpendDesriptionInfo"))?;

        Ok(())
    }

    /// Adds a Sapling address to send funds to.
    pub fn add_sapling_output(
        &mut self,
        ovk: OutgoingViewingKey,
        to: PaymentAddress<Bls12>,
        value: Amount,
        memo: Option<Memo>,
    ) -> Result<()> {
        let output = output::SaplingOutput::new(&mut self.rng, ovk, to, value, memo)?;

        self.mtx.value_balance -= value;

        self.outputs
            .as_mut()
            .map(|outputs| outputs.push(output))
            .ok_or(Error::OperationOnNull("Expected Empty Vec to Push SaplingOutput"))
    }

    /// Adds a transparent address to send funds to.
    pub fn add_transparent_output(
        &mut self,
        to: &TransparentAddress,
        value: Amount,
    ) -> Result<()> {
        if value.is_negative() {
            return Err(Error::InvalidAmount);
        }

        self.mtx.vout.push(TxOut {
            value,
            script_pubkey: to.script(),
        });

        Ok(())
    }

    /// Sets the Sapling address to which any change will be sent.
    ///
    /// By default, change is sent to the Sapling address corresponding to the first note
    /// being spent (i.e. the first call to [`Builder::add_sapling_spend`]).
    pub fn send_change_to(&mut self, ovk: OutgoingViewingKey, to: PaymentAddress<Bls12>) {
        self.change_address = Some((ovk, to));
    }

    
    #[inline]
    fn check_and_fix_consistency(&mut self) -> Result<()> {
        let change: Amount = self.calculate_change();
        if change.is_positive() {
            self.send_change_to_address_or_default(change)?;
            Ok(())
        } else if change.is_negative() {
            Err(Error::ChangeIsNegative(change))
        } else { unimplemented!() }
    }

    #[inline]
    fn calculate_change(&self) -> Amount {
        let balance: _ = self.mtx.value_balance;
        let fee: _ = self.fee;
        let sigma: Amount = self.mtx
            .vout
            .iter()
            .map(|output| output.value)
            .sum();
        balance - fee - sigma
    }

    // todo: fix naming issue with 'send_change_to'
    #[inline]
    fn send_change_to_address_or_default(&mut self, change: Amount) -> Result<()> {
        let (ovk, addr): _ = self.change_address
            .take()
            .or_else(|| {
                self.spends
                    .as_ref()
                    .and_then(|spends| {
                        if !spends.is_empty() {
                            let ovk: _ = spends[0].extsk.expsk.ovk;
                            let pay_addr: _ = PaymentAddress {
                                diversifier: spends[0].diversifier,
                                pk_d: spends[0].note.pk_d.clone(),
                            };
                            Some((ovk, pay_addr))
                        } else { None }
                    })
            }).ok_or(Error::NoChangeAddress)?;
        self.add_sapling_output(ovk, addr, change, None)?;

        Ok(())
    }

    fn take_description_and_map_index(&mut self) -> Result<(Vec<(usize, SpendDescriptionInfo)>, Vec<Option<(usize, output::SaplingOutput)>>)> {
        let spends: Vec<(_, SpendDescriptionInfo)> = self.spends
            .take()
            .ok_or(Error::OperationOnNull("Expected Vec<SpendDescriptionInfo> to Map Index"))?
            .into_iter()
            .enumerate()
            .collect();

        let outputs: Vec<_> = self.outputs
            .take()
            .ok_or(Error::OperationOnNull("Expected Vec<SpendDescriptionInfo> to Map Optional Index"))?
            .into_iter()
            .enumerate()
            .map(|(i, o)| Some((i, o)))
            .collect();

        Ok((spends, outputs))
    } 

    /// Builds a transaction from the configured spends and outputs.
    ///
    /// Upon success, returns a tuple containing the final transaction, and the
    /// [`TransactionMetadata`] generated during the build process.
    ///
    /// `consensus_branch_id` must be valid for the block height that this transaction is
    /// targeting. An invalid `consensus_branch_id` will *not* result in an error from
    /// this function, and instead will generate a transaction that will be rejected by
    /// the network.
    pub fn build<T>(mut self, consensus_branch_id: u32, mut prover: T) -> Result<(Transaction, TransactionMetadata)> 
    where
        T: TxProver
    {
        let mut tx_metadata = TransactionMetadata::new();

        self.check_and_fix_consistency()?;

        let (mut spends, mut outputs): _ = self.take_description_and_map_index()?;

        let anchor = self.anchor.expect("anchor was set if spends were added");

        // Pad Sapling outputs
        let orig_outputs_len = outputs.len();
        if !spends.is_empty() {
            while outputs.len() < MIN_SHIELDED_OUTPUTS {
                outputs.push(None);
            }
        }

        // Randomize order of inputs and outputs
        spends.shuffle(&mut self.rng);
        outputs.shuffle(&mut self.rng);
        tx_metadata.spend_indices.resize(spends.len(), 0);
        tx_metadata.output_indices.resize(orig_outputs_len, 0);

        // Create Sapling SpendDescriptions
        for (i, (pos, spend)) in spends.iter().enumerate() {
            let proof_generation_key = spend.extsk.expsk.proof_generation_key(&JUBJUB);

            let mut nullifier = [0u8; 32];
            nullifier.copy_from_slice(&spend.note.nf(
                &proof_generation_key.into_viewing_key(&JUBJUB),
                spend.witness.position,
                &JUBJUB,
            ));

            let tx_order: _ = (                    
                proof_generation_key,
                spend.diversifier,
                spend.note.r,
                spend.alpha,
                spend.note.value,
                anchor,
                spend.witness.clone()
            );

            let (zkproof, cv, rk) = prover.spend_proof(tx_order)
                .ok_or(Error::SpendProof)?
                .into_tuple();

            self.mtx.shielded_spends.push(SpendDescription {
                cv,
                anchor: anchor,
                nullifier,
                rk: rk,
                zkproof,
                spend_auth_sig: None,
            });

            // Record the post-randomized spend location
            tx_metadata.spend_indices[*pos] = i;
        }

        // Create Sapling OutputDescriptions
        for (i, output) in outputs.into_iter().enumerate() {
            let output_desc = if let Some((pos, output)) = output {
                // Record the post-randomized output location
                tx_metadata.output_indices[pos] = i;

                output.build(&mut prover, &mut self.rng)
            } else {
                // This is a dummy output
                let (dummy_to, dummy_note) = {
                    let (diversifier, g_d) = {
                        let mut diversifier;
                        let g_d;
                        loop {
                            let mut d = [0; 11];
                            self.rng.fill_bytes(&mut d);
                            diversifier = Diversifier(d);
                            if let Some(val) = diversifier.g_d::<Bls12>(&JUBJUB) {
                                g_d = val;
                                break;
                            }
                        }
                        (diversifier, g_d)
                    };

                    let pk_d = {
                        let dummy_ivk = Fs::random(&mut self.rng);
                        g_d.mul(dummy_ivk, &JUBJUB)
                    };

                    (
                        PaymentAddress {
                            diversifier,
                            pk_d: pk_d.clone(),
                        },
                        Note {
                            g_d,
                            pk_d,
                            r: Fs::random(&mut self.rng),
                            value: 0,
                        },
                    )
                };

                let esk = generate_esk(&mut self.rng);
                let epk = dummy_note.g_d.mul(esk, &JUBJUB);

                let (zkproof, cv) = prover.output_proof((
                    esk, 
                    dummy_to, 
                    dummy_note.r, 
                    dummy_note.value
                )).ok_or(Error::SpendProof)?.into_tuple();

                let cmu = dummy_note.cm(&JUBJUB);

                let mut enc_ciphertext = [0u8; 580];
                let mut out_ciphertext = [0u8; 80];
                self.rng.fill_bytes(&mut enc_ciphertext[..]);
                self.rng.fill_bytes(&mut out_ciphertext[..]);

                Some(OutputDescription {
                    cv,
                    cmu,
                    ephemeral_key: epk.into(),
                    enc_ciphertext,
                    out_ciphertext,
                    zkproof,
                })
            };

            self.mtx.shielded_outputs.push(output_desc.ok_or(Error::SpendProof)?);
        }

        //
        // Signatures
        //

        let mut sighash = [0u8; 32];
        sighash.copy_from_slice(&signature_hash_data(
            &self.mtx,
            consensus_branch_id,
            SIGHASH_ALL,
            None,
        ));

        // Create Sapling spendAuth and binding signatures
        for (i, (_, spend)) in spends.into_iter().enumerate() {
            self.mtx.shielded_spends[i].spend_auth_sig = Some(spend_sig(
                PrivateKey(spend.extsk.expsk.ask),
                spend.alpha,
                &sighash,
                &mut self.rng,
                &JUBJUB,
            ));
        }
        self.mtx.binding_sig = Some(
            prover.binding_sig(self.mtx.value_balance, &sighash).ok_or(Error::BindingSig)?
        );

        Ok((
            self.mtx.freeze().expect("Transaction should be complete"),
            tx_metadata,
        ))
    }
}

struct SpendDescriptionInfo {
    extsk: ExtendedSpendingKey,
    diversifier: Diversifier,
    note: Note<Bls12>,
    alpha: Fs,
    witness: CommitmentTreeWitness<Node>,
}

/// Metadata about a transaction created by a [`Builder`].
#[derive(Debug, PartialEq)]
pub struct TransactionMetadata {
    spend_indices: Vec<usize>,
    output_indices: Vec<usize>,
}

impl TransactionMetadata {
    fn new() -> Self {
        TransactionMetadata {
            spend_indices: vec![],
            output_indices: vec![],
        }
    }

    /// Returns the index within the transaction of the [`SpendDescription`] corresponding
    /// to the `n`-th call to [`Builder::add_sapling_spend`].
    ///
    /// Note positions are randomized when building transactions for indistinguishability.
    /// This means that the transaction consumer cannot assume that e.g. the first spend
    /// they added (via the first call to [`Builder::add_sapling_spend`]) is the first
    /// [`SpendDescription`] in the transaction.
    pub fn spend_index(&self, n: usize) -> Option<usize> {
        self.spend_indices.get(n).map(|i| *i)
    }

    /// Returns the index within the transaction of the [`OutputDescription`] corresponding
    /// to the `n`-th call to [`Builder::add_sapling_output`].
    ///
    /// Note positions are randomized when building transactions for indistinguishability.
    /// This means that the transaction consumer cannot assume that e.g. the first output
    /// they added (via the first call to [`Builder::add_sapling_output`]) is the first
    /// [`OutputDescription`] in the transaction.
    pub fn output_index(&self, n: usize) -> Option<usize> {
        self.output_indices.get(n).map(|i| *i)
    }
}

#[cfg(test)]
mod tests {
    use ff::{Field, PrimeField};
    use rand::rngs::OsRng;

    use crate::jubjub::fs::Fs;

    use super::{Builder, Error};
    use crate::{
        legacy::TransparentAddress,
        merkle_tree::{CommitmentTree, IncrementalWitness},
        prover::mock::MockTxProver,
        sapling::Node,
        transaction::components::Amount,
        zip32::{ExtendedFullViewingKey, ExtendedSpendingKey},
        JUBJUB,
    };

    #[test]
    fn fails_on_negative_output() {
        let extsk = ExtendedSpendingKey::master(&[]);
        let extfvk = ExtendedFullViewingKey::from(&extsk);
        let ovk = extfvk.fvk.ovk;
        let to = extfvk.default_address().unwrap().1;

        let mut builder = Builder::new(0);
        assert_eq!(
            builder.add_sapling_output(ovk, to, Amount::from_i64(-1).unwrap(), None),
            Err(Error::InvalidAmount)
        );
    }

    #[test]
    fn fails_on_negative_transparent_output() {
        let mut builder = Builder::new(0);
        assert_eq!(
            builder.add_transparent_output(
                &TransparentAddress::PublicKey([0; 20]),
                Amount::from_i64(-1).unwrap(),
            ),
            Err(Error::InvalidAmount)
        );
    }

    #[cfg(test)]
    fn get_rng_and_key() -> (rand::rngs::OsRng, crate::zip32::ExtendedSpendingKey) {
        let rng = OsRng;
        // Just use the master key as the ExtendedSpendingKey for this test
        let extsk = ExtendedSpendingKey::master(&[]);   
        (rng, extsk)
    }

    #[test]
    fn fails_with_no_input_or_output() {
        // Fails with no inputs or outputs
        // 0.0001 t-ZEC fee
        let builder = Builder::new(0);
        assert_eq!(
            builder.build(1, MockTxProver),
            Err(Error::ChangeIsNegative(Amount::from_i64(-10000).unwrap()))
        );   
    }

    #[test]
    fn fail_with_only_single_output() {
        let (_, extsk): _ = get_rng_and_key();
        let extfvk = ExtendedFullViewingKey::from(&extsk);
        let ovk = extfvk.fvk.ovk;
        let to = extfvk.default_address().unwrap().1;

        let mut builder = Builder::new(0);
        builder.add_sapling_output(
                ovk.clone(),
                to.clone(),
                Amount::from_u64(50000).unwrap(),
                None,
            ).unwrap();

        assert_eq!(
            builder.build(1, MockTxProver),
            Err(Error::ChangeIsNegative(Amount::from_i64(-60000).unwrap()))
        );
    }

    #[test]
    fn fail_with_only_transparent_output() {
        let mut builder = Builder::new(0);
        builder.add_transparent_output(
            &TransparentAddress::PublicKey([0; 20]),
            Amount::from_u64(50000).unwrap(),
        ).unwrap();

        assert_eq!(
            builder.build(1, MockTxProver),
            Err(Error::ChangeIsNegative(Amount::from_i64(-60000).unwrap()))
        );
    }

    #[test]
    fn fail_with_insufficient_input() {
        let (mut rng, extsk): _ = get_rng_and_key();

        let extfvk = ExtendedFullViewingKey::from(&extsk);
        let ovk = extfvk.fvk.ovk;
        let to = extfvk.default_address().unwrap().1;

        let note1 = to.create_note(59999, Fs::random(&mut rng), &JUBJUB).unwrap();
        let cm1 = Node::new(note1.cm(&JUBJUB).into_repr());
        let mut tree = CommitmentTree::new();
        tree.append(cm1).unwrap();
        let witness1 = IncrementalWitness::from_tree(&tree);

        let mut builder = Builder::new(0);
        builder.add_sapling_spend(
            extsk.clone(),
            to.diversifier,
            note1.clone(),
            witness1.clone(),
        ).unwrap();

        builder.add_sapling_output(
            ovk.clone(),
            to.clone(),
            Amount::from_u64(30000).unwrap(),
            None,
        ).unwrap();
        builder.add_transparent_output(
            &TransparentAddress::PublicKey([0; 20]),
            Amount::from_u64(20000).unwrap(),
        ).unwrap();

        assert_eq!(
            builder.build(1, MockTxProver),
            Err(Error::ChangeIsNegative(Amount::from_i64(-1).unwrap()))
        );
    }

    #[test]
    // There was a note that the MockTxProver fails because of BindingSig,
    // but the current failure is that the 'consistency check' didn't account
    // for an amount that was neither positive or negative (ie zero):
    // both returned false on a zero evaluation
    fn succeed_with_sufficient_inputs() {
        let (mut rng, extsk): _ = get_rng_and_key();

        let extfvk = ExtendedFullViewingKey::from(&extsk);
        let ovk = extfvk.fvk.ovk;
        let to = extfvk.default_address().unwrap().1;

        let note1 = to.create_note(59999, Fs::random(&mut rng), &JUBJUB).unwrap();
        let cm1 = Node::new(note1.cm(&JUBJUB).into_repr());
        let mut tree = CommitmentTree::new();
        tree.append(cm1).unwrap();
        let mut witness1 = IncrementalWitness::from_tree(&tree);

        let note2 = to.create_note(1, Fs::random(&mut rng), &JUBJUB).unwrap();
        let cm2 = Node::new(note2.cm(&JUBJUB).into_repr());
        tree.append(cm2).unwrap();
        witness1.append(cm2).unwrap();
        let witness2 = IncrementalWitness::from_tree(&tree);

        let mut builder = Builder::new(0);
        builder.add_sapling_spend(extsk.clone(), to.diversifier, note1, witness1).unwrap();
        builder.add_sapling_spend(extsk, to.diversifier, note2, witness2).unwrap();
        builder.add_sapling_output(ovk, to, Amount::from_u64(30000).unwrap(), None).unwrap();
        builder.add_transparent_output(
            &TransparentAddress::PublicKey([0; 20]),
            Amount::from_u64(20000).unwrap(),
        ).unwrap();

        assert_eq!(builder.build(1, MockTxProver), Err(Error::BindingSig))
    }
}
