use crate::{keys, primitives, note_encryption, prover, JUBJUB};
use transaction::{components, builder};
use pairing::bls12_381;
use rand;
use jubjub;
use jubjub::fs;

use ff::Field;

pub struct SaplingOutput {
    ovk: keys::OutgoingViewingKey,
    to: primitives::PaymentAddress<bls12_381::Bls12>,
    note: primitives::Note<bls12_381::Bls12>,
    memo: note_encryption::Memo,
}

impl SaplingOutput {
    pub fn new<R: rand::RngCore + rand::CryptoRng>(
        rng: &mut R,
        ovk: keys::OutgoingViewingKey,
        to: primitives::PaymentAddress<bls12_381::Bls12>,
        value: components::Amount,
        memo: Option<note_encryption::Memo>,
    ) -> builder::Result<Self> {

        let g_d = to.g_d(&JUBJUB).ok_or(builder::Error::InvalidAddress)?;

        if value.is_negative() { return Err(builder::Error::InvalidAmount) }

        let rcm = fs::Fs::random(rng);

        let note = primitives::Note {
            g_d,
            pk_d: to.pk_d.clone(),
            value: value.into(),
            r: rcm,
        };

        Ok(SaplingOutput {
            ovk,
            to,
            note,
            memo: memo.unwrap_or_default(),
        })
    }

    pub fn build<P, R>(self, prover: &mut P, rng: &mut R) -> Option<builder::OutputDescription> 
    where
        P: prover::TxProver, 
        R: rand::RngCore + rand::CryptoRng
    {
        let encryptor = note_encryption::SaplingNoteEncryption::new(
            self.ovk,
            self.note.clone(),
            self.to.clone(),
            self.memo,
            rng,
        );

        let (zkproof, cv) = prover.output_proof((
            encryptor.esk().clone(),
            self.to,
            self.note.r,
            self.note.value,
        ))?.into_tuple();

        let cmu = self.note.cm(&jubjub::JubjubBls12::new());

        let enc_ciphertext = encryptor.encrypt_note_plaintext();
        let out_ciphertext = encryptor.encrypt_outgoing_plaintext(&cv, &cmu);

        let ephemeral_key = encryptor.epk().clone().into();

        Some(builder::OutputDescription {
            cv,
            cmu,
            ephemeral_key,
            enc_ciphertext,
            out_ciphertext,
            zkproof,
        })
    }
}
