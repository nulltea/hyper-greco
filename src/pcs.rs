// use digest::Digest;
// use fffft::FieldFFT;
// use gkr::{
//     ff_ext::{
//         ff::{Field, PrimeField},
//         ExtensionField,
//     },
//     poly::BoxMultilinearPoly,
//     transcript,
// };
// use lcpc_2d::{LcCommit, LcEvalProof, LcRoot, ProverResult, SizedField};
// use lcpc_ligero_pc::LigeroEncodingRho;

// use merlin::Transcript;
// use plonkish_backend::{pcs::multilinear::{MultilinearBrakedown, MultilinearHyrax}, util::{code::BrakedownSpec6, hash::Keccak256}};
// use serde::{de::DeserializeOwned, Serialize};
// use typenum::U1 as TLo;
// use typenum::U4 as THi;

// type Enc<F> = LigeroEncodingRho<F, TLo, THi>;

// pub trait PCS<F: Field, E: ExtensionField<F>> {
//     type ProverData;
//     type Proof: Clone;

//     fn commit(poly: BoxMultilinearPoly<F, E>) -> (Vec<u8>, Self::ProverData);

//     fn open(p: Self::ProverData, z: Vec<E>, transcript: &mut Transcript) -> (Vec<F>, Self::Proof);
// }

// pub type Brakedown<F: PrimeField + Serialize + DeserializeOwned> = MultilinearBrakedown<F, Keccak256, BrakedownSpec6>;

// impl<F: PrimeField + FieldFFT, E: ExtensionField<F> + SizedField + FieldFFT, D: Digest + Clone>
//     PCS<F, E> for Brakedown<F>
// {
//     type ProverData = (Enc<F>, LcCommit<D, Enc<F>>);
//     type Proof = LcEvalProof<D, Enc<F>>;

//     fn commit(poly: BoxMultilinearPoly<F, E>) -> (Vec<u8>, Self::ProverData) {
        
//         (comm.get_root().as_ref().to_vec(), (enc, comm))
//     }

//     fn open(p: Self::ProverData, z: Vec<E>, transcript: &mut Transcript) -> (Vec<F>, Self::Proof) {
//         let (enc, comm) = p;
//         let proof = comm.prove(&z, &enc, transcript).unwrap();

//         let evals = proof.get_p_eval();

//         (evals, proof)
//     }
// }
