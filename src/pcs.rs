use digest::Digest;
use fffft::FieldFFT;
use gkr::{ff_ext::{ff::{Field, PrimeField}, ExtensionField}, poly::BoxMultilinearPoly};
use lcpc_2d::{LcCommit, LcRoot, ProverResult};
use lcpc_ligero_pc::LigeroEncodingRho;

use typenum::U1 as TLo;
use typenum::U4 as THi;

type Enc<F> = LigeroEncodingRho::<F, TLo, THi>;

pub trait PCS<F: Field, E: ExtensionField<F>> {
    fn commit_polys(poly: BoxMultilinearPoly<F, E>) -> Vec<u8>;
}

pub struct Ligero<F: Field, E: ExtensionField<F>, D: Digest> {
    _f: std::marker::PhantomData<F>,
    _e: std::marker::PhantomData<E>,
    _d: std::marker::PhantomData<D>,
}

impl<F: PrimeField + FieldFFT, E: ExtensionField<F>, D: Digest> PCS<F, E> for Ligero<F, E, D> {
    fn commit_polys(poly: BoxMultilinearPoly<F, E>) -> Vec<u8> {
        let enc = Enc::<F>::new(poly.len());
    
        let comm = LcCommit::<D, _>::commit(&poly.to_dense(), &enc).unwrap();
        comm.get_root().as_ref().to_vec()
    }
}
