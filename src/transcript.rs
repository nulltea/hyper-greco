use gkr::{
    util::{
        arithmetic::{ExtensionField, PrimeField},
        RngCore, SeedableRng, StdRng,
    },
    Error,
};
use itertools::Itertools;
use plonkish_backend::util::{
    arithmetic::fe_mod_from_le_bytes,
    hash::{Hash, Keccak256, Output, Update},
    transcript::{
        FieldTranscript, FieldTranscriptRead, FieldTranscriptWrite, Transcript, TranscriptRead,
        TranscriptWrite,
    },
};
use std::{
    fmt::Debug,
    io::{self},
    iter,
};

pub type StdRngTranscript<S> = RngTranscript<S, StdRng>;

#[derive(Debug)]
pub struct RngTranscript<S, P> {
    stream: S,
    rng: P,
}

impl<P> RngTranscript<Vec<u8>, P> {
    pub fn into_proof(self) -> Vec<u8> {
        self.stream
    }
}

impl<'a> RngTranscript<&'a [u8], StdRng> {
    pub fn from_proof(proof: &'a [u8]) -> Self {
        Self::new(proof)
    }
}

impl<S> RngTranscript<S, StdRng> {
    pub fn new(stream: S) -> Self {
        Self {
            stream,
            rng: StdRng::seed_from_u64(0),
        }
    }
}

impl Default for RngTranscript<Vec<u8>, StdRng> {
    fn default() -> Self {
        Self::new(Vec::new())
    }
}

impl<F: PrimeField, E: ExtensionField<F>, S: Debug, P: Debug + RngCore>
    gkr::transcript::Transcript<F, E> for RngTranscript<S, P>
{
    fn squeeze_challenge(&mut self) -> E {
        let bases = iter::repeat_with(|| F::random(&mut self.rng))
            .take(E::DEGREE)
            .collect_vec();
        E::from_bases(&bases)
    }

    fn common_felt(&mut self, _: &F) {}
}

impl<F: PrimeField, E: ExtensionField<F>, R: Debug + io::Read, P: Debug + RngCore>
    gkr::transcript::TranscriptRead<F, E> for RngTranscript<R, P>
{
    fn read_felt(&mut self) -> Result<F, Error> {
        let mut repr = <F as PrimeField>::Repr::default();
        self.stream
            .read_exact(repr.as_mut())
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;
        repr.as_mut().reverse();
        let felt = F::from_repr_vartime(repr).ok_or_else(err_invalid_felt)?;
        Ok(felt)
    }

    fn read_felt_ext(&mut self) -> Result<E, Error> {
        let bases = iter::repeat_with(|| gkr::transcript::TranscriptRead::<F, E>::read_felt(self))
            .take(E::DEGREE)
            .try_collect::<_, Vec<_>, _>()?;
        Ok(E::from_bases(&bases))
    }
}

impl<F: PrimeField, E: ExtensionField<F>, W: Debug + io::Write, P: Debug + RngCore>
    gkr::transcript::TranscriptWrite<F, E> for RngTranscript<W, P>
{
    fn write_felt(&mut self, felt: &F) -> Result<(), Error> {
        let mut repr = felt.to_repr();
        repr.as_mut().reverse();
        self.stream
            .write_all(repr.as_ref())
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))
    }

    fn write_felt_ext(&mut self, felt: &E) -> Result<(), Error> {
        felt.as_bases()
            .iter()
            .try_for_each(|base| gkr::transcript::TranscriptWrite::<F, E>::write_felt(self, base))
    }
}

fn err_invalid_felt() -> Error {
    Error::Transcript(
        io::ErrorKind::Other,
        "Invalid field element read from stream".to_string(),
    )
}

pub type Keccak256Transcript<S> = FiatShamirTranscript<Keccak256, S>;

#[derive(Debug, Default)]
pub struct FiatShamirTranscript<H, S> {
    state: H,
    stream: S,
}

impl Keccak256Transcript<Vec<u8>> {
    pub fn into_proof(self) -> Vec<u8> {
        self.stream
    }
}

impl<'a> Keccak256Transcript<&'a [u8]> {
    pub fn from_proof(proof: &'a [u8]) -> Keccak256Transcript<io::Cursor<&'a[u8]>> {
        Keccak256Transcript::new(io::Cursor::new(proof))
    }
}

impl<S> Keccak256Transcript<S> {
    pub fn new(stream: S) -> Self {
        Self {
            stream,
            state: Keccak256::default(),
        }
    }
}

impl<F: PrimeField, E: ExtensionField<F>, S: Debug> gkr::transcript::Transcript<F, E>
    for Keccak256Transcript<S>
{
    fn squeeze_challenge(&mut self) -> E {
        E::from_bases(&<Self as FieldTranscript<F>>::squeeze_challenges(
            self,
            E::DEGREE,
        ))
    }

    fn common_felt(&mut self, _: &F) {}
}

impl<F: PrimeField, E: ExtensionField<F>, S: Debug + io::Read> gkr::transcript::TranscriptRead<F, E>
    for Keccak256Transcript<S>
{
    fn read_felt(&mut self) -> Result<F, Error> {
        let mut repr = <F as PrimeField>::Repr::default();
        self.stream
            .read_exact(repr.as_mut())
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))?;
        repr.as_mut().reverse();
        let felt = F::from_repr_vartime(repr).ok_or_else(err_invalid_felt)?;
        Ok(felt)
    }

    fn read_felt_ext(&mut self) -> Result<E, Error> {
        let bases = iter::repeat_with(|| gkr::transcript::TranscriptRead::<F, F>::read_felt(self))
            .take(E::DEGREE)
            .try_collect::<_, Vec<_>, _>()?;
        Ok(E::from_bases(&bases))
    }
}

impl<F: PrimeField, E: ExtensionField<F>, W: Debug + io::Write>
    gkr::transcript::TranscriptWrite<F, E> for Keccak256Transcript<W>
{
    fn write_felt(&mut self, felt: &F) -> Result<(), Error> {
        let mut repr = felt.to_repr();
        repr.as_mut().reverse();
        self.stream
            .write_all(repr.as_ref())
            .map_err(|err| Error::Transcript(err.kind(), err.to_string()))
    }

    fn write_felt_ext(&mut self, felt: &E) -> Result<(), Error> {
        felt.as_bases()
            .iter()
            .try_for_each(|base| gkr::transcript::TranscriptWrite::<F, E>::write_felt(self, base))
    }
}

impl<H: Hash, F: PrimeField, S> FieldTranscript<F> for FiatShamirTranscript<H, S> {
    fn squeeze_challenge(&mut self) -> F {
        let hash = self.state.finalize_fixed_reset();
        self.state.update(&hash);
        fe_mod_from_le_bytes(hash)
    }

    fn common_field_element(&mut self, fe: &F) -> Result<(), plonkish_backend::Error> {
        self.state.update_field_element(fe);
        Ok(())
    }
}

impl<H: Hash, F: PrimeField, R: io::Read> FieldTranscriptRead<F> for FiatShamirTranscript<H, R> {
    fn read_field_element(&mut self) -> Result<F, plonkish_backend::Error> {
        let mut repr = <F as PrimeField>::Repr::default();
        self.stream
            .read_exact(repr.as_mut())
            .map_err(|err| plonkish_backend::Error::Transcript(err.kind(), err.to_string()))?;
        repr.as_mut().reverse();
        let fe = F::from_repr_vartime(repr).ok_or_else(|| {
            plonkish_backend::Error::Transcript(
                io::ErrorKind::Other,
                "Invalid field element encoding in proof".to_string(),
            )
        })?;
        self.common_field_element(&fe)?;
        Ok(fe)
    }
}

impl<H: Hash, F: PrimeField, W: io::Write> FieldTranscriptWrite<F> for FiatShamirTranscript<H, W> {
    fn write_field_element(&mut self, fe: &F) -> Result<(), plonkish_backend::Error> {
        self.common_field_element(fe)?;
        let mut repr = fe.to_repr();
        repr.as_mut().reverse();
        self.stream
            .write_all(repr.as_ref())
            .map_err(|err| plonkish_backend::Error::Transcript(err.kind(), err.to_string()))
    }
}

impl<F: PrimeField, S> Transcript<Output<Keccak256>, F> for Keccak256Transcript<S> {
    fn common_commitment(
        &mut self,
        comm: &Output<Keccak256>,
    ) -> Result<(), plonkish_backend::Error> {
        self.state.update(comm);
        Ok(())
    }
}

impl<F: PrimeField, R: io::Read> TranscriptRead<Output<Keccak256>, F> for Keccak256Transcript<R> {
    fn read_commitment(&mut self) -> Result<Output<Keccak256>, plonkish_backend::Error> {
        let mut hash = Output::<Keccak256>::default();
        self.stream
            .read_exact(hash.as_mut())
            .map_err(|err| plonkish_backend::Error::Transcript(err.kind(), err.to_string()))?;
        Ok(hash)
    }
}

impl<F: PrimeField, W: io::Write> TranscriptWrite<Output<Keccak256>, F> for Keccak256Transcript<W> {
    fn write_commitment(
        &mut self,
        hash: &Output<Keccak256>,
    ) -> Result<(), plonkish_backend::Error> {
        self.stream
            .write_all(hash)
            .map_err(|err| plonkish_backend::Error::Transcript(err.kind(), err.to_string()))?;
        Ok(())
    }
}
