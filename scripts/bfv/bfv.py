from .polynomial import PolynomialRing, Polynomial, get_centered_remainder
from .discrete_gauss import DiscreteGaussian
from .crt import CRTModuli, CRTPolynomial
import numpy as np
from decimal import Decimal

class RLWE:
    def __init__(self, n: int, q: int, t: int, distribution: DiscreteGaussian):
        """
        Initialize a RLWE instance starting from the parameters

        Parameters:
        - n: degree of the f(x) which is the denominator of the polynomial ring, must be a power of 2.
        - q: modulus q of the ciphertext space
        - t: modulus t of the plaintext space
        - distribution: Error distribution (e.g. Gaussian).
        """
        # Ensure that the modulus of the plaintext space is smaller than the modulus of the polynomial ring
        if t > q:
            raise ValueError(
                "The modulus of the plaintext space must be smaller than the modulus of the polynomial ring."
            )

        # Ensure that n is a power of 2
        assert n > 0 and (n & (n - 1)) == 0, "n must be a power of 2"

        # Ensure that p and q are greater than 1
        assert q > 1, "modulus q must be > 1"
        assert t > 1, "modulus t must be > 1"

        # Ensure that q and t are coprime
        # assert np.gcd(q, t) == 1, "modulus q and t must be coprime"

        self.n = n
        self.Rq = PolynomialRing(n, q)
        self.Rt = PolynomialRing(n, t)
        self.distribution = distribution

    def SampleFromTernaryDistribution(self) -> Polynomial:
        """
        Sample a polynomial from the χ Ternary distribution.
        Namely, the coefficients are sampled uniformely from the ternary set {-1, 0, 1}. (coefficients are either of them)

        Returns: Sampled polynomial.
        """

        coefficients = np.random.choice([-1, 0, 1], size=self.n)

        coefficients_int = [int(coeff) for coeff in coefficients]

        return Polynomial(coefficients_int)

    def SampleFromErrorDistribution(self) -> Polynomial:
        """
        Sample a polynomial from the χ Error distribution.

        Returns: Sampled polynomial.
        """
        # Sample a polynomial from the Error distribution
        coefficients = self.distribution.sample(self.n)

        coefficients_int = [int(coeff) for coeff in coefficients]

        return Polynomial(coefficients_int)

class BFV:
    def __init__(self, rlwe: RLWE):
        """
        Initialize a BFV instance starting from a rlwe instance

        Parameters:
        - rlwe: RLWE instance.
        """
        self.rlwe = rlwe

    def SecretKeyGen(self) -> Polynomial:
        """
        Randomly generate a secret key.

        Returns: Generated secret key polynomial.
        """

        return self.rlwe.SampleFromTernaryDistribution()

    def PublicKeyGen(
        self, s: Polynomial, e: Polynomial, a: Polynomial
    ) -> tuple[Polynomial, Polynomial]:
        """
        Generate a public key from a given secret key.

        Parameters:
        - s: Secret key.
        - e: polynomial sampled from the distribution χ Error.
        - a: polynomial sampled from the ring Rq.

        Returns: Generated public key.
        """        

        # a * s
        mul = a * s

        # b = a*s + e.
        b = mul + e

        # pk0 is a polynomial in Rq
        pk0 = b
        pk0.reduce_in_ring(self.rlwe.Rq)

        # pk1 = -a.
        pk1 = a 
        pk1 = pk1.scalar_mul(-1)

        public_key = (pk0, pk1)

        return public_key

    def PubKeyEncrypt(
        self,
        public_key: tuple[Polynomial, Polynomial],
        m: Polynomial,
        e0: Polynomial,
        e1: Polynomial,
        u: Polynomial,
    ) -> tuple[Polynomial, Polynomial]:
        """
        Encrypt a given message m with a given public_key .

        Parameters:
        - public_key: Public key. The public key must be a tuple of polynomials living in the ring of self.rlwe.Rq.
        - m: message. This must be a polynomial in Rt.
        - e0: polynomial sampled from the distribution χ Error.
        - e1: polynomial sampled from the distribution χ Error.
        - u: polynomial sampled from the distribution χ Ternary.

        Returns:
        ciphertext: Generated ciphertext.
        """
        
        # scale the message as round(Q*m/t)
        factor = Decimal(self.rlwe.Rq.modulus) / Decimal(self.rlwe.Rt.modulus)
        scaled_message = [round(coeff * factor) for coeff in m.coefficients]

        # pk0 * u
        pk0_u = public_key[0] * u

        # scaled_message + pk0 * u + e0
        ct_0 = Polynomial(scaled_message) + pk0_u + e0

        # ct_0 will be in Rq
        ct_0.reduce_in_ring(self.rlwe.Rq)

        # pk1 * u
        pk1_u = public_key[1] * u

        # pk1 * u + e1
        ct_1 = pk1_u + e1

        # The result will be in Rq
        ct_1.reduce_in_ring(self.rlwe.Rq)

        ciphertext = (ct_0, ct_1)

        return ciphertext

    def SecretKeyEncrypt(
        self,
        secret_key: Polynomial,
        m: Polynomial,
        a: Polynomial,
        e: Polynomial,
    ) -> tuple[Polynomial, Polynomial]:
        """
        Encrypt a given message m with a given secret key .

        Parameters:
        - secret_key: Polynomial sampled from the ternary distribution χ Ternary.
        - m: message. This must be a polynomial in Rt.
        - a: polynomial sampled from the ring Rq.
        - e: polynomial sampled from the distribution χ Error.

        Returns:
        ciphertext: Generated ciphertext.
        """

        # scale the message as round(Q*m/t)
        factor = Decimal(self.rlwe.Rq.modulus) / Decimal(self.rlwe.Rt.modulus)
        scaled_message = [round(coeff * factor) for coeff in m.coefficients]

        # a * s
        mul = a * secret_key

        # b = a*s + e.
        b = mul + e

        # ct_0 = a*s + e + scaled_message
        ct_0 = b + Polynomial(scaled_message)

        # ct_0 will be in Rq
        ct_0.reduce_in_ring(self.rlwe.Rq)

        # ct_1 = -a
        ct_1 = a 
        ct_1 = ct_1.scalar_mul(-1)

        ciphertext = (ct_0, ct_1)

        return (ciphertext)

    def PubKeyEncryptConst(
        self,
        public_key: tuple[Polynomial, Polynomial],
        m: Polynomial,
        u: Polynomial,
    ):
        """
        Encrypt a given message m with a given public_key setting e0 and e1 to 0. This is used for the constant multiplication and addition.

        Parameters:
        - public_key: Public key.
        - m: message.
        - u: polynomial sampled from the distribution χ Ternary.

        Returns:
        ciphertext: Generated ciphertext.
        """

        # scale the message as round(Q*m/t)
        factor = Decimal(self.rlwe.Rq.modulus) / Decimal(self.rlwe.Rt.modulus)
        scaled_message = [round(coeff * factor) for coeff in m.coefficients]

        # pk0 * u
        pk0_u = public_key[0] * u

        # scaled_message + pk0 * u 
        ct_0 = Polynomial(scaled_message) + pk0_u 

        # ct_0 will be in Rq
        ct_0.reduce_in_ring(self.rlwe.Rq)

        # pk1 * u
        pk1_u = public_key[1] * u

        # pk1 * u 
        ct_1 = pk1_u 

        # The result will be in Rq
        ct_1.reduce_in_ring(self.rlwe.Rq)

        ciphertext = (ct_0, ct_1)

        return ciphertext

    def Decrypt(
        self,
        secret_key: Polynomial,
        ciphertext: tuple[Polynomial, Polynomial],
    ):
        """
        Decrypt a given ciphertext (encrypted using public key encryption) with a given secret key.

        Parameters:
        - secret_key: Secret key.
        - ciphertext: Ciphertext.

        Returns: Decrypted message.
        """

        ct0 = ciphertext[0]
        ct1 = ciphertext[1]
        s = secret_key
        t = self.rlwe.Rt.modulus
        q = self.rlwe.Rq.modulus

        ct1_s = ct1 * s

        # ct0 + ct1*s
        numerator = ct0 + ct1_s

        # reduce the numerator in Rq
        numerator.reduce_in_ring(self.rlwe.Rq)

        # scale the num down by t/q 
        num = [coeff * t for coeff in numerator.coefficients]

        quotient = [coeff / q for coeff in num]

        # round the coefficients of quotient to the nearest integer
        quotient = [round(coeff) for coeff in quotient]

        # trim leading zeros
        quotient = np.trim_zeros(quotient, "f")

        quotient_poly = Polynomial(quotient)

        # reduce the quotient in Rt
        quotient_poly.reduce_in_ring(self.rlwe.Rt)

        return quotient_poly

    def EvalAdd(
        self,
        ciphertext1: tuple[Polynomial, Polynomial],
        ciphertext2: tuple[Polynomial, Polynomial],
    ):
        """
        Add two ciphertexts.

        Parameters:
        - ciphertext1: First ciphertext.
        - ciphertext2: Second ciphertext.

        Returns:
        ciphertext_sum: Sum of the two ciphertexts.
        """
        # ct1_0 + ct2_0
        ct0 = ciphertext1[0] + ciphertext2[0]
        ct0.reduce_in_ring(self.rlwe.Rq)

        # ct1_1 + ct2_1
        ct1 = ciphertext1[1] + ciphertext2[1]
        ct1.reduce_in_ring(self.rlwe.Rq)

        return (ct0, ct1)
    
class BFVCrt:
    def __init__(self, crt_moduli: CRTModuli, n: int, t: int, discrete_gauss: DiscreteGaussian):
        """
        Initialize a BFV instance starting from:

        - crt_moduli: CRTModuli instance representing the CRT decomposition of the modulus q of the ciphertext space.
        Assumes that each modulus qi in the CRT decomposition is a prime number so that we can leverage NTT for fast polynomial multiplication.
        - n: degree of the f(x) which is the denominator of the polynomial ring, must be a power of 2.
        - t: modulus t of the plaintext space
        - discrete_gauss: Error distribution (e.g. Gaussian).

        Using the Chinese Remainder Theorem (CRT), an integer x ∈ Zq can be represented
        by its CRT components {xi = x mod qi ∈ Zqi }i, and operations on x in Zq can
        be implemented by applying the same operations to each CRT component xi
        in its own ring Zqi
        """
        self.crt_moduli = crt_moduli
        self.bfv_q = BFV(RLWE(n, crt_moduli.q, t, discrete_gauss))
        self.bfv_qis = []
        for qi in crt_moduli.qis:
            self.bfv_qis.append(BFV(RLWE(n, qi, t, discrete_gauss)))
            
    def SecretKeyGen(self) -> Polynomial:
        """
        Randomly generate a secret key.

        Returns: Generated secret key polynomial.
        """

        return self.bfv_q.SecretKeyGen()


    def PublicKeyGen(
        self, s: Polynomial, e: Polynomial, ais: list[Polynomial]
    ) -> tuple[Polynomial, Polynomial]:
        """
        Generate a set of public keys for each crt basis from a given secret key.

        Parameters:
        - s: Secret key.
        - e: polynomial sampled from the distribution χ Error.
        - ais: list of polynomials sampled from the ring Rqi.

        Returns: Generated public keys in their CRT representation.
        """

        public_keys = []

        for i in range(len(self.crt_moduli.qis)):

            a = ais[i]

            # a * s
            mul = a * s

            # b = a*s + e.
            b = mul + e

            # pk0 is a polynomial in Rqi
            pk0 = b
            pk0.reduce_in_ring(self.bfv_qis[i].rlwe.Rq)

            # pk1 = -a.
            pk1 = a
            pk1 = pk1.scalar_mul(-1)

            public_key = (pk0, pk1)

            public_keys.append(public_key)
        
        return public_keys
    
    def PubKeyEncrypt(
        self,
        public_keys: tuple[Polynomial, Polynomial],
        m: Polynomial,
        e0: Polynomial,
        e1: Polynomial,
        u: Polynomial,
    ) -> list[tuple[Polynomial, Polynomial]]:
        """
        Encrypt a given message m with a given list of public_keys.

        Parameters:
        - public_keys: Public keys. The public key must be a list of tuple of polynomials living in the ring of self.rlwe.Rqi.
        - m: message. This must be a polynomial in Rt.
        - e0: polynomial sampled from the distribution χ Error.
        - e1: polynomial sampled from the distribution χ Error.
        - u: polynomial sampled from the distribution χ Ternary.

        Returns:
        ciphertexts: Generated ciphertext in their CRT representation.
        """
        ciphertexts = []
        # From https://eprint.iacr.org/2018/117 remark 3.1
        partial_scaled_message = m.scalar_mul(self.bfv_q.rlwe.Rq.modulus)
        partial_scaled_message.reduce_coefficients_by_modulus(self.bfv_q.rlwe.Rt.modulus)

        for i, public_key_qi in enumerate(public_keys):

            negative_mod_inverse_t = pow(-1 * self.bfv_q.rlwe.Rt.modulus, -1, self.crt_moduli.qis[i])
            scaled_message = partial_scaled_message.scalar_mul(negative_mod_inverse_t)
            
            # pk0 * u
            pk0_u = public_key_qi[0] * u

            # scaled_message + pk0 * u + e0
            ct_0 = scaled_message + pk0_u + e0

            # ct_0 will be in Rqi
            ct_0.reduce_in_ring(self.bfv_qis[i].rlwe.Rq)

            # pk1 * u
            pk1_u = public_key_qi[1] * u
            pk1_u = Polynomial(pk1_u)

            # pk1 * u + e1
            ct_1 = pk1_u + e1

            # The result will be in Rqi
            ct_1.reduce_in_ring(self.bfv_qis[i].rlwe.Rq)

            ciphertext = (ct_0, ct_1)

            ciphertexts.append(ciphertext)

        return ciphertexts
    
    def SecretKeyEncrypt(
        self,
        s: Polynomial,
        ais: list[Polynomial],
        e: Polynomial,
        m: Polynomial,
    ) -> list[tuple[Polynomial, Polynomial]]:
        """
        Encrypt a given message m with a given secret key .

        Parameters:
        - s: Secret key.
        - ais: list of polynomials sampled from the ring Rqi.
        - e: polynomial sampled from the distribution χ Error.
        - m: message. This must be a polynomial in Rt.

        Returns:
        ciphertexts: Generated ciphertext in their CRT representation.
        """

        ciphertexts = []

        # From https://eprint.iacr.org/2018/117 remark 3.1
        partial_scaled_message = m.scalar_mul(self.bfv_q.rlwe.Rq.modulus)
        partial_scaled_message.reduce_coefficients_by_modulus(self.bfv_q.rlwe.Rt.modulus)

        for i, a in enumerate(ais):
            
            scaling_factor_den = pow(-1 * self.bfv_q.rlwe.Rt.modulus, -1, self.crt_moduli.qis[i])

            scaled_message = partial_scaled_message.scalar_mul(scaling_factor_den)

            # a * s
            mul = a * s

            # b = a*s + e.
            b = mul + e

            # ct_0 = a*s + e + scaled_message
            ct_0 = b + scaled_message

            # ct_0 will be in Rqi
            ct_0.reduce_in_ring(self.bfv_qis[i].rlwe.Rq)

            # ct_1 = -a
            ct_1 = a
            ct_1 = ct_1.scalar_mul(-1)

            ciphertext = (ct_0, ct_1)

            ciphertexts.append(ciphertext)

        return ciphertexts

    
    def Decrypt(
        self,
        s: Polynomial,
        ciphertexts: list[tuple[Polynomial, Polynomial]],
    ) -> Polynomial:
        """
        Decrypts a set of ciphertexts in their CRT representation given a secret key.
        This decryption implements the technique described in https://eprint.iacr.org/2018/117

        Parameters:
        - ciphertexts: Ciphertexts expressed in their CRT representation.
        - s: Secret key.

        Returns: Decrypted message.
        """

        matrix = []

        for i in range(len(self.crt_moduli.qis)):
            inner_product = ciphertexts[i][0] + ciphertexts[i][1] * s
            inner_product.reduce_in_ring(self.bfv_qis[i].rlwe.Rq)
            matrix.append(inner_product.coefficients)

        # assert that the matrix has k rows and n columns
        assert len(matrix) == len(self.crt_moduli.qis)

        # assert that the matrix has n columns
        assert all(len(row) == len(matrix[0]) for row in matrix)
        assert len(matrix[0]) == self.bfv_q.rlwe.n

        # recover each coefficient of m from the matrix. Procedure based on paragraph 2.3 of https://eprint.iacr.org/2018/117
        message = []
        for i in range(self.bfv_q.rlwe.n):
            message_coeff = 0
            for j in range(len(self.crt_moduli.qis)):
                x_i = matrix[j][i]
                qi_star = 1
                for k in range(len(self.crt_moduli.qis)):
                    if k != j:
                        qi_star *= self.crt_moduli.qis[k]
                qi_tilde = pow(
                    qi_star, -1, self.crt_moduli.qis[j]
                )
                scaling_factor = qi_tilde * self.bfv_q.rlwe.Rt.modulus
                scaling_factor = Decimal(scaling_factor) / Decimal(self.crt_moduli.qis[j])
                message_coeff += x_i * scaling_factor

            # round the coefficient to the nearest integer
            message_coeff = round(message_coeff)

            # message coefficient is in Rt
            message_coeff = get_centered_remainder(message_coeff, self.bfv_q.rlwe.Rt.modulus)

            message.append(message_coeff)

        return Polynomial(message)

    def DecryptDummy(
        self,
        s: Polynomial,
        ciphertexts: list[tuple[Polynomial, Polynomial]],
    ) -> Polynomial:
        """
        Decrypts a set of ciphertexts in their CRT representation given a secret key.
        This dummy approach recovers the ciphertext in the ring Rq and then decrypts it using the BFV scheme.

        Parameters:
        - ciphertexts: Ciphertexts expressed in their CRT representation.
        - s: Secret key.

        Returns: Decrypted message.
        """
        
        ct0_rqis = []
        ct1_rqis = []

        for i in range(len(self.crt_moduli.qis)):
            ct0_rqis.append(ciphertexts[i][0])
            ct1_rqis.append(ciphertexts[i][1])
        
        # Recover ciphertext in Rq
        ct0_rq = CRTPolynomial.from_rqi_polynomials_to_rq_polynomial(
            ct0_rqis, self.bfv_q.rlwe.n, self.crt_moduli
        )
        ct1_rq = CRTPolynomial.from_rqi_polynomials_to_rq_polynomial(
            ct1_rqis, self.bfv_q.rlwe.n, self.crt_moduli
        )

        ciphertext = (ct0_rq, ct1_rq)

        # Perform decryption in Rq basis
        dec = self.bfv_q.Decrypt(s, ciphertext)

        return dec
    
