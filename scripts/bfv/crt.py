import math
from typing import List
from .polynomial import Polynomial, PolynomialRing, get_centered_remainder


class CRTModuli:
    def __init__(self, qis: List[int]):
        """
        Initialize a CRTModuli object from the list of small moduli qis such that q = q1 * q2 * ... * qn
        """
        q = 1
        # qis should be same size single precision integers of 60 bits (or less)
        for qi in qis:
            assert isinstance(qi, int)
            assert qi.bit_length() <= 60, "qi is too large"
            q *= qi

        # qis should be coprime
        for i in range(len(qis)):
            for j in range(i + 1, len(qis)):
                assert math.gcd(qis[i], qis[j]) == 1, "qis are not pairwise coprime"

        self.qis = qis
        self.q = q


class CRTInteger:
    def __init__(self, crt_moduli: CRTModuli, xis: List[int]):
        self.crt_moduli = crt_moduli
        self.xis = xis

    @staticmethod
    def from_crt_components(crt_moduli: CRTModuli, xis: List[int]) -> "CRTInteger":
        """
        Initialize a CRTInteger object from the list of CRT components xis, which are integers in [0, qi) for each qi in qis, and the object CRTModuli, which contains the modulus q and the list of moduli qis such that q = q1 * q2 * ... * qn
        """
        assert len(xis) == len(
            crt_moduli.qis
        ), "xis and qis should have the same length"
        for i in range(len(xis)):
            assert xis[i] < crt_moduli.qis[i], "xi should lie in [0, qi)"

        return CRTInteger(crt_moduli, xis)

    @staticmethod
    def from_integer(crt_moduli: CRTModuli, x: int) -> "CRTInteger":
        """
        Initialize a CRTInteger object from the integer x, which is in [0, q), and the object CRTModuli, which contains the modulus q and the list of moduli qis such that q = q1 * q2 * ... * qn
        """
        assert x < crt_moduli.q, "x should lie in [0, q)"
        xis = []
        for qi in crt_moduli.qis:
            xis.append(x % qi)

        return CRTInteger(crt_moduli, xis)

    def recover(self) -> int:
        """
        Recover the integer x from its CRT representation. The integer x is in [0, q).
        """
        x = 0
        for i in range(len(self.crt_moduli.qis)):
            xi = self.xis[i]
            qi_star = self.crt_moduli.q // self.crt_moduli.qis[i]
            qi_tilde = pow(
                qi_star, -1, self.crt_moduli.qis[i]
            )  # inverse of qi_star mod self.crt_moduli.qis[i]
            x += xi * qi_star * qi_tilde

        return x % self.crt_moduli.q

    def recover_with_centered_remainder(self) -> int:
        """
        Recover the integer x from its CRT representation. The integer x is in (-(q-1)/2, (q-1)/2].
        """
        x = 0
        for i in range(len(self.crt_moduli.qis)):
            xi = self.xis[i]
            qi_star = 1
            for j in range(len(self.crt_moduli.qis)):
                if j != i:
                    qi_star *= self.crt_moduli.qis[j]
            qi_tilde = pow(
                qi_star, -1, self.crt_moduli.qis[i]
            )  # inverse of qi_star mod self.crt_moduli.qis[i]
            x += xi * qi_star * qi_tilde

        return get_centered_remainder(x, self.crt_moduli.q)


class CRTPolynomial:
    @staticmethod
    def from_rq_polynomial_to_rqi_polynomials(
        rq_polynomial: Polynomial, n: int, crt_moduli: CRTModuli
    ) -> List[Polynomial]:
        """
        Reduce polynomial `a`, defined in the ring R_q, to its CRT representations, defined in the rings R_qi

        Parameters:
        - rq_polynomial: polynomial in R_q
        - n: degree of the f(x) which is the denominator of the polynomial ring, must be a power of 2.
        - crt_moduli: object of type CRTModuli, which contains the modulus q and the list of moduli qis such that q = q1 * q2 * ... * qn

        Returns: list of polynomials in R_qi
        """

        rqi_polynomials = []
        for qi in crt_moduli.qis:
            rqi = PolynomialRing(n, qi)
            rq_polynomial_clone = Polynomial(rq_polynomial.coefficients)
            rq_polynomial_clone.reduce_in_ring(rqi)
            rqi_polynomials.append(rq_polynomial_clone)

        return rqi_polynomials

    @staticmethod
    def from_rqi_polynomials_to_rq_polynomial(
        rqi_polynomials: List[Polynomial], n: int, crt_moduli: CRTModuli
    ) -> Polynomial:
        """
        Recover polynomial `rqi_polynomials`, defined in the ring R_q, from its CRT representations, defined in the rings R_qi

        Parameters:
        - rqi_polynomials: list of polynomials in R_qi
        - n: degree of the f(x) which is the denominator of the polynomial ring, must be a power of 2.
        - crt_moduli: object of type CRTModuli, which contains the modulus q and the list of moduli qis such that q = q1 * q2 * ... * qn

        Returns: polynomial in R_q
        """
        assert len(rqi_polynomials) == len(crt_moduli.qis)
        rq_coefficients = []
        for i in range(len(rqi_polynomials[0].coefficients)):
            coeff_crt_components = []
            for j in range(len(rqi_polynomials)):
                coeff_crt_components.append(rqi_polynomials[j].coefficients[i])
            coeff_crt_integer = CRTInteger(crt_moduli, coeff_crt_components)
            coeff = coeff_crt_integer.recover_with_centered_remainder()
            rq_coefficients.append(coeff)

        rq = PolynomialRing(n, crt_moduli.q)
        poly = Polynomial(rq_coefficients)
        poly.reduce_in_ring(rq)
        return poly
