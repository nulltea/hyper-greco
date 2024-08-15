import random

class PolynomialRing:
    def __init__(self, n: int, modulus: int) -> None:
        """
        Initialize a polynomial ring R_modulus = Z_modulus[x]/f(x) where f(x)=x^n+1.
        - modulus is a prime number.
        - n is a power of 2.
        """

        assert n > 0 and (n & (n - 1)) == 0, "n must be a power of 2"

        fx = [1] + [0] * (n - 1) + [1]

        self.denominator = fx
        self.modulus = modulus
        self.n = n

    def sample_polynomial(self) -> "Polynomial":
        """
        Sample polynomial from the ring
        """

        # range for random.randint
        lower_bound = - (self.modulus - 1) / 2 # inclusive
        upper_bound = (self.modulus - 1) / 2 # inclusive

        # assert that the bounds float are integers namely the decimal part is 0
        assert lower_bound % 1 == 0 and upper_bound % 1 == 0
    
        # generate n random coefficients in the range [lower_bound, upper_bound]
        coeffs = [random.randint(int(lower_bound), int(upper_bound)) for _ in range(self.n)]

        return Polynomial(coeffs)

    def __eq__(self, other) -> bool:
        if isinstance(other, PolynomialRing):
            return (
                self.denominator == other.denominator and self.modulus == other.modulus
            )
        return False


class Polynomial:
    def __init__(self, coefficients: list[int]):
        """
        Initialize a polynomial with the given coefficients starting from the highest degree coefficient.
        """
        self.coefficients = coefficients

    def reduce_coefficients_by_modulus(self, modulus: int) -> None:
        """
        Reduce the coefficients of the polynomial by the modulus of the polynomial ring.
        """
        for i in range(len(self.coefficients)):
            self.coefficients[i] = get_centered_remainder(self.coefficients[i], modulus)

    def reduce_coefficients_by_cyclo(self, cyclo: list[int]) -> None:
        """
        Reduce the coefficients by dividing it by the cyclotomic polynomial and returning the remainder.
        The cyclotomic polynomial is x^n+1.
        """
        _, remainder = poly_div(self.coefficients, cyclo)

        n = len(cyclo) - 1

        # pad the remainder with zeroes to make it len=n
        remainder = [0] * (n - len(remainder)) + remainder

        assert len(remainder) == n

        self.coefficients = remainder

    def reduce_in_ring(self, ring: PolynomialRing) -> None:
        """
        Reduce the coefficients of the polynomial by the modulus of the polynomial ring and by the denominator of the polynomial ring.
        """
        self.reduce_coefficients_by_cyclo(ring.denominator)
        self.reduce_coefficients_by_modulus(ring.modulus)

    def __add__(self, other) -> "Polynomial":
        return Polynomial(poly_add(self.coefficients, other.coefficients))

    def __mul__(self, other) -> "Polynomial":
        return Polynomial(poly_mul_naive(self.coefficients, other.coefficients))
    
    def evaluate(self, x: int) -> int:
        """
        Evaluate the polynomial at x.
        """
        result = 0
        for coeff in self.coefficients:
            result = result * x + coeff
        return result
    
    def __eq__(self, other) -> bool:
        if isinstance(other, Polynomial):
            return self.coefficients == other.coefficients
        return False
        
    def scalar_mul(self, scalar: int) -> "Polynomial":
        """
        Multiply the polynomial by a scalar.
        """
        return Polynomial([scalar * coeff for coeff in self.coefficients])
    
    def into_centered_coefficients(self, modulus: int) -> "Polynomial":
        """
        Turn the coefficients of the polynomial into centered coefficients with respect to the modulus, namely in the range [-(modulus-1)/2, (modulus+1)/2].
        """
        centered_coeffs = [get_centered_remainder(coeff, modulus) for coeff in self.coefficients]
        return Polynomial(centered_coeffs)
    
    def into_standard_form(self, modulus: int) -> "Polynomial":
        """
        Turn the coefficients of the polynomial into standard form with respect to the modulus, namely in the range [0, modulus-1].
        """
        standard_coeffs = [get_standard_form(coeff, modulus) for coeff in self.coefficients]
        return Polynomial(standard_coeffs)


def poly_div(dividend: list[int], divisor: list[int]) -> tuple[list[int], list[int]]:
    # assert that the leading coefficient of the divisor is not zero
    assert divisor[0] != 0

    # Initialize quotient and remainder
    quotient = [0] * (len(dividend) - len(divisor) + 1)
    remainder = list(dividend)

    # Main division loop
    for i in range(len(quotient)):
        coeff = (
            remainder[i] // divisor[0]
        )  # Calculate the leading coefficient of quotient
        quotient[i] = coeff

        # Subtract the current divisor*coeff from the remainder
        for j in range(len(divisor)):
            rem = remainder[i + j]
            rem -= divisor[j] * coeff
            remainder[i + j] = rem

    # Remove leading zeroes in remainder, if any
    while remainder and remainder[0] == 0:
        remainder.pop(0)

    return quotient, remainder


def poly_add(poly1: list[int], poly2: list[int]) -> list[int]:
    # Find the length of the longer polynomial
    max_length = max(len(poly1), len(poly2))
    
    # Pad the shorter polynomial with zeros at the beginning
    poly1 = [0] * (max_length - len(poly1)) + poly1
    poly2 = [0] * (max_length - len(poly2)) + poly2

    # print("poly1", [hex(n) for n in poly1[-10:]])
    # print("poly2", [hex(n) for n in poly2[-10:]])


    # Add corresponding coefficients
    result = [poly1[i] + poly2[i] for i in range(max_length)]

    # print("poly_add", [hex(n) for n in result[-10:]])
    
    return result

def get_centered_remainder(x, modulus) -> int:
    """
    Returns the centered remainder of x with respect to modulus.
    """
    r = x % modulus
    return r if r <= modulus / 2 else r - modulus

def get_standard_form(x, modulus) -> int:
    """
    Returns the standard form of x with respect to modulus.
    """
    r = x % modulus
    return r if r >= 0 else r + modulus

def poly_mul_naive(poly1: list[int], poly2: list[int]) -> list[int]:
    """
    Naive polynomial multiplication
    """
    product_len = len(poly1) + len(poly2) - 1
    product = [0] * product_len

    # Multiply each term of the first polynomial by each term of the second polynomial
    for i in range(len(poly1)):
        for j in range(len(poly2)):
            product[i + j] += poly1[i] * poly2[j]

    return product
