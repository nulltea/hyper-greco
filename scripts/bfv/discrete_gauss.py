import numpy as np
import scipy.stats as ss


class DiscreteGaussian:
    def __init__(self, sigma: float):
        """
        Initialize the discrete Gaussian distribution.
        from https://stackoverflow.com/a/37412692

        Parameters:
        - sigma: Standard deviation
        """
        self.sigma = sigma

        # The lower and upper bounds of the range of the discrete Gaussian distribution are rounded to the nearest integer (check section 4 https://inferati.azureedge.net/docs/inferati-fhe-bfv.pdf)
        self.z_lower = np.rint(-6 * sigma)
        self.z_upper = np.rint(6 * sigma)

        # Generate the range which includes both ends
        self.x = np.arange(self.z_lower, self.z_upper + 1)

        # Define upper and lower bounds for each integer value
        self.xU, self.xL = self.x + 0.5, self.x - 0.5

        # Calculate the probability of each integer value in the range
        self.prob = ss.norm.cdf(self.xU, scale=self.sigma) - ss.norm.cdf(
            self.xL, scale=self.sigma
        )
        self.prob = (
            self.prob / self.prob.sum()
        )  # Normalize the probabilities so their sum is 1

    def sample(self, size: int) -> np.ndarray:
        """
        Sample from the Discrete Gaussian distribution

        Parameters:
        - size: Number of samples to generate

        Returns: NumPy array of generated samples
        """
        return np.random.choice(self.x, size=size, p=self.prob)
