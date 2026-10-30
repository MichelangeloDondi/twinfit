"""The laser noise: a Gaussian spread of the frequency, of standard deviation sigma, the same for every trace."""
from __future__ import annotations

import numpy as np


def factor(t, T_celsius: float, values: dict, const: dict) -> np.ndarray:
    """The factor exp(-2 (pi sigma t)^2) of C(t): its Fourier transform is a Gaussian of standard deviation sigma."""
    return np.exp(-2.0 * (np.pi * values["laser_noise"] * t) ** 2)


def alone(T_celsius: float, values: dict, const: dict) -> float:
    """The full width at half maximum of the line that this term gives alone, in MHz: that of a Gaussian of standard
    deviation sigma, 2 sqrt(2 ln 2) sigma."""
    return 2.0 * float(np.sqrt(2.0 * np.log(2.0))) * values["laser_noise"]
