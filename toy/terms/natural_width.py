"""The natural width: the Lorentzian width that the lifetime of the upper state gives, the same for every trace."""
from __future__ import annotations

import numpy as np


def width(T_celsius: float, values: dict, const: dict) -> float:
    """The natural width, in MHz."""
    return values["natural_width"]


def true_of(const: dict) -> tuple:
    """The true natural width and its uncertainty, from the lifetime of const.toml: 1 / (2 pi tau)."""
    life = const["lifetime"]
    value = const["units"]["MHz_for_one_over_a_ns"] / (2 * np.pi * life["value_ns"])
    return value, value * life["uncertainty_ns"] / life["value_ns"]


def alone(T_celsius: float, values: dict, const: dict) -> float:
    """The full width at half maximum of the line that this term gives alone, in MHz: its Lorentzian width."""
    return width(T_celsius, values, const)
