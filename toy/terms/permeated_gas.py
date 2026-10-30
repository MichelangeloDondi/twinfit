"""The permeated gas: one more Lorentzian width, the same for every trace."""
from __future__ import annotations


def width(T_celsius: float, values: dict, const: dict) -> float:
    """The width that the gas adds, in MHz."""
    return values["gas_width"]


def alone(T_celsius: float, values: dict, const: dict) -> float:
    """The full width at half maximum of the line that this term gives alone, in MHz: its Lorentzian width."""
    return width(T_celsius, values, const)
