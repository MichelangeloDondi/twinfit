"""Collisions between Rb atoms: a Lorentzian width that grows with the density of the vapour."""
from __future__ import annotations

from .. import vapour


def width(T_celsius: float, values: dict, const: dict) -> float:
    """The width that the collisions add at this temperature, in MHz: the coefficient times the density of the vapour,
    over the density that the coefficient is written for."""
    return values["collisions"] * vapour.density(T_celsius, const) / const["units"]["density_of_the_collision_unit_per_cm3"]


def alone(T_celsius: float, values: dict, const: dict) -> float:
    """The full width at half maximum of the line that this term gives alone, in MHz: its Lorentzian width."""
    return width(T_celsius, values, const)
