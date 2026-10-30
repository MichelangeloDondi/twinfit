"""The transit of the atoms through the beam: an atom is excited only while it crosses the beam, the same law for every
trace, with a rate that follows the thermal speed of the atoms."""
from __future__ import annotations

import numpy as np


def rate(T_celsius: float, values: dict, const: dict) -> float:
    """The rate b of the transit at this temperature, in MHz. `values` holds it at the reference temperature, and it
    follows the thermal speed: the square root of the temperature in kelvin."""
    kelvin = const["units"]["celsius_to_kelvin"]
    return values["transit"] * np.sqrt((T_celsius + kelvin) / (const["conditions"]["reference"][1] + kelvin))


def factor(t, T_celsius: float, values: dict, const: dict) -> np.ndarray:
    """The factor 1 / (1 + (2 pi b t)^2) of C(t). Alone, its Fourier transform is a line that falls as exp(-|nu| / b) on
    each side of its centre."""
    return 1.0 / (1.0 + (2.0 * np.pi * rate(T_celsius, values, const) * t) ** 2)


def true_of(const: dict) -> tuple:
    """The true rate at the reference temperature and its uncertainty, from the radius of the beam of const.toml:
    b = u / (2 pi w), with u the most probable speed of an atom, the square root of 2 k T / M."""
    units, beam = const["units"], const["beam"]
    kelvin = const["conditions"]["reference"][1] + units["celsius_to_kelvin"]
    speed = np.sqrt(2.0 * units["boltzmann_J_per_K"] * kelvin / (const["atom"]["mass_u"] * units["kg_per_atomic_mass_unit"]))          # m / s
    b = speed / (2.0 * np.pi * beam["waist_um"] / units["um_per_m"]) / units["Hz_per_MHz"]
    return float(b), float(b * beam["uncertainty_um"] / beam["waist_um"])


def from_fit(b: float, its_bar: float, const: dict) -> dict:
    """The fixed values of const.toml that give a fitted rate at the reference temperature: the radius of the beam,
    w = u / (2 pi b), with the uncertainty that the one of the rate gives it, each at the decimals of that uncertainty.
    The other way of `true_of`: when the joint fit of the recorded traces prefers another rate, it is the radius of the
    beam that moves (the owner's message 126)."""
    units = const["units"]
    kelvin = const["conditions"]["reference"][1] + units["celsius_to_kelvin"]
    speed = np.sqrt(2.0 * units["boltzmann_J_per_K"] * kelvin / (const["atom"]["mass_u"] * units["kg_per_atomic_mass_unit"]))          # m / s
    waist = speed / (2.0 * np.pi * b * units["Hz_per_MHz"]) * units["um_per_m"]
    bar = waist * its_bar / b
    decimals = max(0, 1 - int(np.floor(np.log10(bar))))
    return {("beam", "waist_um"): round(float(waist), decimals), ("beam", "uncertainty_um"): round(float(bar), decimals)}


def alone(T_celsius: float, values: dict, const: dict) -> float:
    """The full width at half maximum of the line that this term gives alone, in MHz: that of exp(-|nu| / b), 2 b ln 2."""
    return 2.0 * float(np.log(2.0)) * rate(T_celsius, values, const)
