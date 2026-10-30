"""The saturated vapour of rubidium: its density at a temperature, from the law of const.toml. It has a file of its own so
that the model and each of its terms read one and the same law. No number of the physics is written here."""
from __future__ import annotations


def density(T_celsius: float, const: dict) -> float:
    """Rb atoms per cm^3 in the saturated vapour at this temperature."""
    units, law = const["units"], const["density"]
    T = T_celsius + units["celsius_to_kelvin"]
    pressure = units["pascal_per_atmosphere"] * 10.0 ** (law["a"] - law["b"] / T)
    return pressure / (units["boltzmann_J_per_K"] * T) / units["cm3_per_m3"]
