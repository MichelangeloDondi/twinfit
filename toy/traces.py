"""The 8 synthetic traces of a detail of the toy model, and their files.

A trace is the line of the model at one power and one temperature plus noise. The noise of a point has the standard
deviation that was measured on his recorded traces at that condition: variance = a^2 + b h + c h^2 + d s^2, with h the
signal and s its slope. The random numbers of a trace come from the seed of const.toml, the number of the detail, the power
and the temperature, so that two builds give the same bytes.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import model


def conditions(const: dict) -> list:
    """The 8 conditions as (power in mW, temperature in C): the power sweep, then the temperatures that it does not hold."""
    power = [tuple(c) for c in const["conditions"]["power_sweep"]]
    return power + [tuple(c) for c in const["conditions"]["temperature_sweep"] if tuple(c) not in power]


def sweeps(const: dict) -> dict:
    """The three sets of traces that a joint fit is made over: the power sweep, the temperature sweep, and all 8."""
    return {"power": [tuple(c) for c in const["conditions"]["power_sweep"]],
            "temperature": [tuple(c) for c in const["conditions"]["temperature_sweep"]], "all": conditions(const)}


def noise_sd(P_mW: float, T_celsius: float, signal, nu, const: dict) -> np.ndarray:
    """The standard deviation of the noise at each point of a trace, in V, from the law of its condition."""
    law = next(r for r in const["noise"]["trace"] if (r["power_mW"], r["temperature_C"]) == (P_mW, T_celsius))
    slope = np.gradient(signal, nu)
    return np.sqrt(law["a_V"] ** 2 + law["b_V"] * signal + law["c"] * signal ** 2 + const["noise"]["d_MHz2"] * slope ** 2)


def trace(detail: int, P_mW: float, T_celsius: float, const: dict, config: dict, draw: int = 0) -> dict:
    """One synthetic trace of a detail: its frequencies, its signal with no noise, the standard deviation of its noise, and
    the trace itself. The same at every run. `draw` 0 is the trace of the folder: another number gives another draw of
    the noise, for the check of a fit on many draws."""
    nu = model.axis(const)
    signal = model.line(detail, P_mW, T_celsius, model.truth(config), const, config, nu)
    sd = noise_sd(P_mW, T_celsius, signal, nu, const)
    seed = list(const["seed"]["value"]) + [int(detail), int(P_mW), int(T_celsius)] + ([int(draw)] if draw else [])
    return {"P": P_mW, "T": T_celsius, "nu": nu, "signal": signal, "sd": sd, "y": signal + sd * np.random.default_rng(seed).standard_normal(nu.size)}


def name_of(P_mW: float, T_celsius: float) -> str:
    return f"trace_P{int(P_mW):03d}mW_T{int(T_celsius):03d}C.csv"


def write(detail: int, folder: Path, const: dict, config: dict) -> list:
    """The files of a detail in `folder`: `traces/` with one CSV file a trace, the frequency in MHz and the signal in V, and
    `truth/truth.json` with the values that made them, which an analyst is not given."""
    folder = Path(folder)
    (folder / "traces").mkdir(parents=True, exist_ok=True)
    (folder / "truth").mkdir(parents=True, exist_ok=True)
    written, rows = [], []
    for P, T in conditions(const):
        t = trace(detail, P, T, const, config)
        path = folder / "traces" / name_of(P, T)
        path.write_text("frequency_MHz,signal_V\n" + "".join(f"{x:.6f},{y:.9e}\n" for x, y in zip(t["nu"], t["y"])))
        written.append(path)
        rows.append({"file": path.name, "power_mW": P, "temperature_C": T, "height_V": float(t["signal"].max()),
                     "lorentzian_width_MHz": float(model.width(detail, T, model.truth(config), const)), "density_per_cm3": float(model.density(T, const))})
    held = {name: value for name, value in model.truth(config).items() if config["parameter"][name]["enters"] <= detail}
    (folder / "truth" / "truth.json").write_text(json.dumps({"detail": detail, "parameters": held, "traces": rows}, indent=1) + "\n")
    return written + [folder / "truth" / "truth.json"]
