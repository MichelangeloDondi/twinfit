"""The toy model of his line, one factor a detail (the owner's messages 111 to 114).

The line is the Fourier transform of a function C(t) of the time t. A detail of the model is one more factor of C(t), or one
more term of a width inside a factor. `line(detail, P, T, values, const, config)` is the signal of one trace with the first
`detail` details and no noise.

    detail 1   the natural width: a Lorentzian of that width, the factor exp(-pi Gamma |t|) of C(t)
    detail 2   collisions between Rb atoms: Gamma grows with the density of the vapour
    detail 3   the permeated gas: one more constant width in Gamma
    detail 4   the laser noise: a Gaussian factor of C(t)

Each term is a file of `toy/terms/`, and `const.toml` says at which detail it enters. This program names no term: a new
term is a new file there, and nothing here changes.

His messages 111 to 117 set this work: OWNER_WORDS.md.

No number of the physics is written here. The fixed values are in const.toml. The parameters that a fit may free are in
config.toml with their true values, and `values` is a dictionary of those parameters by name: the true ones to make a
trace, the fitted ones inside a fit. A test refuses a number in this file other than 0, 1, 2 and the 10 of the vapour law.
"""
from __future__ import annotations

import importlib
import json
import tomllib
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from . import vapour

HERE = Path(__file__).resolve().parent
density = vapour.density          # the density of the vapour, under the name that the pages and the checks ask the model for


def load(folder: Path = HERE, root: Path | None = None) -> tuple:
    """The two files: the fixed values, and the free parameters with what the builder makes. The noise law of the peak
    that config.toml names is read with them, from the file of the laws of every peak, and put in the fixed values. A
    parameter that follows the fixed values gets its true value and its uncertainty here, from the term that it names."""
    const, config = tomllib.loads((Path(folder) / "const.toml").read_text()), tomllib.loads((Path(folder) / "config.toml").read_text())
    laws = json.loads(((HERE.parent if root is None else Path(root)) / const["noise"]["file"]).read_text())
    mine = [law for law in laws["laws"] if law["peak_nm"] == config["recorded"]["peak"]]
    if not mine:
        raise ValueError(f"no noise law for the peak {config['recorded']['peak']} nm: the file holds {sorted({law['peak_nm'] for law in laws['laws']})}")
    const["noise"] = dict(const["noise"], peak=config["recorded"]["peak"], d_MHz2=laws["d_MHz2"], trace=mine)
    for entry in config["parameter"].values():          # a value that follows the fixed ones is typed nowhere: its term computes it from const.toml
        if "follows" in entry:
            entry["true"], entry["uncertainty"] = importlib.import_module(f"{__package__}.terms.{entry['follows']}").true_of(const)
    return const, config


def truth(config: dict) -> dict:
    """The true value of each parameter of config.toml, by name: what the synthetic traces are made with."""
    return {name: float(entry["true"]) for name, entry in config["parameter"].items()}


def numbers_of(tree):
    """The numbers of a table of one of the two files, without its words: a source, a name and a sentence are no value."""
    if isinstance(tree, dict):
        got = {k: numbers_of(v) for k, v in tree.items()}
        got = {k: v for k, v in got.items() if v is not None}
        return got or None
    if isinstance(tree, (list, tuple)):
        got = [numbers_of(v) for v in tree]
        return got if got and all(v is not None for v in got) else None
    return tree if isinstance(tree, (bool, int, float)) else None


def values_of(detail: int, const: dict, config: dict) -> dict:
    """The values that the subfolder of a detail is made from (the owner's message 119: "a new combination of values of
    the parameters"): every number of const.toml with the noise law of the peak, the peak and the repetition of the
    recorded traces, the numbers of each parameter that has entered the model by this detail, and the parameters that
    each detail up to this one frees. Not a sentence, a name or a source, and nothing of the code: a build from the same
    values after a repair takes the place of the one before it."""
    fixed = dict(const, details={name: d for name, d in const["details"].items() if d <= detail})          # a term that enters later is no value of this subfolder
    got = {"fixed": numbers_of(fixed), "recorded": dict(config["recorded"]),
           "free": {str(d): list(config["detail"][str(d)]["free"]) for d in sorted(int(k) for k in config["detail"]) if d <= detail},
           "parameters": {name: numbers_of(e) for name, e in config["parameter"].items() if e["enters"] <= detail}}
    return json.loads(json.dumps(got, sort_keys=True))          # as the record of a build reads it back


def axis(const: dict) -> np.ndarray:
    """The frequencies of a trace, in MHz from the centre of the line."""
    half = const["axis"]["half_window_MHz"]
    return np.linspace(-half, half, const["axis"]["points"])


def terms(detail: int, const: dict) -> list:
    """The terms of the model with the first `detail` details, in the order in which they enter: for each the program of
    `toy/terms/` that has its name in const.toml."""
    return [importlib.import_module(f"{__package__}.terms.{name}") for name, d in sorted(const["details"].items(), key=lambda item: item[1]) if d <= detail]


def width(detail: int, T_celsius: float, values: dict, const: dict) -> float:
    """The Lorentzian width of the line with the first `detail` details: its full width at half maximum, in MHz. The sum
    of the widths of the terms that give one."""
    total = 0.0
    for term in terms(detail, const):
        if hasattr(term, "width"):
            total = total + term.width(T_celsius, values, const)
    return total


def fwhm(detail: int, T_celsius: float, values: dict, const: dict) -> float:
    """The full width at half maximum of the line with the first `detail` details, in MHz, read on the line itself: for
    a line that is no Lorentzian no formula gives it."""
    from scipy.optimize import brentq
    half = float(shape(detail, 0.0, T_celsius, values, const)) / 2.0
    return 2.0 * float(brentq(lambda x: float(shape(detail, x, T_celsius, values, const)) - half, 0.0, float(np.max(np.abs(axis(const))))))          # the half height is met inside the frequencies of a trace


def correlation(detail: int, T_celsius: float, values: dict, const: dict) -> np.ndarray:
    """C(t) on the time axis of const.toml, from t = 0: the product of the factors of the first `detail` details."""
    t = const["time"]["step_us"] * np.arange(const["time"]["points"])
    c = np.exp(-np.pi * width(detail, T_celsius, values, const) * t)
    for term in terms(detail, const):
        if hasattr(term, "factor"):
            c = c * term.factor(t, T_celsius, values, const)
    return c


def shape(detail: int, nu, T_celsius: float, values: dict, const: dict) -> np.ndarray:
    """The line of area 1 at the frequencies nu, in MHz from its centre: the Fourier transform of C(t), which is even in t.
    The transform is taken once, on the frequencies that the time axis gives, and read between them by a cubic spline."""
    c = correlation(detail, T_celsius, values, const)
    step = const["time"]["step_us"]
    on_grid = step * (2 * np.fft.rfft(c).real - c[0])
    grid = np.fft.rfftfreq(c.size, step)
    nu = np.abs(np.asarray(nu, dtype=float))
    near = grid <= nu.max() + (grid[2] - grid[0])
    return CubicSpline(grid[near], on_grid[near])(nu)


def area(detail: int, P_mW: float, T_celsius: float, const: dict, config: dict) -> float:
    """The area of the line of a trace, in V MHz: it goes as the square of the power and as the density, and it is set so
    that the reference trace, made with the true values, has the height that was recorded."""
    P_ref, T_ref = const["conditions"]["reference"]
    top = float(shape(detail, 0.0, T_ref, truth(config), const))
    return const["reference_height"]["value"] / top * (P_mW / P_ref) ** 2 * density(T_celsius, const) / density(T_ref, const)


def line(detail: int, P_mW: float, T_celsius: float, values: dict, const: dict, config: dict, nu=None) -> np.ndarray:
    """The signal of one trace in V, with no noise, at the frequencies nu or on the axis of const.toml."""
    nu = axis(const) if nu is None else nu
    return area(detail, P_mW, T_celsius, const, config) * shape(detail, nu, T_celsius, values, const)
