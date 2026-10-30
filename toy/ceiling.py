"""What 8 traces can tell of a detail, at first order (page 7 of a subfolder).

With the right model at its true values and no noise drawn, each trace weighted by the noise law at its true line, first
order gives the smallest uncertainty that a fit of those traces can reach for each free parameter, and how far two
parameters are tied. It is computed by the package, for each of the three joint fits. A fit of a set of traces whose
parameters these traces cannot tell apart says so. On all 8 traces the smallest uncertainty is set against the spread
of the own fit over 100 draws (`draws.py`): the two must agree.
"""
from __future__ import annotations

import numpy as np

import twinfit
from twinfit.identify import identify

from . import fits, model, traces

TIED = fits.TIED          # two parameters whose correlation is beyond this are not told apart: one number, in fits.py


def at_truth(detail: int, conditions: list, const: dict, config: dict, traces_of: int | None = None, free: list | None = None, noise=None) -> tuple:
    """The joint problem of some traces with no noise drawn, each weighted by the noise law at its true line, and the
    height, centre, offset and slope that made each trace. `traces_of` names another detail whose traces the model of
    `detail` is set against: the model that lacks a term, against traces that hold it. `free` names its free parameters
    when they are not those of config.toml. `noise` gives another standard deviation for the samples of a trace than
    the noise law, from the trace: it is used to ask what another weighting of the traces would tell. It is made by
    hand and not by `fits.problem_of`, so a value of theory may be free in it on a sweep alone: it holds the true
    values and no data, and gives no value of its own to a trace, so `fits.one_value_for_all` does not ask it (review 15)."""
    made = detail if traces_of is None else traces_of
    free = list(config["detail"][str(detail)]["free"] if free is None else free)
    entries, values = [config["parameter"][name] for name in free], model.truth(config)
    records, own = [], []
    for P, T in conditions:
        t = traces.trace(made, P, T, const, config)
        records.append(twinfit.Record(t["nu"], t["signal"], t["sd"] if noise is None else noise(t), {"P": P, "T": T}, name=traces.name_of(P, T)))
        own.append([float(model.line(made, P, T, values, const, config, nu=np.array([0.0]))[0]), 0.0, 0.0, 0.0])
    problem = twinfit.Problem(records, fits.curve(detail, const, config, free), tuple(free), fits.OWN, [values[name] for name in free], fits.start_of,
                              centre=lambda s, o, d: o[1], name=f"detail {detail} at its true values",
                              bounds_shared=([e["lower"] for e in entries], [e["upper"] for e in entries]))
    return problem, np.array(own)


def of(detail: int, const: dict, config: dict, free: list | None = None) -> dict:
    """For each of the three joint fits: the smallest uncertainty of each free parameter, the correlation of each pair,
    and whether these traces tell the parameters apart. `free` names other free parameters than those of config.toml:
    those of the fit that frees the values of theory as a test."""
    out = {}
    for sweep, conditions in traces.sweeps(const).items():
        problem, own = at_truth(detail, conditions, const, config, free=free)
        free = list(problem.shared)
        try:
            table = identify(twinfit.Linear(problem, twinfit.fit(problem, shared0=problem.start_shared, own0=own)))
            smallest = {name: float(table["parameters"][name]["standard_deviation"]) for name in free}
            ties = {pair: float(c) for pair, c in table["correlations"].items()}
            told = bool(all(np.isfinite(s) and s > 0 for s in smallest.values()) and all(abs(c) < TIED for c in ties.values()))
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            smallest, ties, told = {name: None for name in free}, {}, False
        out[sweep] = {"traces": len(conditions), "told": told, "smallest": smallest if told else {name: None for name in free}, "correlations": ties}
    return out
