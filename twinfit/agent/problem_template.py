"""The problem of this project, written once for the package twinfit. Fill the four places marked TO FILL.

An agent writes this file with the researcher at the intake (agent/AGENT.md, section 1). Nothing else in the project
needs to know the model: every step of the cycle takes the `problem` that this file builds.
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "{{WHERE_THE_PACKAGE_IS}}")      # written by `python -m twinfit init`: the folder that holds the package on this machine
from twinfit import Problem, Record  # noqa: E402

SHARED = ("first", "second")                   # TO FILL: the parameters that every record shares, which the analysis is after
OWN = ("area", "centre", "offset")             # TO FILL: the parameters that belong to one record
HELD = {}                                      # optional: the values that the model reads and does not fit, each a name that an ansatz of the description may hold with its uncertainty


def model(x, shared, own, design):
    """TO FILL: the mean of one record. `x` is its abscissa, `shared` and `own` are arrays in the order of SHARED and OWN,
    `design` holds the covariates of the record. Start with the small model that is understood in every part."""
    first, second = shared
    area, centre, offset = own
    width = first + second * design["covariate"]
    return area * np.exp(-0.5 * ((x - centre) / width) ** 2) / (width * np.sqrt(2.0 * np.pi)) + offset


def load() -> list:
    """TO FILL: the records. One Record a curve: its abscissa, its values, the standard deviation of each value, and its
    design: every covariate under which it was taken, a name for its block, and its place in the order of acquisition."""
    records = []
    # for each curve of the project:
    #     records.append(Record(x, y, sd, {"covariate": ..., "block": "...", "order": ...}, name="..."))
    return records


def start_own(record: Record):
    base = float(np.median(record.y))
    step = float(np.median(np.diff(record.x)))
    return [float(np.sum(record.y - base) * step), float(record.x[np.argmax(record.y)]), base]


def bounds_own(record: Record):
    """TO FILL when some records hold almost no signal: (lower, upper) for the own parameters of a record, in the order of
    OWN. An area that is not negative and a centre inside the abscissa are the usual two. Without bounds the own
    parameters of such a record wander: the fit then says how many it could not fix, and nothing more."""
    return [0.0, float(record.x[0]), -np.inf], [np.inf, float(record.x[-1]), np.inf]


def problem() -> Problem:
    return Problem(load(), model, SHARED, OWN, start_shared=[1.0, 0.0], start_own=start_own,
                   centre=lambda shared, own, design: own[1], name="the simple model", bounds_own=bounds_own)
