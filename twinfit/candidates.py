"""A dictionary of terms that a working model may lack, each as a discrepancy of the mean, record by record.

A candidate is a small set of columns. Each column gives, for every record, what would be added to its mean for one unit
of the term. The dictionary below is generic: it is built from the fitted model itself and from the design, and it knows
nothing of the field of the problem. A researcher or an agent adds candidates that the theory of the problem suggests with
`Candidate(name, family, columns, parity)`, where a column is a list with one array a record.

    family                          what it is
    baseline                        a curvature and a cubic term of the level under the curves, the same in every record, and each its own in each record
    periodic                        a cosine and a sine of one period, for a ladder of periods
    second component                a copy of each curve displaced by a distance, for a ladder of distances on each side
    pedestal                        a copy of each curve made several times wider about its centre
    tail                            the same on one side of the centre alone
    grows with a covariate          a shared parameter that changes along a numeric covariate of the design, and along its square
    one block apart                 the curves of one block of the design wider than the model says

Parity says how a candidate behaves about the centre of a curve: "even", "odd" or "any". A stage that admits models
symmetric about the centre alone keeps the odd candidates as alarms and proposes none of them as a term.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .identify import stretched


@dataclass
class Candidate:
    name: str
    family: str
    columns: list                   # each column: a list with one array a record, in the units of the values (or None)
    parity: str = "any"
    note: str = ""
    key: str = ""                   # what the candidate is whatever the fit: its family and its place in its ladder. The name can hold a width that a fit computes

    def __post_init__(self):
        if not self.key:
            self.key = self.name     # a name that holds nothing of a fit is its own key


ALONG, ALONG_SQUARE = " changes along ", " changes along the square of "       # how a candidate along a covariate is named: the table reads the key back from it


def _numeric_covariates(problem) -> list:
    keys = []
    for key in problem.records[0].design:
        vals = [rec.design.get(key) for rec in problem.records]
        if all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in vals) and len(set(vals)) > 2:
            keys.append(key)
    return keys


def covariates_of(problem, covariates=None) -> list:
    """The keys of the design that the dictionary reads as covariates. Left to itself it takes every key whose value is a
    number in every record and takes more than two values. The number of a repeat and the place in the order of
    acquisition are such keys, whether they are meant or not: each adds two candidates for every shared parameter, and
    they count in the level of the family. `covariates` names the keys that are meant."""
    numeric = _numeric_covariates(problem)
    if covariates is None:
        return numeric
    for key in covariates:
        if key not in numeric:
            raise ValueError(f"'{key}' is not a covariate that can be read: a covariate is a key of the design of every record whose value is a number and takes more than "
                             f"two values. Those of this problem: {', '.join(numeric) or 'none'}")
    return list(covariates)


def dictionary(lin, periods=None, displacements=(0.35, 0.5, 0.71, 1.0, 1.41, 2.0, 2.83, 4.0), stretches=(2.0, 4.0, 8.0), block: str | None = None,
               covariates=None) -> list:
    """The generic candidates for a problem at a fit. Distances are counted in the median full width of the fitted curves.

    Every key of the design whose value is a number with more than two values is read as a covariate, an index of
    repeats too: `covariates_of(problem)` names them, `covariates=[...]` keeps the ones that are meant, and the table
    of `propose` says which were read. The abscissa of every record rises: a record refuses one that does not."""
    pr, n_rec = lin.problem, len(lin.problem.records)
    centre = [lin.centre_of(k) for k in range(n_rec)]
    unit = float(np.median([lin.width_of(k) for k in range(n_rec)]))
    span = float(np.median([rec.x[-1] - rec.x[0] for rec in pr.records]))
    step = float(np.median([np.median(np.diff(rec.x)) for rec in pr.records]))
    here = lambda k, x=None: pr.predict(k, lin.fit.shared, lin.fit.own[k], x=x)          # noqa: E731
    out = []
    for power, label in ((1, "a slope of the level"), (2, "a curvature of the level"), (3, "a cubic term of the level")):          # the slope first: a model with a constant baseline in each record lacks it before any other, and a drift of the eclipses of a cold test could only be named a cubic (row M169)
        out.append(Candidate(label, "baseline", [[((rec.x - centre[k]) / span) ** power for k, rec in enumerate(pr.records)]], "even" if power == 2 else "odd"))
        if n_rec > 1:          # and its own in each record, as a background that each record carries at its own size: the bench of the toy plants one that grows with the height, which the shared one did not name
            out.append(Candidate(f"{label}, its own in each record", "baseline",
                                 [[((rec.x - centre[k]) / span) ** power if k == j else None for k, rec in enumerate(pr.records)] for j in range(n_rec)], "even" if power == 2 else "odd"))
    if periods is None:                                           # a ladder of multiples of the width: the multiple is the key, and the period in the name moves with the fit
        ladder = [(p, f"a periodic term of period {m:g} widths") for m, p in zip((0.5, 1.0, 2.0, 4.0, 8.0), unit * np.array([0.5, 1.0, 2.0, 4.0, 8.0])) if 4.0 * step <= p <= span]
    else:
        ladder = [(P, "") for P in periods]                       # periods that the caller gives are what they are in every fit
    for P, key in ladder:
        out.append(Candidate(f"a periodic term of period {P:.3g}", "periodic",
                             [[np.cos(2.0 * np.pi * (rec.x - centre[k]) / P) for k, rec in enumerate(pr.records)],
                              [np.sin(2.0 * np.pi * (rec.x - centre[k]) / P) for k, rec in enumerate(pr.records)]], "any", key=key))
    for m in displacements:
        for side, sign in (("above", 1.0), ("below", -1.0)):
            out.append(Candidate(f"a second component {m:g} widths {side} the centre", "second component",
                                 [[here(k, rec.x - sign * m * unit) for k, rec in enumerate(pr.records)]], "any"))
        out.append(Candidate(f"a second component {m:g} widths away on both sides", "second component",
                             [[here(k, rec.x - m * unit) + here(k, rec.x + m * unit) for k, rec in enumerate(pr.records)]], "even"))
    for kappa in stretches:
        wide = [here(k, centre[k] + (rec.x - centre[k]) / kappa) for k, rec in enumerate(pr.records)]
        out.append(Candidate(f"a pedestal {kappa:g} times wider than the curve", "pedestal", [wide], "even"))
        for side, sign in (("upper", 1.0), ("lower", -1.0)):
            out.append(Candidate(f"a tail on the {side} side, {kappa:g} times wider", "tail",
                                 [[w * (sign * (rec.x - centre[k]) > 0) for k, (rec, w) in enumerate(zip(pr.records, wide))]], "odd"))
    for key in covariates_of(pr, covariates):
        z = np.array([rec.design[key] for rec in pr.records], dtype=float)
        zc = (z - z.mean()) / z.std()
        for j, name in enumerate(pr.shared):
            out.append(Candidate(f"{name}{ALONG}{key}", "grows with a covariate",
                                 [[zc[k] * lin.shared_columns[k][:, j] for k in range(n_rec)]], "any"))
            sq = zc ** 2 - np.mean(zc ** 2)
            out.append(Candidate(f"{name}{ALONG_SQUARE}{key}", "grows with a covariate",
                                 [[sq[k] * lin.shared_columns[k][:, j] for k in range(n_rec)]], "any"))
    if block is not None:
        for value, members in pr.blocks(block).items():
            out.append(Candidate(f"the curves of the block {value} are wider", "one block apart",
                                 [[stretched(lin, k, 0.01) if k in members else None for k in range(n_rec)]], "even"))
    return out
