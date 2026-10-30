"""The form of a candidate term, as a discrepancy of the mean, record by record, that a researcher writes for a term of
his own and that `propose` ranks on the residuals of a fit.

A candidate is a small set of columns. Each column gives, for every record, what would be added to its mean for one unit
of the term: `Candidate(name, family, columns, parity)`, where a column is a list with one array a record, and
`from_mechanism` makes one from a term of mechanisms.py. Parity says how a candidate behaves about the centre of a curve:
"even", "odd" or "any".

The package holds no list of candidates (the owner's messages 337 and 340: "fresh chat has to be completely blind: only
data driven and physics driven without any bias or hint from us"). Its report reads the residuals alone, and the term
that fills what it reads is the researcher's.
"""
from __future__ import annotations


from dataclasses import dataclass

import numpy as np



@dataclass
class Candidate:
    name: str
    family: str
    columns: list                   # each column: a list with one array a record, in the units of the values (or None)
    parity: str = "any"
    note: str = ""
    key: str = ""                   # what the candidate is whatever the fit: its family and its place in its ladder. The name can hold a width that a fit computes
    sign: int = 0                   # +1 where the size of the one column cannot be negative, as a width beside the line: after a fit that holds the line's own width a small Gaussian and a small Lorentzian width leave one direction with opposite signs, and their scores tie, so the sign alone tells which kernel is missing (the owner's message 269)

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
    """The keys of the design that are read as covariates. Left to itself it takes every key whose value is a
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


def from_mechanism(lin, mechanism, size: float, parity: str = "even") -> "Candidate":
    """A candidate from a term of mechanisms.py: what the term adds to each record at `size`, read from the fit, so that
    the shape of the term stands among the generic candidates and every family can be matched against it (the owner's
    message 262). The size is the value of the term's parameter, in its own units."""
    pr, ft = lin.problem, lin.fit
    larger = mechanism.extend(pr, ft)
    cols = [np.asarray(larger.predict(k, list(ft.shared) + [size], ft.own[k]), dtype=float) - np.asarray(pr.predict(k, ft.shared, ft.own[k]), dtype=float) for k in range(len(pr.records))]
    return Candidate(mechanism.name, "width" if mechanism.family == "kernel" else mechanism.family, [cols], parity)


