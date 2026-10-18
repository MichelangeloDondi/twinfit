"""Standard checks of the residuals of a fit, to stand beside the readings of the package.

Each is what an analyst runs with no method of this package: the correlation of neighbouring residuals, the runs test on
their signs, the chi square of the residuals summed in equal bins of the abscissa, and the score test against a named
term, which at first order is the likelihood ratio between the model and the model with that term added, its size free.
A table of the package sets its readings against these, the standard technique at its strongest, and not against the
chi square of the fit alone.

The residuals are those of `Linear`: each in the standard deviation of its sample, with what the own parameters of its
record absorb taken out. The levels of the first three are those of a right model with independent errors of the given
standard deviations, by the normal law or by the chi square law. The correlation and the runs test take no account of
what the fit absorbed, which is small on many samples a record. A table that uses a check gives beside it how often the
check reports with nothing missing.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import chi2, norm

from .candidates import Candidate
from .propose import _one


def _in_order(lin) -> list:
    """The residuals of each record in the order of its abscissa, with that abscissa."""
    out = []
    for k, rec in enumerate(lin.problem.records):
        x = np.asarray(rec.x, dtype=float)
        order = np.argsort(x, kind="stable")
        out.append((x[order], np.asarray(lin.r[lin.rows_of(k)], dtype=float)[order]))
    return out


def serial(lin, lag: int = 1) -> dict:
    """The correlation of residuals `lag` samples apart, within each record, over all the records. Under a right model the
    sum of the products of such pairs has mean zero and a variance of one a pair: `z` is that sum over the root of the
    number of pairs, positive when neighbours lie on the same side, and the level is on both sides. No pair joins two
    records."""
    if int(lag) != lag or lag < 1:
        raise ValueError("the lag is a whole number of samples, one or more")
    products = squares = 0.0
    pairs = 0
    for _, r in _in_order(lin):
        if r.size > lag:
            products += float(r[:-lag] @ r[lag:])
            squares += float(r @ r)
            pairs += r.size - lag
    if pairs == 0:
        raise ValueError(f"no record has more than {lag} samples: there is no pair to correlate")
    z = products / np.sqrt(pairs)
    return {"lag": int(lag), "pairs": pairs, "correlation": products / squares if squares > 0.0 else 0.0, "z": float(z), "level": float(2.0 * norm.sf(abs(z)))}


def runs(lin) -> dict:
    """The runs test of Wald and Wolfowitz on the signs of the residuals, in the order of the abscissa of each record, over
    all the records. A run is a stretch of residuals of one sign. With a of one sign and b of the other in a record of n,
    a right model gives 1 + 2ab/n runs with a variance of 2ab(2ab - n) / (n^2 (n - 1)). `z` is the sum over the records of
    the runs less that mean, over the root of the sum of the variances: negative when the residuals keep their sign for
    too long. The level is on both sides. A residual that is exactly zero has no sign and is left out."""
    found = 0
    excess = variance = 0.0
    for _, r in _in_order(lin):
        sign = r[r != 0.0] > 0.0
        n, a = sign.size, int(sign.sum())
        b = n - a
        if a == 0 or b == 0:
            continue          # a record of one sign has one run whatever the model: given its signs it says nothing
        count = 1 + int(np.count_nonzero(sign[1:] != sign[:-1]))
        found += count
        excess += count - (1.0 + 2.0 * a * b / n)
        variance += 2.0 * a * b * (2.0 * a * b - n) / (n * n * (n - 1.0))
    if variance <= 0.0:
        raise ValueError("no record has residuals of both signs: the runs test has nothing to count")
    z = excess / np.sqrt(variance)
    return {"runs": found, "runs_less_their_mean": float(excess), "z": float(z), "level": float(2.0 * norm.sf(abs(z)))}


def binned(lin, bins: int = 32) -> dict:
    """The chi square of the residuals summed in `bins` equal bins of the abscissa of each record. Each sum, over the root of
    the number of its samples, is a standard normal value under a right model, but for what the fit absorbed. The degrees
    of freedom are the bins that hold a sample less the parameters that the fit estimated, as an analyst counts them. The
    level is on the upper side: a misfit that keeps its sign over a bin raises the sum."""
    if int(bins) != bins or bins < 2:
        raise ValueError("the bins of a record are a whole number, two or more")
    statistic = 0.0
    held = 0
    for x, r in _in_order(lin):
        lo, hi = float(x[0]), float(x[-1])
        which = np.zeros(x.size, dtype=int) if hi == lo else np.minimum((bins * (x - lo) / (hi - lo)).astype(int), bins - 1)
        sums = np.bincount(which, weights=r, minlength=bins)
        n = np.bincount(which, minlength=bins)
        full = n > 0
        statistic += float(np.sum(sums[full] ** 2 / n[full]))
        held += int(full.sum())
    fitted = int(lin.n - lin.dof)
    dof = held - fitted
    if dof < 1:
        raise ValueError(f"{held} bins for {fitted} parameters of the fit: no degree of freedom is left")
    return {"bins_a_record": int(bins), "bins_that_hold_a_sample": held, "parameters_of_the_fit": fitted, "degrees_of_freedom": dof,
            "chi_square": statistic, "level": float(chi2.sf(statistic, dof))}


def against(lin, term, name: str = "a named term") -> dict:
    """The score test against one named term. `term` is a list with one array a record: what one unit of the term adds to the
    mean of that record, in the units of the values, or None for a record that it does not touch. The score is the square
    of the part of the residuals along what the fit leaves visible of the term. Under a right model it is a chi square
    with one degree of freedom, at first order, and at first order it is the fall of the chi square of the fit when the
    term is added to the model with its size free: the likelihood ratio of the two models. `size` is the size of the term
    that the residuals ask for, with its standard deviation, in the unit in which the term was written. A term written
    with a fitted value, such as a part of the fitted height, is in a unit that the missing term itself has moved: its
    size then falls short of the true one when the term is not small. The score does not depend on that unit. A term
    that the fit absorbs whole has no test."""
    entry = _one(lin, Candidate(name, "a named term", [list(term)], "any"))
    if entry.get("W") is None:
        return {"name": name, "seen": False, "why": entry["status"], "score": None, "level": None}
    return {"name": name, "seen": True, "score": float(entry["score"]), "degrees_of_freedom": int(entry["degrees_of_freedom"]), "level": float(entry["level_on_its_own"]),
            "size": float(entry["amplitude"][0]), "its_standard_deviation": float(entry["amplitude_standard_deviation"][0]),
            "visible_share": float(entry["visible_share_of_each_column"][0])}
