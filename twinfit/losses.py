"""Chi functions: the fit of a problem under a cost of a residual that is not its square.

The weighted fit asks that the sum of the squares of the residuals be least. A chi function asks the same of another
cost of a residual, taken in a robust scale of the residuals:

    "square"          r^2, the weighted fit itself
    a number p        |r|^p, the family of powers, from 1.5 to 3. Under 2 a large residual counts less, over 2 it counts more
    "huber"           r^2 near zero and |r| beyond a corner: the canonical robust cost
    "logarithm"       log(1 + r^2): a large residual counts next to nothing
    "student"         the likelihood of errors of a Student law, its degrees of freedom and its scale read from the
                      residuals: the canonical technique at its strongest under heavy tails

    got = chi_fit(problem, base, "student")
    got.fit           the estimate: the values of a fit of the records under the last weights. Its own bars are those of
                      that weighted problem and are not the uncertainty of the estimate
    got.factor        the first order standard deviation of this estimate over the one that the weighted fit prints
    got.sd            the standard deviation of each shared parameter: the factor times the bar of the weighted fit
    got.moved         how far the last pass moved the shared parameters, in standard deviations of the weighted fit
    got.converged     whether the passes ended: a search that did not end is not an estimate

The estimate is found by least squares with weights that are set again from the residuals, until it stops moving. The
factor is the first order law of such an estimate: the scale times the root of the mean of the squared score, over the
mean of the derivative of the score. Under normal errors of the stated size every cost but the square has a factor over
1: a chi function buys nothing there. Under errors with heavy tails the factor falls under the real scatter of the
weighted fit, which is the standard deviation of the residuals. A chi function changes the scatter of an estimate. It
does not change, at first order, what a term that the model lacks does to it: for that the families of families.py are
the tool.

Windows along both axes, for each chi function. `grid` lays the residuals of a fit on the plane of position and
residual value and counts them in cells: intervals of position on each side of the centre of a curve, and bands of
value. Under a right model the value of a residual does not depend on where it is, so the table of counts has no
structure, and the reading is the canonical test of independence of a table. `crossed` makes that reading on the
residuals of the fit under each cost, each in the scale of its cost.

The factor is the value of the long run. On one set it scatters about that value, and under normal errors it prints
under 1 on a good share of the sets of a thousand residuals. A power under 1.5 is refused: the derivative of its score
is large near zero, where the fit itself pulls the residuals, and the printed factor was measured wrong by a tenth and
more at the power 1.2, in both directions, on sets of a thousand and of six thousand residuals (review 65). A power over
2 is searched with shorter steps, since a pass of the reweighted fit overshoots there, and a power over 3 is refused: its
search did not end on half of the sets with heavy tails (review 66). A search that did not end has no standard deviation.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.stats import t as student

from .fit import Fit, fit
from .problem import Record

HUBER = 1.345          # the corner of Huber's cost, in robust standard deviations: 95 per cent efficient under a normal law
SMALL = 1e-3           # under this many scales a residual is taken as this far from zero, in the weight of a power under 2
STEP = 0.1             # the step of the difference that reads the derivative of the score of a power, in scales


@dataclass
class ChiFit:
    fit: Fit
    factor: float
    moved: float
    scale: float
    passes: int
    converged: bool
    sd: np.ndarray                                    # of the shared parameters: the factor times the bar of the weighted fit
    degrees_of_freedom: float | None = None          # of the Student law, when that is the cost


def _scale(pooled: np.ndarray) -> float:
    """The robust scale of the residuals: their median absolute deviation, as a standard deviation under a normal law."""
    s = 1.4826 * float(np.median(np.abs(pooled - np.median(pooled))))
    if not s > 0.0:
        raise ValueError("the residuals have no spread: more than half of them are equal, and a robust scale cannot be read")
    return s


def _score(cost, u: np.ndarray, nu: float) -> tuple:
    """The weight of a residual in the next pass, its score and the derivative of its score, for a residual u in scales."""
    a = np.maximum(np.abs(u), SMALL)
    if cost == "square":
        return np.ones_like(u), u, np.ones_like(u)
    if cost == "huber":
        inside = np.abs(u) <= HUBER
        return np.minimum(1.0, HUBER / a), np.clip(u, -HUBER, HUBER), inside.astype(float)
    if cost == "logarithm":
        return 1.0 / (1.0 + u * u), u / (1.0 + u * u), (1.0 - u * u) / (1.0 + u * u) ** 2
    if cost == "student":
        return (nu + 1.0) / (nu + u * u), (nu + 1.0) * u / (nu + u * u), (nu + 1.0) * (nu - u * u) / (nu + u * u) ** 2
    p = float(cost)
    if p < 1.5:
        raise ValueError("a power of the residual is 1.5 or more here: under it the factor that this module prints is not the scatter of the estimate on sets of the usual size")
    if p > 3.0:
        raise ValueError("a power of the residual is 3 or less here: over it the search of this module does not end on errors with heavy tails (review 66)")

    def psi(z):
        return np.sign(z) * np.abs(z) ** (p - 1.0)
    return a ** (p - 2.0), psi(u), (psi(u + STEP) - psi(u - STEP)) / (2.0 * STEP)          # the derivative as a difference over a step: near zero the exact one is not read well from a sample


POSITIONS = (-4.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 4.0)          # the edges of the intervals of position, in widths of the curve from its centre
VALUES = (-1.15, -0.67, -0.32, 0.0, 0.32, 0.67, 1.15)                    # the edges of the bands of value, in scales: eight bands of equal count under a normal law


def grid(problem, got: Fit, scale: float = 1.0, positions=POSITIONS, values=VALUES) -> dict:
    """The residuals of a fit counted on the plane of position and value, all the records together. A cell is an interval
    of position and a band of value. The statistic is the sum over the cells of the squared count less the count that
    the margins give, over that count: the test of independence of a table, a chi square with the rows less one times
    the columns less one degrees of freedom when the value of a residual does not depend on its position. An interval
    of position in which the margins give a cell under five residuals is joined to the interval beside it, and then a
    band of value in the same way. `sparse` says that a cell under five is left all the same, with two intervals or two
    bands remaining: the law of the statistic is then approximate. The scale in which the bands are cut should be the
    spread of the residuals: with a scale twice too large the outer bands are nearly empty."""
    from scipy.stats import chi2
    from .linear import Linear
    if not got.converged:
        raise ValueError("the fit did not end: its residuals are not those of a minimum")
    lin = Linear(problem, got)
    unit = float(np.median([lin.width_of(k) for k in range(len(problem.records))]))
    xe = np.concatenate([[-np.inf], np.asarray(positions, dtype=float), [np.inf]])
    ye = np.concatenate([[-np.inf], np.asarray(values, dtype=float), [np.inf]])
    counts = np.zeros((xe.size - 1, ye.size - 1))
    for k, rec in enumerate(problem.records):
        u = (rec.x - lin.centre_of(k)) / unit
        e = problem.whitened(k, got.shared, got.own[k]) / scale
        counts += np.histogram2d(u, e, bins=(xe, ye))[0]
    def margins(table):
        return table.sum(axis=1, keepdims=True) * table.sum(axis=0, keepdims=True) / table.sum()
    expected = margins(counts)
    def rows_joined(table):          # an interval of position with too few residuals is joined to the one beside it
        table = table.copy()
        while table.shape[0] > 2 and margins(table).min() < 5.0:
            i = int(np.argmin(margins(table).min(axis=1)))
            k = i + 1 if i < table.shape[0] - 1 else i - 1
            table[k] += table[i]
            table = np.delete(table, i, axis=0)
        return table

    def bands_joined(table, scarce_only: bool):          # a band of value is joined to the one beside it: one that holds under half of an even share, or with `scarce_only` False the smallest one
        table = table.copy()
        while table.shape[1] > 2 and margins(table).min() < 5.0:
            band = table.sum(axis=0)
            j = int(np.argmin(band))
            if scarce_only and band[j] * table.shape[1] >= 0.5 * table.sum():
                break
            k = j + 1 if j < table.shape[1] - 1 else j - 1
            table[:, k] += table[:, j]
            table = np.delete(table, j, axis=1)
        return table
    first = rows_joined(counts)
    if margins(first).min() >= 5.0:
        counts = first
    else:          # joining intervals alone leaves cells under five: a band of value is nearly or wholly empty (review 67, moderate 13). The scarce bands are joined first, then the intervals, then whatever band is smallest
        counts = bands_joined(rows_joined(bands_joined(counts, True)), False)
    expected = margins(counts)
    sparse = bool(expected.min() < 5.0)
    if expected.min() <= 0.0:          # an interval or a band holds no residual at all, with two left: no table to test
        return {"chi_square": 0.0, "degrees_of_freedom": 0, "level": 1.0, "intervals_of_position": int(counts.shape[0]), "bands_of_value": int(counts.shape[1]), "sparse": True,
                "smallest_expected_count": float(expected.min()), "standardised_cells": np.zeros_like(counts).tolist()}
    stat = float(np.sum((counts - expected) ** 2 / expected))
    df = int((counts.shape[0] - 1) * (counts.shape[1] - 1))
    return {"chi_square": stat, "degrees_of_freedom": df, "level": float(chi2.sf(stat, df)) if df > 0 else 1.0, "intervals_of_position": int(counts.shape[0]), "bands_of_value": int(counts.shape[1]),
            "sparse": sparse, "smallest_expected_count": float(expected.min()), "standardised_cells": ((counts - expected) / np.sqrt(expected)).tolist()}


def crossed(problem, base: Fit | None = None, costs=("square", 1.5, "huber", "logarithm", "student")) -> dict:
    """The grid of the residuals under each cost: for each chi function, the fit under it and the reading of its residuals
    on the plane of position and value, in the scale that the cost reads. A cost whose search did not end has no reading."""
    base = fit(problem) if base is None else base
    out = {}
    for cost in costs:
        got = chi_fit(problem, base, cost)
        out[str(cost)] = dict(grid(problem, got.fit, got.scale), factor=got.factor, shared=got.fit.shared.tolist()) if got.converged and got.fit.converged else {"its_search_did_not_end": True}
    return out


def chi_fit(problem, base: Fit | None = None, cost="huber", passes: int = 30, tol: float = 1e-3) -> ChiFit:
    """The fit of the problem under a cost of the residual. `base`: the weighted fit, where the search starts. It is made
    when none is given."""
    if not (cost in ("square", "huber", "logarithm", "student") or (isinstance(cost, (int, float)) and not isinstance(cost, bool) and cost == cost)):
        raise ValueError(f"no cost is called {cost!r}: square, huber, logarithm, student, or a power from 1.5 to 3")
    if not len(problem.shared) and cost != "square":
        raise ValueError("a chi function gives another estimate of the shared parameters, and this problem has none")
    base = fit(problem) if base is None else base
    if not base.converged:
        raise ValueError("the weighted fit did not end: a chi function starts from a fit that did")
    n_rec = len(problem.records)
    slow = max(1.0, float(cost) - 1.0) if isinstance(cost, (int, float)) else 1.0          # a pass under a power p over 2 is p - 1 times the step that ends the search
    res, moved, nu, s, done, last = base, np.inf, 5.0, 1.0, 0, True
    for done in range(1, passes + 1):
        e = [problem.whitened(k, res.shared, res.own[k]) for k in range(n_rec)]
        pooled = np.concatenate(e)
        if cost == "student":
            nu, _, s = student.fit(pooled, floc=0.0)
            nu = float(min(max(nu, 2.1), 200.0))
        else:
            s = 1.0 if cost == "square" else _scale(pooled)
        recs = [Record(rec.x, rec.y, rec.sd / np.sqrt(_score(cost, ek / s, nu)[0]), rec.design, rec.name) for rec, ek in zip(problem.records, e)]
        new = fit(problem.with_records(recs), res.shared, res.own)
        last = bool(new.converged)
        if slow > 1.0:
            new = replace(new, shared=res.shared + (new.shared - res.shared) / slow, own=res.own + (new.own - res.own) / slow)
        moved = float(np.max(np.abs(new.shared - res.shared) / base.sd, initial=0.0)) * slow
        res = new
        if moved < tol or cost == "square":
            break
    u = np.concatenate([problem.whitened(k, res.shared, res.own[k]) for k in range(n_rec)]) / s
    _, psi, dpsi = _score(cost, u, nu)
    factor = s * float(np.sqrt(np.mean(psi ** 2)) / np.mean(dpsi))
    ended = bool(last and (moved < tol or cost == "square"))
    return ChiFit(res, factor, moved, float(s), done, ended, factor * base.sd if ended else np.full(base.sd.shape, np.nan), nu if cost == "student" else None)
