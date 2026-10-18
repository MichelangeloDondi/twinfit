"""The weights of a fit when the researcher holds a law of the noise and not the noise of each sample.

A record of the package carries the standard deviation of each of its values. A researcher often holds less: a law,
measured on the instrument, that gives the variance of a sample from the height of the signal there. The height is
known only once the curve is fitted, so the weights and the fit are found in turn:

    flat(x, y)                         a first standard deviation for a record, one number, read from its neighbouring samples
    problem, result, moved = settle(problem, variance)
                                       fit, set every standard deviation from the law at the fitted mean, fit again, and once more
    check(problem, result)             the residuals against the weights, by height of the fitted mean: is the law right
    allowance(problem, result, variants)   how far one bar of each part of the law moves the mean of the chi square
    fitted_law(problem, result, level) the two numbers of a law  variance = a + b * level  that the residuals ask for

`variance(mean, record, own)` is the researcher's law: it receives the fitted mean of a record, the record and the
record's own fitted parameters, and returns the variance of each sample. It is theirs because only they know what in
the mean is signal and what is an offset.

Why the weights are set between fits and not inside one. A sum of squares whose weights move with the parameters is
not a likelihood: a fit can lower it by raising the fitted level, and the estimates are then biased. With the weights
held during each fit, every fit is an ordinary weighted least squares, and its chi square has the law that the
criteria and the thresholds of this package assume, to the extent that the fit before it is close. Three passes are the
default: the first has weights that are not those of the law, the second moves from them, and the third says by how
little it moved that the weights have settled. The weights still come from a fit of these records: on a sample the error of that is of second order, and
`check` is how a wrong law, or a first fit that was far, shows.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from .fit import fit
from .problem import Record


def flat(x, y) -> float:
    """One standard deviation for a record, from the differences of neighbouring samples: for noise that is independent
    from sample to sample, on a curve that changes little over one step. A first weight, to be replaced by the law."""
    d = np.diff(np.asarray(y, dtype=float))
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2.0))


def with_law(problem, result, variance: Callable):
    """The same problem with the standard deviation of every sample taken from the law at the fitted mean."""
    recs = []
    for k, rec in enumerate(problem.records):
        var = np.asarray(variance(problem.predict(k, result.shared, result.own[k]), rec, result.own[k]), dtype=float)
        if var.shape != rec.x.shape or not np.all(np.isfinite(var)) or not np.all(var > 0):
            raise ValueError(f"record {rec.name!r}: the law must give one positive variance for each sample")
        recs.append(Record(rec.x, rec.y, np.sqrt(var), rec.design, rec.name))
    return problem.with_records(recs)


def settle(problem, variance: Callable, passes: int = 3, max_iter: int = 80):
    """Fit, set the weights from the law at the fit, fit again, `passes` times in all. Returns the last problem and its
    fit, and how far the shared estimates moved in the last pass, in bars: under a tenth, the weights have settled. With
    two passes that number is the move from the first weights to those of the law, which is not small."""
    if passes < 2:
        raise ValueError("at least two passes: the first fit has weights that are not those of the law")
    result = fit(problem, max_iter=max_iter)
    moved = None
    for _ in range(passes - 1):
        problem = with_law(problem, result, variance)
        again = fit(problem, result.shared, result.own, max_iter=max_iter)
        free = np.isfinite(again.sd)                              # a shared parameter held at a bound has no bar, and does not move
        moved = float(np.max(np.abs(again.shared - result.shared)[free] / again.sd[free])) if free.any() else 0.0
        result = again
    return problem, result, moved


def check(problem, result, level: Callable | None = None, bands: int = 6) -> list:
    """The mean of the squared whitened residuals in bands of the height of the fitted mean. Under a right model with a
    right law each is 1, within the standard error given. A trend with the height says that the law is wrong in the
    part that grows with the signal, and one value for all bands says that its floor is. `level(mean, record, own)`
    gives the height that the bands are cut on: by default the fitted mean less its lowest value in the record."""
    h, r2 = [], []
    for k, rec in enumerate(problem.records):
        mean = problem.predict(k, result.shared, result.own[k])
        h.append(np.asarray(level(mean, rec, result.own[k]), dtype=float) if level else mean - mean.min())
        r2.append(problem.whitened(k, result.shared, result.own[k]) ** 2)
    h, r2 = np.concatenate(h), np.concatenate(r2)
    edges = np.quantile(h, np.linspace(0.0, 1.0, bands + 1))
    rows = []
    for lo, hi in zip(edges, edges[1:]):
        inside = (h >= lo) & ((h < hi) if hi < edges[-1] else (h <= hi))
        if inside.sum() > 1:
            rows.append({"height_from": float(lo), "to": float(hi), "samples": int(inside.sum()), "mean_squared_residual": float(r2[inside].mean()),
                         "standard_error": float(np.sqrt(2.0 / inside.sum()))})
    return rows


def fitted_law(problem, result, level: Callable) -> dict:
    """The a and b of  variance = a + b * level  that the squared residuals ask for, by weighted least squares, with
    their bars. `level(mean, record, own)` is the height of the signal at each sample. A law fitted on residuals is
    low by the share of the variance that the fit took: about the number of parameters over the number of samples."""
    h, e2, w = [], [], []
    for k, rec in enumerate(problem.records):
        mean = problem.predict(k, result.shared, result.own[k])
        h.append(np.clip(np.asarray(level(mean, rec, result.own[k]), dtype=float), 0.0, None))
        e2.append((rec.y - mean) ** 2)
        w.append(1.0 / (2.0 * rec.sd ** 4))                          # the variance of a squared normal error is twice its variance squared
    h, e2, w = np.concatenate(h), np.concatenate(e2), np.concatenate(w)
    X = np.column_stack([np.ones_like(h), h])
    cov = np.linalg.inv(X.T @ (X * w[:, None]))
    a, b = cov @ (X.T @ (w * e2))
    return {"a": float(a), "b": float(b), "bar_of_a": float(np.sqrt(cov[0, 0])), "bar_of_b": float(np.sqrt(cov[1, 1])), "samples": int(h.size)}


def allowance(problem, result, variants: dict) -> dict:
    """How far the mean of the chi square moves when the law of the noise moves by one bar of one of its parts.

    The weights of a fit are those of the researcher's law. If the noise follows another law, the mean of the chi
    square is no longer its degrees of freedom: each sample adds the variance that it has over the variance that it
    was given. `variants` maps the name of a part of the law to the law with that part raised by one bar, each a
    function like the `variance` of `settle`. What comes back, for each part, is the move of that mean in standard
    deviations of the chi square, and their sum in quadrature for parts whose bars are independent.

    What also comes back is where the chi square of the fit sits in standard deviations of its own scatter and those
    bars together, sqrt(1 + sum^2). No verdict is given: a chi square near its mean on that scale is what a right model
    gives, and it is also what a model with a missing term gives when the bars of the law are wide, since a law known
    to a few per cent can move the mean of the chi square of many samples by more than a missing term does. The share
    of the degrees of freedom of a sample is taken as the same for all."""
    n = sum(rec.x.size for rec in problem.records)
    sd = float(np.sqrt(2.0 * result.dof))
    parts = {}
    for name, law in variants.items():
        move = 0.0
        for k, rec in enumerate(problem.records):
            var = np.asarray(law(problem.predict(k, result.shared, result.own[k]), rec, result.own[k]), dtype=float)
            if var.shape != rec.x.shape or not np.all(np.isfinite(var)) or not np.all(var > 0):
                raise ValueError(f"record {rec.name!r}: the law {name!r} must give one positive variance for each sample")
            move += float(np.sum(var / rec.sd ** 2 - 1.0))
        move *= result.dof / n
        parts[name] = {"move_of_the_mean": move, "in_standard_deviations_of_the_chi_square": move / sd}
    together = float(np.sqrt(sum(v["in_standard_deviations_of_the_chi_square"] ** 2 for v in parts.values())))
    here = float((result.chi2 - result.dof) / sd)
    return {"parts": parts, "together_in_standard_deviations": together, "the_chi_square_sits_at": here,
            "in_standard_deviations_of_its_scatter_and_the_bars_together": float(here / np.sqrt(1.0 + together ** 2))}

