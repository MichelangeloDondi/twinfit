"""What the records tell apart, row M587 (the owner's message 88, read again for his 532: "techniques to identifiability
(for instance in our case about the error bars to give to the lorentzian component and Gaussian component) where
standard techniques struggle with degeneracies of the model", and his 49 and 55 on the bar of each component).

Where shared parameters are correlated beyond `NEAR`, the records hold some combinations of them much better than any
one of them, and others badly, and the bar of the fit, read at first order from the curvature of the chi square, can be
wrong for both. For each group of such parameters, three readings, the first two of mathematics and the third of the
twin:

1. its combinations, the eigenvectors of their covariance (`Ledger.combinations`): of the logarithms of the
   parameters where all are positive, so that a combination is a product of their powers and its bar a share of its
   value, else of the parameters themselves, a sum. Each is normalised so that its largest coefficient is 1;
2. for each parameter of the group, its interval refitted along it: the parameter held at values away from the fit,
   every other parameter refitted, out to where the chi square rises by 1 over its least, the interval of one
   standard deviation where the chi square is that of a normal law of its parameters (Wilks), beside its bar at first
   order;
3. on sets of the twin drawn from the fit, each combination's scatter over its bar and the share of sets whose own
   interval of 1.96 bars holds its value, and each parameter's scatter over its bar at first order and over the half
   width of its interval refitted along it: which bar is honest.

Under a right model and normal noise no reading shrinks the bar that the records allow: these say what the records
hold of a group and which bar to quote. Blind: no term is named.
"""
from __future__ import annotations

import math

import numpy as np

from .fit import fit
from .problem import Problem
from .seeds import seed_of
from .twin import replicas

NEAR = 0.9          # two parameters correlated beyond this are read as a group
TOLD_APART = 0.999          # as best.TOLD_APART: beyond it a parameter gets no bar at first order (rule 5)
HELD = 1e-10          # a combination whose information is under this share of the largest is not held by the records


def groups(corr, names) -> list:
    """The parameters linked by a correlation beyond NEAR, as groups of two or more, each in the order of `names`."""
    n, seen, out = len(names), set(), []
    for start in range(n):
        if start in seen:
            continue
        group, todo = [], [start]
        while todo:
            i = todo.pop()
            if i in group:
                continue
            group.append(i)
            todo += [j for j in range(n) if j != i and j not in group and abs(corr[i, j]) > NEAR]
        seen.update(group)
        if len(group) > 1:
            out.append([names[i] for i in sorted(group)])
    return out


def held_at(problem, j: int, value: float):
    """The problem with its shared parameter j held at `value`: the model reads it there and the fit does not move it."""
    keep = [k for k in range(len(problem.shared)) if k != j]

    def full(shared):
        return np.insert(np.asarray(shared, dtype=float), j, value)

    def model(x, shared, own, design):
        return problem.model(x, full(shared), own, design)
    centre = None if problem.centre is None else (lambda shared, own, design: problem.centre(full(shared), own, design))
    lo, hi = problem.bounds_shared
    bounds = None if problem._bounds_shared is None else ([lo[k] for k in keep], [hi[k] for k in keep])
    return Problem(problem.records, model, [problem.shared[k] for k in keep], problem.own, problem.start_shared[keep], problem._start_own, centre,
                   problem.name, bounds, problem._bounds_own)


def along(problem, result, j: int, most: float = 8.0, max_iter: int = 80) -> dict:
    """The interval of shared parameter j where the chi square of the fit, refitted with j held, rises by 1 over its least:
    each side found by doubling a step of one bar at first order, up to `most` bars or a bound, then by four halvings and
    a line between the last two points. A side not found by then is None."""
    bar, value = float(math.sqrt(result.cov[j, j])), float(result.shared[j])
    lo, hi = problem.bounds_shared[0][j], problem.bounds_shared[1][j]
    keep = [k for k in range(len(problem.shared)) if k != j]

    def rise(v):
        got = fit(held_at(problem, j, v), result.shared[keep], result.own, max_iter=max_iter)
        return float(got.chi2) - float(result.chi2) if got.converged else math.nan
    ends = {}
    for side in (-1.0, 1.0):
        inner, at_inner, step, ends[side] = 0.0, 0.0, bar, None
        while step <= most * bar and lo < value + side * step < hi:
            r = rise(value + side * step)
            if not math.isfinite(r):
                break
            if r >= 1.0:
                a, low_rise, b, high_rise = inner, at_inner, step, r
                for _ in range(4):
                    m = 0.5 * (a + b)
                    rm = rise(value + side * m)
                    if not math.isfinite(rm):
                        break
                    a, low_rise, b, high_rise = (a, low_rise, m, rm) if rm >= 1.0 else (m, rm, b, high_rise)
                ends[side] = a + (1.0 - low_rise) * (b - a) / (high_rise - low_rise) if high_rise > low_rise else b          # names that read as no element (the test of his line)
                break
            inner, at_inner, step = step, r, 2.0 * step
    low = None if ends[-1.0] is None else value - ends[-1.0]
    high = None if ends[1.0] is None else value + ends[1.0]
    half = None if low is None or high is None else 0.5 * (high - low)
    return {"value": value, "bar_at_first_order": bar, "low": low, "high": high, "half_width": half,
            "half_width_over_the_bar": None if half is None else half / bar}


def _combination(coef, values, cov, logs: bool) -> tuple:
    """The value of a combination of the parameters and its bar from their covariance."""
    if logs:
        return float(np.sum(coef * np.log(values))), float(math.sqrt(max(coef @ (cov / np.outer(values, values)) @ coef, 0.0)))
    return float(coef @ values), float(math.sqrt(max(coef @ cov @ coef, 0.0)))


def held_apart(problem, result, sets: int = 30, seed: int | None = None, max_iter: int = 80) -> dict:
    """For each group of shared parameters correlated beyond NEAR: its combinations with their bars and their reading on
    `sets` sets of the twin, and each parameter's interval refitted along it with its reading on the same sets."""
    seed = seed_of(problem, 53) if seed is None else seed          # held_apart: drawn from the records, as every Monte Carlo of the package (row M582)
    names = list(problem.shared)
    sd = np.sqrt(np.clip(np.diag(result.cov), 0.0, None))
    ok = np.isfinite(sd) & (sd > 0)
    corr = np.eye(len(names))
    if ok.any():
        corr[np.ix_(ok, ok)] = result.cov[np.ix_(ok, ok)] / np.outer(sd[ok], sd[ok])
    found = groups(corr, names)
    out = {"near": NEAR, "groups": []}
    if not found:
        return out
    fitted = []
    for one in replicas(problem, result, sets, seed=seed):
        got = fit(one, result.shared, result.own, max_iter=max_iter)
        if got.converged:
            fitted.append((np.asarray(got.shared, dtype=float), np.asarray(got.cov, dtype=float)))
    out["sets"] = len(fitted)
    for group in found:
        idx = [names.index(n) for n in group]
        values, cov = np.asarray(result.shared, dtype=float)[idx], np.asarray(result.cov, dtype=float)[np.ix_(idx, idx)]
        logs = bool(np.all(values > 0))
        told = max(abs(corr[a, b]) for a in idx for b in idx if a != b) <= TOLD_APART
        info = np.linalg.pinv(cov / np.outer(values, values) if logs else cov)          # the information, where the covariance of an exactly singular fit, a pseudo-inverse, put no variance along what the records do not hold at all (review 231, M5)
        w, vec = np.linalg.eigh(info)
        combos = []
        for i in range(len(idx) if told else 0):          # a group not told apart gives no combination: read at first order, it rests on the point where the fit stopped (review 232, M1)
            coef = vec[:, i] / vec[int(np.argmax(np.abs(vec[:, i]))), i]
            value, _ = _combination(coef, values, cov, logs)
            row = {"of": {n: float(c) for n, c in zip(group, coef)}, "kind": "a product of their powers" if logs else "a sum",
                   "value": float(math.exp(value)) if logs else value, "bar_is": "a share of its value" if logs else "in its units"}
            if w[i] <= HELD * float(w.max()):
                combos.append(dict(row, held=False, bar=None, why="the records hold nothing along it: the curvature of the fit is singular there"))
                continue
            bar = float(1.0 / math.sqrt(w[i]) / abs(vec[int(np.argmax(np.abs(vec[:, i]))), i]))
            read = [_combination(coef, s[idx], c[np.ix_(idx, idx)], logs) for s, c in fitted if not logs or np.all(s[idx] > 0)]
            hats = np.array([v for v, _ in read])
            barred = [(v, b) for v, b in read if math.isfinite(b) and b > 0]          # a set whose fit gives no bar is no interval: it was counted as a miss (review 232, N2)
            row.update(held=True, bar=bar)
            if hats.size >= 10:
                spread = float(hats.std(ddof=1))
                row["scatter_over_its_bar"] = float(spread / bar) if bar > 0 else None
                row["bias_over_its_scatter"] = float((hats.mean() - value) / spread) if spread > 0 else None
                row["holds_the_truth"] = float(np.mean([abs(v - value) <= 1.96 * b for v, b in barred])) if barred else None          # each set's interval from its own bar
                row["sets_with_a_bar"] = len(barred)
                band = 2.0 / math.sqrt(2.0 * (hats.size - 1))          # twice the relative standard error of a scatter read on that many sets
                held = 2.0 * math.sqrt(0.95 * 0.05 / max(len(barred), 1))          # twice the standard error of a share of 0.95
                unbiased = row["bias_over_its_scatter"] is not None and abs(row["bias_over_its_scatter"]) <= 2.0 / math.sqrt(hats.size)
                if unbiased and abs(row["scatter_over_its_bar"] - 1.0) <= band and row["holds_the_truth"] is not None and abs(row["holds_the_truth"] - 0.95) <= held:
                    row["bar_to_quote"], row["quote_is"] = bar, "its bar"
                elif unbiased:
                    row["bar_to_quote"], row["quote_is"] = spread, "its scatter on the twin"
                else:          # biased on the twin beyond two standard errors of its mean: no bar holds it (review 232, M1)
                    row["bar_to_quote"], row["quote_is"] = None, f"none: biased on the twin by {row['bias_over_its_scatter']:.2g} of its scatter"
            combos.append(row)
        intervals = {}
        for n, j in zip(group, idx):
            cell = along(problem, result, j, max_iter=max_iter)
            spread = np.array([s[j] for s, _ in fitted])
            if spread.size >= 10:
                cell["scatter_on_the_twin"] = float(spread.std(ddof=1))
                cell["scatter_over_the_bar"] = cell["scatter_on_the_twin"] / cell["bar_at_first_order"]
                cell["scatter_over_the_half_width"] = None if not cell["half_width"] else cell["scatter_on_the_twin"] / cell["half_width"]
            if not told:          # correlated beyond TOLD_APART: no bar at first order, rule 5, as best.md says of it (review 231, M4)
                cell.update(bar_at_first_order=None, half_width_over_the_bar=None, scatter_over_the_bar=None)
            intervals[n] = cell
        out["groups"].append({"parameters": group, "largest_correlation": float(max(abs(corr[a, b]) for a in idx for b in idx if a != b)), "told_apart": told,
                              "combinations": sorted(combos, key=lambda c: (not c["held"], c["bar"] or 0.0)), "intervals": intervals})
    return out
