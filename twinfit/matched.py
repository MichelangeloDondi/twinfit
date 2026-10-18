"""Statistics that are formed on the noisy values of a record, and matched to synthetic records of the model.

A statistic of families.py is a sum of the values of a record under a weight that the fit places. A statistic of this
module is formed on the values themselves:

    the abscissa against the value   where a straight line through the samples near a level of the height meets that
                                     level, for each level and each side
    bands of value                   how many samples have their value in a band, and where they lie, for each side
    the logarithm of the record      sums of the logarithm of the values under kernels of the abscissa

The noise then enters such a statistic twice, in a value and in where the value is taken, and the mean of the statistic
is not the statistic of the mean curve: it has a bias from the noise. So it is not set against the curve of the model.
It is set against its own mean over synthetic records of the model, drawn with the standard deviations of the record,
which is the canonical simulated method of moments, or indirect inference. The estimate is the value of the parameters
at which the statistics of the records agree with the mean of the statistics of the synthetic ones, in the covariance of
the synthetic ones at the fit. The same draws are used at every value of the parameters, so that the mean is a smooth
function of them, and each statistic is made smooth for the same reason: a ratio is softened where its denominator is
no larger than its noise, and a sample near the edge of a band belongs to both bands in part.

    lin = Linear(problem, fit)
    got = match(lin, "the abscissa against the value")
    got.shared, got.own, got.converged

What places a statistic, the centre of a curve, its lowest value, its height and its unit of width, is taken from the
fit and held. The standard deviations are those that the records state. A fit that did not end is refused.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

LEVELS = tuple(float(v) for v in np.arange(1, 10) / 10.0)          # the levels of the abscissa against the value, as parts of the height
REACH = 0.2                                                         # the half width of the interval about a crossing on which the straight line is fitted, in widths
CUTS = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)        # the edges of the bands of value, as parts of the height
EDGE = 0.02                                                         # the width of the edge of a band, as a part of the height
FLOOR = 0.004                                                       # the least value, as a part of the height, of which the logarithm is taken
FAR = 0.75                                                          # beyond this part of the reach of a record on each side: two means that hold its lowest value


def frame(lin, k: int) -> dict:
    """What places the statistics of record k, from the fit and held from then on."""
    pr, ft = lin.problem, lin.fit
    rec = pr.records[k]
    unit = float(np.median([lin.width_of(i) for i in range(len(pr.records))]))
    u = (rec.x - lin.centre_of(k)) / unit
    f = pr.predict(k, ft.shared, ft.own[k])
    low, height = float(f.min()), float(f.max() - f.min())
    reach = max(REACH, 2.5 * float(np.median(np.abs(np.diff(u)))))
    lines = []
    for side in (1.0, -1.0):
        order = np.argsort(side * u)
        d, w = (side * u)[order], ((f - low) / height)[order]
        d, w = d[d >= 0.0], w[d >= 0.0]
        for level in LEVELS:
            under = np.flatnonzero(w < level)
            if not under.size or under[0] == 0:
                continue
            j = int(under[0])
            at = d[j - 1] + (w[j - 1] - level) * (d[j] - d[j - 1]) / (w[j - 1] - w[j])
            idx = np.flatnonzero(np.abs(u - side * at) < reach)
            if idx.size < 3:
                continue
            pinv = np.linalg.pinv(np.column_stack([np.ones(idx.size), u[idx] - side * at]))
            lines.append((idx, pinv, low + level * height, float(np.sqrt(np.sum(pinv[1] ** 2 * rec.sd[idx] ** 2)))))
    far = (u > FAR * u.max(), u < FAR * u.min())
    kernels = np.column_stack([np.ones(u.size), u / np.abs(u).max()] + [1.0 / (u * u + a * a) for a in (0.5, 1.0, 2.0, 4.0)] + [np.log1p(u * u / (a * a)) for a in (1.0, 4.0)])
    return {"u": u, "low": low, "height": height, "lines": lines, "right": u >= 0.0, "far": far, "above": (f - low) > height / 50.0, "kernels": kernels}


def _far(fr: dict, Y: np.ndarray) -> list:
    return [Y[:, side].mean(axis=1) for side in fr["far"] if side.any()]


def crossings(fr: dict, Y: np.ndarray) -> np.ndarray:
    """For each level and side: how far from the place that the fit gave it a straight line through the nearby samples
    meets the level, in widths. The distance is the level less the line, over its slope, and the division is softened:
    the slope over its square plus the square of its standard deviation. Then the mean of the samples far on each side."""
    out = []
    for idx, pinv, target, noise in fr["lines"]:
        ab = Y[:, idx] @ pinv.T
        out.append((target - ab[:, 0]) * ab[:, 1] / (ab[:, 1] ** 2 + noise ** 2))
    return np.column_stack(out + _far(fr, Y))


def bands(fr: dict, Y: np.ndarray) -> np.ndarray:
    """For each side and each band of value: how many samples have their value in the band, and the sum of their distances
    from the centre. For the lowest band, the sum of the distances and the sum of the values. The value is the noisy one,
    and a sample near an edge is shared between two bands over the width of the edge."""
    z = (Y - fr["low"]) / fr["height"]
    cuts = np.asarray(CUTS)
    over = np.clip((z[None, :, :] - cuts[:, None, None]) / EDGE + 0.5, 0.0, 1.0)
    inside = np.concatenate([1.0 - over[:1], over[:-1] - over[1:], over[-1:]])
    out = []
    for side in (fr["right"].astype(float), (~fr["right"]).astype(float)):
        away = np.abs(fr["u"]) * side
        out += [inside[0] @ away, (inside[0] * Y) @ side]
        counts, places = inside[1:] @ side, inside[1:] @ away
        for j in range(counts.shape[0]):
            out += [counts[j], places[j]]
    return np.column_stack(out)


def logarithm(fr: dict, Y: np.ndarray) -> np.ndarray:
    """The sums of the logarithm of the record, over its height and above its lowest fitted value, under eight kernels of
    the abscissa, on the samples at which the fitted curve is above a fiftieth of its height. Then the two far means."""
    above = fr["above"]
    value = np.log(np.clip((Y[:, above] - fr["low"]) / fr["height"], FLOOR, None))
    return np.column_stack([value @ fr["kernels"][above]] + _far(fr, Y))


STATISTICS = {"the abscissa against the value": crossings, "bands of value": bands, "the logarithm of the record": logarithm}


@dataclass
class Matched:
    shared: np.ndarray
    own: np.ndarray
    converged: bool
    evaluations: int
    mismatch: float                 # the squared length of what is left, each statistic in the covariance of the synthetic ones
    statistics: int


def match(lin, statistic: str, draws: int = 48, for_the_covariance: int = 400, seed: int = 1) -> Matched:
    """The estimate that one statistic gives: the parameters at which its values on the records match its mean on synthetic
    records of the model, weighted by the covariance of the synthetic values at the fit. The search starts at the fit."""
    if statistic not in STATISTICS:
        raise ValueError(f"no statistic is called {statistic!r}: {', '.join(STATISTICS)}")
    pr, ft = lin.problem, lin.fit
    if not ft.converged:
        raise ValueError("the fit did not end: a statistic is placed at a minimum of the fit, and this is not one")
    stat, R, s = STATISTICS[statistic], len(pr.records), len(pr.shared)
    rng = np.random.default_rng([97, seed])
    frames = [frame(lin, k) for k in range(R)]
    seen, whiten, noise = [], [], []
    for k, rec in enumerate(pr.records):
        seen.append(stat(frames[k], rec.y[None, :])[0])
        cloud = stat(frames[k], pr.predict(k, ft.shared, ft.own[k]) + rec.sd * rng.standard_normal((for_the_covariance, rec.x.size)))
        w, U = np.linalg.eigh(np.atleast_2d(np.cov(cloud.T)))
        keep = w > 1e-10 * w.max()
        whiten.append((U[:, keep] / np.sqrt(w[keep])).T)
        noise.append(rec.sd * rng.standard_normal((draws, rec.x.size)))          # the same draws at every value of the parameters
    shape = np.shape(ft.own)

    def left(z):
        shared, own = z[:s], z[s:].reshape(shape)
        return np.concatenate([whiten[k] @ (seen[k] - stat(frames[k], pr.predict(k, shared, own[k]) + noise[k]).mean(axis=0)) for k in range(R)])
    start = np.concatenate([np.asarray(ft.shared, dtype=float), np.asarray(ft.own, dtype=float).ravel()])
    if left(start).size < start.size:
        raise ValueError("this statistic gives fewer numbers than the model has parameters on these records")
    got = least_squares(left, start, x_scale="jac", method="trf", diff_step=1e-4, max_nfev=60 * start.size)
    return Matched(got.x[:s], got.x[s:].reshape(shape), bool(got.status > 0), int(got.nfev), float(2.0 * got.cost), int(got.fun.size))
