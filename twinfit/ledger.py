"""The first order ledger of a least squares fit, for any model: what a term that the model lacks does to the estimates, what
it leaves in the residuals, and what no statistic of the same records can tell.

Nothing is fitted here and no noise is drawn. This module holds canonical algebra: the first order shift of a fit under a
discrepancy of the mean, the measures of reliability of the detection of a model error that geodesy uses, the contrast of
two fits, and partialling out. The procedure that uses it is the protocol of this package (agent/AGENT.md).

The setting. The records are n samples y = f(theta) + delta + noise. The working model f is fitted by weighted least squares
and delta is what the working model lacks. J is the n by p matrix of the derivatives of f at the fitted values, each row
divided by the standard deviation of its sample, and d is delta divided in the same way. At first order

    the estimate moves by            S d,  with  S = (J'J)^-1 J'     the absorbed part, which is the bias
    the whitened residual keeps      (I - P) d,  with  P = J S       the visible part, which is all that an alarm can read

The two parts are orthogonal. A discrepancy inside the span of J, d = J b, moves the estimate by b and leaves no residual. The
fit on a subset of the samples and the fit under another loss move by the same b, so no window and no loss sees it either:
only another design or knowledge from outside the records does.

    led = Ledger(J, names)                     J whitened, one column a parameter, shared and nuisance alike
    led.sd, led.correlation()                  the first order bounds of the fit and their correlation
    led.combinations(targets, values)          the well and the badly determined combinations of some parameters
    led.shift(d), led.shift_in_sd(d)           the absorbed part of a discrepancy
    led.visible(d), led.visibility(d)          the visible part, and its share of the norm
    led.harm(d, targets)                       the shift in the plane of the targets, in joint standard deviations
    led.extended(B, targets)                   what freeing the columns B costs each target, and whether B is told from J
    contrast(led, keep, d, targets)            the fit on a subset of the samples against the fit on all of them
    size_at(noncentrality, dof, ...)           the size of a term at which an alarm fires with a stated probability
    silent_zone(led, d, targets, ...)          the harm at that size: the bias that the alarm lets through
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2, ncx2, norm


class Ledger:
    """The first order ledger of one design: J is the whitened n by p matrix of derivatives of the working model."""

    def __init__(self, J, names=None):
        J = np.asarray(J, dtype=float)
        if J.ndim != 2 or J.shape[0] <= J.shape[1]:
            raise ValueError(f"the design must hold more samples than parameters, and it holds {J.shape}")
        scale = np.sqrt((J * J).sum(axis=0))
        if np.any(scale == 0):
            raise np.linalg.LinAlgError(f"the records hold nothing about parameter(s) {np.where(scale == 0)[0].tolist()}")
        Q, R = np.linalg.qr(J / scale)                     # the columns scaled to unit norm: the units of the parameters do not condition the inverse
        diag = np.abs(np.diag(R))
        if diag.size and diag.min() < 1e-10 * diag.max():          # a design with no shared parameter fixes all that it has
            raise np.linalg.LinAlgError("the design does not fix every parameter: a column is a combination of the others")
        self.J, self.Q = J, Q
        self.n, self.p = J.shape
        self._solve = np.linalg.inv(R) / scale[:, None]    # S = _solve Q'
        self.C = self._solve @ self._solve.T               # (J'J)^-1
        self.names = tuple(names) if names is not None else tuple(f"p{i}" for i in range(self.p))
        if len(self.names) != self.p:
            raise ValueError(f"{self.p} parameters and {len(self.names)} names")

    # ---- the design alone
    @property
    def sd(self) -> np.ndarray:
        return np.sqrt(np.diag(self.C))

    def correlation(self) -> np.ndarray:
        s = self.sd
        return self.C / np.outer(s, s)

    def index(self, targets) -> list:
        if targets is None:
            return list(range(self.p))
        return [self.names.index(t) if isinstance(t, str) else int(t) for t in targets]

    def combinations(self, targets, values=None) -> list:
        """The combinations of the targets that the records fix best and worst: the eigenvectors of their covariance, in the
        logarithm of the targets when `values` gives their values. Each entry: (standard deviation, coefficients)."""
        idx = self.index(targets)
        sub = self.C[np.ix_(idx, idx)]
        if values is not None:
            v = np.asarray(values, dtype=float)
            sub = sub / np.outer(v, v)
        w, vec = np.linalg.eigh(sub)
        return [(float(math.sqrt(max(wi, 0.0))), vec[:, i].copy()) for i, wi in enumerate(w)]

    # ---- a discrepancy of the mean
    def shift(self, d) -> np.ndarray:
        """The first order move of every estimate under the whitened discrepancy d: the absorbed part."""
        return self._solve @ (self.Q.T @ np.asarray(d, dtype=float))

    def shift_in_sd(self, d) -> np.ndarray:
        return self.shift(d) / self.sd

    def visible(self, d) -> np.ndarray:
        """What the whitened residual of the fit keeps of d: its part outside the span of the design."""
        d = np.asarray(d, dtype=float)
        return d - self.Q @ (self.Q.T @ d)

    def noncentrality(self, d) -> float:
        e = self.visible(d)
        return float(e @ e)

    def visibility(self, d) -> float:
        """The share of the norm of d that is visible, between 0 (absorbed whole) and 1 (no bias): the sine of the angle
        between d and the span of the design."""
        d = np.asarray(d, dtype=float)
        total = float(d @ d)
        if total == 0.0:
            raise ValueError("a discrepancy of zero has no visibility")
        return math.sqrt(max(self.noncentrality(d), 0.0) / total)

    def harm(self, d, targets=None) -> float:
        """The shift in the plane of the targets, in their joint standard deviations."""
        idx = self.index(targets)
        step = self.shift(d)[idx]
        return float(math.sqrt(step @ np.linalg.solve(self.C[np.ix_(idx, idx)], step)))

    def directed(self, v, d=None) -> tuple:
        """A statistic of the residuals aimed at one direction: v'e, e the whitened residual. Returns (its standard deviation
        under the working model, its mean under the discrepancy d or None). A direction inside the span of the design has no
        variance left and reads nothing: it is refused."""
        v = np.asarray(v, dtype=float)
        left = self.visible(v)
        var = float(left @ left)
        if var < 1e-12 * float(v @ v):
            raise np.linalg.LinAlgError("the direction lies in the span of the design: the fit absorbs it and no residual statistic reads it")
        return math.sqrt(var), (None if d is None else float(left @ np.asarray(d, dtype=float)))

    def extended(self, B, targets=None) -> dict:
        """The working model with the columns B added and free. For each target: the factor by which its standard deviation
        grows. And for each column of B, how much of it is outside the span of the working design (its visibility): a column
        with none cannot be told from the parameters, and the extended model is refused."""
        B = np.asarray(B, dtype=float)
        B = B[:, None] if B.ndim == 1 else B
        vis = [self.visibility(B[:, j]) for j in range(B.shape[1])]
        wide = Ledger(np.column_stack([self.J, B]), self.names + tuple(f"added {j}" for j in range(B.shape[1])))
        idx = self.index(targets)
        return {"visibility_of_each_column": vis, "sd_factor": {self.names[i]: float(wide.sd[i] / self.sd[i]) for i in idx},
                "sd_of_added": [float(v) for v in wide.sd[self.p:]], "ledger": wide}


def contrast(led: Ledger, keep, d=None, targets=None) -> dict:
    """The fit on the samples `keep` (a mask) against the fit on all of them. Under the working model the difference of the
    two estimates has covariance C_keep - C. Returns that covariance on the targets, and under the discrepancy d the difference
    of the two shifts and its size in standard deviations of the difference. A discrepancy inside the span of the design gives
    a difference of zero: the contrast reads the visible part alone."""
    keep = np.asarray(keep, dtype=bool)
    part = Ledger(led.J[keep], led.names)
    idx = led.index(targets)
    cov = part.C[np.ix_(idx, idx)] - led.C[np.ix_(idx, idx)]
    out = {"covariance": cov, "sd_factor": {led.names[i]: float(part.sd[i] / led.sd[i]) for i in idx}, "ledger": part}
    if d is not None:
        d = np.asarray(d, dtype=float)
        diff = part.shift(d[keep])[idx] - led.shift(d)[idx]
        out["difference"] = diff
        out["z"] = float(math.sqrt(max(diff @ np.linalg.solve(cov, diff), 0.0)))
    return out


def alarm_power(noncentrality: float, dof: int, alpha: float = 0.05) -> float:
    """The probability that a chi square statistic with `dof` degrees of freedom and the given noncentrality passes its
    threshold at the level given. The exact noncentral law at every number of degrees of freedom: a normal form above 200
    gave 0.0446 at no noncentrality where 0.05 was asked, and a step at the join (review 32, minor 2)."""
    threshold = chi2.isf(alpha, dof)
    return float(ncx2.sf(threshold, dof, noncentrality)) if noncentrality > 0 else alpha


def size_at(noncentrality: float, dof: int, alpha: float = 0.05, power: float = 0.5) -> float:
    """The factor by which a term must be scaled for a chi square statistic with `dof` degrees of freedom, on which the term
    at its present size has the given noncentrality, to pass its threshold at the level given with the stated probability. With
    dof = 1 it is a two sided test of one directed statistic."""
    if noncentrality <= 0.0:
        return math.inf
    if not alpha < power < 1.0:
        raise ValueError("the probability asked of the alarm must lie between its level and one")

    def gap(log_s):
        return alarm_power(math.exp(2.0 * log_s) * noncentrality, dof, alpha) - power
    lo = hi = -0.5 * math.log(noncentrality)                # the size at which the noncentrality is one
    while gap(lo) > 0.0:
        lo -= 2.0
    while gap(hi) < 0.0:
        hi += 1.0
        if hi > 60.0:
            return math.inf
    return float(math.exp(brentq(gap, lo, hi)))


def silent_zone(led: Ledger, d, targets=None, alpha: float = 0.05, power: float = 0.5, direction=None) -> dict:
    """The size at which an alarm sees the term with the stated probability, and the harm that the term does at that size:
    the bias, in joint standard deviations of the targets, that the alarm lets through. The alarm is the chi square of the
    fit, with n - p degrees of freedom, or the statistic aimed at `direction` when one is given. A harm above one says that
    the term biases the targets by more than their uncertainty before the alarm is likely to fire."""
    d = np.asarray(d, dtype=float)
    if direction is None:
        seen = led.noncentrality(d)
        size = math.inf if seen <= 1e-12 * float(d @ d) else size_at(seen, led.n - led.p, alpha, power)   # under rounding, nothing is visible
    else:
        sd, mean = led.directed(direction, d)
        size = size_at((mean / sd) ** 2, 1, alpha, power)
    harm = led.harm(d, targets)
    return {"size": size, "harm_at_that_size": size * harm, "harm_at_size_one": harm, "size_of_one_joint_sd": (1.0 / harm if harm > 0 else math.inf)}


def two_sided_power(z: float, alpha: float = 0.05) -> float:
    """The probability that a statistic with mean z, in its own standard deviations, passes the two sided threshold at the level given."""
    c = norm.isf(alpha / 2.0)
    return float(norm.sf(c - z) + norm.cdf(-c - z))
