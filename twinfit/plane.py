"""The combination: the plane of position and residual value, cut on both axes, as one estimator.

The weighted fit asks one thing of the residuals: that the sum of their squares be least. This estimator asks more. The
residuals of a candidate, each over its standard deviation, are laid on the plane of position and value. The plane is cut
in cells: slices of position within each record, and bands of value with smooth edges. Under a right model the value of a
residual does not depend on its position, so the table of the cells has no structure, and the estimate is the root of the
weighted equations that say so. The weight is the covariance that the table has when its errors are drawn again from the
residuals at hand, known in closed form: no synthetic curve is drawn for it. Two more equations fix what a table cannot:
the mean of the residuals is zero and their mean square is one. The stack adds, in each slice, three moments of the
residuals with the joint covariance of all of it: their sum, the sum of their sizes, the sum of the roots of their sizes.
Nothing is added as if it were independent: a sum of such terms taken as independent prints an uncertainty that is 0.58
to 0.62 of the real scatter (method/PACKAGE.md, section 10).

    plane_fit(problem, start, kind="cells")     the estimate, from the result of the weighted fit as its start
    kind "cells": the table of the cells. kind "stack": the cells and the three moments by slice
    kind "hybrid": the equations of the weighted fit and those of the cells, as one system with its joint covariance
    columns(problem, start, sets=100)           the three columns for every shared parameter, each with the uncertainty that it
                                                prints and with what synthetic sets of the fitted model say of it

The uncertainty that an estimator prints is a first order one. Where two parameters trade against each other and the noise is
high it can be half of the real scatter, for the weighted fit as for the others (method/PACKAGE.md, section 10: at ten times
the reference noise the printed intervals of the weighted fit held the value put in on 0.66 to 0.75 of 200 sets, not on
0.95). `columns` therefore reads every column on synthetic sets of the fitted model and gives, beside what is printed, the
bias, the spread and the share of the printed intervals that hold the value put in. The spread on those sets is the
uncertainty to quote.

The hybrid. The weighted fit solves its own equations: the sum, over the samples, of each residual times the derivative of the
mean by a parameter is zero. The cells of a slice add up functions of the same residuals, so the two sets of equations are
correlated, and the correlation is known in closed form from the residuals at hand. The hybrid takes the part of the fit's
equations that the table and the mean of the residuals do not already hold, with the covariance of that part, and solves all
of it at once.

What was measured, and where (method/PACKAGE.md, section 10, and reviews/REVIEW_58.md). The experiments of the source
repository measured the cells with a covariance drawn from synthetic curves. Here that covariance is in closed form, which is
another estimator, measured for the first time in this package, and the hybrid is new here. Against the weighted fit, in the
root of the mean squared error of two widths that trade against each other:

    normal errors                    the hybrid loses 3 per cent on 4 records of 365 samples and under 1 per cent on records of
                                     1177, on a thousand sets each. The cells alone lose 10 to 26 per cent: they are not an
                                     estimate to quote there
    errors of a Student law, 3 d.f.  the hybrid has 0.71 to 0.77 of the fit's error. The Student likelihood, a canonical
                                     technique, has 0.67 to 0.69: the hybrid goes under no canonical technique at its
                                     strongest. Like that likelihood, it follows the law of the errors without being told it
    a term that the model lacks      no protection: it is biased as the fit is

Its printed uncertainty is a first order one, and it is short of the real scatter by an amount that falls as the records get
longer: about 3 to 9 per cent on four records of 365 samples (slides/model_steps/decisive_trial.json of the source repository),
and on a linear problem with a known bound a reader of review 59 measured 24, 15 and 6 per cent at 100, 200 and 400 samples a
record, with 95 per cent intervals that held the value on 0.85, 0.90 and 0.93 of a thousand sets (reviews/REVIEW_59.md,
moderate 1). On one record of one curve it was 14 per cent short at 365 samples. So at the floor of 200 samples the printed
intervals hold about 0.90 and not 0.95. At a high noise, where parameters trade against each other, the printed uncertainty
can be half of the scatter, as the weighted fit's is. `columns` measures the uncertainty on synthetic sets and is what to quote. The kind "stack" is not usable: its search ended on 11 of 60 sets
under normal errors, and `plane_fit` warns when it is asked for.

Here the standard deviations of the records are taken as given, up to one common factor that the estimator fits, unless a law
of the noise is given to `plane_fit`: then they follow the mean of the candidate.

What it does not do. It reads the same model as the weighted fit: a term that the model lacks biases it as it biases the
fit. What it adds there is the table itself, `distance`, which grows when the plane holds a structure.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
from scipy.special import expit
from scipy.stats import norm

from .fit import STEP, derivatives, fit
from .pool import Lazy, in_order
from .problem import Record

BANDS = 8                                                # bands of the value of a residual, of equal probability under a normal law
TAU = 0.5                                                # the width of a band's edge, in standard deviations of a residual. 0.25 until review 58: the mean
                                                         # derivative of a band, read from the residuals at hand, is rough with a narrow edge. One step was
                                                         # taken, to 0.5, on runs whose sets were not kept: the figures that stood here are in no file and are
                                                         # withdrawn (review 59, moderate 10). What a file holds: with 0.5 the hybrid loses 2.8 and 2.9 per cent
                                                         # to the weighted fit on 4 records of 365 samples, a thousand sets, every search ended
                                                         # (slides/model_steps/decisive_trial.json of the source repository). It is not an optimum: the
                                                         # reviewer's reader found 0.75 and 1.0 better on two other designs, at the cost of searches that end
EDGES = norm.ppf(np.arange(1, BANDS) / BANDS)
SOFT = 0.05                                              # the size under which the root of |e| is rounded
SLICES = 26                                              # slices of position in each record
FLOOR = 200                                              # the fewest samples of a record that the plane is read on. Under it the table has too many cells for its
                                                         # samples: at 40 samples a record it printed an uncertainty of 16 where the fit's was 0.020, with no word
                                                         # (review 58, moderate 6). The floor is where the table stops failing, not where the printed uncertainty
                                                         # becomes honest: at 200 samples a record it is still about 15 per cent short (the docstring above)
STEP_DONE = 1e-4                                         # the search has ended when the square of the step left, in standard deviations, is under this


@dataclass
class PlaneFit:
    shared: np.ndarray
    own: np.ndarray                 # one row a record
    cov: np.ndarray                 # of the shared parameters
    converged: bool
    iterations: int
    distance: float                 # the sum of squares of the whitened table at the root: how far the plane is from one with no structure
    equations: int                  # the number of those squares
    scale: float                    # the common factor on the standard deviations that the estimator fitted
    kind: str = "cells"

    @property
    def sd(self) -> np.ndarray:
        return np.sqrt(np.diag(self.cov))


def bands(e):
    """The first BANDS - 1 smooth bands of the value of a residual and their derivatives. The bands sum to 1, so the last one
    holds nothing that the others do not."""
    s = expit((np.asarray(e, dtype=float)[:, None] - EDGES[None, :]) / TAU)
    ds = s * (1.0 - s) / TAU
    upper = np.hstack([np.ones((s.shape[0], 1)), s[:, :-1]])
    dupper = np.hstack([np.zeros((s.shape[0], 1)), ds[:, :-1]])
    return upper - s, dupper - ds


def functions(e, kind: str):
    """Values and derivatives of the functions of a residual: the bands, and for the stack e, the rounded |e| and the rounded
    root of |e|. The size stands where the square would: under a law with heavy tails the square has no variance."""
    val, dval = bands(e)
    if kind == "stack":
        q2 = e * e + SOFT * SOFT
        val = np.column_stack([val, e, q2 ** 0.5, q2 ** 0.25])
        dval = np.column_stack([dval, np.ones_like(e), e * q2 ** -0.5, 0.5 * e * q2 ** -0.75])
    return val, dval


def free_columns(kind: str) -> list:
    """The functions whose mean under the law of the errors is not known and is fitted: the bands, and for the stack |e| and its root."""
    return list(range(BANDS - 1)) + ([BANDS, BANDS + 1] if kind == "stack" else [])          # cells and hybrid: the bands alone


CUT_OF_V = 1e-6                                          # the cut for the covariance of what the hybrid takes of the fit, and for that matrix alone


def half_inverse(Cm, rcond: float = 1e-10):
    """H with H' H the generalised inverse of the covariance Cm. A covariance has no negative eigenvalue: one that is
    past rounding says that Cm is not a covariance, and it is told, not dropped in silence (review 58).

    A direction whose eigenvalue is under `rcond` times the largest is dropped. Two matrices are inverted here and they do
    not take the same cut. The covariance of the band functions keeps the cut of a ten billionth: its smallest real
    eigenvalue is 1.2 millionths of the largest (review 60), and a cut of a millionth sat on it. The covariance of what the
    hybrid takes of the fit, V, has one direction that is null by construction only when the standard deviation is constant
    in each record. When it varies inside a record that eigenvalue is small and not zero, about 3.5e-8 of the largest in the
    case of review 59, and with a ten billionth it was kept and weighted by thousands: the search of the hybrid ended on 1
    set of 20. V is cut at `CUT_OF_V`. A first form of that repair raised the cut of both matrices to a millionth: the
    search then stopped ending when the given standard deviations were off by a common factor, on 1 set of 12 at a factor
    of 2 (review 60, severe 1). The cut of V is still a fraction of its largest eigenvalue, so it depends on the units of
    the parameters, and it was measured for a weak dependence of the noise on the signal only (review 60, moderate 1)."""
    w, V = np.linalg.eigh(Cm)
    if w.min() < -1e-6 * w.max():
        warnings.warn(f"a matrix taken as a covariance has an eigenvalue of {w.min() / w.max():.2g} times its largest: the direction is dropped, and the estimate that uses it is not to be trusted", stacklevel=2)
    keep = w > rcond * w.max()
    return (V[:, keep] / np.sqrt(w[keep])).T


def bounded_step(H, g, theta, lower, upper):
    """The solution of H step = g with every parameter that sits at a bound and would leave the domain held there."""
    free, step = np.ones(theta.size, dtype=bool), np.zeros(theta.size)
    for _ in range(theta.size):
        idx = np.where(free)[0]
        step = np.zeros(theta.size)
        if idx.size:
            step[idx] = np.linalg.solve(H[np.ix_(idx, idx)], g[idx])
        out = free & (((theta <= lower) & (step < 0)) | ((theta >= upper) & (step > 0)))
        if not out.any():
            break
        free &= ~out
    return step


def law_derivatives(problem, k: int, shared: np.ndarray, own: np.ndarray, variance) -> tuple:
    """As `fit.derivatives`, with the standard deviation of every sample taken from the law of the noise at the mean of the
    candidate, `variance(mean, record, own)`: the whitened residual, the whitened derivatives of the mean by the shared and
    by the own parameters, and the derivatives of the logarithm of the standard deviation by the same parameters."""
    rec = problem.records[k]

    def at(sh, ow):
        mean = problem.predict(k, sh, ow)
        var = np.asarray(variance(mean, rec, ow), dtype=float)
        if var.shape != rec.x.shape or not np.all(np.isfinite(var)) or not np.all(var > 0):
            raise ValueError(f"record {rec.name!r}: the law must give one positive variance for each sample")
        return mean, var
    f0, v0 = at(shared, own)
    sd0 = np.sqrt(v0)

    def column(vec: np.ndarray, j: int, is_shared: bool) -> tuple:
        h = STEP * max(abs(vec[j]), 1e-3 * (np.max(np.abs(vec)) or 1.0), 1e-9)
        up, down = vec.copy(), vec.copy()
        up[j] += h
        down[j] -= h
        (hi, vh), (lo, vl) = (at(up, own), at(down, own)) if is_shared else (at(shared, up), at(shared, down))
        return (hi - lo) / (2.0 * h) / sd0, (np.log(vh) - np.log(vl)) / (4.0 * h)          # the logarithm of a standard deviation is half that of its variance
    n = rec.x.size
    cols_a = [column(shared, j, True) for j in range(shared.size)]
    cols_b = [column(own, j, False) for j in range(own.size)]
    stack = lambda cols, i: np.column_stack([c[i] for c in cols]) if cols else np.zeros((n, 0))          # noqa: E731
    return (rec.y - f0) / sd0, stack(cols_a, 0), stack(cols_b, 0), stack(cols_a, 1), stack(cols_b, 1)


class Plane:
    """The table of the cells against a table with no structure, as a sum of squares. The parameters are the shared ones, the
    own ones of every record, the logarithm of one factor on the standard deviations, and the means that are not known."""

    def __init__(self, problem, kind: str = "cells", slices: int = SLICES, variance=None):
        if kind not in ("cells", "stack", "hybrid"):
            raise ValueError("kind is cells, stack or hybrid")
        self.problem, self.kind, self.variance = problem, kind, variance
        self.p, self.q, self.R = len(problem.shared), len(problem.own), len(problem.records)
        starts, at = [], 0
        for rec in problem.records:
            n = rec.x.size
            k = max(1, min(slices, n // 8))               # at least 8 samples a slice
            starts += [at + int(i) for i in np.linspace(0, n, k, endpoint=False)]
            at += n
        self.starts, self.N = np.array(starts), at
        self.n = np.diff(np.append(self.starts, self.N)).astype(float)
        self.root = np.sqrt(self.n)
        self.cols = free_columns(kind)
        self.K = BANDS - 1 + (3 if kind == "stack" else 0)
        self.m = self.p + self.R * self.q + 1             # the parameters of the model and the factor on the standard deviations

    def split(self, par):
        own = par[self.p:self.p + self.R * self.q].reshape(self.R, self.q)
        return par[:self.p], own, float(par[self.m - 1]), par[self.m:]

    def pieces(self, par):
        """The residuals over their standard deviations, and for every slice the sum of minus their derivatives by the parameters
        of the model: a record's own parameters reach its own slices alone."""
        shared, own, log_s, _ = self.split(par)
        s = np.exp(log_s)
        e, sums, at, row = np.empty(self.N), np.zeros((self.starts.size, self.m)), 0, 0
        w = self.m - 1                                     # the parameters of the model: those of the equations of the weighted fit
        self.DtD, self.Dte = np.zeros((w, w)), np.zeros(w)
        # for every slice, the sum of what a residual moves by, as a share of itself, when a parameter moves: 1 for the
        # factor on the standard deviations, and with a law of the noise the derivative of the logarithm of its standard deviation
        self.ks = np.zeros((self.starts.size, self.m))
        self.ks[:, self.m - 1] = self.n
        for k in range(self.R):
            if self.variance is None:
                (r, A, B), GA, GB = derivatives(self.problem, k, shared, own[k]), None, None
            else:
                r, A, B, GA, GB = law_derivatives(self.problem, k, shared, own[k], self.variance)
            n = r.size
            A, B, ek = A / s, B / s, r / s
            e[at:at + n] = ek
            local = self.starts[(self.starts >= at) & (self.starts < at + n)] - at
            rows, mine = slice(row, row + local.size), slice(self.p + k * self.q, self.p + (k + 1) * self.q)
            if self.p:
                sums[rows, :self.p] = np.add.reduceat(A, local, axis=0)
                self.DtD[:self.p, :self.p] += A.T @ A
                self.Dte[:self.p] += A.T @ ek
            if self.q:
                sums[rows, mine] = np.add.reduceat(B, local, axis=0)
                self.DtD[mine, mine] += B.T @ B
                self.Dte[mine] += B.T @ ek
                if self.p:
                    self.DtD[:self.p, mine] += A.T @ B
                    self.DtD[mine, :self.p] += B.T @ A
            if GA is not None:
                if self.p:
                    self.ks[rows, :self.p] = np.add.reduceat(GA, local, axis=0)
                if self.q:
                    self.ks[rows, mine] = np.add.reduceat(GB, local, axis=0)
            at, row = at + n, row + local.size
        return e, sums

    def means(self, par):
        g = np.zeros(self.K)
        g[self.cols] = par[self.m:]
        return g

    def start_means(self, par):
        val, _ = functions(self.pieces(par)[0], self.kind)
        return val.mean(axis=0)[self.cols]

    def freeze(self, par):
        """What a step holds fixed, read from the residuals at hand: the whitening of the functions, the variance of e squared,
        and the means of the functions' derivatives and of e times them."""
        e = self.pieces(par)[0]
        val, dval = functions(e, self.kind)
        half = half_inverse(np.cov(val, rowvar=False))
        n = e.size
        with_e = ((e[:, None] * val).mean(axis=0) - e.mean() * val.mean(axis=0)) * n / (n - 1.0)          # the covariance of a residual with each function of it, counted as np.cov counts
        return half, float(max(np.var(e * e), 1e-6)), dval.mean(axis=0), (e[:, None] * dval).mean(axis=0), half @ with_e, float(np.mean(e ** 3)), float(np.var(e, ddof=1))

    def remainder(self, state, ms) -> tuple:
        """The covariance of what the hybrid takes of the fit's equations, with the pieces that build it: (V, M, total, hc).

        What the others hold of the fit's equations (review 58, moderate 2). The cells are read with their means free, so the
        part of the table that is one occupancy for all the slices is spent on those means: what is left of the table holds
        hc.hc of each fit equation less its pooled part, M'M - total total'. The mean of the residuals holds the pooled part
        whole, total total'. The mean square holds nothing more: its covariance with the fit's equations is that of the mean
        times the skewness, and given the mean it adds none. A first form subtracted the three as if they were disjoint and
        counted the pooled part twice: its matrix had an eigenvalue near minus the largest. Every covariance here is that of
        the residuals at hand, counted one way, so that hc.hc cannot pass their variance and V is a covariance."""
        hc, var = state[4], state[6]                       # var: the variance of the residuals at hand, near 1
        w, N = self.m - 1, float(self.N)
        M = ms[:, :w] / np.sqrt(N * self.n)[:, None]       # for each slice: the covariance of the fit's equations with a whitened cell, over hc
        total = ms[:, :w].sum(axis=0) / N                  # the covariance of the fit's equations with the mean of the residuals
        share = float(hc @ hc)
        return var * self.DtD / N - share * (M.T @ M) - (var - share) * np.outer(total, total), M, total, hc

    def fit_part(self, state, ms, r_cells, tails, J_cells=None, J_tails=None):
        """The equations of the weighted fit, less what the table and the two last equations already hold of them, over the
        root of their own covariance. With the Jacobians of the others given, the same for the Jacobian."""
        V, M, total, hc = self.remainder(state, ms)
        w, N = self.m - 1, float(self.N)
        root = half_inverse(V, CUT_OF_V)                   # the one direction that the others hold whole, the pooled offset, drops out here
        g = self.Dte / np.sqrt(N) - M.T @ (r_cells.reshape(self.starts.size, -1) @ hc) - total * tails[0]
        if J_cells is None:
            return root @ g, None
        J = np.zeros((w, J_cells.shape[1]))
        J[:, :w] = self.DtD / np.sqrt(N)
        J -= M.T @ np.einsum("k,skp->sp", hc, J_cells.reshape(self.starts.size, hc.size, -1)) + np.outer(total, J_tails[0])
        return root @ g, root @ J

    def residual(self, par, state):
        half, kappa = state[:2]
        e, ms = self.pieces(par)
        val, _ = functions(e, self.kind)
        T = (np.add.reduceat(val, self.starts, axis=0) - self.n[:, None] * self.means(par)[None, :]) / self.root[:, None]
        r = (T @ half.T).ravel()
        tails = [e.sum() / np.sqrt(e.size), (e * e - 1.0).sum() / np.sqrt(e.size * kappa)]
        if self.kind == "hybrid":
            return np.concatenate([r, tails, self.fit_part(state, ms, r, tails)[0]])
        return np.concatenate([r, tails if self.kind == "cells" else tails[1:]])

    def jacobian(self, par, state):
        """The expected derivative of the residual, with the sign of the solver: what the table moves by, on average over the
        errors, when a parameter moves. The derivative of the table of one set is rough at the scale of a band's edge, and
        its squares would count that roughness as information, so the solver and the uncertainty use the expected one."""
        half, kappa, a, b = state[:4]
        e, ms = self.pieces(par)
        ks = self.ks                                       # the factor on the standard deviations, and the law of the noise where one is given
        ha, hb = half @ a, half @ b
        Jphi = ((ms[:, None, :] * ha[None, :, None] + ks[:, None, :] * hb[None, :, None]) / self.root[:, None, None]).reshape(-1, self.m)
        Jeta = (self.root[:, None, None] * half[None, :, self.cols]).reshape(-1, len(self.cols))
        J = np.hstack([Jphi, Jeta])
        tail = np.zeros((2, J.shape[1]))
        tail[0, :self.m] = ms.sum(axis=0) / np.sqrt(e.size)
        tail[1, :self.m] = 2.0 * ks.sum(axis=0) / np.sqrt(e.size * kappa)
        if self.kind == "hybrid":
            val, _ = functions(e, self.kind)
            T = (np.add.reduceat(val, self.starts, axis=0) - self.n[:, None] * self.means(par)[None, :]) / self.root[:, None]
            tails = [e.sum() / np.sqrt(e.size), (e * e - 1.0).sum() / np.sqrt(e.size * kappa)]
            return np.vstack([J, tail, self.fit_part(state, ms, (T @ half.T).ravel(), tails, J, tail)[1]])
        return np.vstack([J, tail if self.kind == "cells" else tail[1:]])


def root_solve(prob, par0, lower, upper, max_iter: int = 80):
    """The root of the estimating equation J' r = 0, J the expected derivative of the residual r, by the steps (J' J)^-1 J' r,
    each halved until the step left gets shorter. J is not the gradient of r' r, so a search that accepts a step only when
    r' r falls stalls. Returns the root, (J' J)^-1 there, r' r there, the step left and the number of iterations."""
    def parts(par):
        state = prob.freeze(par)
        r, J = prob.residual(par, state), prob.jacobian(par, state)
        H = J.T @ J
        U = J.T @ r
        d = np.clip(np.diag(H), 1e-300, None)
        step = bounded_step(H + 1e-12 * np.diag(d), U, par, lower, upper)
        return r, H, step, float(max(step @ U, 0.0))

    par = np.clip(np.asarray(par0, dtype=float), lower, upper)
    r, H, step, left = parts(par)
    it = 0
    for it in range(1, max_iter + 1):
        if left < STEP_DONE:
            break
        scale, got = 1.0, None
        while scale > 1e-3:
            new = np.clip(par + scale * step, lower, upper)
            trial = parts(new)
            if np.isfinite(trial[3]) and trial[3] < left:
                got = trial
                break
            scale *= 0.5
        if got is None:
            break
        par = new
        r, H, step, left = got
    d = np.sqrt(np.clip(np.diag(H), 1e-300, None))
    cov = np.linalg.pinv(H / np.outer(d, d), rcond=1e-12, hermitian=True) / np.outer(d, d)
    return par, cov, float(r @ r), left, it, r.size


def plane_fit(problem, start, kind: str = "cells", slices: int = SLICES, max_iter: int = 80, variance=None) -> PlaneFit:
    """The root of the weighted equations of the cells, or of the cells and the moments by slice, from the result `start` of
    the weighted fit of the same problem.

    `variance(mean, record, own)`: the law of the noise, as `noise.settle` takes it. With it the standard deviation of every
    sample follows the mean of the candidate, where without it the standard deviations of the records are held. Where the
    variance of the noise follows the signal, the size of the residuals then says something of the mean, which the weighted
    fit does not use: at first order nothing at a signal to noise ratio of tens, and much where a record barely shows its
    curve (method/PACKAGE.md, section 10)."""
    if kind == "stack":
        warnings.warn("the kind 'stack' is not usable yet: its search ended on 11 of 60 sets under normal errors (method/PACKAGE.md, section 10)", stacklevel=2)
    start_own = np.asarray(start.own, dtype=float)
    if np.asarray(start.shared).shape != (len(problem.shared),) or start_own.shape != (len(problem.records), len(problem.own)):
        raise ValueError(f"the start does not belong to this problem: it has {np.asarray(start.shared).size} shared values and own values of shape {start_own.shape}, "
                         f"and the problem has {len(problem.shared)} shared parameters and {len(problem.records)} records of {len(problem.own)} own ones")
    if not (np.all(np.isfinite(np.asarray(start.shared, dtype=float))) and np.all(np.isfinite(start_own))):
        raise ValueError("the start holds a value that is not a number")
    short = [rec.name or str(k) for k, rec in enumerate(problem.records) if rec.x.size < FLOOR]
    if short:
        raise ValueError(f"the plane of position and residual value needs at least {FLOOR} samples in every record, and {len(short)} of the {len(problem.records)} records have fewer "
                         f"(the first: {short[0]}). With fewer the table has too many cells for its samples: quote the weighted fit, or a robust fit")
    if slices < 4:
        raise ValueError("at least 4 slices a record: with fewer the table says nothing of where along a record a residual lies")
    prob = Plane(problem, kind, slices, variance)
    R, q, p = prob.R, prob.q, prob.p
    lo_s, hi_s = problem.bounds_shared
    lo_o = np.concatenate([np.asarray(problem.own_bounds(k)[0], dtype=float) if problem.bounded else np.full(q, -np.inf) for k in range(R)]) if q else np.zeros(0)
    hi_o = np.concatenate([np.asarray(problem.own_bounds(k)[1], dtype=float) if problem.bounded else np.full(q, np.inf) for k in range(R)]) if q else np.zeros(0)
    more = len(prob.cols) - (BANDS - 1)
    lower = np.concatenate([np.asarray(lo_s, dtype=float), lo_o, [-3.0], np.full(BANDS - 1, 1e-6), np.full(more, 0.05)])
    upper = np.concatenate([np.asarray(hi_s, dtype=float), hi_o, [3.0], np.full(BANDS - 1, 1.0), np.full(more, 5.0)])
    phi = np.concatenate([np.asarray(start.shared, dtype=float), np.asarray(start.own, dtype=float).ravel(), [0.0]])
    phi = np.clip(phi, lower[:prob.m], upper[:prob.m])
    eta = prob.start_means(np.concatenate([phi, np.zeros(len(prob.cols))]))
    par, cov, dist, left, it, rows = root_solve(prob, np.concatenate([phi, eta]), lower, upper, max_iter)
    shared, own, log_s, _ = prob.split(par)
    ended = bool(left < STEP_DONE)
    if not ended:                     # said, and not left to a field that a caller may not read (review 59, moderate 2)
        warnings.warn(f"the search of the kind '{kind}' did not end in {it} steps: its estimate and its uncertainty are not to be quoted. Quote the weighted fit", stacklevel=2)
    return PlaneFit(shared.copy(), own.copy(), cov[:p, :p], ended, it, dist, rows, float(np.exp(log_s)), kind)


KINDS = ("fit", "cells", "hybrid")


def _columns_of_one_set(args):
    """The three estimates of one synthetic set and the uncertainty that each prints. A column whose search fails or does
    not end is None for that set."""
    twin, start, slices, max_iter = args
    out = dict.fromkeys(KINDS)
    try:
        base = fit(twin, start.shared, start.own, max_iter=max_iter)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return out
    if base.converged:
        out["fit"] = (base.shared, base.sd)
    for kind in KINDS[1:]:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)          # a search that does not end is counted in `ended_on`, set by set
                got = plane_fit(twin, base, kind, slices, max_iter)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            continue
        if got.converged:
            out[kind] = (got.shared, got.sd)
    return out


def columns(problem, start, sets: int = 100, seed: int = 1, workers: int = 1, errors: str = "residuals", slices: int = SLICES, max_iter: int = 80,
            functions: dict | None = None) -> dict:
    """The three columns for every shared parameter: the weighted fit `start`, the cells and the hybrid. For each, the estimate,
    the uncertainty that it prints, and what `sets` synthetic sets of the fitted model say: the share of searches that ended,
    the bias, the spread, and the share of the printed 95 per cent intervals that hold the value put in. The spread is the
    uncertainty to quote: it is `honest` in the table.

    The synthetic sets hold the mean of the weighted fit and the standard deviations of the records. `errors="residuals"`:
    their errors are the residuals of the weighted fit over their standard deviations, drawn again with replacement, so
    that a law of the errors that is not normal is kept without being named. `errors="normal"`: they are drawn from a
    normal law. Every column is read on the same sets.

    `functions`: name -> a function of the vector of the shared parameters, for a quantity that is a combination of them. It
    gets the same three columns. Its printed uncertainty is the first order one, from the covariance of the shared parameters.

    At a high noise, where the parameters of each record are themselves poorly known, the spread on these sets is still too
    small: on a test at ten times the reference noise it was two thirds of the real spread, and its intervals held the value on
    0.85 and 0.77 of 48 sets where the printed ones held it on 0.75 and 0.69 (method/PACKAGE.md, section 10, step 77). One
    reason that could account for it, NOT TESTED: the sets are drawn at fitted values, and a model built on estimates may tell
    its parts apart better than the true one does (review 59, moderate 4: this docstring gave it as the cause).

    What it does not show. The sets are those of the fitted model: a term that the model lacks is not in them, and the bias
    that it gives is in no column. With residuals drawn again the position of a structure in the residuals is lost, and
    its size is kept as if it were noise."""
    if errors not in ("residuals", "normal"):
        raise ValueError("errors is residuals or normal")
    if sets < 20:
        raise ValueError("at least 20 synthetic sets: a spread from fewer is not known to one part in six")
    found = {"fit": start, "cells": plane_fit(problem, start, "cells", slices, max_iter), "hybrid": plane_fit(problem, start, "hybrid", slices, max_iter)}
    means = [problem.predict(k, start.shared, start.own[k]) for k in range(len(problem.records))]
    pool = np.concatenate([(rec.y - m) / rec.sd for rec, m in zip(problem.records, means)])
    pool = pool - pool.mean()

    def one(rep: int):
        rng = np.random.default_rng([seed, rep])
        recs = []
        for rec, mean in zip(problem.records, means):
            e = rng.choice(pool, rec.x.size) if errors == "residuals" else rng.standard_normal(rec.x.size)
            recs.append(Record(rec.x, mean + rec.sd * e, rec.sd, rec.design, rec.name))
        return (problem.with_records(recs), start, slices, max_iter)

    read = in_order(_columns_of_one_set, Lazy(one, sets), workers)
    truth = np.asarray(start.shared, dtype=float)
    out = {"sets": sets, "errors_of_the_synthetic_sets": errors, "ended_on_the_records": {k: bool(r.converged) for k, r in found.items()}, "parameters": {}}
    # every quantity as a function of the shared vector: a parameter is the function that picks it
    quantities = {name: (lambda v, j=j: float(v[j])) for j, name in enumerate(problem.shared)}
    quantities.update(functions or {})

    def printed_of(fn, shared, cov) -> float:
        """The first order uncertainty of a function of the shared parameters."""
        shared = np.asarray(shared, dtype=float)
        grad = np.zeros(shared.size)
        for i in range(shared.size):
            h = 1e-6 * max(abs(shared[i]), 1e-3)
            up, down = shared.copy(), shared.copy()
            up[i] += h
            down[i] -= h
            grad[i] = (fn(up) - fn(down)) / (2.0 * h)
        return float(np.sqrt(max(grad @ cov @ grad, 0.0)))
    names = list(problem.shared)
    for j, (name, fn) in enumerate(quantities.items()):
        row = {}
        own_bar = j < len(names)                                # a shared parameter: the uncertainty that each set printed is kept. A function: the first order one of the records
        for kind in KINDS:
            got = [r[kind] for r in read if r[kind] is not None]
            cell = {"estimate": fn(found[kind].shared), "printed": printed_of(fn, found[kind].shared, found[kind].cov), "ended_on": len(got) / sets}
            if len(got) >= 20:
                est = np.array([fn(g[0]) for g in got])
                bar = np.array([g[1][j] for g in got]) if own_bar else np.full(est.size, cell["printed"])
                spread = float(est.std(ddof=1))
                cell.update({"bias": float(est.mean() - fn(truth)), "honest": spread, "printed_over_honest": float(np.median(bar) / spread) if spread > 0 else None,
                             "printed_intervals_that_hold_the_value_put_in": float(np.mean(np.abs(est - fn(truth)) < 1.96 * bar))})
            row[kind] = cell
        # is a column better than the weighted fit: the ratio of the two spreads over the sets on which both ended, with the
        # standard error that drawing those sets again gives. The columns are read on the same sets, so the ratio is known
        # better than either spread
        best, least = "fit" if "honest" in row["fit"] else None, 1.0
        for kind in KINDS[1:]:
            both = [(fn(r["fit"][0]), fn(r[kind][0])) for r in read if r["fit"] is not None and r[kind] is not None]
            if len(both) < 20 or best is None:
                continue
            a, b = np.array(both).T
            rng = np.random.default_rng([seed, j, KINDS.index(kind)])
            draws = rng.integers(0, a.size, (400, a.size))
            ratios = b[draws].std(axis=1, ddof=1) / a[draws].std(axis=1, ddof=1)
            ratio, se = float(b.std(ddof=1) / a.std(ddof=1)), float(ratios.std(ddof=1))
            row[kind].update({"spread_over_the_weighted_fit": ratio, "its_standard_error": se})
            if found[kind].converged and row[kind]["ended_on"] >= 0.9 and ratio + 2.0 * se < least:
                best, least = kind, ratio + 2.0 * se
        row["smallest_honest_uncertainty"] = best
        out["parameters"][name] = row
    out["rule_of_the_smallest"] = ("the weighted fit, unless the spread of another column over the spread of the weighted fit, on the sets on which both ended, "
                                   "is under 1 by more than 2 of its standard errors, on a column whose searches ended on the records and on 0.9 of the sets")
    return out


def lines(table: dict) -> str:
    """The table of `columns` as text: one block a parameter, one line a column."""
    out = [f"{table['sets']} synthetic sets of the fitted model, errors: {table['errors_of_the_synthetic_sets']}"]
    names = {"fit": "weighted fit", "cells": "combination", "hybrid": "hybrid"}
    for name, row in table["parameters"].items():
        out.append(f"{name}")
        for kind in KINDS:
            c = row[kind]
            text = f"   {names[kind]:13s} {c['estimate']:.6g}   printed {c['printed']:.3g}"
            if "honest" in c:
                text += f"   honest {c['honest']:.3g}   bias {c['bias']:+.2g}   printed intervals hold the value on {c['printed_intervals_that_hold_the_value_put_in']:.2f}"
            if "spread_over_the_weighted_fit" in c:
                text += f"   spread over the weighted fit {c['spread_over_the_weighted_fit']:.2f} +- {c['its_standard_error']:.2f}"
            text += f"   searches ended on {c['ended_on']:.2f}"
            out.append(text + ("   <- the smallest honest uncertainty" if row["smallest_honest_uncertainty"] == kind else ""))
    return "\n".join(out)

