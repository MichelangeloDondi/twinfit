"""The canonical fit: weighted least squares of every record at once, the shared parameters common to all.

The method is the separable form of the Gauss and Newton iteration with the damping of Levenberg and Marquardt. For each
record the columns of its own parameters are taken out of its whitened residual and of the columns of the shared
parameters, the shared step is solved on what is left, and each record then takes its own step. The derivatives of the
model are central finite differences, so a model is any function of its parameters.

    result = fit(problem)                       from the starts of the problem
    result = fit(problem, shared0, own0)        from other starts, for instance the values that a researcher brings

The covariance of the shared parameters is the inverse of the reduced normal matrix at the end: the first order bound.
A fit that did not end says so, and nothing downstream takes it as one that did.

Two things keep a record with almost no signal from stopping the fit of all the others. A problem may give bounds: a
step is then cut at them, and an own parameter that sits at a bound and is asked to go further out is held, a number and
no longer a parameter, until its record asks it back in. And the own columns of a record that the record does not tell
apart, such as the centre of a curve whose area is zero, are left out of its basis. With no bounds and no such record
the search is what it was.

Three things that a fit now says of itself, since a search that ends is not all that a reader needs. `held`: how many
own parameters sit at a bound. `undetermined`: how many own parameters their record does not fix, with or without
bounds: such a parameter stays where the search left it and means nothing. `at_bound`: which shared parameters are
held at a bound at the end, that is, sit at it while the chi square falls only beyond it. Such a parameter has no
bar, its row of the covariance is not a number, and the bars of the other parameters are those of a fit with it
held. The degrees of freedom count the parameters that are free: not an own parameter that is held or that its record
does not fix, and not a shared one that is held at a bound. The linear picture of such a fit (linear.py) still treats
a held shared parameter as free: its tables then say what a move of that parameter could take, which the bound
forbids.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

STEP = 1e-5              # the relative step of a central difference
RANK = 1e-10             # an own column under this share of the largest is not told apart by its record


@dataclass
class Fit:
    shared: np.ndarray
    own: np.ndarray                 # one row a record
    cov: np.ndarray                 # of the shared parameters
    chi2: float
    dof: int
    converged: bool
    iterations: int
    held: int = 0                   # own parameters that sit at a bound of the problem at the end
    undetermined: int = 0           # own parameters that their record does not fix: columns left out of a basis
    at_bound: tuple = ()            # the places of the shared parameters that sit at a bound at the end
    finished_without_bounds: bool = False   # the search with bounds stopped, and this is the end of the search with none, inside the bounds
    ended_by_the_arithmetic: bool = False   # the search stopped short of its test of an end, where a whole step promised a fall of chi square under `ARITHMETIC` of it

    @property
    def sd(self) -> np.ndarray:
        return np.sqrt(np.diag(self.cov))


def derivatives(problem, k: int, shared: np.ndarray, own: np.ndarray) -> tuple:
    """(whitened residual, whitened derivatives of the mean by the shared parameters, by the own parameters) of record k."""
    rec = problem.records[k]
    f0 = problem.predict(k, shared, own)

    def column(vec: np.ndarray, j: int, is_shared: bool) -> np.ndarray:
        h = STEP * max(abs(vec[j]), 1e-3 * (np.max(np.abs(vec)) or 1.0), 1e-9)
        up, down = vec.copy(), vec.copy()
        up[j] += h
        down[j] -= h
        hi = problem.predict(k, up, own) if is_shared else problem.predict(k, shared, up)
        lo = problem.predict(k, down, own) if is_shared else problem.predict(k, shared, down)
        return (hi - lo) / (2.0 * h) / rec.sd
    A = np.column_stack([column(shared, j, True) for j in range(shared.size)]) if shared.size else np.zeros((rec.x.size, 0))
    B = np.column_stack([column(own, j, False) for j in range(own.size)]) if own.size else np.zeros((rec.x.size, 0))
    return (rec.y - f0) / rec.sd, A, B


def cost(problem, shared: np.ndarray, own: np.ndarray) -> float:
    return float(sum(np.sum(problem.whitened(k, shared, own[k]) ** 2) for k in range(len(problem.records))))


def basis(B: np.ndarray) -> np.ndarray:
    """An orthonormal basis of what the columns of B span. A column that is zero, or a combination of the others to
    rounding, adds nothing to it."""
    if not B.shape[1]:
        return np.zeros((B.shape[0], 0))
    Q, R = np.linalg.qr(B)
    norms = np.linalg.norm(B, axis=0)
    if norms.min() > 0.0 and np.all(np.abs(np.diag(R)) > RANK * norms):      # each column against its own size: the units of the parameters do not enter
        return Q
    scaled = B[:, norms > 0.0] / norms[norms > 0.0]
    if not scaled.shape[1]:
        return np.zeros((B.shape[0], 0))
    U, s, _ = np.linalg.svd(scaled, full_matrices=False)
    return U[:, s > RANK * s.max()]


def reduced(problem, shared: np.ndarray, own: np.ndarray, release: bool = False) -> list:
    """For every record: its whitened residual, its shared columns, the columns of its own parameters that are free, an
    orthonormal basis of those, and which own parameters are held. An own parameter is held when it sits at a bound of
    the problem. With `release`, one that its record asks to move back inside is free again: the search uses that."""
    parts = []
    for k in range(len(problem.records)):
        r, A, B = derivatives(problem, k, shared, own[k])
        hold = np.zeros(B.shape[1], dtype=bool)
        if problem.bounded and B.shape[1]:
            lo, hi = problem.own_bounds(k)
            low, high = own[k] <= lo, own[k] >= hi
            if (low | high).any():
                if release:
                    step = np.linalg.lstsq(B, r, rcond=None)[0]
                    hold = (low & (step <= 0.0)) | (high & (step >= 0.0))
                else:
                    hold = low | high
                B = B[:, ~hold]
        parts.append((r, A, B, basis(B), hold))
    return parts


def _asked_of(B: np.ndarray, r: np.ndarray) -> float:
    """The largest step that a record asks of one of its free own parameters, in standard deviations of that parameter."""
    M = B.T @ B
    try:
        Minv = np.linalg.inv(M)
        d = np.diag(Minv)
        if not (np.all(np.isfinite(Minv)) and np.all(d > 0.0)):
            raise np.linalg.LinAlgError
    except np.linalg.LinAlgError:                       # the record does not tell its own parameters apart: ask of those that it fixes
        Minv = np.linalg.pinv(M, rcond=RANK, hermitian=True)
        d = np.diag(Minv)
    ok = d > 0.0
    return float(np.max(np.abs(Minv @ (B.T @ r))[ok] / np.sqrt(d[ok]))) if ok.any() else 0.0


ARITHMETIC = 1e-12          # the share of chi square under which the fall that a whole step promises is set by the arithmetic and not by the records


def _fall_of(B: np.ndarray, r: np.ndarray) -> float:
    """The fall of chi square that a whole step of a record's own parameters promises at first order: r B (B'B)^-1 B' r."""
    if not B.shape[1]:
        return 0.0
    g = B.T @ r
    return float(g @ np.linalg.lstsq(B.T @ B, g, rcond=RANK)[0])


def fall_promised(problem, shared: np.ndarray, own: np.ndarray) -> float:
    """The fall of chi square that a whole step of every parameter promises at first order, from a point of the search:
    of the shared parameters g' H^-1 g, and of each record's own parameters that of its own step. Where it is a share of
    chi square that the arithmetic sets, no step can lower chi square by the records any more."""
    parts = reduced(problem, shared, own, release=True)
    p = shared.size
    H, g = np.zeros((p, p)), np.zeros(p)
    for r, A, B, Q, hold in parts:
        At = A - Q @ (Q.T @ A)
        H += At.T @ At
        g += At.T @ (r - Q @ (Q.T @ r))
    fall = float(g @ np.linalg.lstsq(H, g, rcond=RANK)[0]) if p else 0.0
    for r, A, B, Q, hold in parts:
        fall += _fall_of(B, r)
    return fall


def _inverse(M: np.ndarray) -> np.ndarray:
    """The inverse of a curvature, or its pseudo-inverse where it is exactly singular: two parameters that move the
    curves alike then get a correlation of 1, which says that the records cannot tell them apart (rule 5). On Linux the
    first fits of the toy raised numpy's "Singular matrix" here, where the library of the Mac gave a near inverse (the
    red team of the cloud session, row M213)."""
    try:
        return np.linalg.inv(M)
    except np.linalg.LinAlgError:
        return np.linalg.pinv(M)


def _own_step(B: np.ndarray, rhs: np.ndarray, lam: float) -> np.ndarray:
    M = B.T @ B
    damped = M + lam * np.diag(np.diag(M))
    try:
        step = np.linalg.solve(damped, rhs)
        if np.all(np.isfinite(step)):
            return step
    except np.linalg.LinAlgError:
        pass
    return np.linalg.lstsq(damped, rhs, rcond=RANK)[0]  # a parameter that the record does not fix does not move


def fit(problem, shared0=None, own0=None, max_iter: int = 80, tol: float = 1e-4) -> Fit:
    """The weighted least squares fit of every record at once. The search has ended when the step that first order asks of
    every parameter, shared or own, is under `tol` of its own standard deviation."""
    n_rec = len(problem.records)
    shared = np.array(problem.start_shared if shared0 is None else shared0, dtype=float)
    own = np.array([problem.own_start(k) for k in range(n_rec)] if own0 is None else own0, dtype=float).reshape(n_rec, len(problem.own))
    bounded = problem.bounded
    lo_s, hi_s = problem.bounds_shared
    box = [problem.own_bounds(k) for k in range(n_rec)] if bounded else None
    if bounded:                                                   # a start outside the bounds is brought to them
        shared = np.clip(shared, lo_s, hi_s)
        own = np.array([np.clip(own[k], *box[k]) for k in range(n_rec)])
    chi2 = cost(problem, shared, own)
    if not np.isfinite(chi2):
        raise FloatingPointError("the model does not give a number at the start")
    lam, converged, it = 1e-3, False, 0
    p = shared.size
    for it in range(1, max_iter + 1):
        parts = reduced(problem, shared, own, release=True)
        H, g = np.zeros((p, p)), np.zeros(p)
        for r, A, B, Q, hold in parts:
            At = A - Q @ (Q.T @ A)
            H += At.T @ At
            g += At.T @ (r - Q @ (Q.T @ r))
        dead = ~(np.diag(H) > 0) if p else np.zeros(0, dtype=bool)
        if dead.any() and it == 1:
            raise np.linalg.LinAlgError("the records hold nothing about a shared parameter" + (f", {', '.join(n for n, h in zip(problem.shared, np.diag(H)) if not h > 0)}" if len(problem.shared) == p else "")
                                        + ": at its start a change of it moves no curve. Look at its start and its bounds")
        cov = _inverse(H) if p else np.zeros((0, 0))
        want = cov @ g if p else np.zeros(0)
        sds = np.sqrt(np.diag(cov)) if p else np.zeros(0)
        fixed = _held(shared, g, lo_s, hi_s) if bounded and p else np.zeros(p, dtype=bool)
        fixed = fixed | dead          # a parameter that moves no curve after a step, since another stands at a bound, as the spacing of a second line whose share the first step took to zero: held for this step and read again at the next. The whole fit was refused so, said to move no curve at its start, where a start near its truth was the start that failed (round 18 of the cold tests)
        if fixed.any():
            # A shared parameter at a bound whose chi square falls only beyond the bound is held for this step, and the
            # others are solved without it: left in the solve, its pull is cut at the bound while the others still move
            # as if it had moved, and the search creeps. The test is the slope of the chi square at this point and not
            # the step that the whole solve asks: with parameters that trade, the step of one can point out of its
            # bound while its slope points in, and to hold it then ends the search at a corner that is not the best
            # point (review 41, moderate 1). The set is read again at every step, so a parameter is released when its
            # slope turns.
            free = ~fixed
            want, sds = np.zeros(p), np.ones(p)
            if free.any():
                part = _inverse(H[np.ix_(free, free)])
                want[free], sds[free] = part @ g[free], np.sqrt(np.diag(part))
        asked = float(np.max(np.abs(want) / sds)) if p else 0.0
        for r, A, B, Q, hold in parts:                            # the step that each record asks of its own parameters
            if B.shape[1]:
                asked = max(asked, _asked_of(B, r))
        if asked < tol:
            converged = True
            break
        moved = False
        while lam < 1e12:
            if fixed.any():
                ds, free = np.zeros(p), ~fixed
                if free.any():
                    Hf = H[np.ix_(free, free)]
                    ds[free] = np.linalg.solve(Hf + lam * np.diag(np.diag(Hf)), g[free])
            else:
                ds = np.linalg.solve(H + lam * np.diag(np.diag(H)), g) if p else np.zeros(0)
            trial_shared = np.clip(shared + ds, lo_s, hi_s) if bounded else shared + ds
            moved_by = trial_shared - shared if bounded else ds   # with no bounds the step itself: its digits are those of the search before bounds
            trial_own = own.copy()
            for k, (r, A, B, Q, hold) in enumerate(parts):
                if B.shape[1]:
                    step = _own_step(B, B.T @ (r - A @ moved_by), lam)
                    if bounded:
                        full = np.zeros(own.shape[1])
                        full[~hold] = step
                        trial_own[k] = np.clip(own[k] + full, *box[k])
                    else:
                        trial_own[k] = own[k] + step
            trial = cost(problem, trial_shared, trial_own)
            if np.isfinite(trial) and trial <= chi2:
                shared, own, chi2 = trial_shared, trial_own, trial
                lam = max(lam / 5.0, 1e-12)
                moved = True
                break
            lam *= 4.0
        if not moved:
            break
    by_the_arithmetic = False
    if not converged and np.isfinite(chi2) and fall_promised(problem, shared, own) < ARITHMETIC * max(chi2, 1e-300):
        converged = by_the_arithmetic = True          # the second test of an end (his message 156): on records much cleaner than their noise, the step that the first test asks shrinks with the noise and the arithmetic does not, so a fit of a model that lacks a term stopped at its minimum and said that it had not ended
    parts = reduced(problem, shared, own)
    H, g = np.zeros((p, p)), np.zeros(p)
    for r, A, B, Q, hold in parts:
        At = A - Q @ (Q.T @ A)
        H += At.T @ At
        g += At.T @ (r - Q @ (Q.T @ r))
    fixed = _held(shared, g, lo_s, hi_s) if bounded and p else np.zeros(p, dtype=bool)
    if fixed.any():                                                            # the bars of the others are those of a fit with these held: a parameter that
        cov = np.full((p, p), np.nan)                                          # cannot move takes no part in their scatter, and has no bar of its own
        if (~fixed).any():
            cov[np.ix_(~fixed, ~fixed)] = _inverse(H[np.ix_(~fixed, ~fixed)])
    else:
        cov = _inverse(H) if p else np.zeros((0, 0))
    if p and not np.all(np.diag(H) > 0):          # a parameter that moves no curve where the search ended has no bar: the pseudo-inverse gave it a bar of zero
        gone = ~(np.diag(H) > 0)
        cov[gone, :], cov[:, gone] = np.nan, np.nan
    free = sum(Q.shape[1] for _, _, _, Q, _ in parts)                         # the own parameters that are free and that their record fixes
    held = int(sum(int(hold.sum()) for *_, hold in parts))
    at_bound = tuple(int(j) for j in np.where(fixed)[0])
    out = Fit(shared, own, cov, chi2, int(problem.samples - p - free + int(fixed.sum())), bool(converged), it, held, int(own.size - held - free), at_bound, ended_by_the_arithmetic=by_the_arithmetic)
    finite = bounded and (np.any(np.isfinite(lo_s)) or np.any(np.isfinite(hi_s)) or any(np.any(np.isfinite(box[k][0])) or np.any(np.isfinite(box[k][1])) for k in range(n_rec)))
    if finite and not converged:
        # A search with bounds that did not end is tried once with no bounds, from the same start. If that search ends at a
        # lower chi square with every parameter inside the bounds, its point is the fit, and the result says so in
        # `finished_without_bounds`. This finishes some stalls and does not trace them: the joint step of the shared and
        # the own parameters can point uphill at every damping, with a bound active or with none (review 56). Found on
        # traces of four lines with a kernel in the model, where the bounded search stopped 12 in chi square above that
        # point. Where the best point lies on a bound nothing changes and the search is reported as not ended, as before.
        # The second search calls the model outside the bounds, where a model may refuse: any refusal of the model or of
        # the arithmetic there leaves the first result as it is. With no finite bound the second search would repeat the
        # first, and is not made.
        from .problem import Problem
        loose = Problem(problem.records, problem.model, problem.shared, problem.own, problem.start_shared, problem._start_own, problem.centre, problem.name)
        try:
            other = fit(loose, shared0, own0, max_iter, tol)
        except (ValueError, ArithmeticError):
            return out
        inside = bool(np.all(other.shared >= lo_s) and np.all(other.shared <= hi_s) and all(np.all(other.own[k] >= box[k][0]) and np.all(other.own[k] <= box[k][1]) for k in range(n_rec)))
        if other.converged and inside and other.chi2 < chi2:
            other.finished_without_bounds = True
            return other
    return out


def _held(shared, g, lo, hi) -> np.ndarray:
    """Which shared parameters are held at a bound at this point: those that sit at a bound and whose chi square falls only
    beyond it. `g` is half of minus the slope of the chi square: positive where a larger value of the parameter lowers it."""
    return ((shared <= lo) & (g < 0.0)) | ((shared >= hi) & (g > 0.0))
