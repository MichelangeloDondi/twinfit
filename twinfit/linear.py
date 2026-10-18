"""The linear picture of a problem at a fit: everything that first order needs, with the own parameters taken out.

At the fitted values the whitened residual of each record and the whitened columns of the shared parameters are freed of
the columns of the record's own parameters. What is left is one design matrix for the shared parameters alone, to which
the ledger applies, and a rule that turns any discrepancy of the mean, given record by record in the units of the values,
into the same picture. No noise is drawn.
"""
from __future__ import annotations

import numpy as np

from .fit import reduced
from .ledger import Ledger


class Linear:
    def __init__(self, problem, result):
        self.problem, self.fit = problem, result
        self.parts, rows, res, self.edges, self.shared_columns = [], [], [], [0], []
        for k, (r, A, B, Q, hold) in enumerate(reduced(problem, result.shared, result.own)):
            rows.append(A - Q @ (Q.T @ A))
            res.append(r - Q @ (Q.T @ r))
            self.parts.append((problem.records[k].sd, Q))
            self.shared_columns.append(A * problem.records[k].sd[:, None])        # the derivatives of the mean by each shared parameter, in the units of the values
            self.edges.append(self.edges[-1] + r.size)
        self.J = np.vstack(rows)
        self.r = np.concatenate(res)
        self.ledger = Ledger(self.J, problem.shared)
        self.n = self.r.size
        self.dof = self.n - len(problem.shared) - sum(Q.shape[1] for _, Q in self.parts) + len(getattr(result, "at_bound", ()))     # a shared parameter held at a bound takes none

    def rows_of(self, k: int) -> slice:
        return slice(self.edges[k], self.edges[k + 1])

    def centre_of(self, k: int) -> float:
        """The place about which the curve of record k is read: the centre that the problem gives, or the middle of its abscissa."""
        c = self.problem.centre_of(k, self.fit.shared, self.fit.own[k])
        x = self.problem.records[k].x
        return float(0.5 * (x[0] + x[-1])) if c is None else c

    def width_of(self, k: int) -> float:
        """The full width of the fitted curve of record k at half of its height above its median value: the unit in which
        slices and displacements are counted when none is given. Where the wings of a curve fill its window it is narrower
        than the width at half the peak: 19.8 against 24 on Lorentzian bands (round 7 of the cold tests). It said above its
        lowest value, which the code does not do."""
        x = self.problem.records[k].x
        f = self.problem.predict(k, self.fit.shared, self.fit.own[k])
        up = np.abs(f - np.median(f))
        inside = x[up >= 0.5 * up.max()]
        return float(max(inside.max() - inside.min(), 2.0 * np.median(np.diff(x))))

    def reduce(self, per_record) -> np.ndarray:
        """A discrepancy of the mean, one array a record in the units of the values (None for a record that has none),
        whitened and freed of the own parameters of each record."""
        out = np.zeros(self.n)
        for k, d in enumerate(per_record):
            if d is None:
                continue
            sd, Q = self.parts[k]
            w = np.asarray(d, dtype=float) / sd
            out[self.rows_of(k)] = w - Q @ (Q.T @ w)
        return out

    def norm2(self, per_record) -> float:
        """The squared norm of a discrepancy in units of the noise, before anything is taken out."""
        return float(sum(np.sum((np.asarray(d, dtype=float) / self.parts[k][0]) ** 2) for k, d in enumerate(per_record) if d is not None))

    def absorbed(self, per_record) -> np.ndarray:
        """The first order move of each shared parameter under a discrepancy: the bias that it causes."""
        return self.ledger.shift(self.reduce(per_record))

    def visible(self, per_record) -> np.ndarray:
        """What the residuals keep of a discrepancy: all that a statistic of these records can read of it."""
        return self.ledger.visible(self.reduce(per_record))

    def visibility(self, per_record) -> float:
        """The share of the norm of a discrepancy that the residuals keep, between 0 (absorbed whole) and 1."""
        total = self.norm2(per_record)
        if total == 0.0:
            raise ValueError("a discrepancy of zero has no visibility")
        v = self.visible(per_record)
        return float(np.sqrt(max(float(v @ v), 0.0) / total))
