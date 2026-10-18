"""Readings of the residuals by position: pairs of slices about the centre of each curve, read as candidates.

A pair of slices is the set of samples whose distance from the centre of their curve lies between two edges, on both
sides. Its even reading is a discrepancy of one standard deviation of every sample of the pair. Its odd reading is the same
with the sign of the side. Each is a candidate of one column, so the table of `propose` gives it a score that is a standard
chi square with one degree of freedom under a right model, at first order, whatever the fit absorbed, and it gives the bias
that a misfit of that shape would mean. The edges are a geometric progression of the width of the curve by default: 0 to
1/4, 1/4 to 1/2, 1/2 to 1, 1 to 2, 2 to 4 and 4 to 8 widths.
"""
from __future__ import annotations

import numpy as np

from .candidates import Candidate

GEOMETRIC = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)


def slices(lin, edges=None, by: str | None = None, odd: bool = True) -> list:
    """The pairs of slices as candidates, for all the records together or for each block of a covariate of the design."""
    pr = lin.problem
    n_rec = len(pr.records)
    unit = float(np.median([lin.width_of(k) for k in range(n_rec)]))
    widths = list(GEOMETRIC) if edges is None else None          # the edges as multiples of the width, which no fit changes, when they are the package's own
    edges = [e * unit for e in GEOMETRIC] if edges is None else list(edges)
    centre = [lin.centre_of(k) for k in range(n_rec)]
    groups = {"all the records": list(range(n_rec))} if by is None else {f"{by} {v}": m for v, m in pr.blocks(by).items()}
    out = []
    for label, members in groups.items():
        for i, (lo, hi) in enumerate(zip(edges, edges[1:])):
            key = f"{label}: the samples {widths[i]:g} to {widths[i + 1]:g} widths from the centre" if widths else ""
            even, lean = [], []
            for k, rec in enumerate(pr.records):
                if k not in members:
                    even.append(None)
                    lean.append(None)
                    continue
                u = rec.x - centre[k]
                inside = (np.abs(u) >= lo) & (np.abs(u) < hi)
                even.append(rec.sd * inside)
                lean.append(rec.sd * inside * np.sign(u))
            if not any(e is not None and e.any() for e in even):
                continue
            out.append(Candidate(f"{label}: the samples {lo:.3g} to {hi:.3g} from the centre, both sides together", "slice", [even], "even", key=key and key + ", both sides together"))
            if odd:
                out.append(Candidate(f"{label}: the samples {lo:.3g} to {hi:.3g} from the centre, one side less the other", "slice", [lean], "odd", key=key and key + ", one side less the other"))
    return out
