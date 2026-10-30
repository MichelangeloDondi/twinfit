"""The package's own list of generic terms, its dictionary: second components, pedestals, tails, widths beside the line,
the dispersive partner, a block apart, the levels, the periodic terms and a parameter along a covariate, each as the
shape that it adds at a fit. Read only where the researcher asks it with `--dictionary` (the owner's messages 337, 338,
341 and 343: "the dictionary should be present, but used only if specifically required by the user with a flag"), and
never in a blind run: the package is for the term that nobody suspects, and with no flag it reads the residuals alone
and the researcher's own list, empty at the start (342).
"""
from __future__ import annotations

import dataclasses

import numpy as np

from .candidates import ALONG, ALONG_SQUARE, Candidate, covariates_of, from_mechanism
from .identify import stretched


def candidates(lin, periods=None, displacements=(0.35, 0.5, 0.71, 1.0, 1.41, 2.0, 2.83, 4.0), stretches=(2.0, 4.0, 8.0), block: str | None = None,
               covariates=None, widths=(0.03, 0.1, 0.3)) -> list:
    """The generic candidates for a problem at a fit. Distances are counted in the median full width of the fitted curves.

    Every key of the design whose value is a number with more than two values is read as a covariate, an index of
    repeats too: `covariates_of(problem)` names them, `covariates=[...]` keeps the ones that are meant, and the table
    of `propose` says which were read. The abscissa of every record rises: a record refuses one that does not."""
    pr, n_rec = lin.problem, len(lin.problem.records)
    centre = [lin.centre_of(k) for k in range(n_rec)]
    unit = float(np.median([lin.width_of(k) for k in range(n_rec)]))
    span = float(np.median([rec.x[-1] - rec.x[0] for rec in pr.records]))
    step = float(np.median([np.median(np.diff(rec.x)) for rec in pr.records]))
    here = lambda k, x=None: pr.predict(k, lin.fit.shared, lin.fit.own[k], x=x)          # noqa: E731
    out = []
    for power, label in ((1, "a slope of the level"), (2, "a curvature of the level"), (3, "a cubic term of the level")):          # the slope first: a model with a constant baseline in each record lacks it before any other, and a drift of the eclipses of a cold test could only be named a cubic (row M169)
        out.append(Candidate(label, "baseline", [[((rec.x - centre[k]) / span) ** power for k, rec in enumerate(pr.records)]], "even" if power == 2 else "odd"))
        if n_rec > 1:          # and its own in each record, as a background that each record carries at its own size: the bench of the toy plants one that grows with the height, which the shared one did not name
            out.append(Candidate(f"{label}, its own in each record", "baseline",
                                 [[((rec.x - centre[k]) / span) ** power if k == j else None for k, rec in enumerate(pr.records)] for j in range(n_rec)], "even" if power == 2 else "odd"))
    if periods is None:                                           # a ladder of multiples of the width: the multiple is the key, and the period in the name moves with the fit
        ladder = [(p, f"a periodic term of period {m:g} widths") for m, p in zip((0.5, 1.0, 2.0, 4.0, 8.0), unit * np.array([0.5, 1.0, 2.0, 4.0, 8.0])) if 4.0 * step <= p <= span]
    else:
        ladder = [(P, "") for P in periods]                       # periods that the caller gives are what they are in every fit
    for P, key in ladder:
        out.append(Candidate(f"a periodic term of period {P:.3g}", "periodic",
                             [[np.cos(2.0 * np.pi * (rec.x - centre[k]) / P) for k, rec in enumerate(pr.records)],
                              [np.sin(2.0 * np.pi * (rec.x - centre[k]) / P) for k, rec in enumerate(pr.records)]], "any", key=key))
    for m in displacements:
        for side, sign in (("above", 1.0), ("below", -1.0)):
            out.append(Candidate(f"a second component {m:g} widths {side} the centre", "second component",
                                 [[here(k, rec.x - sign * m * unit) for k, rec in enumerate(pr.records)]], "any"))
        out.append(Candidate(f"a second component {m:g} widths away on both sides", "second component",
                             [[here(k, rec.x - m * unit) + here(k, rec.x + m * unit) for k, rec in enumerate(pr.records)]], "even"))
    for kappa in stretches:
        wide = [here(k, centre[k] + (rec.x - centre[k]) / kappa) for k, rec in enumerate(pr.records)]
        out.append(Candidate(f"a pedestal {kappa:g} times wider than the curve", "pedestal", [wide], "even"))
        for side, sign in (("upper", 1.0), ("lower", -1.0)):
            out.append(Candidate(f"a tail on the {side} side, {kappa:g} times wider", "tail",
                                 [[w * (sign * (rec.x - centre[k]) > 0) for k, (rec, w) in enumerate(zip(pr.records, wide))]], "odd"))
    for key in covariates_of(pr, covariates):
        z = np.array([rec.design[key] for rec in pr.records], dtype=float)
        zc = (z - z.mean()) / z.std()
        for j, name in enumerate(pr.shared):
            out.append(Candidate(f"{name}{ALONG}{key}", "grows with a covariate",
                                 [[zc[k] * lin.shared_columns[k][:, j] for k in range(n_rec)]], "any"))
            sq = zc ** 2 - np.mean(zc ** 2)
            out.append(Candidate(f"{name}{ALONG_SQUARE}{key}", "grows with a covariate",
                                 [[sq[k] * lin.shared_columns[k][:, j] for k in range(n_rec)]], "any"))
    even_steps = all(np.allclose(np.diff(rec.x), rec.x[1] - rec.x[0], rtol=1e-6, atol=0.0) for rec in pr.records)
    rising = all(np.all(np.diff(rec.x) > 0.0) for rec in pr.records)          # the widths read any abscissa that rises: on steps that are not equal their kernels go through a grid of equal steps (his messages 294 to 296)
    if rising:          # the width families on a ladder of sizes in widths of the curve (the owner's message 262): at a third of the width a Gaussian reshapes the core past its first order, and a pedestal led it
        from .mechanisms import exponential_width, gaussian_width, lorentzian_width
        for m in widths:
            for kind, width in (("Gaussian", gaussian_width), ("Lorentzian", lorentzian_width), ("exponential", exponential_width)):
                mechanism = width(name=f"a{'n' if kind == 'exponential' else ''} {kind} width of {m:g} widths beside the line")
                if mechanism.parameter in pr.shared:          # a model that holds this kernel already: a second of its kind is the same kernel wider, and its parameter would collide. In 0.1.0rc2 the first Gaussian width kept stopped the next run of the loop (his message 280)
                    continue
                out.append(dataclasses.replace(from_mechanism(lin, mechanism, m * unit), sign=1))
    if even_steps:          # the dispersive partner of each curve about its baseline (the owner's message 266): an odd residual of this shape is a phase, which the asymmetry windows see and do not name
        from scipy.signal import hilbert
        out.append(Candidate("a dispersive admixture of the curve", "dispersive", [[np.imag(hilbert(here(k) - np.median(here(k)))) for k in range(n_rec)]], "odd"))
    if block is not None:
        for value, members in pr.blocks(block).items():
            out.append(Candidate(f"the curves of the block {value} are wider", "one block apart",
                                 [[stretched(lin, k, 0.01) if k in members else None for k in range(n_rec)]], "even"))
    return out
