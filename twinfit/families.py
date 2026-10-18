"""Families of statistics of a curve, and what each gives by itself.

A statistic of this module is a sum of the values of a record under a weight function. A family is a rule that gives a
record its list of weight functions. At first order every statistic that a curve offers is such a sum:

    a generalised moment             the weight is a kernel of the abscissa: a power, a fractional power of the distance,
                                     an inverse or a logarithmic kernel, a cosine or a sine
    a window along the abscissa      the weight is 1 on a pair of intervals at equal distances from the centre, the two
                                     sides added. With one side less the other it is the difference f(c + u) - f(c - u),
                                     summed over the samples of the interval on one side, the other side read at the
                                     mirrored places by a straight line between its two nearest samples
    a window along the value         the weight is 1 where the curve lies between two levels of its height
    the abscissa against the value   the place where the curve crosses a level moves as the sum of the values on a short
                                     interval about that place
    a derivative or an integral      a moment of the derivative of the record by its abscissa, or of its running integral
    a ratio of two moments           at first order the upper moment less the lower one times the fitted ratio
    a transformation of the value    a moment of the logarithm of the curve is a sum of the values under the kernel over
                                     the fitted curve, baseline and all, and a moment of its reciprocal under the kernel
                                     over the square of the fitted curve
    the window of a fit              every sample within a distance of the centre, each with a weight of its own

    lin = Linear(problem, fit)
    fam = Family(lin, moments(range(9)))
    fam.bounds()       the first order standard deviation of each shared parameter from the family alone, over the fit's
    fam.reading()      what the family reads in the residuals of the fit: a chi square, its degrees of freedom, its level
    fam.shift(d)       the first order move of the family's own estimate under a discrepancy of the mean
    fam.power(d)       the chance that the reading reports that discrepancy, at first order
    fam.estimate()     the estimate that the family gives by itself
    fam.contrast()     that estimate less the fit's, as a chi square: the test of two estimates of which one is efficient

The statistics of a family are used with the weight that is best for them, the inverse of their covariance, so that a
family is judged at its strongest: its estimate is the weighted least squares fit of the model to its statistics. Under
a right model and a noise of known variance no family has a smaller first order bound than the fit of every sample,
which is the last family of the table with no limit of distance. That holds when the family estimates every own
parameter. A family that holds one at the fit borrows it from the fit, and its bound leaves out the scatter of what it
borrows: it can then read under 1. What a family can give is a reading that stands on few
numbers, a bias under a missing term that is not the fit's, and so a contrast with the fit.

The weights are placed where the fit puts each curve: its centre, one unit of width for all the records, its fitted
values. They are held there while the estimate of the family is searched. A weight that is placed by the values of a
noisy record is another statistic, with a bias of its own from the noise: it is formed and matched on synthetic
records, and it is not in this module. The bounds of a problem are not applied in the search of a family.

An own parameter that a family does not see is held where the fit put it: `Family(lin, rule, hold=("centre",))`. Weights
that are even about the centre see next to nothing of the centre, and an estimate that is left to find it wanders:
`fam.seen()` gives, for each own parameter, the largest part of its column that the family sees in any record. What the
error of a held parameter does to the estimate is not in `shift`. It is nothing at first order when the family is blind
to that parameter, which is the case that holding is for: a hold of a parameter that the family sees in some record at
more than HOLD_SEEN of its column is refused.

An own parameter that a family sees by rounding alone, under BLIND of its column in every record, is not estimated by it
and costs it nothing: the derivative of a record does not see a constant that is added to it. Such a parameter stays
where the fit put it, as a held one does, and `fam.blind` names it. Before this was written the rounding of such a
column counted as a direction that the family had to give up, one a record, and the bound of the moments of the
derivative was too large by a factor of 1.9 to 2.5 on the design of four records of the tests and on the design of 40
(found on 2026-10-12 by a rebuild of the bound from the definitions, which the tests now hold).

The choice is made record by record: `fam.free_at` says, for each record, which own parameters the family estimates
there, and `fam.blind` names those that it estimates in no record. A rule that sees a parameter in some records and by
rounding alone in the others estimates it in the first and leaves it at the fit in the second. With one choice for
all the records the rounding counted again in the second kind (review 68). On the design of 40 records the two sides
subtracted see between 5 parts in 100 million and 2 parts in 100 000 of the two heights: BLIND is not a line between
rounding and a small real part. It need not be one: a parameter that is seen under BLIND and left at the fit moves a
statistic by under a millionth of its standard deviation.

The contrast of a family with the fit has its chi square law at first order only. Two estimates from the same records
also differ by terms of the square of the noise, and for a family that is nearly as good as the fit, or one that holds
little of the shared parameters, those terms are larger than the first order difference: the contrast then reports with
nothing missing far more often than its level says. `contrast()` gives a first order level only inside the range where
that law was found to hold, and outside it asks for `draws`: the level is then read on synthetic records of the fit.

A fit that did not end is refused: its residuals are not those of a minimum, and every reading of them would be wrong.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2, ncx2

from .fit import basis, derivatives

GEOMETRIC = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)          # the edges of readings.py, in widths of the curve
TENTHS = tuple(float(v) for v in np.arange(1, 10) / 10.0)
LARGEST = 2000          # a family of more statistics than this is the window of a fit: it is not read as one chi square
RANK = 1e-8
READING = 80.0            # a family whose bound is this many times the fit's or more, in some combination of the shared parameters, is a reading and not an estimator. At 100 the cut fell inside the values that one family takes from set to set of one design, 98 to 170, so that it was searched on one set in 40: at 80 it lies in the gap between that family and every other, which are under 67 or over 160 there (review 69). For one set of records the verdict is made at its fit and follows the level of the noise. A table over several terms takes it once, at the fit without a term (`planted_table`)
HOLD_SEEN = 0.05          # a family may hold at the fit an own parameter of which it sees no more than this part of the column, in any record
BLIND = 1e-6              # an own parameter of which a family sees less than this part of the column in every record is seen by rounding alone: it is not estimated by the family
GAP_FLOOR = 0.05          # the first order law of a contrast is used when the family's variance exceeds the fit's by this part of it in every direction,
BOUND_CEILING = 3.0       # and when its bound is under this many times the fit's. Both are read from one small design and one of 40 records (review 65)


def _span(M: np.ndarray) -> np.ndarray:
    """An orthonormal basis of what the columns of M span, also when M has fewer rows than columns: a family of few
    statistics sees the own parameters of a record through fewer numbers than there are parameters."""
    if M.shape[0] >= M.shape[1]:
        return basis(M)
    norms = np.linalg.norm(M, axis=0)
    if not (norms > 0.0).any():
        return np.zeros((M.shape[0], 0))
    U, sv, _ = np.linalg.svd(M[:, norms > 0.0] / norms[norms > 0.0], full_matrices=False)
    return U[:, sv > RANK * sv.max()]


def moments(orders=(0, 1, 2, 3, 4), kernel: str = "power", scales=(1.0,)):
    """Generalised moments. The abscissa u is counted from the centre in widths of the curve. `kernel`: "power", u to each
    whole order. "distance", the distance to each order, which may be fractional. "inverse", 1 / (u^2 + a^2) for each scale
    a. "logarithm", log(1 + u^2 / a^2) for each scale. "cosine" and "sine", of each order as a frequency in cycles a width."""
    def rule(u, v):
        reach = float(np.abs(u).max()) or 1.0          # a power is taken of the abscissa over its largest distance, so that the columns keep a size
        if kernel == "power":
            return [(u / reach) ** int(p) for p in orders]
        if kernel == "distance":
            return [np.abs(u / reach) ** float(p) for p in orders]
        if kernel == "inverse":
            return [1.0 / (u * u + a * a) for a in scales]
        if kernel == "logarithm":
            return [np.log1p(u * u / (a * a)) for a in scales]
        if kernel in ("cosine", "sine"):
            wave = np.cos if kernel == "cosine" else np.sin
            return [wave(2.0 * np.pi * float(f) * u) for f in orders]
        raise ValueError(f"no kernel is called {kernel!r}: power, distance, inverse, logarithm, cosine or sine")
    return rule


def derivative(orders=range(5), kernel: str = "power", scales=(1.0,)):
    """Moments of the derivative of the record by its abscissa: each kernel summed against the difference of the two
    neighbours of a sample over their distance. Summed by parts it is a moment of the record under the derivative of the
    kernel, with a term at each end: a family of its own only through those ends and through what the differences do
    to the noise."""
    inner = moments(orders, kernel, scales)

    def rule(u, v):
        if np.any(np.diff(u) <= 0.0):
            raise ValueError("a derivative needs the abscissa of a record in increasing order")
        out = []
        for k in inner(u, v):
            w, g = np.zeros(u.size), k[1:-1] / (u[2:] - u[:-2])
            w[2:] += g
            w[:-2] -= g
            out.append(w)
        return out
    return rule


def integral(orders=range(5), kernel: str = "power", scales=(1.0,)):
    """Moments of the running integral of the record from its first sample: each kernel summed against the running sum
    of the values times their spacing. A value then carries the sum of the kernel over every sample after it."""
    inner = moments(orders, kernel, scales)

    def rule(u, v):
        if np.any(np.diff(u) <= 0.0):
            raise ValueError("a running integral needs the abscissa of a record in increasing order")
        step = np.gradient(u)
        return [step * np.cumsum(k[::-1])[::-1] for k in inner(u, v)]
    return rule


def ratios(orders=(1.0, 2.0, 3.0, 4.0), over: float = 0.0, kernel: str = "distance"):
    """Ratios of two moments of the record, each moment of an order over the moment of the order `over`. At first order
    a ratio is the moment of its order less the moment below it times the fitted ratio, so it does not move when the
    record is multiplied by a number: it is blind to the height of a curve that has no baseline."""
    top, bottom = moments(orders, kernel), moments((over,), kernel)

    def rule(u, v, f):
        below = bottom(u, v)[0]
        m0 = float(below @ f)
        if m0 == 0.0:
            raise ValueError("the moment under the ratios is zero on the fitted curve")
        return [k / m0 - float(k @ f) / m0 ** 2 * below for k in top(u, v)]
    rule.curve = True
    return rule


def slices(edges=GEOMETRIC, sides: str = "both", beyond: bool = True, mirror: bool = False):
    """Windows along the abscissa: pairs of intervals at equal distances from the centre, their edges in widths. `sides`:
    "sum", the two sides added. "difference", one side less the other. "both". With `beyond`, a last pair holds every
    sample past the last edge. The difference is always taken at mirrored places. With `mirror` the sum is too: it is
    then f(c + u) + f(c - u) over the samples of one side, which sees nothing of the centre of a symmetric curve, where
    a plain sum over the samples of both sides sees what their uneven placing about the centre lets through."""
    if sides not in ("sum", "difference", "both"):
        raise ValueError("sides is sum, difference or both")

    def mirrored(u, inside, sign):
        """The weights of the sum, over the samples of an interval on the right of the centre, of the value there less
        (sign -1) or plus (sign +1) the value at the mirrored place on the left, which a straight line between its two
        nearest samples gives. A sample whose mirrored place is outside the record is left out. The abscissa is taken
        in increasing order."""
        w = np.zeros(u.size)
        for i in np.nonzero(inside & (u > 0.0))[0]:
            j = int(np.searchsorted(u, -u[i], side="right"))          # a mirrored place that is the first sample itself is inside the record
            if j <= 0 or j >= u.size:
                continue
            part = (-u[i] - u[j - 1]) / (u[j] - u[j - 1])
            w[i] += 1.0
            w[j - 1] += sign * (1.0 - part)
            w[j] += sign * part
        return w

    def rule(u, v):
        a = np.abs(u)
        e = list(edges) + ([np.inf] if beyond else [])
        inside = [(a >= lo) & (a < hi) for lo, hi in zip(e, e[1:])]
        inside = [m for m in inside if m.any()]
        pairs = [m.astype(float) for m in inside]
        if sides == "sum" and not mirror:
            return pairs
        if np.any(np.diff(u) <= 0.0):
            raise ValueError("the two sides at mirrored places need the abscissa of a record in increasing order")
        if mirror:
            pairs = [w for w in (mirrored(u, m, 1.0) for m in inside) if w.any()]
        if sides == "sum":
            return pairs
        lean = [w for w in (mirrored(u, m, -1.0) for m in inside) if w.any()]
        return (pairs if sides == "both" else []) + lean
    return rule


def bands(levels=(0.05,) + TENTHS):
    """Windows along the value: the samples whose fitted value lies between two levels, given as parts of the height of the
    curve above its lowest value, on each side of the centre."""
    def rule(u, v):
        e = [-np.inf] + sorted(levels) + [np.inf]
        out = []
        for lo, hi in zip(e, e[1:]):
            for side in (u >= 0.0, u < 0.0):
                col = ((v >= lo) & (v < hi) & side).astype(float)
                if col.any():
                    out.append(col)
        return out
    return rule


def crossings(levels=TENTHS, half: float = 0.05):
    """The abscissa against the value: for each level and each side, the samples within `half` widths of the place where
    the fitted curve crosses that level, or within a sample and a half when the samples are further apart than that. The
    curve is taken to fall from its top on each side: the first crossing from the centre is the one that is used."""
    def rule(u, v):
        reach = max(half, 1.5 * float(np.median(np.diff(np.sort(u)))))
        out = []
        for side in (1.0, -1.0):
            order = np.argsort(side * u)
            d, w = (side * u)[order], v[order]
            d, w = d[d >= 0.0], w[d >= 0.0]
            for level in levels:
                under = np.nonzero(w < level)[0]
                if not under.size or under[0] == 0:
                    continue
                j = int(under[0])
                at = d[j - 1] + (w[j - 1] - level) * (d[j] - d[j - 1]) / (w[j - 1] - w[j])
                col = (np.abs(side * u - at) < reach).astype(float)
                if col.any():
                    out.append(col)
        return out
    return rule


def transformed(how: str = "logarithm", orders=range(5), above: float = 0.1):
    """Moments of a transformation of the value, at first order: power kernels over the fitted curve for its logarithm,
    over the square of the fitted curve for its reciprocal, on the samples whose fitted value is above a part of the
    height over the lowest fitted value. Below that part the transformation of a noisy value has no first order. The
    curve is the fitted mean as the record holds it, its baseline included: the logarithm that an analyst takes."""
    if how not in ("logarithm", "reciprocal"):
        raise ValueError("the transformation is logarithm or reciprocal")

    def rule(u, v, f):
        if np.any(f[v > above] <= 0.0):
            raise ValueError(f"the {how} of a curve needs fitted values over zero where it is taken")
        reach = float(np.abs(u).max()) or 1.0
        keep = v > above
        safe = np.where(keep, f, 1.0) ** (1 if how == "logarithm" else 2)
        return [keep * (u / reach) ** int(p) / safe for p in orders]
    rule.curve = True          # this rule is given the fitted curve itself as a third argument
    return rule


def window(within: float = np.inf):
    """Every sample within a distance of the centre, in widths, each with a weight of its own: the window of a fit. With no
    limit it is the fit of every sample."""
    def rule(u, v):
        return np.abs(u) <= within
    return rule


def union(*rules):
    """The weights of several rules as one family."""
    def rule(u, v, f):
        out = []
        for one in rules:
            got = one(u, v, f) if getattr(one, "curve", False) else one(u, v)
            if isinstance(got, np.ndarray) and got.dtype == bool:
                raise ValueError("the window of a fit is not joined to other weights: it holds every weight of its samples")
            out += list(got)
        return out
    rule.curve = True
    return rule


def catalogue(centre: str | None = None) -> dict:
    """The families of the package by name, one rule each. Orders of a cosine or a sine are cycles a width, scales of an
    inverse or a logarithmic kernel are widths. `centre`: the name of the own parameter that places a curve, when the
    problem has one. The family of even weights holds it, and is left out when no name is given."""
    low = moments((0, 1))
    out = {
        "moments: powers 0 to 4": moments(range(5)),
        "moments: powers 0 to 8": moments(range(9)),
        "moments: fractional powers of the distance": union(low, moments((0.5, 1.5, 2.0, 2.5, 3.0), "distance")),
        "moments: inverse and logarithmic kernels": union(low, moments(kernel="inverse", scales=(0.5, 1.0, 2.0, 4.0)), moments(kernel="logarithm", scales=(1.0, 4.0))),
        "moments: cosines and sines": union(moments((0,)), moments((0.1, 0.25, 0.5, 1.0, 2.0), "cosine"), moments((0.1, 0.25, 0.5), "sine")),
        "windows along the abscissa: pairs of intervals, the two sides added and subtracted": slices(),
        "the two sides added: pairs of intervals, the centre held": (slices(sides="sum", mirror=True), ("centre",)),
        "the two sides subtracted: pairs of intervals": slices(sides="difference"),
        "windows along the value: bands between levels of the height, each side": bands(),
        "the abscissa against the value: the crossings of nine levels, each side": crossings(),
        "moments of the derivative: powers 0 to 8": derivative(range(9)),
        "moments of the running integral: powers 0 to 8": integral(range(9)),
        "ratios of moments of the distance: orders 0.5 to 4 over order 0": union(moments((0, 1)), ratios((0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0))),
        "the logarithm of the value: power kernels 0 to 8": transformed("logarithm", range(9)),
        "the reciprocal of the value: power kernels 0 to 8": transformed("reciprocal", range(9)),
        "every sample: the fit": window(),
    }
    even = "the two sides added: pairs of intervals, the centre held"
    if centre is None:
        del out[even]
    else:
        out[even] = (out[even][0], (centre,))
    return out


def family(lin, name: str, centre: str | None = None) -> "Family":
    """A family of the catalogue at a fit."""
    entry = catalogue(centre)[name]
    rule, hold = entry if isinstance(entry, tuple) else (entry, ())
    return Family(lin, rule, name, hold=hold)


def table(lin, terms: dict | None = None, centre: str | None = None, level: float = 0.05) -> dict:
    """Every family of the catalogue at a fit, at first order and with no noise drawn. For each: its number of statistics,
    whether it fixes every shared parameter, its bounds over the fit's, the largest part of an own parameter that it sees
    in any record, and for each named term that the model lacks, given as a discrepancy of the mean record by record: the
    bias of the family's own estimate in standard deviations of the fit, the chance that its reading reports the term,
    and the chance that its contrast with the fit does, where the first order law of the contrast holds. A family that
    the package refuses on these records has the reason in place of its row."""
    fit_sd = np.sqrt(np.diag(np.linalg.inv(lin.J.T @ lin.J)))

    def finite(values):
        return [float(v) if np.isfinite(v) else None for v in values]
    out = {}
    for name in catalogue(centre):
        try:
            fam = family(lin, name, centre)
        except ValueError as err:
            out[name] = {"refused": str(err)}
            continue
        b = fam.bounds()
        seen = fam.seen()
        why = None if name == "every sample: the fit" else fam.first_order_holds()
        row = {"statistics": fam.size, "fixes_every_shared_parameter": bool(b["fixed"]), "is_an_estimator": bool(b["estimator"]),
               "standard_deviation_over_the_fit": finite(b["over_the_fit"]) if b["estimator"] else None,
               "largest_standard_deviation_over_the_fit_in_any_direction": finite([b["largest_over_the_fit_in_any_direction"]])[0],          # the figure that the verdict reads, also for a reading
               "largest_part_seen_of_each_own_parameter": seen, "held_at_the_fit": list(fam.hold), "seen_by_rounding_alone_and_not_estimated": list(fam.blind),
               "the_contrast_with_the_fit_has_its_first_order_law": why is None and name != "every sample: the fit", "why_not": why, "terms": {}}
        for label, d in (terms or {}).items():
            row["terms"][label] = {"bias_in_standard_deviations_of_the_fit": finite(fam.shift(d) / fit_sd) if b["estimator"] else None,
                                   "chance_that_the_reading_reports_it": fam.power(d, level),
                                   "chance_that_the_contrast_with_the_fit_reports_it": fam.contrast_power(d, level) if why is None and name != "every sample: the fit" else None}
        out[name] = row
    return out


def planted_table(lin, terms: dict, centre: str | None = None, level: float = 0.05) -> dict:
    """The table of `table` with nothing of first order in its biases. For each named term the records are made again with
    no noise, the fitted mean plus the term, and fitted. Each family is placed on that fit, where an analyst would place
    it, and its own estimate is searched. A bias is an estimate less the values of the fit of the records without the
    term, in standard deviations of that fit. The chance of a reading is that of a chi square whose noncentrality is the
    statistic of the residuals with no noise. A family that is a reading, at READING times the bound of the fit, is not
    searched. That verdict is made once, at the fit of the records without a term, and kept for every term: made at
    each fit that holds a term, a family that lies near the cut was searched under one term and not under the next
    (review 70: two families of sixteen on a design of 40 records, by a part in 160)."""
    from .fit import fit as fit_
    from .linear import Linear
    from .problem import Record
    pr, base = lin.problem, lin.fit
    fit_sd = np.sqrt(np.diag(np.linalg.inv(lin.J.T @ lin.J)))

    def finite(values):
        return [float(v) if np.isfinite(v) else None for v in values]
    searched = {}
    for name in catalogue(centre):
        try:
            searched[name] = bool(family(lin, name, centre).bounds()["estimator"])
        except ValueError:
            searched[name] = False
    out = {}
    for label, d in terms.items():
        recs = [Record(r.x, pr.predict(k, base.shared, base.own[k]) + (0.0 if d[k] is None else np.asarray(d[k], dtype=float)), r.sd, r.design, r.name) for k, r in enumerate(pr.records)]
        planted = pr.with_records(recs)
        got = fit_(planted, base.shared, base.own)
        if not got.converged:
            out[label] = {"the_fit_did_not_end": True}
            continue
        at = Linear(planted, got)
        row = {}
        for name in catalogue(centre):
            try:
                fam = family(at, name, centre)
            except ValueError as err:
                row[name] = {"refused": str(err)}
                continue
            lam, df = fam._statistic(at.r)
            entry = {"chance_that_the_reading_reports_it": float(ncx2.sf(chi2.isf(level, df), df, lam)) if df else float(level),
                     "largest_standard_deviation_over_the_fit_in_any_direction": finite([fam.bounds()["largest_over_the_fit_in_any_direction"]])[0]}          # at the fit that holds the term: the figure that a verdict made at this fit would read
            if name == "every sample: the fit":
                entry["bias_in_standard_deviations_of_the_fit"] = finite((got.shared - base.shared) / fit_sd)
            elif searched[name]:
                try:
                    est = fam.estimate(verdict=True)
                except ValueError:          # the family leaves a direction free at the fit that holds the term
                    entry["its_search_ended"], entry["bias_in_standard_deviations_of_the_fit"] = False, None
                else:
                    entry["its_search_ended"] = bool(est.converged)
                    entry["bias_in_standard_deviations_of_the_fit"] = finite((est.shared - base.shared) / fit_sd) if est.converged else None
            row[name] = entry
        out[label] = row
    return out


def influence(fam: "Family") -> np.ndarray:
    """The first order map from the whitened errors of the records to the error of the family's estimate of the shared
    parameters: one row a shared parameter, one column a sample, the records one after the other."""
    lin = fam.lin
    J, _, parts = fam._picture(lin.fit.shared, lin.fit.own)
    back = np.linalg.pinv(J)
    out, at = [], 0
    for k, q in enumerate(fam.Q):
        Bq = parts[k][3]
        m = Bq.shape[0]
        block = back[:, at:at + m] - (back[:, at:at + m] @ Bq) @ Bq.T
        if q.dtype == bool:
            whole = np.zeros((block.shape[0], q.size))
            whole[:, q] = block
        else:
            whole = block @ q.T
        out.append(whole)
        at += m
    return np.hstack(out)


def free_of(lin, names, terms: dict | None = None, centre: str | None = None, parameter: int = 0) -> dict:
    """The combination of several estimates of one shared parameter that is free of the bias of each named term, at first
    order, and has the least variance among those. The estimates are those of the families that are named, the fit of
    every sample among them if it is named. Their errors are linear in the errors of the records at first order, so their
    covariance is known with no draw, and the weights are those of least variance under two kinds of condition: they add
    to one, and the first order biases of the estimates under each term, weighted, add to zero. With no term the answer
    is the fit, when the fit is among the estimates: no combination has a smaller variance. The freedom from bias is first
    order: at a size where the bias of a term is not linear in its size the combination keeps a part of it."""
    fams = []
    for name in names:
        fam = family(lin, name, centre)
        if not fam.bounds()["estimator"]:
            raise ValueError(f"{name!r} holds next to nothing of the shared parameters: it is a reading, and gives no estimate to combine")
        fams.append(fam)
    rows = np.array([influence(f)[parameter] for f in fams])
    C = rows @ rows.T
    labels = list(terms or {})
    B = np.array([[f.shift(terms[label])[parameter] for label in labels] for f in fams]).reshape(len(fams), len(labels))
    G = np.column_stack([np.ones(len(fams)), B])
    if len(fams) <= len(labels):
        raise ValueError("there are as many conditions as estimates, or more: name more families than terms")
    inner = np.linalg.pinv(C) @ G
    w = inner @ np.linalg.pinv(G.T @ inner)[:, 0]
    fit_sd = float(np.sqrt(np.linalg.inv(lin.J.T @ lin.J)[parameter, parameter]))
    return {"weights": dict(zip(names, (float(v) for v in w))), "standard_deviation_over_the_fit": float(np.sqrt(max(w @ C @ w, 0.0))) / fit_sd,
            "first_order_bias_under_each_term": {label: float(w @ B[:, i]) for i, label in enumerate(labels)},
            "bias_of_each_estimate_in_standard_deviations_of_the_fit": {name: {label: float(B[k, i]) / fit_sd for i, label in enumerate(labels)} for k, name in enumerate(names)}}


def least_variance(C: np.ndarray) -> np.ndarray:
    """The weights that add to one and give the least variance to a weighted sum of estimates whose covariance is C. When
    some weighted sum of the estimates has no variance at all, as two estimates whose errors are opposite, the weights are
    those of that sum."""
    C = np.atleast_2d(np.asarray(C, dtype=float))
    w, V = np.linalg.eigh(C)
    null = V[:, w <= 1e-12 * max(float(w.max()), 0.0)]
    ones = np.ones(C.shape[0])
    still = null @ (null.T @ ones)          # what the directions of no variance hold of the condition that the weights add to one
    if float(ones @ still) > 1e-9:
        return still / float(ones @ still)
    inv = np.linalg.pinv(C)
    return inv.sum(axis=1) / inv.sum()


def together(lin, names, centre: str | None = None, draws: int = 40, check: int = 40, seed: int = 1) -> dict:
    """One estimate of each shared parameter from several, weighted by how they vary together on synthetic records of the
    fit. The estimates are those of the named families, and the fit of every sample if it is named. On `draws` synthetic
    sets of records of the fit, each fitted and read as the records are, the errors of the estimates give their
    covariance, with no first order law assumed: terms of the square of the noise are in it. The weights of a shared
    parameter add to one and give the least variance under that covariance, which is the canonical combination of
    estimates by generalised least squares. The standard deviation that is returned is read on `check` further sets with
    the weights held: the one that the first sets give with their own weights is too small, since the weights were
    chosen on those sets, and it is returned beside it under its own name. The mean error on the further sets is the
    bias that the noise gives the combination, with its standard error: `corrected` is the estimate less that bias.

    What it cannot do. Under a right model and normal errors of the stated size no combination has a smaller variance
    than the fit at first order, so the gain there is none and the weights cost a little. The synthetic records hold
    normal errors of the stated size and the fitted model: a noise of another law, or a term that the model lacks, is
    not in them, and the standard deviation returned is then not that of the estimate. For a named term see `free_of`."""
    from . import twin
    from .fit import fit as fit_
    from .linear import Linear
    names = list(names)
    known = catalogue(centre)
    unknown = [n for n in names if n not in known]
    if unknown:
        raise ValueError(f"no family is called {unknown[0]!r}")
    if len(names) < 2:
        raise ValueError("a combination is of two estimates or more")
    if draws < 4 * len(names):
        raise ValueError(f"the covariance of {len(names)} estimates is not read on fewer than four times as many sets: {draws} are asked")
    pr, ft = lin.problem, lin.fit
    at = np.asarray(ft.shared, dtype=float)

    def estimates(l):
        rows = []
        for name in names:
            if name == "every sample: the fit":
                rows.append(np.asarray(l.fit.shared, dtype=float))
                continue
            fam = family(l, name, centre)
            if not fam.bounds()["estimator"]:
                return name
            got = fam.estimate()
            if not got.converged:
                return None
            rows.append(got.shared)
        return np.array(rows)
    here = estimates(lin)
    if isinstance(here, str):
        raise ValueError(f"{here!r} holds next to nothing of the shared parameters: it is a reading, and gives no estimate to combine")
    if here is None:
        raise ValueError("the search of a named family did not end on the records: there is no estimate to combine")

    def cloud(count: int, seed_: int) -> tuple:
        errors, of_the_fit = [], []
        for replica in twin.replicas(pr, ft, count, level=1.0, seed=seed_):
            base = fit_(replica, ft.shared, ft.own)
            if not base.converged:
                continue
            got = estimates(Linear(replica, base))
            if got is None or isinstance(got, str):
                continue
            errors.append(got - at)
            of_the_fit.append(np.asarray(base.shared, dtype=float) - at)
        return np.array(errors), np.array(of_the_fit)
    first, _ = cloud(draws, seed)
    if first.shape[0] < 4 * len(names):
        raise ValueError(f"only {first.shape[0]} of the {draws} synthetic sets ended: too few for the covariance of {len(names)} estimates")
    weights = np.array([least_variance(np.cov(first[:, :, j].T)) for j in range(at.size)])          # one row a shared parameter
    on_the_first = np.array([float(np.sqrt(max(w @ np.cov(first[:, :, j].T) @ w, 0.0))) for j, w in enumerate(weights)])
    second, fit_errors = cloud(check, seed + 1)
    if second.shape[0] < 8:
        raise ValueError(f"only {second.shape[0]} of the {check} further sets ended: too few to read a standard deviation")
    errors = np.einsum("jn,dnj->dj", weights, second)
    sd, bias = errors.std(axis=0, ddof=1), errors.mean(axis=0)
    value = np.einsum("jn,nj->j", weights, here)
    return {"estimate": value, "corrected": value - bias, "standard_deviation": sd, "bias_on_synthetic_records": bias, "its_standard_error": sd / np.sqrt(errors.shape[0]),
            "standard_deviation_of_the_fit_on_the_same_sets": fit_errors.std(axis=0, ddof=1), "over_the_fit": sd / fit_errors.std(axis=0, ddof=1),
            "standard_deviation_on_the_sets_of_the_weights": on_the_first, "weights": {name: [float(v) for v in weights[:, k]] for k, name in enumerate(names)},
            "sets_for_the_weights": int(first.shape[0]), "further_sets": int(errors.shape[0])}


OBJECTS = ("the residuals of a fit", "the curve itself", "its two sides", "the abscissa against the value", "its reciprocal and its logarithm",
           "its derivative and its running integral", "its cosine and sine transform")
TECHNIQUES = ("windows along the abscissa", "windows along the value, and the grid", "the cost of a residual", "the two sides", "moments")
RUN, FIRST, SPAN, NONE = "run with noise", "computed with no noise", "in the span of another cell", "not built"


def cells() -> list:
    """The 35 cells of objects and techniques, each with the function of the package that holds it and how far it was read,
    or the reason why it is empty. `state` is one of four words: run with noise on synthetic records, computed with no
    noise, in the span of another cell (its statistics are combinations of those of the cell that is named), not built.
    The arrangement in seven objects and five techniques is one way to lay out the list, not a law."""
    held = {
        (0, 0): (RUN, "readings.slices", "pairs of slices of the residuals, the two sides together"),
        (0, 1): (RUN, "plane.plane_fit, losses.grid", "the plane of position and residual value"),
        (0, 2): (RUN, "losses.crossed", "the grid of the residuals of the fit under each cost"),
        (0, 3): (RUN, "readings.slices", "pairs of slices of the residuals, one side less the other"),
        (0, 4): (RUN, "families.Family.reading", "the residuals summed under the kernels of a family of moments"),
        (1, 0): (RUN, "families.slices", "pairs of intervals of the curve"),
        (1, 1): (RUN, "families.bands, matched.bands", "bands between levels of the value, each side, placed by the fitted curve or by the values of the record"),
        (1, 2): (RUN, "losses.chi_fit", "the fit of the curve under a cost that is not the square"),
        (1, 3): (SPAN, "families.slices", "the cell of the two sides under windows along the abscissa"),
        (1, 4): (RUN, "families.moments, families.ratios", "moments under powers, fractional powers, inverse and logarithmic kernels. The ratios were computed with no noise only"),
        (2, 0): (RUN, "families.slices", "the two sides added and subtracted at mirrored places, over pairs of intervals"),
        (2, 1): (SPAN, "families.bands", "each side has its own bands, so the sum and the difference of the two sides in a band are combinations of them"),
        (2, 3): (SPAN, "families.slices", "the two sides are the object and the technique at once: the cell under windows along the abscissa"),
        (2, 4): (SPAN, "families.moments", "the odd and the even orders of the moments of the curve"),
        (3, 1): (RUN, "families.crossings, matched.crossings", "the places where the curve crosses nine levels of its height, each side, on the fitted curve or on the values of the record"),
        (4, 4): (RUN, "families.transformed, matched.logarithm", "power kernels over the fitted curve and over its square, and the logarithm of the values of the record under kernels"),
        (5, 4): (RUN, "families.derivative, families.integral", "moments of the derivative and of the running integral"),          # run with noise on 100 sets a level since 2026-10-12, and the map still said computed with no noise (review 72)
        (6, 4): (RUN, "families.moments", "cosine and sine kernels"),
    }
    why = {2: "a cost is a property of a fit of the curve or of its residuals, and this object is not fitted by itself", 3: "the two sides of this object were not built",
           0: "this object has no windows along the abscissa of its own in the package", 1: "this object has no windows along the value of its own in the package", 4: "moments of this object were not built"}
    out = []
    for i, obj in enumerate(OBJECTS):
        for j, tech in enumerate(TECHNIQUES):
            state, where, what = held.get((i, j), (NONE, "", why[j]))
            out.append({"object": obj, "technique": tech, "state": state, "in_the_package": where, "what": what, "row": i, "column": j})
    return out


@dataclass
class FamilyFit:
    shared: np.ndarray
    own: np.ndarray                 # one row a record
    cov: np.ndarray                 # of the shared parameters, at first order
    cost: float                     # the sum of the squares of the statistics, each in its standard deviation
    dof: int
    converged: bool
    iterations: int

    @property
    def sd(self) -> np.ndarray:
        return np.sqrt(np.diag(self.cov))


class Family:
    def __init__(self, lin, rule, name: str = "", hold=()):
        self.lin, self.rule, self.name = lin, rule, name
        pr, ft = lin.problem, lin.fit
        if not ft.converged:
            raise ValueError("the fit did not end: a family is read at a minimum of the fit, and this is not one")
        unknown = [h for h in hold if h not in pr.own]
        if unknown:
            raise ValueError(f"no own parameter is called {unknown[0]!r}: {', '.join(pr.own)}")
        self.hold, self.blind = tuple(hold), ()
        self.unit = float(np.median([lin.width_of(k) for k in range(len(pr.records))]))
        self.Q = []          # for each record: an orthonormal basis of its whitened weights, or the mask of the samples of a window
        for k, rec in enumerate(pr.records):
            u = (rec.x - lin.centre_of(k)) / self.unit
            f = pr.predict(k, ft.shared, ft.own[k])
            v = (f - f.min()) / (float(f.max() - f.min()) or 1.0)
            got = rule(u, v, f) if getattr(rule, "curve", False) else rule(u, v)
            if isinstance(got, np.ndarray) and got.dtype == bool:
                self.Q.append(got)
            else:
                W = np.column_stack([np.asarray(c, dtype=float) for c in got]) if len(got) else np.zeros((rec.x.size, 0))
                self.Q.append(_span(W * rec.sd[:, None]))          # a sum of values under a weight w has the whitened weight w times the standard deviation
        self.size = sum(int(q.sum()) if q.dtype == bool else q.shape[1] for q in self.Q)
        self._form = None          # the covariance of the statistics of the residuals, made once
        parts, _ = self._parts()
        seen = dict(zip(pr.own, (float(v) for v in parts.max(axis=0))))
        wide = [h for h in hold if seen[h] > HOLD_SEEN]
        self.free_at = (parts >= BLIND) & np.array([n not in hold for n in pr.own], dtype=bool)          # record by record: the own parameters that the family estimates there
        self.free = self.free_at.any(axis=0)          # estimated in some record
        self.blind = tuple(n for j, n in enumerate(pr.own) if n not in hold and not self.free[j])          # seen by rounding alone in every record: not a parameter of this family
        if wide:
            raise ValueError(f"this family sees {seen[wide[0]]:.2f} of the column of {wide[0]!r} in some record: a parameter that a family sees is estimated by it, not held (over {HOLD_SEEN} of the column)")

    @staticmethod
    def _see(q, M):
        """What the family sees of the rows M of one record."""
        return M[q] if q.dtype == bool else q.T @ M

    def _picture(self, shared, own) -> tuple:
        """At given values: the columns of the shared parameters and the residual as the family sees them, freed of the
        own parameters that the family estimates, and for each record what it sees before that."""
        pr = self.lin.problem
        rows, res, parts = [], [], []
        for k in range(len(pr.records)):
            r, A, B = derivatives(pr, k, np.asarray(shared, dtype=float), np.asarray(own[k], dtype=float))
            rF, AF, BF = self._see(self.Q[k], r), self._see(self.Q[k], A), self._see(self.Q[k], B[:, self.free_at[k]])
            Bq = _span(BF)
            rows.append(AF - Bq @ (Bq.T @ AF))
            res.append(rF - Bq @ (Bq.T @ rF))
            parts.append((rF, AF, BF, Bq))
        return np.vstack(rows), np.concatenate(res), parts

    def _parts(self) -> tuple:
        """For each record and each own parameter: the part of its whitened column that the family sees, between 0 and 1,
        and whether the column is other than zero in that record."""
        pr, ft = self.lin.problem, self.lin.fit
        parts, alive = np.zeros((len(pr.records), len(pr.own))), np.zeros((len(pr.records), len(pr.own)), dtype=bool)
        for k in range(len(pr.records)):
            _, _, B = derivatives(pr, k, ft.shared, ft.own[k])
            norms = np.linalg.norm(B, axis=0)
            alive[k] = norms > 0.0
            parts[k] = np.linalg.norm(np.atleast_2d(self._see(self.Q[k], B)), axis=0) / np.where(alive[k], norms, 1.0)
        return parts, alive

    def seen(self) -> dict:
        """For each own parameter, the largest part of its whitened column that the family sees in any record, between 0
        and 1. A parameter that is seen at a few parts in a hundred in every record is one to hold."""
        return dict(zip(self.lin.problem.own, (float(v) for v in self._parts()[0].max(axis=0))))

    def bounds(self) -> dict:
        """The first order standard deviation of each shared parameter from the family alone, that deviation over the
        fit's, and the correlations. `fixed` is False when the statistics of the family leave a direction of the shared
        parameters free: the family is then a reading and not an estimator, and its deviations are infinite."""
        lin = self.lin
        s = len(lin.problem.shared)
        if s == 0:          # a model with no shared parameter: the family has nothing to estimate, and is a reading
            return {"fixed": True, "estimator": False, "statistics": self.size, "sd": np.zeros(0), "over_the_fit": np.zeros(0), "largest_over_the_fit_in_any_direction": float("nan"), "correlation": np.zeros((0, 0))}
        J, _, _ = self._picture(lin.fit.shared, lin.fit.own)
        scale = np.linalg.norm(lin.J, axis=0)
        sv = np.linalg.svd(J / scale, compute_uv=False) if J.shape[0] >= s else np.zeros(s)
        fit_sd = np.sqrt(np.diag(np.linalg.inv(lin.J.T @ lin.J)))
        if sv.size < s or sv.min() <= RANK * max(sv.max(), 1.0):
            return {"fixed": False, "estimator": False, "statistics": self.size, "sd": np.full(s, np.inf), "over_the_fit": np.full(s, np.inf), "largest_over_the_fit_in_any_direction": float("inf"),
                    "correlation": np.full((s, s), np.nan)}
        cov = np.linalg.inv(J.T @ J)
        sd = np.sqrt(np.diag(cov))
        over = sd / fit_sd
        low = np.linalg.cholesky(np.linalg.inv(lin.J.T @ lin.J))          # the largest ratio of the two standard deviations over every combination of the shared parameters: it does not change when they are mixed or written in other units, as the ratio of one parameter does (review 67, moderate 2)
        half = np.linalg.solve(low, np.linalg.solve(low, cov).T).T
        worst = float(np.sqrt(np.linalg.eigvalsh(0.5 * (half + half.T)).max()))
        # a family at READING times the bound of the fit in some direction holds next to nothing there: what it holds is rounding or the placing of the samples, and it is a reading and not an estimator
        return {"fixed": True, "estimator": bool(worst < READING), "statistics": self.size, "sd": sd, "over_the_fit": over, "largest_over_the_fit_in_any_direction": worst, "correlation": cov / np.outer(sd, sd)}

    def _statistic(self, r: np.ndarray) -> tuple:
        """The squared length of the statistics of a whitened residual, in the covariance that the statistics of the
        residuals of the fit have under a right model and at first order: what is left when the directions of the model
        are taken out of the family. With it, the degrees of freedom: the directions of the family that the model does
        not hold."""
        lin = self.lin
        if all(q.dtype == bool and q.all() for q in self.Q):
            return float(r @ r), int(lin.dof)          # the family of every sample: the chi square of the fit
        if self.size > LARGEST:
            raise ValueError(f"a family of {self.size} statistics is not read as one chi square: read the contrast of its estimate with the fit's")
        if self._form is None:
            U = _span(lin.J)
            own_blocks = [self._see(q, lin.parts[k][1]) for k, q in enumerate(self.Q)]
            E = np.vstack([self._see(q, U[lin.rows_of(k)]) for k, q in enumerate(self.Q)])
            G = np.eye(E.shape[0]) - E @ E.T
            at = 0
            for D in own_blocks:
                m = D.shape[0]
                G[at:at + m, at:at + m] -= D @ D.T
                at += m
            w, V = np.linalg.eigh(G)
            keep = w > 1e-9
            self._form = (V[:, keep], w[keep])
        V, w = self._form
        t = np.concatenate([self._see(q, r[lin.rows_of(k)]) for k, q in enumerate(self.Q)])
        return float(np.sum((V.T @ t) ** 2 / w)), int(w.size)

    def reading(self) -> dict:
        """What the family reads in the residuals of the fit: a chi square with as many degrees of freedom as the family
        has directions that the model does not hold, and its level."""
        stat, df = self._statistic(self.lin.r)
        return {"chi_square": stat, "degrees_of_freedom": df, "level": float(chi2.sf(stat, df)) if df else 1.0, "statistics": self.size}

    def power(self, per_record, level: float = 0.05) -> float:
        """The chance that the reading reports, at the level, a discrepancy of the mean that the model lacks. First order,
        no noise drawn: the reading of what the residuals keep of the discrepancy is the noncentrality of its law."""
        lam, df = self._statistic(self.lin.visible(per_record))
        return float(ncx2.sf(chi2.isf(level, df), df, lam)) if df else float(level)

    def _gap(self) -> tuple:
        """The covariance of the family's estimate less the fit's, under a right model and a noise of known variance: the
        family's own less the fit's, since the fit is the estimate of least variance. With it, its directions that are not
        null and their sizes, in units of the variance of the fit."""
        lin = self.lin
        J, _, _ = self._picture(lin.fit.shared, lin.fit.own)
        of_the_fit = np.linalg.inv(lin.J.T @ lin.J)
        C = np.linalg.inv(J.T @ J) - of_the_fit
        L = np.linalg.cholesky(of_the_fit)          # the sizes are the generalised eigenvalues of the gap against the covariance of the fit: they do not change when the shared parameters are mixed (review 66)
        half = np.linalg.solve(L, np.linalg.solve(L, C).T).T
        w, V = np.linalg.eigh(0.5 * (half + half.T))
        keep = w > 1e-9
        return np.linalg.solve(L.T, V[:, keep]), w[keep]

    def first_order_holds(self) -> str | None:
        """None when the contrast of this family may be given its first order level, else the reason why not. The law of
        the difference of two estimates is first order: it was found to keep its level for a family whose variance
        exceeds the fit's by at least GAP_FLOOR of it in every direction and whose bound is under BOUND_CEILING times
        the fit's, and to report several times too often outside that range. The range is read from two designs at the
        noise that their records state, and it does not know the level of the noise: the terms of the square of the noise
        grow with it. At ten times that noise two families inside the range reported on 10 of 28 sets with nothing
        missing (review 67). At any other level of noise than the one of those designs the level is read with `draws`."""
        b = self.bounds()
        if not b["estimator"]:
            return "this family holds next to nothing of the shared parameters: it is a reading, and has no estimate to set against the fit's"
        if self.hold:
            return "this family holds an own parameter at the fit: its estimate borrows from the fit, and the difference of the two has not the law that this module computes"
        if b["largest_over_the_fit_in_any_direction"] > BOUND_CEILING:
            return f"this family holds little of the shared parameters, {b['largest_over_the_fit_in_any_direction']:.1f} times the bound of the fit in some direction: its estimate is not in the range of first order"
        _, w = self._gap()
        if w.size < len(self.lin.problem.shared) or float(w.min()) < GAP_FLOOR:
            return (f"this family is nearly as good as the fit in some direction, {float(w.min()) if w.size else 0.0:.3f} of its variance apart: terms of the square of the noise are "
                    "larger than the first order difference")
        return None

    def _contrast_of(self, got: "FamilyFit") -> tuple:
        V, w = self._gap()
        z = V.T @ (got.shared - self.lin.fit.shared)
        return float(np.sum(z * z / w)), int(w.size)

    def contrast(self, estimate: "FamilyFit | None" = None, draws: int = 0, seed: int = 1) -> dict:
        """The estimate of the family less the fit's, read in the covariance that the difference has under a right model,
        which is the canonical specification test of two estimates of which one is efficient. A family whose bias under
        a term is not the fit's sees that term here, with the few numbers of the shared parameters. A level to trust is
        read with `draws`: on that many synthetic sets of records of the fit, each fitted and read as the records are,
        it is one more than the count of sets with a statistic as large, over one more than the sets whose searches
        ended. With 19 draws it is 0.05 at its smallest, so a test at 5 in 100 reads a level at or under 0.05. With no
        draws the level is the first order one, a chi square, and it comes with a flag: `first_order_verified` is False
        outside the range where that law was found to hold on two designs, and the range itself is not a rule: it was
        read at one level of noise and says nothing of a higher one. A family
        that is a reading and an estimate whose search did not end are refused."""
        if not self.bounds()["estimator"]:
            raise ValueError("this family holds next to nothing of the shared parameters, or leaves a direction of them free: it has no estimate to set against the fit's")
        got = self.estimate() if estimate is None else estimate
        if not got.converged:
            raise ValueError("the search of the family did not end: its estimate is not set against the fit's")
        stat, df = self._contrast_of(got)
        out = {"chi_square": stat, "degrees_of_freedom": df, "difference": got.shared - self.lin.fit.shared, "ended": True}
        if not draws:
            why = self.first_order_holds()
            return dict(out, level=float(chi2.sf(stat, df)) if df else 1.0, level_by="first order", first_order_verified=why is None, why_not=why,
                        verified_at="the noise that the records of two designs state, and no higher one")
        from . import twin
        from .fit import fit as fit_
        from .linear import Linear
        over, ended = 0, 0
        for replica in twin.replicas(self.lin.problem, self.lin.fit, draws, level=1.0, seed=seed):
            base = fit_(replica, self.lin.fit.shared, self.lin.fit.own)
            if not base.converged:
                continue
            again = Family(Linear(replica, base), self.rule, self.name, hold=self.hold)
            est = again.estimate()
            if not est.converged:
                continue
            ended += 1
            over += again._contrast_of(est)[0] >= stat
        return dict(out, level=(1.0 + over) / (1.0 + ended), level_by="synthetic records of the fit", draws=int(draws), draws_that_ended=int(ended))

    def contrast_power(self, per_record, level: float = 0.05) -> float:
        """The chance that the contrast with the fit reports a discrepancy of the mean that the model lacks, at first order.
        It is refused outside the range where the first order law of the contrast holds."""
        why = self.first_order_holds()
        if why:
            raise ValueError(why)
        V, w = self._gap()
        z = V.T @ (self.shift(per_record) - self.lin.absorbed(per_record))
        return float(ncx2.sf(chi2.isf(level, w.size), w.size, float(np.sum(z * z / w)))) if w.size else float(level)

    def shift(self, per_record) -> np.ndarray:
        """The first order move of the family's estimate of each shared parameter when the records hold a discrepancy of
        the mean that the model lacks: one array a record in the units of the values, None for a record that has none."""
        J, _, parts = self._picture(self.lin.fit.shared, self.lin.fit.own)
        rhs = []
        for k, d in enumerate(per_record):
            Bq = parts[k][3]
            if d is None:
                rhs.append(np.zeros(Bq.shape[0]))
                continue
            w = self._see(self.Q[k], np.asarray(d, dtype=float) / self.lin.problem.records[k].sd)
            rhs.append(w - Bq @ (Bq.T @ w))
        return np.linalg.lstsq(J, np.concatenate(rhs), rcond=None)[0]

    def cost(self, shared, own) -> float:
        pr = self.lin.problem
        return float(sum(np.sum(self._see(q, pr.whitened(k, np.asarray(shared, dtype=float), np.asarray(own[k], dtype=float))) ** 2) for k, q in enumerate(self.Q)))

    def _whole(self, shared, own) -> tuple:
        """The residual as the family sees it and its derivatives by every parameter that it estimates, the shared ones
        first and then the own ones record by record."""
        pr = self.lin.problem
        s = len(pr.shared)
        at = s + np.concatenate([[0], np.cumsum(self.free_at.sum(axis=1))])          # where the own parameters of each record begin: a record gives those that the family estimates in it
        res, rows = [], []
        for k in range(len(pr.records)):
            r, A, B = derivatives(pr, k, shared, own[k])
            block = np.zeros((self._see(self.Q[k], r).size, int(at[-1])))
            block[:, :s] = self._see(self.Q[k], A)
            block[:, at[k]:at[k + 1]] = self._see(self.Q[k], B[:, self.free_at[k]])
            res.append(self._see(self.Q[k], r))
            rows.append(block)
        return np.concatenate(res), np.vstack(rows)

    def estimate(self, shared0=None, own0=None, tol: float = 1e-3, verdict: bool | None = None) -> FamilyFit:
        """The estimate that the family gives by itself: the values at which its statistics of the records agree best with
        its statistics of the model, each statistic in its standard deviation. The search is a least squares search over
        every parameter, each in its own scale, so that an own parameter that the family hardly sees does not carry it
        away. It starts at the fit, or where it is told. It is taken as ended when one more whole step would move the
        shared parameters by less than `tol` standard deviations of the family: on a distance, not on a search that slows.
        `verdict`: whether the family is an estimator, when that was decided at another fit of the same design, as a
        table over several terms does. A family that leaves a direction free at this fit is refused all the same."""
        here = self.bounds()
        if not here["fixed"] or not (here["estimator"] if verdict is None else verdict):
            raise ValueError("this family leaves a direction of the shared parameters free, or holds next to nothing of them: it is a reading, not an estimator")
        pr = self.lin.problem
        s = len(pr.shared)
        held = np.array(self.lin.fit.own if own0 is None else own0, dtype=float)          # the own parameters that the family does not estimate stay at these values
        start = np.concatenate([np.array(self.lin.fit.shared if shared0 is None else shared0, dtype=float), held[self.free_at]])

        def split(z):
            own = held.copy()
            own[self.free_at] = z[s:]          # record by record, in the order of _whole
            return z[:s], own
        got = least_squares(lambda z: -self._whole(*split(z))[0], start, jac=lambda z: self._whole(*split(z))[1], x_scale="jac", method="trf", xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=400)
        shared, own = split(got.x)
        J, r, parts = self._picture(shared, own)
        left = float(np.linalg.norm(J @ np.linalg.lstsq(J, r, rcond=None)[0]))          # what one more whole step would move, in standard deviations of the family
        dof = self.size - s - sum(p[3].shape[1] for p in parts)
        return FamilyFit(shared, own, np.linalg.inv(J.T @ J), float(2.0 * got.cost), int(dof), bool(got.status > 0 and left < tol), int(got.nfev))
