"""Where to search a new term of the model: one report, read from the whole machinery (the owner's messages 139 to 142).

A fit of a working model leaves residuals. The report asks of them the questions of a senior researcher, and gives each
answer its level, counted on synthetic sets of the fitted model that are fitted again as the records were:

    where    the core of the curves, within half a width of their centre, their flanks, from half a width to two, their
             wings, beyond two widths, or their asymmetry
    shape    how the residuals are shaped, read from them alone: a slow trend over the span of each curve, by the
             Legendre polynomials of orders 1 to 3, a periodic structure by a dense scan of periods about each centre,
             fixed to the abscissa and on each side, a part in quadrature with the line by the Hilbert transform of the
             fitted curve, and a structure in one place and at one scale by a scan of wavelets. Each names the shapes
             that these records cannot tell from it
    what     glitches, a noise that grows with the signal, or a noise correlated from sample to sample
    changes along  which place of the residuals changes along each covariate of the design whose value changes from
             record to record. One that takes one value in every record is named among those that cannot be read
    seen by  each family of statistics of the package: the generalized moments, of the curve, of its derivative and of
             its running integral, their ratios, the windows along the abscissa and along the value, the crossings, the
             transformations of the value, and every sample, which is the chi square of the fit, and beside them the
             eight chi functions of `CHI_FUNCTIONS`, the plane of position and residual value, `residual_plane`, the
             chi functions in windows of position and in bands of the value, `windowed_chi` and `value_windowed_chi`,
             in cells of position and of the curve's level or of the residual's value, `gridded_chi` and
             `value_gridded_chi`, the fit made again under seven chi functions, `chi_fits`, the largest excursion of each
             curve's running sum, `running_sum`, and each curve predicted by the fit of the others, `left_out`. The
             reading that sees most says by which statistics the missing term is read: a periodic function by the
             cosines, glitches by a chi function that weighs large residuals more than their square, the exponential of
             half the residual first on the planted glitches of the tests

    report = where(problem, result, sets=100, centre="position")      the answers, each with its statistic, its level,
                                                                       its excess, and the median and 95th percentile of
                                                                       the twin from which the excess is counted
    statistics(problem, result, centre="position")                    the statistics alone, on the records as they are
    null = twin_sets(problem, result, sets=100)                        the statistics on the sets of the twin, once, for
    where(other, its_fit, null=null)                                   records of the same design read against them, as
                                                                       a count of many draws reads one twin of its case

`centre` names the own parameter that places a curve, as `families.catalogue` takes it. A family that the package
refuses on the records, or on more than a few sets of the twin, has no answer.

Each statistic is larger when the residuals hold more of what it reads. Those of the mean grow as the square of the
residuals, and so do the four powers of the chi functions once their mean is mapped by `ON_THE_SQUARE`. The logarithm of one and
the square grows a little slower, Huber's cost as the square inside its bound and as the first power beyond it, and the
exponentials about as the square for residuals near the normal law, faster for skewed or heavy tailed residuals, or for a
few large ones that do not yet hold their mean, until the largest residual bounds them: twice the logarithm of a mean
of exp(r/2) is the sum of the cumulants of the residuals, k1 + k2/4 + k3/24 + k4/192 and on. The readings of the shape of the noise and the plane of
position and value are free of the scale of the residuals: records much cleaner than the twin's put every answer of the
mean beyond every set, where only the excesses order them, and a statistic that grows faster than the square leads by arithmetic, as the raw fourth power or exponential of
the residuals did before their map: on a clean ripple of the toy, at a thousandth of its noise, the two statistics of
the law of the noise once stood at 1e12 and 2.5e10 against 2.3e6 for the best shape.

Its level is the share of the synthetic sets on which it is as large, and its excess says how far beyond them it stands: the statistic less the median of the sets, over
the span from that median to their 95th percentile, so that an excess of one is the 95th percentile of the twin. Where
several answers lie beyond every set their levels are all the smallest, and their excess orders them.

Many answers are read at once, so some answer stands beyond its own 95th percentile on most sets of a right model.
The twin says how large the largest excess is when nothing is missing: `largest_excess_on_the_twin` is the 95th
percentile, over the sets, of the largest excess of each set against the others. An answer beyond it is beyond the
twin, which some answer is on 5 of 100 of the twin's own sets. On other records of a right model it can be more often:
the twin is made from the one fit of these records, so a count on draws of a right model says how often, as the
problem's own twin makes them. The first answer of each question is the first of it that stands beyond the twin, or
none. The glitches and the law of the noise are read once the shape of the residuals that leads is taken out of them: a
shape that the model lacks is not a law of the noise.

The report is blind (the owner's messages 337 and 340): it holds no list of terms and names none, and every answer is a
reading of the residuals. The package names a place, a shape and a covariate, never a cause. Which term of physical origin fills them is the
choice of the researcher, or of the agent that stands in the researcher's place. Nothing of any one field is read here.
"""
from __future__ import annotations

import hashlib

import numpy as np

from .losses import HUBER

from .candidates import Candidate, covariates_of, from_mechanism
from .families import catalogue, family
from .fit import fit
from .linear import Linear
from .propose import ABSORBED, RANK, _one
from .report import by_position, by_position_and_record
from .twin import replicas

CORE, WINGS = 0.5, 2.0          # the core is within half a width of the centre of a curve, the flanks from there to two widths, the wings beyond: a term that acts between the core and the wings, had no place of its own (the owner's messages 184 and 187)
PLACES = (("the core", lambda r: r["to"] <= CORE), ("the flanks", lambda r: r["from"] >= CORE and r["to"] <= WINGS), ("the wings", lambda r: r["from"] >= WINGS))          # each band of report.BANDS in one place
SHAPES_READ = (("trend", "a slow trend over the span of each curve"), ("periodic", "a periodic structure along the abscissa"),
               ("periodic fixed", "a periodic structure fixed to the abscissa, in phase from curve to curve"),
               ("periodic local", "a periodic structure on one side of the curves"),
               ("quadrature", "a part in quadrature with the line, the Hilbert transform of the fitted curve"),
               ("localized", "a structure in one place and at one scale, a wavelet of the residuals"))          # the four after the first two since rows M589 and M590, his message 533: the Fourier transform in two frames and on each side, and the Hilbert transform          # the answers of the question shape, each a reading of the residuals
TREND_CYCLES = 2.0          # a structure of fewer cycles over the span of a curve than this is a slow trend, of more a periodic one: the Legendre polynomials of orders 1 to 3 of `slow_trend` cross zero up to three times
SIMILAR = 0.9        # two classes whose leading shapes, as the records see them, have a cosine beyond this are not told apart
TIE = 1.96          # two places whose excesses differ by less than this many standard deviations of their difference on the records are named together, his decision of message 388
MISFIT_SETS = 30          # the sets of the twin that carries the misfit, which reads the scatter of an excess on the records, his message 389
SMOOTH = 9          # the samples over which a residual is averaged into the misfit that the twin carries, as the noise of neighbouring samples is read


def statistics(problem, result, centre: str | None = None, terms=(), dictionary: bool = False) -> dict:
    """The statistics of the report on the records as they are: for each answer a number that grows with what the
    residuals hold of it. A whitened residual is a residual over its standard deviation. The glitches and the law of
    the noise are read on the residuals once the shape that leads is taken out of them. `terms`: the researcher's list
    of terms, Mechanisms, empty by default (the owner's message 342). `dictionary`: the package's own list beside it,
    only where the researcher asks it (337 to 343)."""
    return _read(problem, result, centre, terms, dictionary)[0]


def changes(problem, result) -> dict:
    """How the structure of what the fit leaves changes along each setting of the design, the owner's message 170: "the
    whole point of the package is to look at the structure (and how it's varying with the parameters!) of the residuals".
    For each covariate that takes more than two values across the records, and each place, the core, the flanks, the wings and the
    asymmetry, the readings of `report.by_position_and_record` of each band of that place are set against the
    covariate, centred and of unit length over the records: the sum over the records of covariate times reading, a
    standard normal number at first order where nothing changes along it. The statistic of the place is the sum of the
    squares of those numbers over its bands, so that it grows as the square of the residuals, as every statistic of
    the mean of the report, and the four powers of the chi functions once mapped by `ON_THE_SQUARE`. A covariate with one value, or two, says nothing of a trend and is not read."""
    rows = by_position_and_record(problem, result)
    out = {}
    for key in covariates_of(problem):
        z = np.array([float(rec.design[key]) for rec in problem.records])
        z = z - z.mean()
        if len(set(np.round(z, 12))) <= 2 or not np.any(z):
            continue
        z = z / np.sqrt(z @ z)
        for place, pick, side in [(name, pick, "both_sides") for name, pick in PLACES] + [("the asymmetry", lambda r: True, "one_side_less_the_other")]:
            trend = [float(np.nansum(z * np.asarray(r[side], dtype=float))) for r in rows if pick(r)]
            if trend:
                out[("changes along", f"{place} along {key}")] = float(sum(t * t for t in trend))
    return out


def slow_trend(lin) -> Candidate:
    """The Legendre polynomials of orders 1 to 3 over the span of each curve, as the columns of one reading: how much
    of the residuals is a slow trend over the span. A reading of the residuals and no term."""
    cols = []
    for order in (1, 2, 3):
        col = []
        for rec in lin.problem.records:
            x = np.asarray(rec.x, dtype=float)
            u = 2.0 * (x - x.min()) / max(float(np.ptp(x)), 1e-300) - 1.0
            col.append(np.polynomial.legendre.legval(u, [0.0] * order + [1.0]))
        cols.append(col)
    return Candidate("a slow trend over the span of each curve", "trend", cols, "any")


def matched_terms(lin, terms, centre: str | None = None, dictionary: bool = False) -> dict:
    """The researcher's list of terms set against the residuals: for each term, the part of what every family reads
    that its shape accounts for, `matched_families`, as the statistic of the question matched. A term asked at a sign
    that it cannot have scores zero. The package holds no such list, and an empty one reads nothing (the owner's message
    342: "For the list of terms, allow the user to define it, but keep it empty at the beginning"). With `dictionary`,
    the package's generic list stands beside it, only where the researcher asks it with a flag (337 to 343)."""
    if not terms and not dictionary:
        return {}
    cands = [from_mechanism(lin, m, float(m.start) if m.start else 0.1) for m in terms]
    if dictionary:
        from .dictionary import candidates as generic
        cands += generic(lin, covariates=covariates_of(lin.problem))
    try:
        got = matched_families(lin, cands, centre)
    except ValueError:
        return {}          # no family reads these records, so nothing is matched
    out = {("matched", c.get("key") or c["name"]): float(c["score"]) for c in got["candidates"]}          # by its key: a period in the name moves with each fit, and a twin's sets would read another name than the records'
    out.update({("matched", r.get("key") or r["name"]): 0.0 for r in got["of_the_sign_they_cannot_have"]})
    return out


CHI_FUNCTIONS = (          # each cost of a whitened residual, averaged over the records and mapped by ON_THE_SQUARE, the exponentials with the residual clipped at 60 so that they stay finite: his messages 1 and 4, "the whole family of residuals obtained with chi^p - and maybe even exp(chi), log(chi)"
    ("the chi function: the residual to the power 1", lambda r: np.abs(r)),
    ("the chi function: the residual to the power 1.5", lambda r: np.abs(r) ** 1.5),
    ("the chi function: the residual to the power 3", lambda r: np.abs(r) ** 3),
    ("the chi function: the residual to the power 4", lambda r: r ** 4),
    ("the chi function: the logarithm of one and the square", lambda r: np.log1p(r * r)),
    ("the chi function of Huber", lambda r: np.where(np.abs(r) <= HUBER, 0.5 * r * r, HUBER * (np.abs(r) - 0.5 * HUBER))),
    ("the chi function: the exponential of half the residual, the side above", lambda r: np.exp(0.5 * np.clip(r, -60.0, 60.0))),
    ("the chi function: the exponential of half the residual, the side below", lambda r: np.exp(-0.5 * np.clip(r, -60.0, 60.0))),
)
ON_THE_SQUARE = {          # the mean m of each chi function mapped so that near the normal law the four powers grow as the square of the residuals and the others about as fast or slower: a monotone map, so its level on the twin is that of the mean. The logarithm of one and the square grows a little slower, Huber's cost as the first power beyond its bound, and the exponentials faster than the square for skewed or heavy tailed residuals, until the largest residual bounds them (the stamp's fifth reading) (before any map, the clock of the README read the exponential at 368590 against 28522 for the place along the power)
    "the chi function: the residual to the power 1": lambda m: m ** 2,
    "the chi function: the residual to the power 1.5": lambda m: m ** (4.0 / 3.0),
    "the chi function: the residual to the power 3": lambda m: m ** (2.0 / 3.0),
    "the chi function: the residual to the power 4": lambda m: m ** 0.5,
    "the chi function: the logarithm of one and the square": lambda m: np.expm1(m),
    "the chi function of Huber": lambda m: m,          # its mean grows as the square of the residuals inside its bound and slower beyond it: squared it grew as their fourth power near the noise
    "the chi function: the exponential of half the residual, the side above": lambda m: 2.0 * np.log(m),          # for residuals near the normal law, log m already grows as their square: squared again it grew as their fourth power (the stamp's fourth reading)
    "the chi function: the exponential of half the residual, the side below": lambda m: 2.0 * np.log(m),
}
PLANE = "the plane of position and residual value: intervals of position by bands of value"


def residual_plane(problem, result, lin) -> dict:
    """The readings of the value of the residuals, his messages 1, 4, 20, 39, 48, 69, 76, 88 and 98 read again for his
    message 527: each chi function of `CHI_FUNCTIONS` averaged over the whitened residuals and mapped
    by `ON_THE_SQUARE`, times their number, and the windows on the value of the
    residuals, `losses.grid`, the residuals counted on the plane of position and value in their robust scale, whose
    statistic is the test of independence of that table. Each is read against the twin as every answer of the report."""
    from .losses import _scale, grid
    r = np.asarray(lin.r, dtype=float)
    out = {said: float(r.size * ON_THE_SQUARE[said](float(np.mean(cost(r))))) for said, cost in CHI_FUNCTIONS}          # n times the mean of the cost, mapped so that near the normal law none grows much faster than the square
    try:
        out[PLANE] = float(grid(problem, result, _scale(r))["chi_square"])
    except ValueError:
        pass          # residuals with no spread, or a fit that did not end, have no table of position and value
    return out


TUKEY, WELSCH = 4.685, 2.985          # the constants of the biweight of Tukey and of the cost of Welsch, each 95 per cent efficient under a normal law
WINDOW_COSTS = CHI_FUNCTIONS + (          # the costs that the scans of windows read, wider than the eight answers of the report: his 536 and 537, review 235, N5
    ("the chi function: the residual to the power 0.5", lambda r: np.abs(r) ** 0.5),
    ("the chi function: the residual to the power 6", lambda r: np.asarray(r, dtype=float) ** 6),
    ("the chi function of Tukey, the biweight", lambda r: np.where(np.abs(r) <= TUKEY, TUKEY ** 2 / 6.0 * (1.0 - (1.0 - (np.asarray(r, dtype=float) / TUKEY) ** 2) ** 3), TUKEY ** 2 / 6.0)),
    ("the chi function of Welsch", lambda r: 0.5 * WELSCH ** 2 * (1.0 - np.exp(-(np.asarray(r, dtype=float) / WELSCH) ** 2))),
)
WINDOWED = "the chi functions in windows of position: every cost of the residual in every smooth window, about the centres or fixed to the abscissa"          # row M593, the owner's message 536
WINDOW_SCALES = (0.25, 0.5, 1.0, 2.0)          # the scales of the windows of `windowed_chi`, in widths of the curves, as the wavelet's


def windowed_chi(lin, scales=WINDOW_SCALES) -> dict | None:
    """Generalised windowing of the generalised chi functions, row M593 (the owner's messages 536, 537 and 539): each
    cost of `WINDOW_COSTS`, the eight of `CHI_FUNCTIONS` and four more, and the square, in every smooth window of position,
    a Gaussian of a scale of `scales` about places half a scale apart, in two frames: counted in widths from the centre
    of each curve, the same for every curve, and fixed to the abscissa, in widths of the median curve. Under each window
    the cost of each whitened residual less its mean over all the residuals, summed with the window's weights and
    divided by the root of the sum of their squares times the cost's variance: a standardised excess of the cost in that
    window, about normal on a right model. The largest square over every cost and every window, with its cost, place and
    scale. A cost gathered in one place, which a sum over every residual spreads thin, stands out here, and a level of
    the cost over all the residuals does not, which the chi functions read. None where a curve holds fewer than eight
    samples or no cost has a spread."""
    pr = lin.problem
    r = np.asarray(lin.r, dtype=float)
    about, fixed = np.empty(r.size), np.empty(r.size)
    unit = float(np.median([lin.width_of(k) for k in range(len(pr.records))]))
    for k, rec in enumerate(pr.records):
        x = np.asarray(rec.x, dtype=float)
        if x.size < 8:
            return None
        about[lin.rows_of(k)] = (x - lin.centre_of(k)) / lin.width_of(k)
        fixed[lin.rows_of(k)] = x / unit          # the frame of the abscissa, in widths of the median curve, his 539
    costs = [("the square of the residual", r * r)] + [(said, cost(r)) for said, cost in WINDOW_COSTS]          # the wider family, row M593
    best = None
    for frame, u in (("about each centre", about), ("fixed to the abscissa", fixed)):
        for scale in scales:
            places = np.arange(np.floor(u.min() / (0.5 * scale)), np.ceil(u.max() / (0.5 * scale)) + 1.0) * 0.5 * scale
            W = np.exp(-0.5 * ((u[:, None] - places[None, :]) / scale) ** 2)
            root = np.sqrt((W * W).sum(axis=0))
            for said, c in costs:
                z = c - c.mean()
                v = float(z @ z) / z.size
                if not v > 0.0:
                    continue
                stat = (W.T @ z) / (root * np.sqrt(v))
                j = int(np.argmax(stat * stat))
                if best is None or stat[j] ** 2 > best["largest"]:
                    best = {"largest": float(stat[j] ** 2), "cost": said, "place": float(places[j]), "scale": float(scale), "sign": "more" if stat[j] > 0 else "less",
                            "frame": frame, "unit": unit if frame == "fixed to the abscissa" else None}          # the place in widths of each curve from its centre, or in widths of the median curve along the abscissa
    return best


WINDOWED_VALUE = "the chi functions in windows of the residual's value: every cost of the residual in every band of its value"          # row M594, the owner's messages 76 and 537
VALUE_BANDS = (-2.0, -1.0, 0.0, 1.0, 2.0)          # the edges of the bands of the whitened residual's value of `value_windowed_chi`, open at both ends


def _normal_parts(cost, low: float, high: float) -> tuple:
    """The mean and the variance, under a normal law of unit deviation, of `cost` of a residual times the indicator of its
    value lying between `low` and `high`: by the trapezoid rule on a grid of 4001 points inside the band, cut at 12, so
    that no edge of the band falls between two points of the grid."""
    lo, hi = max(low, -12.0), min(high, 12.0)
    if not hi > lo:
        return 0.0, 0.0
    z = np.linspace(lo, hi, 4001)
    density = np.exp(-0.5 * z * z) / np.sqrt(2.0 * np.pi)
    c = np.asarray(cost(z), dtype=float)
    first, second = float(np.trapezoid(c * density, z)), float(np.trapezoid(c * c * density, z))
    return first, second - first * first


def value_windowed_chi(lin, edges=VALUE_BANDS) -> dict | None:
    """The chi functions in windows of the residual's value, row M594 (the owner's messages 76 and 537): each cost of
    `WINDOW_COSTS`, the eight of `CHI_FUNCTIONS` and four more, and the count, summed over the whitened residuals whose value lies in one band of `edges`, less n
    times its mean there under a normal law of unit deviation, divided by the root of n times its variance there. The
    largest square over every cost and every band, with its cost, band and sign. The residuals are read in their own
    robust scale, 1.4826 times their median absolute deviation, so that a spread that a missing term widens is not read
    as a law of the noise: a ripple over every curve led the report's families at 171 when they were not. A law of the
    noise wrong in one band of its values, a tail heavier than the normal one, a tail on one side or a core too narrow,
    stands out, and the twin of the report sets it against what a right model's residuals give. None where fewer than
    30 residuals are read or their scale is nought."""
    r = np.asarray(lin.r, dtype=float)
    if r.size < 30:
        return None
    scale = 1.4826 * float(np.median(np.abs(r - np.median(r))))          # the robust scale of the residuals: a shape of the law and not its width
    if not scale > 0.0:
        return None
    r = (r - np.median(r)) / scale
    bounds = [-np.inf] + list(edges) + [np.inf]
    costs = [("the count of the residuals", lambda z: np.ones_like(np.asarray(z, dtype=float)))] + list(WINDOW_COSTS)          # the wider family, row M593
    best = None
    for low, high in zip(bounds[:-1], bounds[1:]):
        inside = (r > low) & (r <= high)
        for said, cost in costs:
            mean, var = _normal_parts(cost, low, high)
            if not var > 0.0:
                continue
            stat = (float(np.sum(np.asarray(cost(r[inside]), dtype=float))) - r.size * mean) / np.sqrt(r.size * var)
            if best is None or stat * stat > best["largest"]:
                best = {"largest": float(stat * stat), "cost": said, "band": (float(low), float(high)), "sign": "more" if stat > 0 else "less"}
    return best


GRIDDED_LEVEL = "the chi functions in cells of position and of the curve's level: every cost in every place where the curve stands at a height"          # row M588, the owner's message 572
GRIDDED_VALUE = "the chi functions in cells of position and of the residual's value: every cost in every band of the residual's value, place by place"          # row M588, the owner's message 572
GRID_SCALES = (0.5, 1.0, 2.0)          # the scales of the windows of position of the two grids, in widths of the curves
LEVEL_SCALES = (0.125, 0.25, 0.5)          # the scales of the windows of the curve's level, in fractions of its height above its median
VALUE_SCALES = (0.5, 1.0)          # the scales of the windows of the residual's value, in its robust scale
GRID_LEAST = 8          # the fewest samples within one scale of a cell's centre, in each of its directions, that a cell of a grid reads: a lone sample beyond the others, which the tails of a smooth window weigh in, read 523 in the cells of the value counted by the weights alone, against 7 counted in this way, on 2026-11-12


def _places_and_levels(lin):
    """The place of each whitened residual about the centre of its curve, in the curve's widths, and the level of its
    curve there, the fitted value's distance from the curve's median over the largest such distance, from 0 at the
    median to 1 at the top of a peak or the floor of a dip. None where a curve holds fewer than eight samples."""
    pr = lin.problem
    n = np.asarray(lin.r).size
    u, h = np.empty(n), np.empty(n)
    for k, rec in enumerate(pr.records):
        x = np.asarray(rec.x, dtype=float)
        if x.size < 8:
            return None
        rows = lin.rows_of(k)
        u[rows] = (x - lin.centre_of(k)) / lin.width_of(k)
        f = np.asarray(pr.predict(k, lin.fit.shared, lin.fit.own[k]), dtype=float)
        up = np.abs(f - np.median(f))
        h[rows] = up / up.max() if up.max() > 0.0 else 0.0
    return u, h


def _windows(v, scale: float, low: float | None = None, high: float | None = None):
    """Gaussian windows of `scale` along the values `v`, about places half a scale apart from `low` to `high`, the range
    of `v` where none is given: the matrix of weights, a column a window, and the places."""
    low = float(np.min(v)) if low is None else low
    high = float(np.max(v)) if high is None else high
    places = np.arange(np.floor(low / (0.5 * scale)), np.ceil(high / (0.5 * scale)) + 1.0) * 0.5 * scale
    return np.exp(-0.5 * ((np.asarray(v, dtype=float)[:, None] - places[None, :]) / scale) ** 2), places


def gridded_chi(lin, scales=GRID_SCALES, level_scales=LEVEL_SCALES) -> dict | None:
    """Generalised griding of the data, row M588 (the owner's message 572: "generalised windowing of the data and of the
    residuals either in the horizontal direction and in the vertical direction"): each cost of `WINDOW_COSTS`, and the
    square, in every cell of a smooth window of position about each centre, a Gaussian of a scale of `scales`, times a
    smooth window of the curve's level, a Gaussian of a scale of `level_scales` on the fitted value's height. Under each
    cell the cost of each whitened residual less its mean over all the residuals, summed with the cell's weights and
    divided by the root of the sum of their squares times the cost's variance: a standardised excess of the cost there,
    about normal on a right model. The largest square over every cost and every cell that holds `GRID_LEAST` samples or
    more within one scale of its centre in each direction, with its cost, place, level, scales and sign. Where the curves hold several lines, the level tells a miss at
    the top of the tallest from one at the same place of a smaller, which a window of position alone joins. None where a
    curve holds fewer than eight samples or no cost has a spread."""
    got = _places_and_levels(lin)
    if got is None:
        return None
    u, h = got
    r = np.asarray(lin.r, dtype=float)
    costs = [("the square of the residual", r * r)] + [(said, np.asarray(cost(r), dtype=float)) for said, cost in WINDOW_COSTS]
    best = None
    for scale in scales:
        Wx, places = _windows(u, scale)
        for level_scale in level_scales:
            Wh, levels = _windows(h, level_scale, 0.0, 1.0)
            root = np.sqrt((Wx * Wx).T @ (Wh * Wh))
            held = (np.abs(u[:, None] - places[None, :]) < scale).astype(float).T @ (np.abs(h[:, None] - levels[None, :]) < level_scale).astype(float) >= GRID_LEAST
            for said, c in costs:
                z = c - c.mean()
                v = float(z @ z) / z.size
                if not v > 0.0:
                    continue
                stat = np.where(held, (Wx.T @ (z[:, None] * Wh)) / np.where(held, root, 1.0) / np.sqrt(v), 0.0)
                i, j = np.unravel_index(int(np.argmax(stat * stat)), stat.shape)
                if held[i, j] and (best is None or stat[i, j] ** 2 > best["largest"]):
                    best = {"largest": float(stat[i, j] ** 2), "cost": said, "place": float(places[i]), "scale": float(scale),
                            "level": float(levels[j]), "level_scale": float(level_scale), "sign": "more" if stat[i, j] > 0 else "less"}
    return best


def _normal_windowed(cost, centre: float, scale: float) -> tuple:
    """The mean and the variance, under a normal law of unit deviation, of `cost` of a residual times a Gaussian window of
    `scale` about `centre` on its value: by the trapezoid rule on 8001 points from -12 to 12."""
    z = np.linspace(-12.0, 12.0, 8001)
    density = np.exp(-0.5 * z * z) / np.sqrt(2.0 * np.pi)
    c = np.asarray(cost(z), dtype=float) * np.exp(-0.5 * ((z - centre) / scale) ** 2)
    first, second = float(np.trapezoid(c * density, z)), float(np.trapezoid(c * c * density, z))
    return first, second - first * first


def value_gridded_chi(lin, scales=GRID_SCALES, value_scales=VALUE_SCALES) -> dict | None:
    """Generalised griding of the residuals, row M588 (the owner's message 572): each cost of `WINDOW_COSTS`, and the
    count, in every cell of a smooth window of position about each centre, a Gaussian of a scale of `scales`, times a
    smooth window of the residual's value, a Gaussian of a scale of `value_scales` about values from -3 to 3, the
    residuals read in their robust scale as `value_windowed_chi` reads them. Under each cell the cost times the value's
    window, less its mean under a normal law of unit deviation, summed with the position's weights and divided by the
    root of the sum of their squares times the variance there: a standardised excess, about normal on a right model. The
    largest square over every cost and every cell whose window of position holds `GRID_LEAST` samples or more within one scale of its centre, with its
    cost, place, value, scales and sign. A law of the noise wrong in one place alone, glitches gathered in the wings or
    a tail on one side near the centre, which the bands of the value over every place spread thin, stands out here.
    None where fewer than 30 residuals are read, a curve holds fewer than eight samples, or their scale is nought."""
    r = np.asarray(lin.r, dtype=float)
    if r.size < 30:
        return None
    got = _places_and_levels(lin)
    if got is None:
        return None
    u, _ = got
    spread = 1.4826 * float(np.median(np.abs(r - np.median(r))))
    if not spread > 0.0:
        return None
    r = (r - np.median(r)) / spread
    costs = [("the count of the residuals", lambda z: np.ones_like(np.asarray(z, dtype=float)))] + list(WINDOW_COSTS)
    best = None
    for scale in scales:
        Wx, places = _windows(u, scale)
        held = (np.abs(u[:, None] - places[None, :]) < scale).sum(axis=0) >= GRID_LEAST
        root = np.sqrt((Wx * Wx).sum(axis=0))
        for value_scale in value_scales:
            values = np.arange(-3.0, 3.0 + 1e-9, 0.5 * value_scale)
            G = np.exp(-0.5 * ((r[:, None] - values[None, :]) / value_scale) ** 2)
            for said, cost in costs:
                parts = [_normal_windowed(cost, b, value_scale) for b in values]
                mean = np.array([p[0] for p in parts])
                var = np.array([p[1] for p in parts])
                good = var > 0.0
                if not good.any():
                    continue
                C = np.asarray(cost(r), dtype=float)[:, None] * G - mean[None, :]
                stat = np.where(held[:, None] & good[None, :], (Wx.T @ C) / (root[:, None] * np.sqrt(np.where(good, var, 1.0))[None, :]), 0.0)
                i, j = np.unravel_index(int(np.argmax(stat * stat)), stat.shape)
                if held[i] and good[j] and (best is None or stat[i, j] ** 2 > best["largest"]):
                    best = {"largest": float(stat[i, j] ** 2), "cost": said, "place": float(places[i]), "scale": float(scale),
                            "value": float(values[j]), "value_scale": float(value_scale), "sign": "more" if stat[i, j] > 0 else "less"}
    return best


CHI_FITS = "the fit under every chi function: how far each shared parameter moves from the fit by the square, in the spread that a right model gives"          # row M579, the owner's messages 4, 25 and 573
CHI_FIT_COSTS = (("the chi function of Huber", "huber"), ("the logarithm of one plus the residual's square", "logarithm"),
                 ("the residual to the power 1.5", 1.5), ("the residual to the power 3", 3.0), ("the residual to the power 4", 4.0),
                 ("the exponential of half the residual", "exponential up"), ("the exponential of minus half the residual", "exponential down"))          # the owner's messages 4 and 25, "chi^4, o exp(chi)": each with its full fit, `losses.chi_fit` or `losses.direct_fit`, to check the first order
CHI_FIT_FLOOR = 0.05          # the least square of the spread of a move, over the square of the bar of the fit: a cost as efficient as the square moves by nothing that a right model would not


def chi_fits(lin, costs=CHI_FIT_COSTS) -> dict | None:
    """Generalised chi functions, the owner's messages 4, 25, 573 and 574: "not only the residuals obtaining the min of
    chi^2 but ... the whole family of residuals obtained with chi^p - and maybe even exp(chi), log(chi)", and "non è
    nemmeno necessario che si lavori con chi^2 per calcolare la MLE. Si può per esempio usare chi^4, o exp(chi) o altre
    funzioni, che creano diverse strutture di residui". The model fitted again under each cost of `costs` by
    `losses.direct_fit`, from the fit by the square and in its robust scale, and each shared parameter's move. Fitted
    again and not at first order: one Newton step from the fit by the square was within a few per cent of the full fit
    for the four milder costs, and four to five times off for the fourth power and the exponentials on two glitches,
    so it is not used for any of them. On a right model every cost estimates the same values,
    and a cost moves only by the spread of its own estimate less that of the square's, the cost's scatter factor squared
    less one, times the square of the robust scale times the fit's bar (Hausman 1978). A term that the model lacks pulls the fit by the square,
    and a cost that weighs the large residuals less, or more, does not follow it. Each move in that spread and in the
    fit's bars, for every shared parameter and every cost whose fit ended, with the largest square, its cost and its
    parameter, and the costs whose fit did not end. None where the problem shares no parameter, fewer than 30 residuals
    are read, or their scale is nought."""
    from .losses import cost_parts, direct_fit
    r = np.asarray(lin.r, dtype=float)
    J = np.asarray(lin.J, dtype=float)
    if J.ndim != 2 or J.shape[1] == 0 or r.size < 30:
        return None
    try:
        sd = np.sqrt(np.diag(np.linalg.inv(J.T @ J)))
    except np.linalg.LinAlgError:
        return None
    s = 1.4826 * float(np.median(np.abs(r - np.median(r))))
    if not s > 0.0:
        return None
    u = r / s
    names = list(lin.problem.shared)
    base = np.asarray(lin.fit.shared, dtype=float)
    moves, best, unended = {}, None, []
    for said, cost in costs:
        _, psi, dpsi = cost_parts(cost, u)
        if not float(np.mean(dpsi)) > 0.0:
            continue
        got = direct_fit(lin.problem, lin.fit, cost, scale=s)
        if not got["converged"]:
            unended.append(said)
            continue
        delta = np.asarray(got["shared"], dtype=float) - base
        factor = float(np.sqrt(np.mean(psi * psi)) / np.mean(dpsi))          # the cost's scatter over the square's in the robust scale: its estimate's bar is s times this times the fit's bar, and the square's is s times the fit's bar
        z = delta / (sd * s * np.sqrt(max(factor * factor - 1.0, CHI_FIT_FLOOR)))          # s outside the root: written inside it, a record whose deviations were stated twice too large read its moves 40 per cent apart (2026-11-12)
        moves[said] = {name: {"bars": float(d / b), "spreads": float(zz)} for name, d, b, zz in zip(names, delta, sd, z)}
        j = int(np.argmax(z * z))
        if best is None or z[j] ** 2 > best["largest"]:
            best = {"largest": float(z[j] ** 2), "cost": said, "parameter": names[j], "bars": float(delta[j] / sd[j]), "spreads": float(z[j])}
    if best is None:
        return None
    return dict(best, moves=moves, unended=unended)


FOURIER_CHI = "the chi functions in windows of frequency of the residuals' Fourier transform: every cost of the coefficients, band by band"          # the owner's messages 533, 574 and 575
DIFFERENCE_CHI = "the chi functions in windows of position of the residuals' differences: every cost of each step from sample to sample"          # the owner's message 575
RUNNING_SUM = "the running sum of the residuals: its largest excursion along each curve"          # the owner's message 575
FREQUENCY_SCALES = (0.1, 0.2, 0.4)          # the scales of the windows of frequency, in decades of the frequency over the samples' own


def _robust(lin) -> np.ndarray | None:
    """The whitened residuals in their robust scale, 1.4826 times their median absolute deviation about their median, as
    `value_windowed_chi` reads them: a transform then reads the shape of the residuals and not the deviations as they
    were stated, as `chi_fits` needs too. None where that scale is nought."""
    r = np.asarray(lin.r, dtype=float)
    s = 1.4826 * float(np.median(np.abs(r - np.median(r))))
    return (r - np.median(r)) / s if s > 0.0 else None


def _coefficients(lin):
    """The Fourier coefficients of each record's whitened residuals, a cosine and a sine for every whole number of cycles
    j over its span from 1 to under half its samples, scaled by the root of 2 over its samples so that, on a white noise
    of unit deviation sampled evenly, each is about normal of unit deviation, with the frequency of each as a fraction of
    the record's own sampling, from 1/n to 1/2. None where a record holds fewer than eight samples."""
    pr = lin.problem
    robust = _robust(lin)
    if robust is None:
        return None
    z, nu = [], []
    for k, rec in enumerate(pr.records):
        x = np.asarray(rec.x, dtype=float)
        n = x.size
        if n < 8:
            return None
        r = robust[lin.rows_of(k)]
        span = float(x[-1] - x[0]) + float(np.median(np.diff(x)))
        j = np.arange(1, (n - 1) // 2 + 1)
        phase = 2.0 * np.pi * np.outer(x - x[0], j) / span
        scale = np.sqrt(2.0 / n)
        z += [scale * (r @ np.cos(phase)), scale * (r @ np.sin(phase))]
        nu += [j / n, j / n]
    return np.concatenate(z), np.concatenate(nu)


def fourier_chi(lin, scales=FREQUENCY_SCALES) -> dict | None:
    """Generalised chi functions on the Fourier transform of the residuals, the owner's messages 533, 574 and 575: "Improve
    the statistical machinery also with Fourier transform or other transforms (for instance to detect periodical terms in
    the data and/or in the residuals)". Each cost of `WINDOW_COSTS`, and the square, of the residuals' Fourier
    coefficients, `_coefficients`, in every smooth window of frequency, a Gaussian of a scale of `scales` in decades about
    places half a scale apart, all the records together: the cost of each coefficient less its mean under a normal law of
    unit deviation, summed with the window's weights and divided by the root of the sum of their squares times its
    variance there. The scan of periods reads one frequency by the square of its two coefficients: these read a band, a
    noise whose power rises or falls with the frequency, several lines, or a chirp that no one period holds, and every
    cost tells a band of large coefficients from one of many moderate ones. A window is read where it holds `GRID_LEAST`
    coefficients within one scale of its centre, which binds only on a few records of few samples, at the lowest
    frequencies. The largest square, with its cost, its band as a fraction of the sampling and in cycles over the median
    curve's width, its scale and sign. None where a record holds fewer than eight samples."""
    got = _coefficients(lin)
    if got is None:
        return None
    z, nu = got
    lg = np.log10(nu)
    unit = float(np.median([lin.width_of(k) for k in range(len(lin.problem.records))]))
    step = float(np.median([np.median(np.diff(np.asarray(rec.x, dtype=float))) for rec in lin.problem.records]))
    costs = [("the square of the coefficient", lambda v: np.asarray(v, dtype=float) ** 2)] + list(WINDOW_COSTS)
    best = None
    for said, cost in costs:
        mean, var = _normal_parts(cost, -np.inf, np.inf)
        if not var > 0.0:
            continue
        c = np.asarray(cost(z), dtype=float) - mean
        for scale in scales:
            W, places = _windows(lg, scale)
            held = (np.abs(lg[:, None] - places[None, :]) < scale).sum(axis=0) >= GRID_LEAST
            stat = np.where(held, (W.T @ c) / (np.sqrt((W * W).sum(axis=0)) * np.sqrt(var)), 0.0)
            j = int(np.argmax(stat * stat))
            if held[j] and (best is None or stat[j] ** 2 > best["largest"]):
                band = float(10.0 ** places[j])
                best = {"largest": float(stat[j] ** 2), "cost": said, "frequency": band, "cycles_per_width": band / step * unit, "scale": float(scale), "sign": "more" if stat[j] > 0 else "less"}
    return best


def difference_chi(lin, scales=WINDOW_SCALES) -> dict | None:
    """Generalised chi functions on the residuals' differences from sample to sample, the owner's message 575: each cost
    of `WINDOW_COSTS`, and the square, of the differences of the whitened residuals over the root of two, in every smooth
    window of position about the centre of each curve, standardised as `windowed_chi` standardises the residuals. A step
    in a curve is one large difference where the residuals hold a level on each side, and a noise that is rougher in one
    place than the model says reads there, where the residuals themselves spread both over a width. The largest square,
    with its cost, place, scale and sign. None where a curve holds fewer than eight samples or no cost has a spread."""
    pr = lin.problem
    robust = _robust(lin)
    if robust is None:
        return None
    d, u = [], []
    for k, rec in enumerate(pr.records):
        x = np.asarray(rec.x, dtype=float)
        if x.size < 8:
            return None
        r = robust[lin.rows_of(k)]
        d.append(np.diff(r) / np.sqrt(2.0))
        u.append((0.5 * (x[1:] + x[:-1]) - lin.centre_of(k)) / lin.width_of(k))
    d, u = np.concatenate(d), np.concatenate(u)
    costs = [("the square of the difference", d * d)] + [(said, np.asarray(cost(d), dtype=float)) for said, cost in WINDOW_COSTS]
    best = None
    for scale in scales:
        W, places = _windows(u, scale)
        held = (np.abs(u[:, None] - places[None, :]) < scale).sum(axis=0) >= GRID_LEAST
        root = np.sqrt((W * W).sum(axis=0))
        for said, c in costs:
            z = c - c.mean()
            v = float(z @ z) / z.size
            if not v > 0.0:
                continue
            stat = np.where(held, (W.T @ z) / (root * np.sqrt(v)), 0.0)
            j = int(np.argmax(stat * stat))
            if held[j] and (best is None or stat[j] ** 2 > best["largest"]):
                best = {"largest": float(stat[j] ** 2), "cost": said, "place": float(places[j]), "scale": float(scale), "sign": "more" if stat[j] > 0 else "less"}
    return best


def running_sum(lin) -> dict | None:
    """The running sum of each record's whitened residuals over the root of its samples, the owner's message 575: on a
    right model a bridge of a random walk whose largest excursion stays near one, the statistic of Kolmogorov and
    Smirnov's test. A level that shifts partway along a curve, or a drift that the model does not hold, bends it far
    from nought where the residuals change, which a sum over windows reads only at the window that happens to fit. The
    square of the largest excursion over every record, with its record and its place in the curve's widths from the
    centre. None where a record holds fewer than eight samples."""
    pr = lin.problem
    robust = _robust(lin)
    if robust is None:
        return None
    best = None
    for k, rec in enumerate(pr.records):
        x = np.asarray(rec.x, dtype=float)
        if x.size < 8:
            return None
        r = robust[lin.rows_of(k)]
        s = np.cumsum(r - r.mean()) / np.sqrt(r.size)
        i = int(np.argmax(np.abs(s)))
        if best is None or s[i] ** 2 > best["largest"]:
            best = {"largest": float(s[i] ** 2), "record": k, "place": float((x[i] - lin.centre_of(k)) / lin.width_of(k)), "sign": "more" if s[i] > 0 else "less"}
    return best


LEFT_OUT = "each curve predicted by the fit of the others"          # the answer of the cross validation, row M590


def left_out(lin) -> dict:
    """Each curve predicted by the fit of the others, row M590 (the owner's message 534): the cross validation of a joint
    fit, at first order. Without curve k the shared parameters move by d = -(A - Jk'Jk)^-1 Jk' rk, A the curvature of the
    whole fit, Jk and rk the reduced columns and the whitened residuals of curve k, its own parameters freed as the fit
    frees them. The rise of its chi square there, |rk - Jk d|^2 - |rk|^2, is what the others mispredict of it. The sum
    over the curves, which the report sets against the twin, and the rise of each curve. Nothing where the model shares
    no parameter or holds one curve."""
    pr = lin.problem
    p = len(pr.shared)
    if not p or len(pr.records) < 2:
        return {}
    J, r = np.asarray(lin.J, dtype=float), np.asarray(lin.r, dtype=float)
    A = J.T @ J
    rises = []
    for k in range(len(pr.records)):
        rows = lin.rows_of(k)
        Jk, rk = J[rows], r[rows]
        try:
            d = -np.linalg.solve(A - Jk.T @ Jk, Jk.T @ rk)
        except np.linalg.LinAlgError:          # the others alone do not fix the shared parameters: no prediction of this curve
            return {}
        e = rk - Jk @ d
        rises.append(float(e @ e - rk @ rk))
    return {"sum": float(sum(rises)), "each": rises}


def _read(problem, result, centre: str | None = None, terms=(), dictionary: bool = False) -> tuple:
    """The statistics, each a reading of the residuals alone, and the leading direction of each shape as the records
    see it, of length one. No list of terms of the package is read (the owner's messages 337 and 340): `terms` is the
    researcher's own, empty by default (342)."""
    lin = Linear(problem, result)
    out = {}
    rows = by_position(problem, result)
    for place, pick in PLACES:
        out[("where", place)] = float(sum(r["both_sides"] ** 2 for r in rows if pick(r)))
    out[("where", "the asymmetry")] = float(sum(r["one_side_less_the_other"] ** 2 for r in rows))
    out.update(changes(problem, result))
    read, leading = {}, None
    for kind, said in SHAPES_READ:
        side = one_side(lin) if kind == "periodic local" else None          # one side beyond the other, scanned once
        cand = side[0] if side is not None else (None if kind == "periodic local" else shape_candidate(lin, kind))
        entry = _one(lin, cand) if cand is not None else {"status": "none"}
        if entry["status"] == "ranked":
            if side is not None:
                entry = dict(entry, score=side[1])
            out[("shape", said)] = float(entry["score"])
            read[said] = entry["unit"]
            if leading is None or entry["score"] > leading["score"]:
                leading = entry
    r = np.asarray(lin.r, dtype=float)
    if leading is not None:
        r = r - leading["W"] @ (leading["W"].T @ r)          # the noise is read once the shape that leads is taken out: a shape the model lacks is not a law of the noise
    for said, value in noise_reader(lin)(r).items():          # the skewness and the tails of the residuals, and the glitches, the noise that grows and the correlated noise (the owner's message 266)
        out[("noise" if said in NOISE_SHAPES else "what", said)] = value
    for name in catalogue(centre):
        try:
            out[("seen by", name)] = float(family(lin, name, centre).reading()["chi_square"])
        except (ValueError, np.linalg.LinAlgError):
            pass          # a family that the package refuses on these records reads nothing here
    for said, value in beyond_reader(lin)(lin.r).items():          # the records together, which no family gives (the owner's message 269)
        out[("seen by", said)] = value
    loo = left_out(lin)          # the cross validation of the joint fit, row M590
    if loo:
        out[("seen by", LEFT_OUT)] = loo["sum"]
    for said, value in residual_plane(problem, result, lin).items():          # the chi functions and the windows on the residual's value (527)
        out[("seen by", said)] = value
    windowed = windowed_chi(lin)          # the chi functions in windows of position, row M593
    if windowed is not None:
        out[("seen by", WINDOWED)] = windowed["largest"]
    banded = value_windowed_chi(lin)          # the chi functions in windows of the residual's value, row M594
    if banded is not None:
        out[("seen by", WINDOWED_VALUE)] = banded["largest"]
    gridded = gridded_chi(lin)          # generalised griding of the data, row M588, the owner's message 572
    if gridded is not None:
        out[("seen by", GRIDDED_LEVEL)] = gridded["largest"]
    value_gridded = value_gridded_chi(lin)          # generalised griding of the residuals, row M588, the owner's message 572
    if value_gridded is not None:
        out[("seen by", GRIDDED_VALUE)] = value_gridded["largest"]
    refitted = chi_fits(lin)          # the fit under every chi function, the owner's messages 4, 25 and 573
    if refitted is not None:
        out[("seen by", CHI_FITS)] = refitted["largest"]
    summed = running_sum(lin)          # the residuals' running sum, the owner's message 575: the one transform that saw a smaller miss than the report's readings, a step at a quarter of its size where they needed half (runs/reports/transforms-596/README.md). The Fourier transform tied at best and the differences saw nothing, so they stand beside the report and not in its family, as `families.wide` says a reading joins
    if summed is not None:
        out[("seen by", RUNNING_SUM)] = summed["largest"]
    out.update(matched_terms(lin, terms, centre, dictionary))
    return out, read


from .families import NOISE_SHAPES, beyond_reader, catalogue, family, noise_reader, noise_shape  # noqa: E402,F401  the readings of the noise live with the families, so that no kept reading depends on this file through them (review 64, minor 2)


def matched_families(lin, candidates, centre: str | None = None, names=None) -> dict:
    """The owner's message 262: the clues of every family of statistics as the pillars of the decision of a turn. For each
    family of the catalogue but the fit of every sample, the fingerprint of the residuals, `Family.fingerprint`, and for
    each candidate the fingerprint of what the residuals keep of it, `Linear.visible`, both standardised in the
    covariance of the family's statistics under a right model. The fingerprints are joined over the families, and a
    candidate is scored by the part of the residuals' fingerprint that its own accounts for, the square of their product
    over the squared length of its own, with the amplitude that gives, and what each family reads once the candidate is
    taken out of the residuals at that amplitude. A family that the package refuses on the records is left out and
    named. A candidate that no family sees is left out. The score is a chi square of one degree of freedom under a right
    model only where the families' statistics are independent, which they are not: the twin gives its level, as it
    gives every level of the report. Returns the families read, the reading of each before any term, and the candidates
    ranked by their score, each with its amplitude and the reading of each family after it."""
    fams, left = {}, []
    for name in (list(names) if names is not None else [n for n in catalogue(centre) if n != "every sample: the fit"]):
        try:
            fams[name] = family(lin, name, centre)
        except (ValueError, np.linalg.LinAlgError) as err:
            left.append({"family": name, "why": str(err)})
    if not fams:
        raise ValueError("no family of the catalogue can be read on these records: " + "; ".join(f"{x['family']}: {x['why']}" for x in left))
    prints = {}
    for name, fam in list(fams.items()):          # a family read on its construction can still refuse its fingerprint, as one of more than LARGEST statistics on a hundred records: it is left out and named, and the match goes on (the code review of 0.1.0rc2, moderate 2)
        try:
            prints[name] = fam.fingerprint(lin.r)
        except (ValueError, np.linalg.LinAlgError) as err:
            left.append({"family": name, "why": str(err)})
            del fams[name]
    if not fams:
        raise ValueError("no family of the catalogue can be read on these records: " + "; ".join(f"{x['family']}: {x['why']}" for x in left))
    Z = np.concatenate(list(prints.values()))
    refused = []
    before = {name: float(p @ p) for name, p in prints.items()}
    ranked = []
    for cand in candidates:
        D = np.column_stack([lin.reduce(col) for col in cand.columns])
        V = np.column_stack([lin.ledger.visible(D[:, j]) for j in range(D.shape[1])])
        total = np.array([lin.norm2(col) for col in cand.columns])
        if total.max() == 0.0 or float(np.max(np.sqrt(np.maximum((V * V).sum(axis=0), 0.0) / np.maximum(total, 1e-300)))) < ABSORBED:          # absorbed whole: no statistic of these records sees it, as `propose._one` reads it
            continue
        T = np.column_stack([np.concatenate([fams[name].fingerprint(V[:, j]) for name in fams]) for j in range(V.shape[1])])
        G = T.T @ T
        w, U = np.linalg.eigh(G)
        keep = w > 1e-9 * max(float(w.max()), 1e-300)
        if not keep.any():
            continue
        amp = (U[:, keep] / w[keep]) @ (U[:, keep].T @ (T.T @ Z))          # the amplitude of each column that matches the fingerprint best
        if getattr(cand, "sign", 0) and amp.size == 1 and float(amp[0]) * cand.sign < 0.0:          # a width asked to be negative is no width, and its score ties with the kernel of the other kind (the owner's message 269)
            refused.append({"name": cand.name, "key": cand.key or cand.name, "amplitude": float(amp[0])})
            continue
        fitted = T @ amp
        score = float(fitted @ fitted)
        at, after = 0, {}
        for name, p in prints.items():
            m = p.size
            left_over = p - fitted[at:at + m]
            after[name] = float(left_over @ left_over)
            at += m
        ranked.append({"name": cand.name, "key": cand.key or cand.name, "family": cand.family, "parity": cand.parity, "score": score, "amplitude": amp.tolist(), "after": after})
    ranked.sort(key=lambda e: -e["score"])
    return {"families": list(fams), "left_out": left, "readings": before, "candidates": ranked, "of_the_sign_they_cannot_have": refused}


def shape_candidate(lin, kind: str):
    """The candidate of a shape of `SHAPES_READ` by its kind, or None where the records hold none: the report scores it,
    and the proposals measure it (row M589)."""
    if kind == "trend":
        return slow_trend(lin)
    if kind == "quadrature":
        return quadrature(lin)
    if kind == "localized":          # the wavelet scan, row M590
        got = wavelet_best(lin)
        return got[0] if got is not None else None
    if kind == "periodic local":
        got = one_side(lin)
        return got[0] if got is not None else None
    got = periodic_best(lin, frame="abscissa" if kind == "periodic fixed" else "centre") if kind in ("periodic", "periodic fixed") else None
    return got[0] if got is not None else None


def one_side(lin):
    """The periodic scan of each side of the centres, the left and the right, row M589: the candidate of the side that
    scores more, and its score less the best of the other side, so that a structure over both sides, which the scan of
    the whole span reads, is not read on one side. A ripple over every curve and a pedestal in both wings were named "a
    periodic structure on one side of the curves" (the plant trials of 18:37 of 2026-11-08). None where neither side
    reads one."""
    sides = sorted((got for got in (periodic_best(lin, side=side) for side in ("left", "right")) if got is not None), key=lambda got: got[1], reverse=True)
    if not sides:
        return None
    return sides[0][0], sides[0][1] - (sides[1][1] if len(sides) > 1 else 0.0)


def quadrature(lin):
    """A part in quadrature with the line, row M589 (the owner's message 533: "other transforms"): the Hilbert transform
    of each fitted curve, its level taken off, which turns a peak into an odd shape about it, as a phase of the
    detection adds to a line. Named by the transform and not by a term of a field, as the report is blind (rows M455
    and M466). Read on an even grid over the span of each curve and put back on its samples. One size for
    all the curves, each in proportion to its own fitted curve. None where a curve has fewer than eight samples."""
    from scipy.signal import hilbert
    pr = lin.problem
    columns = []
    for k, rec in enumerate(pr.records):
        x = np.asarray(rec.x, dtype=float)
        if x.size < 8:
            return None
        f = np.asarray(pr.predict(k, lin.fit.shared, lin.fit.own[k]), dtype=float)
        order = np.argsort(x)
        grid = np.linspace(x.min(), x.max(), x.size)
        g = np.interp(grid, x[order], f[order])
        g = g - np.median(g)          # the level off: the transform of a constant is no shape of the line
        columns.append(np.interp(x, grid, np.imag(hilbert(g))))
    return Candidate("a part in quadrature with the line, the Hilbert transform of the fitted curve", "quadrature", [columns], "any")


WAVELETS = (0.25, 0.5, 1.0, 2.0)          # the scales of the wavelet scan, in widths of the curves


def wavelet_best(lin, scales=WAVELETS, chunk: int = 64):
    """A structure in one place and at one scale, row M590 (the owner's messages 533 and 534): Morlet wavelets, a cosine
    and a sine of one cycle over two of its scales under a Gaussian of that scale, at every place half a scale apart over
    the span about the centre of each curve, one size for all the curves. Each place and scale is scored as
    `periodic_best` scores a period, on what the fit cannot absorb, all at once and with nothing else. The candidate and
    its score, or None where no wavelet is seen."""
    pr = lin.problem
    unit = float(np.median([lin.width_of(k) for k in range(len(pr.records))]))
    centre = [lin.centre_of(k) for k in range(len(pr.records))]
    lo = min(float(np.min(rec.x)) - centre[k] for k, rec in enumerate(pr.records))
    hi = max(float(np.max(rec.x)) - centre[k] for k, rec in enumerate(pr.records))
    grid = [(s, p) for s in scales for p in np.arange(lo, hi + 1e-9, 0.5 * s * unit)]
    best, at = 0.0, None
    for start in range(0, len(grid), chunk):
        part = grid[start:start + chunk]
        m = len(part)
        W, total = np.zeros((lin.n, 2 * m)), np.zeros(2 * m)
        for k, rec in enumerate(pr.records):
            sd, Q = lin.parts[k]
            u = np.asarray(rec.x, dtype=float)[:, None] - centre[k] - np.array([p for _, p in part])[None, :]
            sigma = unit * np.array([s for s, _ in part])[None, :]
            envelope = np.exp(-0.5 * (u / sigma) ** 2)
            block = np.hstack([envelope * np.cos(np.pi * u / sigma), envelope * np.sin(np.pi * u / sigma)]) / sd[:, None]
            total += (block * block).sum(axis=0)
            W[lin.rows_of(k)] = block - Q @ (Q.T @ block)          # freed of the own parameters, as `periodic_best` does
        W -= lin.ledger.Q @ (lin.ledger.Q.T @ W)          # and of the shared ones
        b = W.T @ lin.r
        cc, ss, cs = (W[:, :m] ** 2).sum(axis=0), (W[:, m:] ** 2).sum(axis=0), (W[:, :m] * W[:, m:]).sum(axis=0)
        seen = np.sqrt(np.maximum(np.stack([cc, ss]), 0.0) / np.maximum(np.stack([total[:m], total[m:]]), 1e-300)).max(axis=0) >= ABSORBED
        for j in np.flatnonzero(seen):
            w, U = np.linalg.eigh(np.array([[cc[j], cs[j]], [cs[j], ss[j]]]))
            keep = w > RANK * w.max()
            c = U[:, keep].T @ np.array([b[j], b[m + j]])
            score = float(np.sum(c * c / w[keep]))
            if score > best:
                best, at = score, part[j]
    if at is None:
        return None
    s, p = at
    cols = []
    for part_ in (np.cos, np.sin):
        cols.append([np.exp(-0.5 * ((np.asarray(rec.x, dtype=float) - centre[k] - p) / (s * unit)) ** 2) * part_(np.pi * (np.asarray(rec.x, dtype=float) - centre[k] - p) / (s * unit))
                     for k, rec in enumerate(pr.records)])
    return Candidate(f"a structure of {s:g} widths at {p / unit:+.3g} widths from the centre", "localized", cols, "any"), best


def periodic_best(lin, chunk: int = 64, frame: str = "centre", side: str | None = None):
    """The periodic candidate whose score is largest over a dense ladder: a cosine and a sine for every period from half
    the span of the records, `TREND_CYCLES`, to four of their steps, half a Fourier step apart, so that the nearest keeps a correlation of
    0.9 or more with any sine over the span. Each period is scored as `propose._one` scores it, the residuals on
    the part of its columns that the fit cannot absorb, all at once and with nothing else. The candidate and its score by
    the scan, or None when no period fits. `frame`: the phases about the centre of each curve, or about the abscissa
    itself, the same at equal abscissa in every curve, as an instrument's fringes are. `side`: the columns on the left or
    the right of each centre alone, a ripple on one side (row M589, the owner's message 533)."""
    pr = lin.problem
    centre = [lin.centre_of(k) if frame == "centre" else 0.0 for k in range(len(pr.records))]
    middle = [lin.centre_of(k) for k in range(len(pr.records))] if side is not None else None
    span = float(np.median([rec.x[-1] - rec.x[0] for rec in pr.records]))
    step = float(np.median([np.median(np.diff(rec.x)) for rec in pr.records]))
    cycles = np.arange(TREND_CYCLES, span / (4.0 * step) + 1e-9, 0.5)          # fewer cycles are a slow trend, which `slow_trend` reads: from one, a polynomial of the level was named a periodic structure of 1.5 cycles
    best, at = 0.0, None
    for lo in range(0, cycles.size, chunk):
        f = cycles[lo:lo + chunk] / span
        m = f.size
        W, total = np.zeros((lin.n, 2 * m)), np.zeros(2 * m)
        for k, rec in enumerate(pr.records):
            sd, Q = lin.parts[k]
            phase = 2.0 * np.pi * np.outer(np.asarray(rec.x, dtype=float) - centre[k], f)
            block = np.hstack([np.cos(phase), np.sin(phase)]) / sd[:, None]
            if side is not None:          # one side of the centre alone
                block = block * (((np.asarray(rec.x, dtype=float) < middle[k]) if side == "left" else (np.asarray(rec.x, dtype=float) > middle[k])).astype(float))[:, None]
            total += (block * block).sum(axis=0)
            W[lin.rows_of(k)] = block - Q @ (Q.T @ block)          # freed of the own parameters, as `Linear.reduce` does
        W -= lin.ledger.Q @ (lin.ledger.Q.T @ W)          # and of the shared ones, as `Ledger.visible` does: what the residuals keep
        b = W.T @ lin.r
        cc, ss, cs = (W[:, :m] ** 2).sum(axis=0), (W[:, m:] ** 2).sum(axis=0), (W[:, :m] * W[:, m:]).sum(axis=0)
        seen = np.sqrt(np.maximum(np.stack([cc, ss]), 0.0) / np.maximum(np.stack([total[:m], total[m:]]), 1e-300)).max(axis=0) >= ABSORBED
        for j in np.flatnonzero(seen):
            w, U = np.linalg.eigh(np.array([[cc[j], cs[j]], [cs[j], ss[j]]]))
            keep = w > RANK * w.max()
            c = U[:, keep].T @ np.array([b[j], b[m + j]])
            score = float(np.sum(c * c / w[keep]))
            if score > best:
                best, at = score, float(cycles[lo + j])
    if at is None:
        return None
    P = span / at
    keep = [np.ones(len(rec.x)) if side is None else ((np.asarray(rec.x, dtype=float) < middle[k]) if side == "left" else (np.asarray(rec.x, dtype=float) > middle[k])).astype(float)
            for k, rec in enumerate(pr.records)]
    where_said = {None: "", "left": " on the left of the curves", "right": " on the right of the curves"}[side] + (", fixed to the abscissa" if frame == "abscissa" else "")
    return Candidate(f"a periodic structure of {at:g} cycles over the span{where_said}", "periodic",
                     [[keep[k] * np.cos(2.0 * np.pi * (rec.x - centre[k]) / P) for k, rec in enumerate(pr.records)],
                      [keep[k] * np.sin(2.0 * np.pi * (rec.x - centre[k]) / P) for k, rec in enumerate(pr.records)]], "any"), best


def _one_set(args):
    twin, result, max_iter, centre = args[:4]
    terms = args[4] if len(args) > 4 else ()
    dictionary = args[5] if len(args) > 5 else False
    try:
        got = fit(twin, result.shared, result.own, max_iter=max_iter)
        return (statistics(twin, got, centre, terms, dictionary) if terms or dictionary else statistics(twin, got, centre)) if got.converged else None          # the list passed only where there is one
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return None


SEED = 23          # beside the digest of the records, the seed of the twin of a report


def seed_of(problem, salt: int = SEED) -> int:
    """The seed of the twin of a report, drawn from the records themselves: the digest of each record's abscissa, values
    and standard deviations, with `salt` beside it. The same records give the same twin, so a report is the same each time
    it is read, and two sets of records give two twins. A fixed seed gave every set of records the same vectors of noise, so a
    level read over many sets was the level of that one set of vectors: on the toy's right models the report alarmed on 0.188 of
    250 draws at 30 sets with the seed 23 and on 0.068 with a seed drawn for each set of records (his message 526, rows M564 to
    M576 of the analysis repository)."""
    from .seeds import seed_of as drawn          # one digest for every Monte Carlo of the package, row M582
    return drawn(problem, salt)


def twin_sets(problem, result, sets: int = 100, seed: int | None = None, centre: str | None = None, max_iter: int = 80, terms=(), dictionary: bool = False, planted=None) -> list:
    """The statistics of the report on `sets` synthetic sets of the fitted model, each fitted again as the records were.
    The sets whose fit did not end are left out. `planted`: one array a record added to the mean of every set. `seed`: the
    seed of the sets, by default drawn from the records, `seed_of`."""
    seed = seed_of(problem) if seed is None else seed
    return [s for s in (_one_set((twin, result, max_iter, centre, terms, dictionary)) for twin in replicas(problem, result, sets, seed=seed, planted=planted)) if s is not None]


def misfit(problem, result, span: int = SMOOTH) -> list:
    """What the fit leaves in each record, averaged over `span` neighbouring samples: the misfit that the records hold,
    with the noise averaged down, which the twin that carries the misfit adds to the mean of each of its sets."""
    out = []
    for k, rec in enumerate(problem.records):
        r = np.asarray(rec.y, dtype=float) - np.asarray(problem.predict(k, result.shared, result.own[k]), dtype=float)
        n = max(1, min(span, r.size)); n -= 1 - n % 2
        out.append(np.convolve(np.pad(r, n // 2, mode="edge"), np.ones(n) / n, mode="valid"))
    return out


def on_the_records(problem, result, answers: list, seed: int, centre, max_iter: int) -> None:
    """The scatter of the excess of each place beyond the twin as the records hold it, `scatter_on_the_records`, read on
    MISFIT_SETS sets of the twin that carries the misfit, each fitted again, in the units of the excess. A twin with no
    misfit scatters an excess by 0.3 to 0.6, while on draws of the toy with a term hidden the same excess scattered 2.5
    to 11 times more, and a tie read on the first named the flanks of the fifth hidden term with its core on 8 draws of 30 where the
    records do not tell them apart (his message 389). Read only where two places stand beyond the twin."""
    places = [a for a in answers if a["question"] == "where" and a.get("beyond_the_twin")]
    if len(places) < 2:
        return
    held = twin_sets(problem, result, MISFIT_SETS, seed + 1, centre, max_iter, planted=misfit(problem, result))
    for a in places:
        on = [t[("where", a["answer"])] for t in held if ("where", a["answer"]) in t]
        span = a["on_the_twin"]["95th percentile"] - a["on_the_twin"]["median"]
        a["scatter_on_the_records"] = float(np.std(on) / max(span, 1e-12)) if len(on) >= 10 else None


def where(problem, result, sets: int = 100, seed: int | None = None, fixed=(), centre: str | None = None, max_iter: int = 80, null: list | None = None, terms=(), dictionary: bool = False) -> dict:
    """The report of where to search: each answer with its statistic on the records, its level and its excess on
    `sets` synthetic sets of the fitted model, and whether it stands beyond the twin, first the smallest level and then
    the largest excess. `first_of_each_question` names the answer that leads each question. `fixed`: the covariates that
    the description says take one value in every record, named as what cannot be read. `null`: the statistics of the
    sets of a twin already read, by `twin_sets`, in place of new ones. `terms`: the researcher's list of terms, which
    adds the question matched, empty by default. Fewer than 20 sets that end are refused: a level and an excess need
    them."""
    seen, leading_shapes = _read(problem, result, centre, terms, dictionary)
    seed = seed_of(problem) if seed is None else seed          # the twin drawn from the records, `seed_of`
    null = twin_sets(problem, result, sets, seed, centre, max_iter, terms, dictionary) if null is None else null
    if len(null) < 20:
        raise ValueError(f"only {len(null)} synthetic sets of the fitted model ended: a level needs 20 at least")
    answers = []
    read = [key for key in seen if sum(key in t for t in null) >= 0.9 * len(null)]          # a family refused on more than a tenth of the sets of the twin is read on none
    for key in read:
        value, on_twin = seen[key], [t[key] for t in null if key in t]
        typical, high = np.quantile(on_twin, [0.5, 0.95])
        answers.append({"question": key[0], "answer": key[1], "statistic": value, "level": float((1 + sum(v >= value for v in on_twin)) / (len(on_twin) + 1)),
                        "excess": _excess(value, on_twin), "scatter": _scatter(on_twin), "on_the_twin": {"median": float(typical), "95th percentile": float(high)}})
        if key[1] in leading_shapes:
            answers[-1]["not_told_apart_from"] = [other for other, unit in leading_shapes.items() if other != key[1] and abs(float(unit @ leading_shapes[key[1]])) >= SIMILAR]
    largest = [max(_excess(t[key], [u[key] for u in null if u is not t and key in u]) for key in read if key in t) for t in null]
    beyond = float(np.quantile(largest, 0.95))
    for a in answers:
        a["beyond_the_twin"] = a["excess"] > beyond
    on_the_records(problem, result, answers, seed, centre, max_iter)
    answers.sort(key=lambda a: (a["question"] == "noise", a["level"], -a["excess"]))          # the shape of the noise after the answers on the mean
    first = {q: next((a for a in answers if a["question"] == q and a["beyond_the_twin"]), None) for q in ("where", "shape", "what", "changes along", "seen by", "noise") + (("matched",) if terms or dictionary else ())}
    along = None
    if first["where"] is None and first["changes along"] is not None:          # no place beyond the twin alone, and one that changes along a setting: its place, read along it. The report named no place on 26 draws of 30 with a width that grows along a setting hidden, which it read along that setting (row M485)
        along = first["changes along"]["answer"]
        place = along.split(" along ", 1)[0]
        first["where"] = next((a for a in answers if a["question"] == "where" and a["answer"] == place), {"question": "where", "answer": place})
    varies = {key for key in covariates_of(problem) if any(k[0] == "changes along" and k[1].endswith(f" along {key}") for k in seen)}
    return {"answers": answers, "first_of_each_question": {q: None if a is None else a["answer"] for q, a in first.items()}, "where_read_along": along,
            "tied_with_the_first": {"where": tied_with(first["where"], answers)},
            "sets": len(null), "largest_excess_on_the_twin": beyond, "mode": "default" if not terms and not dictionary else ", ".join(([f"the researcher's list of {len(terms)} term{'s' if len(terms) != 1 else ''}"] if terms else []) + (["the package's dictionary, asked by its flag"] if dictionary else [])),
            "cannot_be_read": [f"along {name}: it takes one value in every record" for name in fixed]
            + [f"along {name}: no place of the residuals is read along it" for name in covariates_of(problem) if name not in varies and name not in fixed]}


def _scatter(on_twin) -> float:
    """The standard deviation of an excess on the sets of the twin, in the units of the excess: one in about 1.6 for a
    statistic near a normal law."""
    typical, high = np.quantile(on_twin, [0.5, 0.95])
    return float(np.std(on_twin) / max(high - typical, 1e-12))


def tied_with(lead, answers: list) -> list:
    """The places that the records do not tell from the first, `lead`: beyond the twin at its level, their excess within
    TIE standard deviations of the difference of the two as the records scatter them, `on_the_records`. On traces with no
    noise the misfit of a hidden term stands at 21.5 in the flanks and 21.1 in the core, and the report put the core first on 24
    draws of 30 where the flanks are the place of its physics: the two are named together, his decisions of messages 388
    and 389. A place whose scatter on the records was not read ties with none."""
    if not lead or lead.get("scatter_on_the_records") is None:
        return []
    return [a["answer"] for a in answers if a is not lead and a["question"] == lead["question"] and a.get("beyond_the_twin") and a["level"] == lead["level"]
            and a.get("scatter_on_the_records") is not None and abs(a["excess"] - lead["excess"]) <= TIE * float(np.hypot(a["scatter_on_the_records"], lead["scatter_on_the_records"]))]


def _excess(value: float, on_twin) -> float:
    """How far beyond the sets of the twin a statistic stands: one at their 95th percentile, zero at their median."""
    typical, high = np.quantile(on_twin, [0.5, 0.95])
    return float((value - typical) / max(high - typical, 1e-12))
