"""Terms that models of recorded curves often lack, each as the working model with one shared parameter more.

A candidate (candidates.py) is a discrepancy of the mean at the fitted values: it is ranked with no fit,
at first order. A mechanism is the same guess made into a model. The working model is wrapped, whatever its formula, and
gains one shared parameter, so the guess can be fitted, compared by a criterion and put into synthetic records. First
order is not enough for a term that is not small: on a curve of width 3.7, a kernel of width 1.4 moved an estimate by
1.3 bars at first order and by 8.7 once fitted.

    widened(size)            the curve convolved with a two sided exponential kernel of scale size(design) * sqrt(a), size a
                             number or a function of the settings of a curve, as lambda design: 0.5
    displaced(reach)         a mixture of copies of the curve displaced by reach(design) * sqrt(a) * v, v from 0 to 1 of
                             density 2 v, reach a number or a function of the settings, as for widened. About its mean,
                             its long side points toward smaller abscissa where reach is positive, toward larger where
                             it is negative: the opposite of second, whose positive distance is after the curve
    along(parameter, key, power=1)
                             a shared parameter, named as SHARED names it, that changes along a setting of the description,
                             named as its [covariates] names it: p (1 + a z), z the setting less its mean over the
                             records, over its standard deviation over them, taken with their number, not one less.
                             power=2 takes z squared less its mean
    gaussian_width(along=None)
                             the curve convolved with a Gaussian of standard deviation a, in the unit of the abscissa, on
                             the steps of the records, which need equal steps: a width that the line lacks, as the noise
                             of a source or the response of an instrument. along, the name of a setting or a function of the settings of a curve, scales a in it
    lorentzian_width(along=None)
                             the same with a Lorentzian of full width a at half its height: a width of a lifetime, and
                             along="concentration" makes it in proportion to that setting over its mean over the
                             records, a the width at that mean, and names it
    pedestal(stretch)        a copy of the curve `stretch` times wider under it, of relative area a, 0 or more unless
                             `lower` says otherwise, stretch a number
    bent()                   the abscissa is not linear: x + a (x - middle)^2 / half its span
    second(distance)         a second copy of the curve, `distance` away along the abscissa, in its unit, after the curve
                             where it is positive and before it where negative, of relative area a, 0 or more unless
                             `lower` says otherwise
    tail(stretch, side)      a copy of the curve `stretch` times wider on one side of its centre alone, of relative area a,
                             side "upper" or "lower". A wider copy, not the exponential tailing of a column or of a
                             detector, which a Mechanism of your own writes
    level(power)             a power of the abscissa under every curve, of one size for all: a ((x - middle) / half)^power,
                             the middle and the half span those of the abscissas of all the records together, power 1, 2
                             or 3 for a slope, a curvature or a cubic term of the baseline
    periodic(period, part)   a cosine or a sine of the abscissa of one period, in its unit, part "cosine" by default or
                             "sine": a cos or a sin(2 pi (x - centre) / period), about the centre of each curve. A ripple
                             fixed to the abscissa, the same on every curve whatever its centre, as the fringe of an
                             etalon, is a Mechanism of your own with x in the place of x - centre:
                             section 11 of the tutorial writes one

Glitches and a change of the noise are no term of the mean: they are read in the weights of the fit, the noise of each
sample that the description gives. Each term here is made of copies of the curve, wider, displaced or on one side: a shape
that the physics names otherwise, as the exponential tail of a detector or of a column, is a Mechanism of your own. The report
gives the place, the shape and the setting of a misfit, and never the name of a term: which term of physical origin fills it
is the researcher's to say. A Mechanism is Mechanism(name, kind, parameter, start, wrap,
lower=-inf, upper=inf): kind is a word that the package keeps with the term and reads in no result, lower and upper bound
the new parameter, and wrap(model, problem, n) returns any curve, the line replaced as well as a term added to it.

    m = widened(size=lambda design: 1.03, ground="the theory of the crossing time")
    larger = m.extend(problem, result)        the next working model, started at the fit, the new parameter at m.start
    held = m.hold(problem)                    the working model with the term held at m.start: no parameter more

`size` and `reach` say how the term scales across the design, as the theory of the problem or an educated guess gives
it. For these two terms the new parameter a is the variance that the term adds over the variance that was given: 1 means
the size that was given, 0 means no term, and the curve changes in proportion to a near zero. A parameter that
multiplied the scale itself would leave the curve unchanged at first order near zero, and a search on records without the
term could not end: the first form of this file did that, and its first run on such records showed it. A mechanism with a `ground`, a sentence of
the theory that it rests on, is added whatever a criterion says of it: refine.py holds that rule.

What a wrapped model costs: `displaced` calls the working model 8 times, `pedestal` twice, the others once, and `widened`
needs an abscissa of equal steps.

What each does to the offset and the slope of a record, which are part of the model that it wraps. `widened` leaves
both as they are. `displaced` leaves both as they are, since its mixture is taken about its mean. `second`, `pedestal`
and `tail` take out of their copy the straight line that the model draws far from its curve, ten spans of the record
away, so that they add the shape of the curve and no share of its offset or slope. Their first forms carried that share:
the record's own offset took it and read low by the size of the term, 8 per cent under a second component of 0.08, and
on a curve that stands on a level of 1, as a normalised flux, the size and the offset moved the curve alike and the
search stopped (round 3 of the cold tests). A baseline that bends is drawn by that line as a chord, and what the chord
misses is a constant under the copy, which the record's own offset takes. `bent` bends the slope into a parabola of
curvature a times the slope over half the span, which no own parameter takes: on a record with a slope, part of what
`bent` gains is that.

Limits. The eight nodes of `displaced` give the mixture to one part in 1e13 of the height of a curve while the reach is
under half its width, to 3 parts in 1e5 at twice its width and to 2 per cent at five times. The parameter is bounded
at ten times the `reach` that the caller gives: a ratio, so it keeps a search inside those widths only when the reach
is under about a fifth of a width of the curve. Give a reach of that size, or read the fitted reach in widths. The cost of `widened` grows with its scale, so its parameter is bounded too.
`bent` folds the abscissa when its parameter passes a half, and is bounded at a quarter.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .problem import Problem

_g, _w = np.polynomial.legendre.leggauss(8)
V, WV = 0.5 * (_g + 1.0), 0.5 * _w * 2.0 * 0.5 * (_g + 1.0)        # the fraction v on (0, 1) and its weight 2 v dv
WV = WV / WV.sum()
REACH = 25.0                                                      # a kernel is carried this many scales to each side
FINE = 64.0                                                       # a kernel is applied on steps of at most this share of the width of a curve
WIDEST = 100.0                                                    # the largest variance, over the one that was given, that a search may give a kernel or a displacement


@dataclass
class Mechanism:
    name: str
    family: str
    parameter: str                       # the name of the shared parameter that it adds
    start: float
    wrap: Callable                       # wrap(model, problem, n) -> the model with the new parameter at place n of the shared ones
    lower: float = -np.inf
    upper: float = np.inf
    ground: str = ""                     # the sentence of the theory that the term rests on. Empty: an educated guess
    theory: tuple | None = None          # (centre, standard deviation) of the new parameter where the theory gives them
    formula: str = ""                    # what the term does to the model f before it, as mathtext, a its parameter: the equation of each turn in the visual log, his message 309

    def extend(self, problem: Problem, result=None) -> Problem:
        """The working model of `problem` with this term: one shared parameter more, started at the fit if one is given."""
        name, k = self.parameter, 1          # his messages 283 and 284: several terms make one kind of width, two mechanisms that broaden a line each with a law of its own, and a further term with a physical ground is not refused a priori: it takes the next free name and the records judge it, or its theory holds it where they cannot
        while name in problem.shared:
            k += 1
            name = f"{self.parameter} {k}"
        n = len(problem.shared)
        lo, hi = problem.bounds_shared
        start = list(problem.start_shared if result is None else result.shared) + [self.start]
        own = problem._start_own if result is None else [np.array(row, dtype=float) for row in result.own]
        centre = None if problem.centre is None else (lambda shared, own_k, design, c=problem.centre: c(shared[:n], own_k, design))
        return Problem(problem.records, self.wrap(problem.model, problem, n), problem.shared + (name,), problem.own, start, own, centre,
                       f"{problem.name} with {self.name}", (list(lo) + [self.lower], list(hi) + [self.upper]), problem._bounds_own)


    def hold(self, problem: Problem, value: float | None = None) -> Problem:
        """The working model with this term held at one value, by default the one that it starts at: no parameter is
        added. For a term that stands on theory and that the records cannot fix, such as a component that the offset of
        every record takes: it enters at the size that the theory gives, and its bar is carried (refine.carried)."""
        n = len(problem.shared)
        v = float(self.start if value is None else value)
        wrapped = self.wrap(problem.model, problem, n)

        def model(x, shared, own, design):
            return wrapped(x, np.append(np.asarray(shared, dtype=float), v), own, design)
        return Problem(problem.records, model, problem.shared, problem.own, problem.start_shared, problem._start_own, problem.centre,
                       f"{problem.name} with {self.name} held at {v:g}", problem._bounds_shared, problem._bounds_own)


def _line_under(model, shared, own, design, c: float, span: float):
    """The straight line that the model draws far from its curve, read ten spans of the record away on each side of
    `c`: the offset and the slope that a copy of the curve would carry. Where the model gives no finite value there, the
    line is read at the ends of the record, `span` apart about `c`."""
    for reach in (10.0, 0.5):
        far = c + reach * span * np.array([-1.0, 1.0])
        v = np.asarray(model(far, shared, own, design), dtype=float)
        if np.all(np.isfinite(v)):
            slope = (v[1] - v[0]) / (far[1] - far[0])
            return lambda t, v0=v[0], f0=far[0], s=slope: v0 + s * (np.asarray(t, dtype=float) - f0)
    raise ValueError("the model gives no finite value at the ends of the record: its baseline cannot be read")


def _rising(problem: Problem) -> None:
    """A kernel needs an abscissa that rises from point to point, and no more: on steps that are not equal the kernels of
    this file are applied on a grid of equal steps and read back at the record's own points. Records of a scan calibrated
    against a ruler stand on a step that drifts by a per cent or so along the scan, two steps from equal ones at the most,
    and a kernel that asked equal steps refused every one of them (the owner's messages 294 to 296)."""
    for rec in problem.records:
        if rec.x.size > 1 and not np.all(np.diff(rec.x) > 0.0):
            raise ValueError(f"record {rec.name!r}: a kernel is applied on an abscissa that rises from point to point, and this one does not: give its points in the order of the abscissa")


def _exponential_kernel(b2: float, step: float) -> np.ndarray:
    """Weights on a grid of spacing `step` for the two sided exponential kernel whose scale squared is b2: the kernel over
    each step, then the three weights about zero moved so that the variance of the weights is 2 b2 exactly. A kernel far
    narrower than the step then still adds its variance, as a second difference, and the weights are a smooth function of
    b2 through zero: with b2 negative they take that variance away, which no kernel does and a derivative at zero needs."""
    if b2 > 0.0:
        b = np.sqrt(b2)
        m = max(1, int(np.ceil(REACH * b / step)))
        edge = np.exp(-(np.arange(0, m + 1) + 0.5) * step / b)                  # exp(-(j + 1/2) step / b) for j = 0 to m
        half = 0.5 * (np.concatenate([[edge[0]], edge[:-1]]) - edge)            # the kernel over the step at distance j, for j from 1
        half[0] = 1.0 - edge[0]                                                # and over the step about zero
        half[0] += 1.0 - (half[0] + 2.0 * half[1:].sum())                      # what was cut beyond the reach goes back to the middle
    else:
        half = np.array([1.0, 0.0])
    lack = (2.0 * b2 - 2.0 * float(np.sum((np.arange(half.size) * step) ** 2 * half))) / (2.0 * step * step)
    half[0] -= 2.0 * lack
    half[1] += lack
    return np.concatenate([half[:0:-1], half])


def _on_a_grid(x, step: float, curve) -> np.ndarray:
    """A curve that a kernel makes on equal steps, read at points `x` that are not on equal steps: `curve` on a grid of
    `step` that covers x with two points of room on each side, then the cubic of `_cubic` at x."""
    first, last = int(np.floor(float(x.min()) / step)) - 2, int(np.ceil(float(x.max()) / step)) + 2
    grid = step * np.arange(first, last + 1)
    return _cubic(step * first, step, np.asarray(curve(grid), dtype=float), x)


def _cubic(first: float, step: float, values: np.ndarray, x) -> np.ndarray:
    """The values at x of a curve given on equal steps from `first`: the cubic through the two points on each side of x
    (Catmull and Rom). The curve and its slope are then continuous in x, which a search along a parameter that moves x
    needs. A line between neighbours has a corner at every point of the grid, and such a search stopped at one."""
    u = (np.asarray(x, dtype=float) - first) / step
    i = np.clip(np.floor(u).astype(int), 1, values.size - 3)
    t = u - i
    p0, p1, p2, p3 = values[i - 1], values[i], values[i + 1], values[i + 2]
    return p1 + 0.5 * t * (p2 - p0 + t * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3 + t * (3.0 * (p1 - p2) + p3 - p0)))


def widened(size: Callable = lambda design: 1.0, name: str = "a kernel with exponential wings", parameter: str = "kernel variance", start: float = 1.0,
            ground: str = "", theory: tuple | None = None) -> Mechanism:
    """The curve convolved with exp(-|u| / b) / (2 b), with b^2 = a * size(design)^2: each curve gains a variance of
    2 b^2 and wings that fall as an exponential. The parameter a is a variance over the variance that `size` gives, so
    a = 1 is the size that was given, the curve changes in proportion to a near zero, and a search can end at a = 0, where
    the term is absent. A scale that comes back as a times `size` squared has the ratio sqrt(a) to that size.
    The model is read on steps of at most a 64th of the width of the tallest curve. A term that bends the abscissa can
    stand outside this one: the curve is then read on a grid of equal steps that does not move with the bend, a
    fraction of the step of the tallest record, and taken back to the bent points by a cubic whose slope is continuous. Until candidate 45
    it was taken back by a line between neighbours on a grid that moved with the bend: the chi square then had a corner
    at every point of that grid, and a search along the bend stopped at one, 23 steps into a fit of the tutorial."""
    if not callable(size):          # a number is taken as the size of every curve: widened(1.0) left a trace of the code (the cold tests of row M169, round 2)
        _size = float(size)
        size = lambda design, _size=_size: _size          # noqa: E731

    def wrap(model, problem, n):
        _rising(problem)
        tall = max(problem.records, key=lambda rec: float(np.max(rec.y) - np.median(rec.y)))      # the record whose curve stands highest gives the width of a curve
        up = tall.y - np.median(tall.y)
        base = float(np.median(np.diff(tall.x)))                                                   # the step of the records, the median where they are not equal: an uneven abscissa is read on a grid of this step
        width = max(float(np.sum(up > 0.5 * up.max()) * base), base)
        fine = max(1, int(np.ceil(FINE * base / width)))                                           # the model is read on steps of at most a FINE-th of that width

        def new(x, shared, own, design):
            b2 = float(shared[n]) * float(size(design)) ** 2
            x = np.asarray(x, dtype=float)
            gaps = np.diff(x)
            even = bool(np.all(np.abs(gaps - gaps[0]) <= 1e-9 * abs(gaps[0])))
            # On equal steps the curve is read on a finer grid of the same points. An abscissa that a term outside this one
            # has bent or stretched unevenly is read on a grid of equal steps that does not depend on it, and the result
            # is taken back to its own points: until candidate 44 the first step was taken for all of them
            if even:
                low, step = float(x[0]), float(gaps[0]) / fine
                kernel = _exponential_kernel(b2, step)
                m = (kernel.size - 1) // 2
                wide = low + step * np.arange(-m, (x.size - 1) * fine + 1 + m)
                return np.convolve(np.asarray(model(wide, shared[:n], own, design), dtype=float), kernel, mode="valid")[::fine]
            step = base / fine                                            # the curve is as smooth on the bent points as on the others: the step of the records reads it
            kernel = _exponential_kernel(b2, step)
            m = (kernel.size - 1) // 2
            first, last = int(np.floor(float(x.min()) / step)) - 2, int(np.ceil(float(x.max()) / step)) + 2       # two points of room each side, for the cubic
            wide = step * np.arange(first - m, last + m + 1)
            smooth = np.convolve(np.asarray(model(wide, shared[:n], own, design), dtype=float), kernel, mode="valid")
            return _cubic(step * first, step, smooth, x)
        return new
    return Mechanism(name, "kernel", parameter, start, wrap, lower=0.0, upper=WIDEST, ground=ground, theory=theory,
                     formula=r"$f \ast \frac{e^{-|u|/b}}{2b},\ \ b^2 = a\,s^2$")


def displaced(reach: Callable = lambda design: 1.0, name: str = "a displacement of linear density", parameter: str = "depth squared", start: float = 1.0,
              ground: str = "", theory: tuple | None = None) -> Mechanism:
    """A mixture of copies of the curve displaced toward larger abscissa by d v, with v between 0 and 1 of density 2 v and
    d^2 = a * reach(design)^2. The mixture is taken about its mean, so the centre of a record stays the mean position of
    its curve and need not move when the term enters: the place where nothing is displaced is 2 d / 3 below that centre.
    What the term adds to a curve is a variance of d^2 / 18, in proportion to a, and a lean toward larger abscissa. The
    parameter a is that variance over the one that `reach` gives: a = 1 is the reach that was given, a = 0 no term."""
    if not callable(reach):          # a number is taken as the reach of every curve: widened(1.0) left a trace of the code (the cold tests of row M169, round 2)
        _reach = float(reach)
        reach = lambda design, _reach=_reach: _reach          # noqa: E731

    def wrap(model, problem, n):
        def new(x, shared, own, design):
            d = np.sqrt(max(float(shared[n]), 0.0)) * float(reach(design))
            x = np.asarray(x, dtype=float)
            if d == 0.0:
                return np.asarray(model(x, shared[:n], own, design), dtype=float)
            return sum(w * np.asarray(model(x - d * (v - 2.0 / 3.0), shared[:n], own, design), dtype=float) for v, w in zip(V, WV))
        return new
    return Mechanism(name, "displacement", parameter, start, wrap, lower=0.0, upper=WIDEST, ground=ground, theory=theory,
                     formula=r"$\int_0^1 2v\, f(x - d\,(v - 2/3))\, dv,\ \ d^2 = a\,r^2$")


def along(parameter: str, key: str, power: int = 1, ground: str = "", theory: tuple | None = None) -> Mechanism:
    """A shared parameter that is not the same in every record: p (1 + a z), z the covariate `key` of the design, centred
    and in units of its standard deviation over the records. With power 2, the square of z less its mean."""
    def wrap(model, problem, n):
        if parameter not in problem.shared[:n]:
            raise ValueError(f"along: {parameter!r} is not a shared parameter of the model, which are {', '.join(problem.shared[:n]) or 'none'}")
        missing = [rec.name for rec in problem.records if key not in rec.design]
        if missing:
            raise ValueError(f"along: the records carry no setting {key!r}: the settings of the first are {', '.join(problem.records[0].design) or 'none'}")
        j = problem.shared.index(parameter)
        z = np.array([rec.design[key] for rec in problem.records], dtype=float)
        mean, sd = float(z.mean()), float(z.std())
        if sd == 0.0:
            raise ValueError(f"the covariate {key} is the same in every record")
        offset = float(np.mean(((z - mean) / sd) ** 2)) if power == 2 else 0.0

        def new(x, shared, own, design):
            zc = ((float(design[key]) - mean) / sd) ** power - offset
            s = np.array(shared[:n], dtype=float)
            s[j] = s[j] * (1.0 + float(shared[n]) * zc)
            return model(x, s, own, design)
        return new
    square = " the square of" if power == 2 else ""
    return Mechanism(f"{parameter} changes along{square} {key}", "grows with a covariate", f"{parameter} along {key}" + ("^2" if power == 2 else ""), 0.0, wrap, ground=ground, theory=theory,
                     formula=r"$p \rightarrow p\,(1 + a\,z)$" if power == 1 else r"$p \rightarrow p\,(1 + a\,(z^2 - \overline{z^2}))$")


def pedestal(stretch: float = 8.0, ground: str = "", theory: tuple | None = None, lower: float = 0.0, upper: float = np.inf) -> Mechanism:
    """A copy of the curve `stretch` times wider about the centre of the record, of area a times the area of the curve.
    `lower` and `upper` bound that area, at zero below unless `lower` says otherwise: a pedestal adds to the curve. Left
    free, a fit returns a negative pedestal where the working model draws its wings too high, and keeps it: on eight
    curves drawn with a kernel of exponential wings and fitted with a Gaussian in its place, a pedestal 3 times wider was
    kept at -0.1347 +- 0.0082 and moved the width of the curve from 12 to 21 of its bars from its truth, where bounded at
    zero it is refused at its bound (the blind red team of 0.1.0rc19, mislead 6). A dip that takes away from the curve
    is written with `lower` below zero."""
    def wrap(model, problem, n):
        if problem.centre is None:
            raise ValueError("a pedestal is placed about the centre of a curve: the problem gives none")

        def new(x, shared, own, design):
            x = np.asarray(x, dtype=float)
            c = float(problem.centre(shared[:n], own, design))
            f = np.asarray(model(x, shared[:n], own, design), dtype=float)
            under = _line_under(model, shared[:n], own, design, c, float(np.ptp(x)) or 1.0)
            t = c + (x - c) / stretch
            return f + float(shared[n]) * (np.asarray(model(t, shared[:n], own, design), dtype=float) - under(t)) / stretch
        return new
    return Mechanism(f"a pedestal {stretch:g} times wider than the curve", "pedestal", f"pedestal {stretch:g}", 0.0, wrap, lower=lower, upper=upper, ground=ground, theory=theory,
                     formula=r"$f(x) + \frac{a}{s}\,[\,f(c + (x - c)/s) - \ell\,],\ \ s = %g$" % stretch)


def _width(shape: str, along=None, name: str | None = None, parameter: str | None = None, ground: str = "", theory: tuple | None = None, start: float | None = None) -> Mechanism:
    """A width of one value for every curve, of the shape `shape`, a Gaussian of standard deviation s, a Lorentzian of
    full width g, or an exponential kernel exp(-|u| / b) / (2 b) of decay length b, convolved with the curve on the steps
    of its records, each step of the kernel its area, so that a width under the step still adds its area right.
    `along(design)` scales the width of each curve, or `along` the name of a setting, which the name of the term then says,
    the width then in proportion to that setting over its mean over the records, so that its parameter is the width at
    that mean, in the units of the abscissa: the raw setting, a density of 3.6e13, took the start of the width to a size
    of the order of 1e12 and refused the term at its bound (the blind red team of 0.1.0rc18, block 1). The records need equal
    steps. `start`, in the units of the abscissa, is where a search begins, 0.3 for a Gaussian and an exponential and 0.1
    for a Lorentzian unless given: the change of the curve under a Gaussian or an exponential kernel grows as the square
    of its width, since their variances are s^2 and 2 b^2, so a search started at zero has no slope and ends there,
    refused at its bound. The owner's message 262: a blind run on a line that lacked a Gaussian width took it up in
    another term for want of the candidate. The exponential kernel is the transform of a coherence that falls as
    1 / (1 + (2 pi b t)^2), as that of an emitter crossing a narrow region does, and the owner's message 336 put it beside
    the other two."""
    from math import erf
    from scipy.signal import fftconvolve
    gauss, expo = shape == "gaussian", shape == "exponential"
    key = along if isinstance(along, str) else None          # a setting by its name, which the name of the term then says: two widths of one kind, along two settings, read as one term (the fresh red team of 0.1.0rc16, mislead 3)

    def wrap(model, problem, n):
        _rising(problem)
        scale = along
        if key is not None:
            missing = [rec.name for rec in problem.records if key not in rec.design]
            if missing:
                raise ValueError(f"the width along {key}: the records carry no setting {key!r}: the settings of the first are {', '.join(problem.records[0].design) or 'none'}")
            mean = float(np.mean([float(rec.design[key]) for rec in problem.records]))
            if mean <= 0.0:
                raise ValueError(f"the width along {key}: the mean of the setting over the records is {mean:.4g}, and a width in proportion to it needs it above zero")
            scale = lambda design: float(design[key]) / mean          # noqa: E731

        def new(x, shared, own, design):
            size, x = float(shared[n]) * (1.0 if scale is None else float(scale(design))), np.asarray(x, dtype=float)
            if size <= 0.0:
                return np.asarray(model(x, shared[:n], own, design), dtype=float)
            gaps = np.diff(x)
            if not np.all(np.abs(gaps - gaps[0]) <= 1e-9 * abs(gaps[0])):          # steps that are not equal: the kernel on a grid of the median step, read back at the record's points by the cubic
                return _on_a_grid(x, float(np.median(gaps)), lambda grid: new(grid, shared, own, design))
            step = float(x[1] - x[0])
            m = min(int(np.ceil((6.0 if gauss else 36.0) * size / step)), x.size) if gauss or expo else x.size          # a Gaussian ends at six deviations, an exponential kernel at 36 decay lengths, under 1e-15 of its area; a Lorentzian's wings carry area over the whole record, so its kernel reaches as far as the record: a reach fixed in the unit of one problem cut them. No kernel reaches past the record: a search that wandered to a Gaussian far wider than it built arrays without end, 6.4 GB on one draw of the count of Mercury on 2026-10-20
            edges = step * (np.arange(-m, m + 2) - 0.5)
            cdf = (np.array([0.5 * (1.0 + erf(e / (size * np.sqrt(2.0)))) for e in edges]) if gauss
                   else np.where(edges < 0.0, 0.5 * np.exp(np.minimum(edges, 0.0) / size), 1.0 - 0.5 * np.exp(-np.maximum(edges, 0.0) / size)) if expo
                   else np.arctan(edges / (0.5 * size)) / np.pi)
            kernel = np.diff(cdf)
            kernel /= kernel.sum()
            wide = float(x[0]) + step * np.arange(-m, x.size + m)
            return fftconvolve(np.asarray(model(wide, shared[:n], own, design), dtype=float), kernel, mode="valid")
        return new
    kind = "Gaussian" if gauss else "exponential" if expo else "Lorentzian"
    return Mechanism(name or f"a{'n' if expo else ''} {kind} width beside the line" + (f", in proportion to {key}" if key else ""), "kernel",
                     parameter or ("gaussian standard deviation" if gauss else "exponential decay length" if expo else "lorentzian full width") + (f" at the mean of {key}" if key else ""),
                     float((0.3 if gauss or expo else 0.1) if start is None else start), wrap, lower=0.0, upper=np.inf, ground=ground, theory=theory,
                     formula=(r"$f \ast G_{\sigma},\ \ \sigma = a$" if gauss else r"$f \ast E_{b},\ \ b = a$" if expo else r"$f \ast L_{\Gamma},\ \ \Gamma = a$").replace("= a$", "= a$" if along is None else r"= a\,q$" if key else r"= a\,w(z)$"))


def gaussian_width(along=None, name: str | None = None, parameter: str | None = None, ground: str = "", theory: tuple | None = None, start: float | None = None) -> Mechanism:
    """A Gaussian width beside the line, its parameter the standard deviation in the units of the abscissa: see `_width`."""
    return _width("gaussian", along, name, parameter, ground, theory, start)


def lorentzian_width(along=None, name: str | None = None, parameter: str | None = None, ground: str = "", theory: tuple | None = None, start: float | None = None) -> Mechanism:
    """A Lorentzian width beside the line, its parameter the full width at half maximum: see `_width`."""
    return _width("lorentzian", along, name, parameter, ground, theory, start)


def exponential_width(along=None, name: str | None = None, parameter: str | None = None, ground: str = "", theory: tuple | None = None, start: float | None = None) -> Mechanism:
    """An exponential width beside the line, the curve convolved with exp(-|u| / b) / (2 b), its parameter the decay length
    b in the units of the abscissa: see `_width`. The transform of a coherence that falls as 1 / (1 + (2 pi b t)^2) (the
    owner's message 336)."""
    return _width("exponential", along, name, parameter, ground, theory, start)


def bent(ground: str = "", theory: tuple | None = None) -> Mechanism:
    """The abscissa is not linear: the model is read at x + a (x - middle)^2 / half, with the middle and the half span of
    the abscissa of the first record. The curve leans, to one side at its low end and to the same side at its high end."""
    def wrap(model, problem, n):
        x0 = problem.records[0].x
        mid, half = 0.5 * float(x0[0] + x0[-1]), 0.5 * float(x0[-1] - x0[0])

        def new(x, shared, own, design):
            x = np.asarray(x, dtype=float)
            return model(x + float(shared[n]) * (x - mid) ** 2 / half, shared[:n], own, design)
        return new
    return Mechanism("a quadratic map of the abscissa", "abscissa", "bend", 0.0, wrap, lower=-0.25, upper=0.25, ground=ground, theory=theory,
                     formula=r"$f(x + a\,(x - m)^2/h)$")


def second(distance: float, ground: str = "", theory: tuple | None = None, lower: float = 0.0, upper: float = np.inf) -> Mechanism:
    """A second copy of the curve, `distance` away along the abscissa, in the units of the abscissa, of a times the size
    of the curve: an impurity that elutes later, a line of another species or isotope beside the first. The copy is the
    curve less the straight line that the model draws far from it, so that it carries no share of the offset or the slope
    of the record. A negative distance places it before. Its size is bounded at zero unless `lower` says otherwise: a
    second component adds to the curve, as a pedestal does (the blind red team of 0.1.0rc19, mislead 6)."""
    def wrap(model, problem, n):
        def new(x, shared, own, design):
            x = np.asarray(x, dtype=float)
            c = float(problem.centre(shared[:n], own, design)) if problem.centre is not None else 0.5 * float(x[0] + x[-1])
            under = _line_under(model, shared[:n], own, design, c, float(np.ptp(x)) or 1.0)
            t = x - distance
            return np.asarray(model(x, shared[:n], own, design), dtype=float) + float(shared[n]) * (np.asarray(model(t, shared[:n], own, design), dtype=float) - under(t))
        return new
    return Mechanism(f"a second component {distance:g} away", "second component", f"second component at {distance:g}", 0.0, wrap, lower=lower, upper=upper, ground=ground, theory=theory,
                     formula=r"$f(x) + a\,[\,f(x - D) - \ell\,],\ \ D = %g$" % distance)


def tail(stretch: float = 2.0, side: str = "upper", ground: str = "", theory: tuple | None = None, lower: float = 0.0, upper: float = np.inf) -> Mechanism:
    """A copy of the curve `stretch` times wider about its centre, on the `side` of larger abscissa, "upper", or of
    smaller, "lower", alone, of relative area a: a tail of one side, as the slow release of a column or the decay of a
    detector leaves. The copy is taken without the baseline that it carries: that line, read on the copy far from the
    curve on the side of the tail, ten spans of the record away, is taken out of it, so that the tail adds the shape of
    the curve and no step of its baseline, which no own parameter could take. Its first form carried that step, and a
    tail of a tenth of the curve was fitted at half its size. The tail rises at the centre over one sample of the record,
    and not at once: an edge that jumps as the centre crosses a sample stops the search of the place, which ended after
    three steps with the size right. Its area is bounded at zero unless `lower` says otherwise: a tail adds to the curve,
    and a tail of minus 1.12 times the area of the curve was kept in a cold test, standing in for a slope (round 2)."""
    if side not in ("upper", "lower"):
        raise ValueError(f"a tail is on the upper or the lower side, not {side!r}")
    sign = 1.0 if side == "upper" else -1.0

    def wrap(model, problem, n):
        if problem.centre is None:
            raise ValueError("a tail is placed about the centre of a curve: the problem gives none")

        def new(x, shared, own, design):
            x = np.asarray(x, dtype=float)
            c = float(problem.centre(shared[:n], own, design))
            f = np.asarray(model(x, shared[:n], own, design), dtype=float)
            wide = np.asarray(model(c + (x - c) / stretch, shared[:n], own, design), dtype=float) / stretch
            far = c + sign * 10.0 * float(np.ptp(x)) * np.array([1.0, 2.0])          # where the copy holds its baseline alone
            line = np.asarray(model(c + (far - c) / stretch, shared[:n], own, design), dtype=float) / stretch
            base = line[0] + (line[1] - line[0]) * (x - far[0]) / (far[1] - far[0])
            step = float(np.median(np.abs(np.diff(x)))) if x.size > 1 else 1.0
            rise = 0.5 * (1.0 + np.tanh(sign * (x - c) / (step or 1.0)))          # the side of the tail, smooth over one sample
            return f + float(shared[n]) * (wide - base) * rise
        return new
    return Mechanism(f"a tail on the {side} side, {stretch:g} times wider", "tail", f"tail {side} {stretch:g}", 0.0, wrap, lower=lower, upper=upper, ground=ground, theory=theory,
                     formula=r"$f(x) + a\,[\,f(c + (x - c)/s)/s - \ell\,]\ \Theta(%s(x - c)),\ \ s = %g$" % ("" if side == "upper" else "-", stretch))


def level(power: int = 2, ground: str = "", theory: tuple | None = None) -> Mechanism:
    """A power of the abscissa under every curve, of one size a for all: a ((x - middle) / half)^power, with the middle
    and the half span of the abscissas of all the records together, so that one size is one slope or one curvature in
    the unit of the abscissa on every record. A curvature of the baseline, as a drift of the detector or a background
    that the straight baseline of each record does not hold. Its first form took the middle and the half span of each
    record, and of whatever abscissa it was called on: a size meant another slope on a record of another span, 6 bars
    apart in a cold test (round 3), and a term that read the model far from its curve read another slope there."""
    if int(power) != power or power < 1:
        raise ValueError(f"level takes a whole power of 1 or more, not {power!r}")
    def wrap(model, problem, n):
        lo = min(float(np.min(r.x)) for r in problem.records)
        hi = max(float(np.max(r.x)) for r in problem.records)
        mid, half = 0.5 * (lo + hi), 0.5 * (hi - lo) or 1.0

        def new(x, shared, own, design):
            x = np.asarray(x, dtype=float)
            return np.asarray(model(x, shared[:n], own, design), dtype=float) + float(shared[n]) * ((x - mid) / half) ** power
        new.middle, new.half = mid, half
        return new
    return Mechanism(f"a power {power} of the abscissa under every curve", "baseline", f"level power {power}", 0.0, wrap, ground=ground, theory=theory,
                     formula=r"$f(x) + a\,((x - m)/h)^{%d}$" % power)


def periodic(period: float, part: str = "cosine", ground: str = "", theory: tuple | None = None) -> Mechanism:
    """A cosine or a sine of the abscissa of one `period`, about the centre of each curve, of one amplitude a for all
    records, in the units of the values: a fringe of an etalon, a ripple of a supply. A periodic function of unknown
    phase is the two, the cosine and the sine, added one after the other, since each term adds one parameter."""
    if part not in ("cosine", "sine"):
        raise ValueError(f"a periodic term is a cosine or a sine, not {part!r}")
    wave = np.cos if part == "cosine" else np.sin

    def wrap(model, problem, n):
        if problem.centre is None:
            raise ValueError("a periodic term is placed about the centre of a curve: the problem gives none")

        def new(x, shared, own, design):
            x = np.asarray(x, dtype=float)
            c = float(problem.centre(shared[:n], own, design))
            return np.asarray(model(x, shared[:n], own, design), dtype=float) + float(shared[n]) * wave(2.0 * np.pi * (x - c) / period)
        return new
    return Mechanism(f"a {part} of period {period:g}", "periodic", f"{part} of period {period:g}", 0.0, wrap, ground=ground, theory=theory,
                     formula=r"$f(x) + a\,\%s(2\pi\,(x - c)/P),\ \ P = %g$" % ("cos" if part == "cosine" else "sin", period))

