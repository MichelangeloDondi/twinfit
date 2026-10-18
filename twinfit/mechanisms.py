"""Terms that models of recorded curves often lack, each as the working model with one shared parameter more.

A candidate of the dictionary (candidates.py) is a discrepancy of the mean at the fitted values: it is ranked with no fit,
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
    pedestal(stretch)        a copy of the curve `stretch` times wider under it, of relative area a, stretch a number
    bent()                   the abscissa is not linear: x + a (x - middle)^2 / half its span
    second(distance)         a second copy of the curve, `distance` away along the abscissa, in its unit, after the curve
                             where it is positive and before it where negative, of relative area a
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

Which term answers which kind of the report of where to search:

    a second component of the curve     second
    a wider pedestal under the curve    pedestal
    a tail on one side                  tail
    a polynomial of the level           level
    a periodic function                 periodic, twice, its cosine and its sine, since each term adds one parameter:
                                        a scan of the cosine alone can find no period where the sine would
    a setting along which it changes    along
    the wings of the curves             widened, a kernel that widens the curve, beside a pedestal
    the flanks of the curves            widened, a width of another shape beside the one of the model, as a Gaussian beside a
                                        Lorentzian, where the misfit stands from half a width to two
    an asymmetry that grows from it     displaced, a mixture of copies displaced to one side
    an abscissa that is not straight    bent

Glitches and a change of the noise are no term of the mean: they are read in the weights of the fit, the noise of each
sample that the description gives. Each term here is made of copies of the curve, wider, displaced or on one side: a shape
that the physics names otherwise, as the exponential tail of a detector or of a column, is a Mechanism of your own, and the
report then gives its place and its kind, and not its form. A Mechanism is Mechanism(name, kind, parameter, start, wrap,
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

    def extend(self, problem: Problem, result=None) -> Problem:
        """The working model of `problem` with this term: one shared parameter more, started at the fit if one is given."""
        if self.parameter in problem.shared:
            raise ValueError(f"the model already holds a parameter named {self.parameter}")
        n = len(problem.shared)
        lo, hi = problem.bounds_shared
        start = list(problem.start_shared if result is None else result.shared) + [self.start]
        own = problem._start_own if result is None else [np.array(row, dtype=float) for row in result.own]
        centre = None if problem.centre is None else (lambda shared, own_k, design, c=problem.centre: c(shared[:n], own_k, design))
        return Problem(problem.records, self.wrap(problem.model, problem, n), problem.shared + (self.parameter,), problem.own, start, own, centre,
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


def _equal_steps(problem: Problem) -> None:
    for rec in problem.records:
        d = np.diff(rec.x)
        if not np.allclose(d, d[0], rtol=1e-6, atol=0.0) or d[0] <= 0:
            flat = np.linspace(rec.x[0], rec.x[-1], rec.x.size)
            off = float(np.max(np.abs(rec.x - flat)) / abs(flat[1] - flat[0])) if rec.x.size > 1 and flat[1] != flat[0] else float("inf")
            raise ValueError(f"record {rec.name!r}: a kernel is applied on an abscissa of equal rising steps, and its points lie up to {off:.2g} of a step from equal steps. "
                             "If that is the rounding of the file that the abscissa was read from, build the records with `twinfit.even(x)`")


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
        _equal_steps(problem)
        tall = max(problem.records, key=lambda rec: float(np.max(rec.y) - np.median(rec.y)))      # the record whose curve stands highest gives the width of a curve
        up = tall.y - np.median(tall.y)
        width = max(float(np.sum(up > 0.5 * up.max()) * (tall.x[1] - tall.x[0])), float(tall.x[1] - tall.x[0]))
        fine = max(1, int(np.ceil(FINE * float(tall.x[1] - tall.x[0]) / width)))                   # the model is read on steps of at most a FINE-th of that width
        base = float(tall.x[1] - tall.x[0])                                                        # the step of the records: an uneven abscissa is read on a grid of this step

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
    return Mechanism(name, "kernel", parameter, start, wrap, lower=0.0, upper=WIDEST, ground=ground, theory=theory)


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
    return Mechanism(name, "displacement", parameter, start, wrap, lower=0.0, upper=WIDEST, ground=ground, theory=theory)


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
    return Mechanism(f"{parameter} changes along{square} {key}", "grows with a covariate", f"{parameter} along {key}" + ("^2" if power == 2 else ""), 0.0, wrap, ground=ground, theory=theory)


def pedestal(stretch: float = 8.0, ground: str = "", theory: tuple | None = None, lower: float = -np.inf, upper: float = np.inf) -> Mechanism:
    """A copy of the curve `stretch` times wider about the centre of the record, of area a times the area of the curve.
    `lower` and `upper` bound that area. Give `lower=0.0` when the ground of the term says that it adds to the curve:
    left free, a fit can return a negative pedestal, which stands for wings that the working model draws too high, and
    a term on theory would then be added with a sign that its own ground forbids."""
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
    return Mechanism(f"a pedestal {stretch:g} times wider than the curve", "pedestal", f"pedestal {stretch:g}", 0.0, wrap, lower=lower, upper=upper, ground=ground, theory=theory)


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
    return Mechanism("a quadratic map of the abscissa", "abscissa", "bend", 0.0, wrap, lower=-0.25, upper=0.25, ground=ground, theory=theory)


def second(distance: float, ground: str = "", theory: tuple | None = None, lower: float = -np.inf, upper: float = np.inf) -> Mechanism:
    """A second copy of the curve, `distance` away along the abscissa, in the units of the abscissa, of a times the size
    of the curve: an impurity that elutes later, a line of another species or isotope beside the first. The copy is the
    curve less the straight line that the model draws far from it, so that it carries no share of the offset or the slope
    of the record. A negative distance places it before."""
    def wrap(model, problem, n):
        def new(x, shared, own, design):
            x = np.asarray(x, dtype=float)
            c = float(problem.centre(shared[:n], own, design)) if problem.centre is not None else 0.5 * float(x[0] + x[-1])
            under = _line_under(model, shared[:n], own, design, c, float(np.ptp(x)) or 1.0)
            t = x - distance
            return np.asarray(model(x, shared[:n], own, design), dtype=float) + float(shared[n]) * (np.asarray(model(t, shared[:n], own, design), dtype=float) - under(t))
        return new
    return Mechanism(f"a second component {distance:g} away", "second component", f"second component at {distance:g}", 0.0, wrap, lower=lower, upper=upper, ground=ground, theory=theory)


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
    return Mechanism(f"a tail on the {side} side, {stretch:g} times wider", "tail", f"tail {side} {stretch:g}", 0.0, wrap, lower=lower, upper=upper, ground=ground, theory=theory)


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
    return Mechanism(f"a power {power} of the abscissa under every curve", "baseline", f"level power {power}", 0.0, wrap, ground=ground, theory=theory)


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
    return Mechanism(f"a {part} of period {period:g}", "periodic", f"{part} of period {period:g}", 0.0, wrap, ground=ground, theory=theory)

