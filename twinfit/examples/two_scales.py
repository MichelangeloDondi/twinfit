"""A worked example on synthetic records: curves with two scales across a design, and a term that the working model lacks.

The curves are the convolution of a Cauchy and a normal density, with an area, a location, an offset and a slope of
their own. The Cauchy scale grows with a covariate of the design: scale = g + k z. The normal scale s is the same for
every curve. The three shared parameters are k, g and s, and the two scales are weakly told apart: this is the design of
the first toy of the repository that this package comes from.

    records(hidden=..., noise=...)      synthetic records, with or without a term that the working model does not hold
    problem(records)                    the working model on them
The hidden terms that can be planted are those of HIDDEN. Each is small beside the curve.
"""
from __future__ import annotations

import numpy as np
from scipy.special import voigt_profile

from ..problem import Problem, Record

TRUTH = {"k": 0.02, "g": 0.6, "s": 1.1}
LEVELS = (10.0, 30.0, 60.0, 100.0)             # the covariate of the design, one block a level
REPEATS = 3                                    # curves a block
HALF, STEP = 30.0, 0.25                        # the abscissa runs from -HALF to HALF
AREA, NOISE = 1.0, 0.004                       # the noise of a sample, about 2 per cent of the height of a curve at its centre
HIDDEN = ("none", "second component", "pedestal", "periodic", "scale grows with the square of the covariate", "one block wider")
# Three further terms, which the generic dictionary of candidates.py does not hold in form or holds off its grid. The five
# above are of its own families, and the displaced copy of "second component" sits on a point of its grid of distances.
OUTSIDE = ("second component off the grid", "step of the level", "narrow dip")
# What `size` means for each hidden term, in its own units (review 35, moderate 3): it is not one size.
SIZES = {"second component": "the area of the displaced copy, as a share of the area of the curve",
         "pedestal": "0.47 times the size is the height of the broad copy at its centre, as a share of the height at its centre of the curve at z = 30: 0.47 per cent for a size of 0.01",
         "periodic": "the amplitude of the cosine, as a share of the height at its centre of the curve at z = 30. The curves of the four blocks are 1.2, 1.0, 0.78 and 0.59 times that height",
         "scale grows with the square of the covariate": "the relative growth of k at z = 100 is 10 times the size: 10 per cent for a size of 0.01",
         "one block wider": "the curves of the last block are wider by 10 times the size: 10 per cent for a size of 0.01",
         "second component off the grid": "the area of the displaced copy, as a share of the area of the curve",
         "step of the level": "the height of the step, as a share of the height of the curve at z = 30",
         "narrow dip": "the depth of the dip, as a share of the height of the curve at z = 30"}


def curve(x, k, g, s, z):
    return voigt_profile(x, s, g + k * z)


def model(x, shared, own, design):
    area, centre, offset, slope = own
    return area * curve(x - centre, shared[0], shared[1], shared[2], design["z"]) + offset + slope * x


def hidden_term(name: str, x, z, block, size: float):
    """What is added to the mean of a record for a hidden term of a given size, as a share of the height of the curve at z = 30 at its centre."""
    peak = AREA * curve(0.0, TRUTH["k"], TRUTH["g"], TRUTH["s"], 30.0)
    if name == "none":
        return np.zeros_like(x)
    if name == "second component":                       # a smaller copy of the curve, displaced
        return size * AREA * curve(x - 6.0, TRUTH["k"], TRUTH["g"], TRUTH["s"], z)
    if name == "pedestal":                               # a broad copy under the curve
        return size * peak * 10.0 * voigt_profile(x, 8.0, 0.5)
    if name == "periodic":
        return size * peak * np.cos(2.0 * np.pi * x / 7.0 + 0.4)
    if name == "scale grows with the square of the covariate":
        wider = AREA * curve(x, TRUTH["k"] * (1.0 + size * 10.0 * z / 100.0), TRUTH["g"], TRUTH["s"], z)
        return wider - AREA * curve(x, TRUTH["k"], TRUTH["g"], TRUTH["s"], z)
    if name == "one block wider":                        # the curves of the last block are wider, for no cause in the model
        if z != LEVELS[-1]:
            return np.zeros_like(x)
        f = 1.0 + size * 10.0
        return AREA * (curve(x / f, TRUTH["k"], TRUTH["g"], TRUTH["s"], z) / f - curve(x, TRUTH["k"], TRUTH["g"], TRUTH["s"], z))
    if name == "second component off the grid":         # displaced by 5.0, between two points of the dictionary's grid of distances
        return size * AREA * curve(x - 5.0, TRUTH["k"], TRUTH["g"], TRUTH["s"], z)
    if name == "step of the level":                      # the level under the curve jumps at the centre: no family of the dictionary
        return size * peak * 0.5 * np.sign(x)
    if name == "narrow dip":                             # a dip narrower than the curve at its centre: the dictionary holds wider copies only
        return -size * peak * np.exp(-0.5 * (x / 0.5) ** 2)
    raise KeyError(name)


def records(hidden: str = "none", size: float = 0.01, noise: float = NOISE, seed: int = 0) -> list:
    rng = np.random.default_rng([20261009, seed])
    x = np.arange(-HALF, HALF + STEP / 2.0, STEP)
    out = []
    for b, z in enumerate(LEVELS):
        for rep in range(REPEATS):
            centre = 0.3 * rng.standard_normal()
            mean = AREA * curve(x - centre, TRUTH["k"], TRUTH["g"], TRUTH["s"], z) + 0.01 + 0.0002 * x + hidden_term(hidden, x - centre, z, b, size)
            y = mean + noise * rng.standard_normal(x.size) if noise > 0 else mean
            out.append(Record(x, y, max(noise, NOISE), {"z": z, "block": f"z = {z:g}", "repeat": rep}, name=f"z {z:g}, curve {rep}"))
    return out


def problem(recs: list) -> Problem:
    def start_own(rec):
        base = np.median(rec.y[:20])
        area = np.sum(rec.y - base) * STEP
        centre = rec.x[np.argmax(rec.y)]
        return [area, centre, base, 0.0]
    return Problem(recs, model, shared=("k", "g", "s"), own=("area", "centre", "offset", "slope"), start_shared=[0.03, 0.4, 0.9],
                   start_own=start_own, centre=lambda shared, own, design: own[1], name="two scales, one growing with the covariate")
