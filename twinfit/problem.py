"""The description of a problem: recorded curves across a design, and a working model with shared and own parameters.

A record is one curve: its abscissa, its values, the standard deviation of each value, and the covariates of the design
under which it was recorded (a block, a level of a factor, an order of acquisition). The working model gives the mean of a
record from two sets of parameters: the shared ones, which are the same for every record and are what the analysis is
after, and the own ones, which belong to one record (an area, a location, an offset) and are nuisances.

    def model(x, shared, own, design):          # arrays in, the mean of the record out
        ...
    problem = Problem(records, model, shared=("k", "g", "s"), own=("area", "centre", "offset", "slope"),
                      start_shared=[...], start_own=lambda record: [...], centre=lambda shared, own, design: own[1])

`centre` is optional: where it is given, slices of the abscissa are taken in pairs about it and a curve can be stretched
about it. Bounds are optional too: `bounds_shared=(lower, upper)` and `bounds_own=(lower, upper)`, or a function of a
record that returns such a pair, keep a search inside the values that have a meaning, for instance an area that is not
negative and a centre inside the abscissa. A record that holds almost no signal needs them: without them its own
parameters wander and the fit of every record fails with it. Nothing here fits or draws anything.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np


def even(x, within: float = 0.05) -> np.ndarray:
    """The same abscissa in equal steps, from its first value to its last. An abscissa that was written to a file is
    rounded, and its steps then differ by the last digit that was kept: the terms that apply a kernel refuse it, since
    they count in steps. Where every point lies within `within` of a step from its place in equal steps the difference
    is that rounding, and the equal steps are returned. A larger difference is refused: such an abscissa is not one of
    equal steps, and the records must be put on one by the researcher, who knows how they were taken."""
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or x.size < 2:
        raise ValueError("an abscissa of at least two values")
    out = np.linspace(x[0], x[-1], x.size)
    step = (x[-1] - x[0]) / (x.size - 1)
    off = float(np.max(np.abs(x - out)) / step) if step > 0 else np.inf
    if not off <= within:
        raise ValueError(f"REFUSED: the abscissa lies up to {off:.2g} of a step from equal steps, more than the {within:g} that the rounding of a file explains")
    return out


@dataclass
class Record:
    """One recorded curve. `sd` is the standard deviation of each value, a number or an array."""
    x: np.ndarray
    y: np.ndarray
    sd: np.ndarray
    design: dict = field(default_factory=dict)
    name: str = ""

    def __post_init__(self):
        self.x = np.asarray(self.x, dtype=float)
        self.y = np.asarray(self.y, dtype=float)
        self.sd = np.broadcast_to(np.asarray(self.sd, dtype=float), self.x.shape).copy()
        if self.x.shape != self.y.shape or self.x.ndim != 1:
            raise ValueError(f"record {self.name!r}: the abscissa and the values must be two arrays of one length")
        if not np.all(np.isfinite(self.x)) or not np.all(np.isfinite(self.y)):
            raise ValueError(f"record {self.name!r}: a value or an abscissa is not a number")
        if not np.all(self.sd > 0):
            raise ValueError(f"record {self.name!r}: every standard deviation must be positive")
        if self.x.size > 1 and not np.all(np.diff(self.x) > 0):
            raise ValueError(f"record {self.name!r}: the abscissa must rise strictly from its first value to its last. Reverse or sort the record before it is read: "
                             "the widths, the periods and the bounds of a centre are counted from the ends and the steps of the abscissa")


class Problem:
    """Records, a working model, and the names and starts of its parameters."""

    def __init__(self, records, model: Callable, shared, own, start_shared, start_own, centre: Callable | None = None, name: str = "the working model",
                 bounds_shared=None, bounds_own=None):
        self.records = list(records)
        if not self.records:
            raise ValueError("a problem needs at least one record")
        self.model = model
        self.shared, self.own = tuple(shared), tuple(own)
        self.start_shared = np.asarray(start_shared, dtype=float)
        if self.start_shared.shape != (len(self.shared),):
            raise ValueError("one start for each shared parameter")
        self._start_own = start_own
        self.centre = centre
        self.name = name
        self._bounds_shared, self._bounds_own = bounds_shared, bounds_own
        self.bounded = bounds_shared is not None or bounds_own is not None
        self.bounds_shared = self._box(bounds_shared, len(self.shared), "the shared parameters")
        if bounds_own is not None and not callable(bounds_own):
            self._box(bounds_own, len(self.own), "the own parameters")      # read once here, so that a pair that cannot be one is refused at the start

    @staticmethod
    def _box(bounds, n: int, what: str) -> tuple:
        if bounds is None:
            return np.full(n, -np.inf), np.full(n, np.inf)
        lo, hi = (np.asarray(b, dtype=float) for b in bounds)
        if lo.shape != (n,) or hi.shape != (n,) or not np.all(lo < hi):
            raise ValueError(f"the bounds of {what}: one lower and one upper value for each, the lower under the upper")
        return lo, hi

    def own_bounds(self, k: int) -> tuple:
        """(lower, upper) for the own parameters of record k. Without bounds: no limit either side."""
        b = self._bounds_own(self.records[k]) if callable(self._bounds_own) else self._bounds_own
        return self._box(b, len(self.own), f"the own parameters of record {k}")

    def own_start(self, k: int) -> np.ndarray:
        start = self._start_own(self.records[k]) if callable(self._start_own) else self._start_own[k]
        start = np.asarray(start, dtype=float)
        if start.shape != (len(self.own),):
            raise ValueError(f"record {k}: one start for each own parameter")
        return start

    def predict(self, k: int, shared, own, x=None) -> np.ndarray:
        rec = self.records[k]
        return np.asarray(self.model(rec.x if x is None else x, np.asarray(shared, dtype=float), np.asarray(own, dtype=float), rec.design), dtype=float)

    def whitened(self, k: int, shared, own) -> np.ndarray:
        rec = self.records[k]
        return (rec.y - self.predict(k, shared, own)) / rec.sd

    def centre_of(self, k: int, shared, own):
        return None if self.centre is None else float(self.centre(np.asarray(shared, dtype=float), np.asarray(own, dtype=float), self.records[k].design))

    def blocks(self, key: str) -> dict:
        """The records grouped by one covariate of the design: value -> the places of its records."""
        out = {}
        for k, rec in enumerate(self.records):
            out.setdefault(rec.design.get(key), []).append(k)
        return out

    def with_model(self, model: Callable, shared, own, start_shared, start_own, centre: Callable | None = None, name: str = "the next model",
                   bounds_shared=None, bounds_own=None) -> "Problem":
        """The same records under another working model: the next rung of a ladder. Its bounds are its own."""
        return Problem(self.records, model, shared, own, start_shared, start_own, self.centre if centre is None else centre, name, bounds_shared, bounds_own)

    def with_records(self, records) -> "Problem":
        """The same working model on other records: synthetic ones, or a subset."""
        return Problem(records, self.model, self.shared, self.own, self.start_shared, self._start_own, self.centre, self.name, self._bounds_shared, self._bounds_own)

    @property
    def samples(self) -> int:
        return int(sum(rec.x.size for rec in self.records))
