"""A term whose number is read from the dictionary HELD of the model file at every call, so that an ansatz on that
number reaches the term as it reaches the model, and its uncertainty becomes a bar of theory of every other value.

    from twinfit import held, mechanisms
    term = held.reading(lambda: mechanisms.second(HELD["distance"]))

A term of the package takes its number when its file is read: `mechanisms.second(HELD["distance"])` keeps the distance
of that moment, so the loop, which moves a held value one bar either way to read what it moves, moved nothing, and the
bar of theory of every value read 0 (round 6 of the cold tests). `reading` builds the term again at each call of the
model, from the value that HELD holds then. Its name and its parameter are those of the term as first built."""
from __future__ import annotations

from typing import Callable

from .mechanisms import Mechanism


def reading(make: Callable[[], Mechanism]) -> Mechanism:
    """The term that `make` builds, built again at each call of the model, so that a number it reads from HELD is read
    at that call. `make` takes nothing and returns a Mechanism."""
    first = make()
    if not isinstance(first, Mechanism):
        raise ValueError(f"held.reading takes a function that returns a term of the package or a Mechanism of your own, and it returned a {type(first).__name__}")

    def wrap(model, problem, n):
        def new(x, shared, own, design):
            return make().wrap(model, problem, n)(x, shared, own, design)
        return new
    return Mechanism(first.name, first.family, first.parameter, first.start, wrap, first.lower, first.upper, first.ground, first.theory)
