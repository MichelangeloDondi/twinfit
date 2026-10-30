"""A stand-in for the queue of the author's machine, which is not public. A count writes into its table, beside its
shares, what the queue spent on its waves, through `spent`. Here no queue runs: `toy/count.py`, or a slice made by
hand, makes the rows, so `spent` says that no wave of a queue was recorded, and a table written here says so. Nothing
else of the author's queue is needed to make, count or check a count."""
from __future__ import annotations

from pathlib import Path


def spent(root: Path, program: str, made_from: str) -> dict:
    """What a queue spent on the units of `program` made from `made_from`: none here."""
    return {"recorded": False, "waves": [], "ended": 0, "failed": 0, "overrun": 0, "untimed": 0, "seconds": {},
            "why": "no queue in this repository: the rows were made by toy/count.py or by hand"}
