"""One draw of the count of where the report points, `toy/pointing.py`, made again from its seed and set beside the row
that the count kept: the first reproduction of a reader of the public repository (his messages 441 and 442, the bus
test). The row is made by the count's own function `one`, against the twin of its case that runs/ keeps, so nothing of
the count is written again here, and the count's digest and its table stay as they are.

    python toy/again.py                list the cases whose kept rows and twin runs/ holds, with their draws
    python toy/again.py CASE DRAW      make draw DRAW of CASE again, a few minutes on one core, and print each field of
                                       the row beside the kept one. Exit 0 when the two are the same row

A case is named as the table names it, in quotes: python toy/again.py "2 hidden under 1" 7
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from toy import model, pointing, power  # noqa: E402


def held(root: Path = ROOT) -> list:
    """The cases whose kept rows and whole twin runs/ holds, as (name, the draws kept)."""
    const, config = model.load(Path(root) / "toy", root)
    return [(c[0], sorted(power.read(pointing.kept(root, c[0])))) for c in pointing.cases(config)
            if pointing.kept(root, c[0]).exists() and pointing.twin_path(root, c[0]).exists()]


def again(name: str, draw: int, root: Path = ROOT) -> tuple:
    """The row made again and the kept row of draw `draw` of case `name`, each as JSON reads it."""
    const, config = model.load(Path(root) / "toy", root)
    case = next((c for c in pointing.cases(config) if c[0] == name), None)
    if case is None:
        raise SystemExit(f"REFUSED: no case {name!r}: `python toy/again.py` lists them")
    rows = power.read(pointing.kept(root, name))
    if draw not in rows:
        raise SystemExit(f"REFUSED: case {name!r} keeps no draw {draw}: it keeps {sorted(rows)[:5]} and {len(rows) - 5} more")
    if not pointing.twin_path(root, name).exists():
        raise SystemExit(f"REFUSED: runs/ holds no twin of {name!r}, {pointing.twin_path(root, name).name}: made again it takes hours")
    null = pointing.null_of(case, const, config, root)
    made = pointing.one((case[0], case[1], case[2], case[3], case[4], draw, const, config, null))[2]
    kept = {k: v for k, v in rows[draw].items() if k != "draw"}
    return json.loads(json.dumps(made)), kept


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        for name, draws in held():
            print(f"{name!r}: draws {draws[0]} to {draws[-1]}, {len(draws)} kept")
        sys.exit(0)
    if len(args) != 2 or not args[1].isdigit():
        print(__doc__)
        sys.exit(1)
    began = time.time()
    made, kept = again(args[0], int(args[1]))
    for key in sorted(set(made) | set(kept)):
        print(f"{key}: {'the same' if made.get(key) == kept.get(key) else 'DIFFERS'}\n  made again: {json.dumps(made.get(key))}\n  kept:       {json.dumps(kept.get(key))}")
    same = made == kept
    print(f"draw {args[1]} of {args[0]!r} made again in {time.time() - began:.0f} s: " + ("the kept row, field for field" if same else "NOT the kept row"))
    sys.exit(0 if same else 1)
