"""Another field: where the report of the package points on the projects of its gallery, row M443 of the plan, the red team
of his message 331, moderate 2 ("The package is model agnostic!", his messages 141 and 331).

The gallery of the package, twinfit/examples, holds eleven projects of other fields, each with a truth that its folder
does not hold. Here each is drawn again on `DRAWS` draws of its noise, its first model is run through the loop once, and
the first answers of its report are read against `KEY`, written from the truth of each project and its physics before
any count: where its misfit shows, of which kind, along which setting, and which shape of the noise. A question that
the truth does not settle is not asked of that project.

    python3 tools/autoqueue.py add toy/gallery_count.py   the count, through the queue (his message 147): each draw of
                                                         each project a unit, kept in runs/ as it ends
    python toy/gallery_count.py table                     writes toy/gallery.json, which `tools/met.py M443` reads
    python toy/gallery_count.py --check                   ends without error when toy/gallery.json holds the shares of
                                                         its rows, made from the files as they are
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from twinfit import loop  # noqa: E402
from twinfit.examples import gallery  # noqa: E402

from toy import power  # noqa: E402

DRAWS = 20          # a share near one half known to 0.11 a project, and to 0.04 pooled
SETS = 30          # the sets of the twin of the first model, those of the loop and of the command line
KEY = {          # written again 2026-10-23 for the blind report, from the truth of each project and its physics, before any count of it (the red team of his message 343, S1): the first answer that each question should have, or the answers that the physics admits where a width moves the core and the flanks together. A question absent is not asked of that project
    "pedestal": {"where": "the wings"},
    "satellite": {"where": "the asymmetry"},
    "saturation": {"changes along": ("the core along power", "the flanks along power")},
    "filter": {"what": "a noise correlated from sample to sample"},
    "fano": {"where": "the asymmetry"},
    "counts": {"what": "a noise that grows with the signal"},
    "clock": {"changes along": ("the asymmetry along power",)},
    "gel": {"where": "the asymmetry"},
    "transit": {"changes along": ("the core along wavelength", "the flanks along wavelength")},
    "fringes": {"shape": "a periodic structure along the abscissa"},
    "returns": {"noise": "a noise with heavy tails"},
}


def written(name: str, folder: Path, draw: int) -> None:
    """The project `name` written into `folder` with the noise of draw `draw`: the seed of the project with the draw
    beside it, so that draw 0 is the project as the gallery draws it."""
    source = (gallery.PROJECTS / f"{name}.py").read_text()
    if draw:
        source, n = re.subn(r"np\.random\.default_rng\((\d+)\)", lambda m: f"np.random.default_rng([{m.group(1)}, {int(draw)}])", source)
        if n != 1:
            raise ValueError(f"the project {name} holds {n} seeds, where one is drawn again")
    folder.mkdir(parents=True)
    space = dict(gallery.drawing(1.0, None), home=folder, __name__=f"twinfit.examples.projects.{name}")
    exec(compile(source, str(gallery.PROJECTS / f"{name}.py"), "exec"), space)


def one(name: str, draw: int, root: Path) -> dict:
    """One row: the first answers of the report of the first model on one draw of a project, against its key."""
    home = Path(root) / "runs" / "gallery_count" / f"{name}_{draw}"
    if home.exists():
        shutil.rmtree(home)
    try:
        written(name, home, draw)
        loop.run(home, n=0, sets=SETS, report=True, workers=1)
        table = json.loads((home / "loop" / "model_0.json").read_text())
        first = table["where"]["first_of_each_question"]
        beyond = [[a["question"], a["answer"]] for a in table["where"]["answers"] if a["beyond_the_twin"]]
        return {"first": {q: first.get(q) for q in ("where", "shape", "what", "changes along", "noise")}, "beyond": beyond,
                "right": {q: bool(beyond) and (first.get(q) in want if isinstance(want, tuple) else first.get(q) == want) for q, want in KEY[name].items()}}
    except (Exception, SystemExit) as why:
        return {"failed": f"{type(why).__name__}: {why}"[:300]}
    finally:
        if home.exists():
            shutil.rmtree(home)          # the folder of a draw is no result: its row is


MAKES = ("DRAWS", "SETS", "KEY", "written", "one")


def code_of(package: Path) -> str:
    """The digest of what the package runs, its code and its data, by name and bytes: its pages run nothing of the gallery,
    and the digest of the whole package made a page written for a reader, as 0.1.0.post1's, a count to make again (his
    messages 441 and 442, the bus test)."""
    h = hashlib.sha256()
    for p in sorted(Path(package).rglob("*")):
        name = p.relative_to(package).as_posix()
        if p.is_file() and "__pycache__" not in p.parts and p.suffix in (".py", ".json", ".toml", ".txt", ".csv") and name != "version.py":
            h.update(name.encode() + b"\0" + p.read_bytes() + b"\0")
    return h.hexdigest()


def digest(root: Path = ROOT) -> str:
    """What a row depends on: the parts of this program that make it, and the package with the projects of its gallery."""
    root = Path(root)
    h = hashlib.sha256(power.sources((root / "toy" / "gallery_count.py").read_text(), MAKES))
    h.update(code_of(root / "twinfit").encode())          # the whole package, its projects among it: the gallery imports the command line, whose imports are made while it runs, so the files that a row imports cannot be read in their text
    return h.hexdigest()


def kept(root: Path, name: str) -> Path:
    return Path(root) / "runs" / f"gallery_count_{name}_{digest(root)[:12]}.jsonl"


def units(root: Path = ROOT) -> list:
    """What the queue asks: the draws of each project not yet kept, a kind a project, since their costs differ."""
    return [{"condition": name, "units": [k for k in range(1, DRAWS + 1) if k not in power.read(kept(root, name))], "after": [], "kind": name} for name in KEY]


def slice_of(condition: str, first: int, last: int, limit: float = float("inf"), root: Path = ROOT) -> dict:
    """The draws of a project from `first` to `last` not yet kept, each kept as it ends, until `limit` seconds."""
    import resource
    root = Path(root)
    began, ids, seconds = time.time(), [], []
    path = kept(root, condition)
    path.parent.mkdir(exist_ok=True)
    for k in range(first, last + 1):
        if k in power.read(path):
            continue
        t = time.time()
        row = one(condition, k, root)
        with path.open("a") as fh:
            fh.write(json.dumps(dict(row, draw=k)) + "\n")
        ids.append(k)
        seconds.append(time.time() - t)
        if time.time() - began > limit:
            break
    return {"done": ids, "seconds": seconds, "rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


def shares(rows_of: dict) -> dict:
    """For each question, and pooled over every question asked of every project, the share of draws on which the report
    puts first the answer of the key, with its standard error. A failed draw is right on nothing and is named."""
    def share(hits):
        n = len(hits)
        p = sum(hits) / n if n else math.nan
        return {"share": p, "standard_error": math.sqrt(p * (1.0 - p) / n) if n else math.nan, "of": n}
    by = {}
    for name, rows in rows_of.items():
        for r in rows:
            for q in KEY[name]:
                by.setdefault(q, []).append(bool(r.get("right", {}).get(q)))
    return {"by_question": {q: share(v) for q, v in by.items()}, "pooled": share([h for v in by.values() for h in v]),
            "failed": sum(1 for rows in rows_of.values() for r in rows if "failed" in r)}


def table(root: Path = ROOT) -> dict:
    """toy/gallery.json from the rows kept under the files as they are, with the version of the package under them.
    Refused while a draw is not kept."""
    from twinfit.version import __version__
    root = Path(root)
    got = {name: list(power.read(kept(root, name)).values()) for name in KEY}
    lacking = {n: DRAWS - len(rows) for n, rows in got.items() if len(rows) < DRAWS}
    if lacking:
        raise SystemExit("REFUSED: draws are not kept under the files as they are: " + ", ".join(f"{n} {k}" for n, k in lacking.items()) + ". Add toy/gallery_count.py to the queue again")
    s = shares(got)
    out = {"digest": digest(root), "package": __version__, "key": KEY, "projects": {n: {"right": shares({n: rows})["pooled"], "rows": rows} for n, rows in got.items()},
           "summary": {"pooled": s["pooled"]}, "by_question": s["by_question"], "failed": s["failed"]}
    (root / "toy" / "gallery.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def check(root: Path = ROOT) -> list:
    path = Path(root) / "toy" / "gallery.json"
    if not path.exists():
        return ["toy/gallery.json is not written: `python3 tools/autoqueue.py add toy/gallery_count.py`, then `python toy/gallery_count.py table`"]
    held = json.loads(path.read_text())
    if held.get("digest") != digest(root):
        return ["toy/gallery.json was made from other files than those on disk now: count it again"]
    return [] if held["summary"]["pooled"] == shares({n: p["rows"] for n, p in held["projects"].items()})["pooled"] else ["toy/gallery.json does not hold the shares of its own rows"]


def progress(root: Path = ROOT) -> str:
    """The preliminary result after a wave: for each project with rows, on how many draws the report puts first the
    answers of its key."""
    said = []
    for name in KEY:
        rows = [r for r in power.read(kept(root, name)).values() if "failed" not in r]
        if rows:
            said.append(f"{name}: all of its key first on {sum(all(r['right'].values()) for r in rows)} of {len(rows)}")
    return "; ".join(said) or "no row yet"


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["units"]:
        print(json.dumps(units()))
    elif args[:1] == ["slice"] and len(args) >= 4:
        if not os.environ.get("TOY_QUEUE"):
            sys.exit("REFUSED: a slice is a wave of the queue: `python3 tools/autoqueue.py add toy/gallery_count.py` (his message 147)")
        print(json.dumps(slice_of(args[1], int(args[2]), int(args[3]), float(args[5]) if args[4:5] == ["--for"] else float("inf"))))
    elif args == ["progress"]:
        print(progress())
    elif args == ["digest"]:
        print(digest())
    elif args == ["table"]:
        print(json.dumps(table()["summary"], indent=1))
    elif args == ["--check"]:
        wrong = check()
        print("\n".join(wrong) or "another field: made from the files as they are")
        sys.exit(1 if wrong else 0)
    else:
        print(__doc__)
