"""The chance that the sign names the kernel of a missing width, row M348 and paper A, `PAPERS.md`: "the chance that the
sign names the right kernel against the size of the width and the noise, with the scores of the tie in a file".

After a fit that fits the line's own width, a small missing Gaussian width and a small missing Lorentzian width leave
residuals along nearly one direction with opposite signs, so their scores tie, and the package refuses a width beside the
line whose size the residuals ask negative (`propose._one`, his message 269). Here each kernel is planted in turn on the
eight traces of the fifth detail of the toy, at each size of SIZES and at each level of NOISES of the toy's law, and the
traces are fitted by the model that lacks it and fits the width of the other kind: the laser noise missing, a Gaussian,
with the permeated gas free, and the permeated gas missing, a Lorentzian, with the laser noise free. The six widths beside
the line of the package's dictionary are read on each draw, with their status, their size and their score: on each rung
of widths the sign names a kernel where one of the two is ranked with a positive size and the other is refused for its
sign, and the name is right where it is the kernel planted.

    python3 tools/autoqueue.py add toy/sign.py    the count through the queue: a draw of a condition is a unit, kept in
                                                  runs/sign_units as it ends
    python toy/sign.py table                       writes the chances and the scores into toy/matrix.json under `sign`
    python toy/sign.py --check                     ends without error when that table was made from the files as they are
"""
from __future__ import annotations

import copy
import dataclasses
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from toy import fits, model, traces  # noqa: E402

DETAIL = 5
KINDS = {"Gaussian": ("laser_noise", "gas_width"), "Lorentzian": ("gas_width", "laser_noise")}          # the kernel planted: the parameter of the toy that carries it, and the width of the other kind that the fit frees
SIZES = (0.05, 0.1, 0.2, 0.4)          # the planted parameter, in MHz: the toy's laser noise is 0.369 and its permeated gas 0.037
NOISES = (0.1, 0.3, 1.0, 3.0)          # the noise, over the toy's law
DRAWS = 400          # the third stage, his message 502 on decision 27. 200 at the second stage of row M529, his message 468
RUNGS = (0.03, 0.1, 0.3)          # the widths beside the line that the toy's bench reads, in widths of the curve: the package holds no list of terms (the owner's message 340), and this study of named widths keeps its own


def widths(lin, rungs=RUNGS) -> list:
    """A Gaussian and a Lorentzian width beside the line on a ladder of sizes, in widths of the curve, each made from a term of
    `twinfit.mechanisms` and of a size that cannot be negative: the ladder that this study reads."""
    from twinfit.candidates import from_mechanism
    from twinfit.mechanisms import gaussian_width, lorentzian_width
    unit = float(np.median([lin.width_of(k) for k in range(len(lin.problem.records))]))
    return [dataclasses.replace(from_mechanism(lin, width(name=f"a {kind} width of {m:g} widths beside the line"), m * unit), sign=1)
            for m in rungs for kind, width in (("Gaussian", gaussian_width), ("Lorentzian", lorentzian_width))]


def conditions() -> list:
    """Every condition, as (its name, the kernel, the size, the noise), the noise first so that a slice holds one level."""
    return [(f"a {kind} width of {size:g} MHz, noise {noise:g}", kind, size, noise) for noise in NOISES for kind in KINDS for size in SIZES]


def named(read: dict, rung: float) -> str:
    """The kernel that the sign names on one rung: the kind whose width is ranked with a positive size while the other is
    refused for its sign, or `none`."""
    got = {kind: read.get(f"a {kind} width of {rung:g} widths beside the line", {}) for kind in KINDS}
    ranked = [k for k, r in got.items() if r.get("status") == "ranked" and (r.get("amplitude") or 0.0) > 0.0]
    refused = [k for k, r in got.items() if r.get("status") == "of the sign that it cannot have"]
    return ranked[0] if len(ranked) == 1 and len(refused) == 1 and ranked[0] != refused[0] else "none"


def unit(condition: str, draw: int, root: Path = ROOT) -> dict:
    """One draw of one condition: the traces with the kernel planted, fitted by the model that lacks it with the noise's
    weights at its level, and each width beside the line read by the package, its status, size and score."""
    import twinfit
    from twinfit import noise as N, propose
    index = [c[0] for c in conditions()].index(condition)
    _, kind, size, noise = conditions()[index]
    planted, free = KINDS[kind]
    const, config = model.load(Path(root) / "toy", root)
    made, fitted = copy.deepcopy(config), copy.deepcopy(config)
    made["parameter"][planted]["true"] = size
    fitted["parameter"][planted]["true"] = 0.0          # the model lacks the kernel: held at nothing
    rng = np.random.default_rng([2026, 10, 21, index, draw])
    rows = []
    for P, T in traces.conditions(const):
        t = traces.trace(DETAIL, P, T, const, made, 0)
        rows.append({**t, "y": t["signal"] + noise * t["sd"] * rng.standard_normal(t["nu"].size)})
    law = fits.variance_of(const)
    try:
        problem, result, _ = N.settle(fits.problem_of(DETAIL, rows, const, fitted, [free]), lambda mean, rec, own: noise * noise * law(mean, rec, own))
        lin = twinfit.Linear(problem, result)
        read = {}
        for c in widths(lin):
            if " width of " in c.name and c.name.endswith("beside the line"):
                e, unsigned = propose._one(lin, c), propose._one(lin, dataclasses.replace(c, sign=0))          # the score of a width refused for its sign is read with its sign lifted: the two scores of the tie
                read[c.name] = {"status": e["status"], "amplitude": float(unsigned["amplitude"][0]) if unsigned.get("amplitude") else None, "score": unsigned.get("score")}
    except (ValueError, np.linalg.LinAlgError, RuntimeError) as e:          # a fit that does not end is said and counted apart
        return {"failed": f"{type(e).__name__}: {e}"[:300]}
    return {"widths": read, "named": {f"{rung:g}": named(read, rung) for rung in RUNGS}, "free_value": float(result.shared[0]), "chi_square_a_degree": float(result.chi2) / max(int(result.dof), 1)}


def unit_path(root: Path, condition: str, draw: int) -> Path:
    index = [c[0] for c in conditions()].index(condition)
    return Path(root) / "runs" / "sign_units" / f"{index:02d}_{draw:03d}.json"


def digest(root: Path = ROOT) -> str:
    """What a unit depends on: this program, the model, the traces, the fits, the two files of values, and the files of
    the package that it imports."""
    from toy import draws
    h = hashlib.sha256()
    names = ["toy/sign.py", "toy/model.py", "toy/fits.py", "toy/traces.py", "toy/const.toml", "toy/config.toml"]
    draws.hashed(h, Path(root), names + draws.package_of(Path(root), ("toy/sign.py",)))
    return h.hexdigest()


def done(root: Path, made_from: str, condition: str) -> dict:
    """The draws of a condition kept under the digest of the files as they are, by their number."""
    out = {}
    for draw in range(1, DRAWS + 1):
        path = unit_path(root, condition, draw)
        if path.exists():
            kept = json.loads(path.read_text())
            if kept.get("digest") == made_from:
                out[draw] = kept["result"]
    return out


def units(root: Path = ROOT) -> list:
    """What the queue asks first: for each condition the draws not yet kept, its kind the level of the noise, so that the
    pre-waves time the levels apart."""
    made_from = digest(root)
    out = []
    for name, kind, size, noise in conditions():
        left = [d for d in range(1, DRAWS + 1) if d not in done(root, made_from, name)]
        if left:
            out.append({"condition": name, "units": left, "after": [], "kind": f"noise {noise:g}"})
    return out


def slice_of(condition: str, first: int, last: int, limit: float = float("inf"), root: Path = ROOT) -> dict:
    """The draws from `first` to `last` of a condition that are not yet kept, each kept as it ends, until `limit` seconds."""
    import resource
    made_from, began, ids, seconds = digest(root), time.time(), [], []
    kept = done(root, made_from, condition)
    for draw in range(first, last + 1):
        if draw in kept:
            continue
        t = time.time()
        result = unit(condition, draw, root)
        path = unit_path(root, condition, draw)
        path.parent.mkdir(parents=True, exist_ok=True)
        aside = path.with_suffix(".part")
        aside.write_text(json.dumps({"digest": made_from, "result": result}) + "\n")
        aside.replace(path)
        ids.append(draw)
        seconds.append(time.time() - t)
        if time.time() - began > limit:
            break
    return {"done": ids, "seconds": seconds, "rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


def summary(results: list, kind: str) -> dict:
    """The draws of one condition read together: on each rung the share that the sign names rightly, wrongly and not at
    all, and the median score of each kind's width, the scores of the tie."""
    ended = [r for r in results if "failed" not in r]
    out = {"draws": len(results), "failed": len(results) - len(ended), "rungs": {}}
    other = next(k for k in KINDS if k != kind)
    for rung in RUNGS:
        names = [r["named"][f"{rung:g}"] for r in ended]
        scores = {k: [r["widths"].get(f"a {k} width of {rung:g} widths beside the line", {}).get("score") for r in ended] for k in KINDS}
        out["rungs"][f"{rung:g}"] = {"right": names.count(kind) / max(len(ended), 1), "wrong": names.count(other) / max(len(ended), 1), "none": names.count("none") / max(len(ended), 1),
                                     "median_score": {k: float(np.median([s for s in v if s is not None])) if any(s is not None for s in v) else None for k, v in scores.items()}}
    return out


def progress(root: Path = ROOT) -> str:
    made_from = digest(root)
    kept = {name: done(root, made_from, name) for name, *_ in conditions()}
    total = sum(len(v) for v in kept.values())
    lines = [f"draws kept: {total} of {len(conditions()) * DRAWS}"]
    for name, kind, size, noise in conditions():
        if kept[name]:
            s = summary(list(kept[name].values()), kind)["rungs"]["0.1"]
            lines.append(f"{name}: {len(kept[name])} draws, at 0.1 widths the sign names it on {s['right']:.2f}, the other kernel on {s['wrong']:.2f}")
    return "\n".join(lines)


def table(root: Path = ROOT, path: Path | None = None) -> dict:
    """The chances into toy/matrix.json, or `path`, under `sign`, beside what that file holds. Refused while a draw is not
    kept."""
    made_from = digest(root)
    out = {"digest": made_from, "draws": DRAWS, "sizes_mhz": list(SIZES), "noises": list(NOISES), "rungs": list(RUNGS), "conditions": {}}
    for name, kind, size, noise in conditions():
        got = done(root, made_from, name)
        if len(got) < DRAWS:
            raise SystemExit(f"REFUSED: {name} holds {len(got)} of {DRAWS} draws under the files as they are: add toy/sign.py to the queue again")
        out["conditions"][name] = {"kind": kind, "size_mhz": size, "noise": noise, **summary([got[d] for d in sorted(got)], kind)}
    path = Path(path or Path(root) / "toy" / "matrix.json")
    held = json.loads(path.read_text()) if path.exists() else {}
    held["sign"] = out
    path.write_text(json.dumps(held, indent=1) + "\n")
    return out


def check(root: Path = ROOT, path: Path | None = None) -> list:
    path = Path(path or Path(root) / "toy" / "matrix.json")
    held = json.loads(path.read_text()).get("sign") if path.exists() else None
    if held is None:
        return ["toy/matrix.json holds no `sign`: `python3 tools/autoqueue.py add toy/sign.py`, then `python toy/sign.py table`"]
    return [] if held.get("digest") == digest(root) else ["the `sign` of toy/matrix.json was made from other files than those on disk now: count it again"]


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["units"]:
        print(json.dumps(units()))
    elif args[:1] == ["slice"] and len(args) >= 4:
        if not os.environ.get("TOY_QUEUE"):
            sys.exit("REFUSED: a slice is a wave of the queue: `python3 tools/autoqueue.py add toy/sign.py` (his message 147)")
        print(json.dumps(slice_of(args[1], int(args[2]), int(args[3]), float(args[5]) if args[4:5] == ["--for"] else float("inf"))))
    elif args == ["progress"]:
        print(progress())
    elif args == ["digest"]:
        print(digest())
    elif args == ["table"]:
        got = table()
        for name, c in got["conditions"].items():
            print(f"{name}: " + ", ".join(f"at {rung} widths right {r['right']:.2f} wrong {r['wrong']:.2f}" for rung, r in c["rungs"].items()))
    elif args == ["--check"]:
        wrong = check()
        print("\n".join(wrong) or "the sign of the kernel: made from the files as they are")
        sys.exit(1 if wrong else 0)
    else:
        print(__doc__)
        sys.exit(1)
