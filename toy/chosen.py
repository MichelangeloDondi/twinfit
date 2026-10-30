"""The loop of the package on the toy, with the choices of the toy's physics written beforehand (row M87 of the plan).

The loop without an agent (twinfit/loop.py) reads a folder: a description, the file of a first model, and the term files
of a researcher. Here the researcher's choices are those of the toy, written before any run: the four terms of details
2 to 5, in their order, each with the sentence of its physics, and last a fifth term that the toy lacks, `ABSENT`, which
the loop should refuse. Every description holds his ansatz on the waist, 43 +- 1 um. The first model of the loop takes
the waist from the constants of the toy and not from its description, so the bar that the ansatz gives is zero on every
run: the count says so, and the output does not yet say that a held value moves nothing (review 22, finding 5). The
folder is written for three sets of traces:

    one draw     the 8 synthetic traces of the fifth detail: does the loop keep the four terms, and reach the fifth model
    his traces   his 8 recorded traces, as toy/description.toml describes them: which terms the loop keeps
    draws        `DRAWS` draws of the fifth detail with the noise of neighbouring samples, and `LAW_DRAWS` with the toy's
                 law: how often the loop keeps each term, refuses the fifth, and how often the interval of one bar of
                 each parameter of its best model holds the value that made the traces

Each run reads only what keeps or refuses a term, with no report of where to search and no proposals (loop.run with
report=False), but for the two runs named with report, on one draw and on his traces with the toy's law, which read
the report and the six proposals of each model. The count of where the report points is toy/pointing.py.

    python3 tools/autoqueue.py add toy/chosen.py    the loop through the queue (his message 147): eight conditions,
                                                    each set of traces with the noise of neighbouring samples and,
                                                    named with law, with the noise of each sample from the toy's law as
                                                    a column of the description. Review 21 found that the noise read
                                                    from neighbouring samples is several times too small at the top of
                                                    a line, which no twin of the loop can see. A run is a unit, kept in
                                                    runs/chosen_units as it ends
    python toy/chosen.py table                      writes toy/chosen.json from the units kept
    python toy/chosen.py --check                    ends without error when toy/chosen.json was made from the files as they are

Two choices of the session: the first model holds the natural width at its value of the literature, since a free
natural width and the width of the permeated gas are one parameter on these traces. The values of theory, the collision
coefficient and the rate of the transit, are held at the centre that theory gives, `HELD`, and their bars reach the
other values, as in the main fits of the folder (rule 10, row M133): until version 0.0.4.dev7 of the package the loop
fitted every term free.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from twinfit import description, loop  # noqa: E402
from twinfit.mechanisms import Mechanism  # noqa: E402

from toy import fits, model, real, traces  # noqa: E402

ORDER = ("collisions", "gas_width", "laser_noise", "transit")          # the terms of details 2 to 5, in the order in which the toy adds them
HELD = ("collisions", "transit")          # the terms whose size theory gives, held at it by the loop and not fitted: rule 10, his messages 133 and 135, row M133
WHY = {"collisions": "collisions between Rb atoms widen the line in proportion to the density of the vapour, which grows with the temperature",
       "gas_width": "a gas that permeated the cell widens the line by one Lorentzian width, the same in every trace",
       "laser_noise": "the noise of the frequency of the laser widens the line by a Gaussian of one width, the same in every trace",
       "transit": "the atoms cross the beam in a time that falls as the root of the temperature rises, which widens the line"}
ABSENT = ("a pedestal four times wider than the line", "atoms that change their velocity in a collision could put a pedestal under the line, wider than it")          # a term that the toy lacks, offered last: the loop should refuse it (review 21)
ABSENT_TERM = '''from twinfit import mechanisms

term = mechanisms.pedestal(4.0, lower=0.0)          # it adds to the line, or it is not there
term.name = {name!r}
why = {why!r}
'''
REASONS = {"its_bound": "its fit ended at the bound of its own parameter", "the_criterion": "the Bayesian criterion "}          # how twinfit/loop.py opens the two reasons for which it refuses a term
DRAWS, LAW_DRAWS = 20, 20          # 20 with the toy's law, the first stage of rows M474 and M481, whose second stage of 100 read the bar of β_self 0.82 at 0.1.0rc16. Until then 100: the draws with the noise of neighbouring samples, which 20 show wrong, and with the toy's law, where a share near 68 per cent is then known to 0.060: 100 would end after the night. With the toy's law 60 until version 0.0.5, and 100 read a share near 80 of 100 within 8 points, row M204
SETS = 30


def design_of(design: dict) -> dict:
    """The covariates of the toy, from those of the description."""
    return {"T": design["temperature"], "P": design["power"]}


def mechanism_of(name: str, const: dict, config: dict) -> Mechanism:
    """A term of the toy as the loop takes it: the model of the detail that the term enters, with the shared parameters of
    the model before and the size of the term free, every other value held at config.toml."""
    entry = config["parameter"][name]

    def wrap(_, problem, n):
        of = fits.curve(entry["enters"], const, config, list(problem.shared) + [name])
        return lambda x, shared, own, design: of(x, shared, own, design_of(design))
    theory = (float(entry["true"]), float(entry["uncertainty"])) if entry.get("kind") == "theory" else None          # a value of theory is held at its centre, with its bar (rule 10, row M133)
    return Mechanism(config["detail"][str(entry["enters"])]["name"], "a term of the toy", name, float(entry["start"]), wrap,
                     lower=float(entry.get("lower", 0.0)), upper=float(entry.get("upper", np.inf)), theory=theory)


MODEL = '''"""The first model of the loop on the toy: the line of the first detail, with the natural width held at its value of the
literature, written by toy/chosen.py."""
import sys
sys.path.insert(0, {root!r})
from toy import chosen, fits, model as M

const, config = M.load()
SHARED, OWN, START_SHARED, CENTRE = (), fits.OWN, [], "centre"
NAME = "the line of the natural width, held at its value of the literature"
_of = fits.curve(1, const, config, [])


def model(x, shared, own, design):
    return _of(x, shared, own, chosen.design_of(design))


def start_own(record):
    return fits.start_of(record)


def centre(shared, own, design):
    return own[1]
'''
TERM = '''import sys
sys.path.insert(0, {root!r})
from toy import chosen, model as M

const, config = M.load()
term = chosen.mechanism_of({name!r}, const, config)
why = {why!r}
hold = {hold!r}
'''


def folder(home: Path, description: str) -> Path:
    """A folder for the loop: the description given, the first model and the four term files."""
    home = Path(home)
    (home / "terms").mkdir(parents=True, exist_ok=True)
    (home / "description.toml").write_text(description)
    (home / "model.py").write_text(MODEL.format(root=str(ROOT)))
    for k, name in enumerate(ORDER, start=1):
        (home / "terms" / f"{k:02d}.py").write_text(TERM.format(root=str(ROOT), name=name, why=WHY[name], hold=name in HELD))
    (home / "terms" / f"{len(ORDER) + 1:02d}.py").write_text(ABSENT_TERM.format(name=ABSENT[0], why=ABSENT[1]))
    return home


def law_sd(row: dict, const: dict) -> np.ndarray:
    """The noise of each sample of a trace by the toy's law of the noise (review 21): the one that a synthetic trace was
    drawn with, and for a recorded trace the law of its condition at the height of the trace, smoothed, above its ends."""
    if "sd" in row:
        return np.asarray(row["sd"], dtype=float)
    y, nu = np.asarray(row["y"], dtype=float), np.asarray(row["nu"], dtype=float)
    smooth = np.convolve(y, np.ones(31) / 31.0, mode="same")
    signal = np.clip(smooth - np.median(np.concatenate([y[:40], y[-40:]])), 0.0, None)
    return traces.noise_sd(row["P"], row["T"], signal, nu, const)


def described(rows: list, home: Path, name: str, const: dict | None = None) -> str:
    """The description of 8 traces of the toy, each written as a file beside it, with the levers of his description.
    With `const`, each file holds a third column, the noise of each sample by the toy's law, and the description reads
    the noise from it. Without, the noise is read from neighbouring samples, one number a trace."""
    (home / "traces").mkdir(parents=True, exist_ok=True)
    text = ['[problem]', f'name = "{name}"', 'what = "the two photon line of rubidium from 5S to 6S, one hyperfine peak, in 8 traces at four temperatures and five powers"',
            'abscissa = "frequency from the top of the line"', 'abscissa_unit = "MHz"', 'value = "signal"', 'value_unit = "V"', '',
            '[covariates.temperature]', 'unit = "C"', '', '[covariates.power]', 'unit = "mW"', '', '[covariates.waist]', 'unit = "um"', '']
    for k, row in enumerate(rows):
        columns = [row["nu"], row["y"]] + ([law_sd(row, const)] if const is not None else [])
        (home / "traces" / f"{k}.csv").write_text(("nu,y,sd\n" if const is not None else "nu,y\n") + "".join(",".join(repr(float(v)) for v in values) + "\n" for values in zip(*columns)))
        text += ['[[trace]]', f'file = "traces/{k}.csv"', f'temperature = {row["T"]}', f'power = {row["P"]}', 'waist = 43.0', '']
    if const is not None:
        text += ['[noise]', 'from = "column"', 'column = "sd"', '']
    his_text = (ROOT / "toy" / "description.toml").read_text()
    ansatz = "[addenda.ansatz." + his_text.split("[addenda.ansatz.", 1)[1].split("\n\n", 1)[0]          # his ansatz on the waist, the same dignity
    levers = his_text.split("[addenda.lever.", 1)[1]          # the levers of his description, the same experiment
    setup = ("[setup." + his_text.split("[setup.", 1)[1].split("\n[noise]", 1)[0].rstrip()) if "[setup." in his_text else ""          # the setup that his levers of the setup name: left out until review 28, and the package refused every run of the count, which kept the refusals as its rows
    return "\n".join(text + ([setup, ''] if setup else []) + ['[addenda]', 'model = "model.py"', '', ansatz, '', '[addenda.lever.' + levers])


def his() -> str:
    """The description of his 8 recorded traces, with the paths of their files made whole and the model of the loop."""
    text = (ROOT / "toy" / "description.toml").read_text().replace('file = "../', f'file = "{ROOT}/')
    return text.replace("[addenda.lever.repetitions]", '[addenda]\nmodel = "model.py"\n\n[addenda.lever.repetitions]', 1)


def read(out: dict, config: dict) -> dict:
    """What a run of the loop says: the terms kept and refused, and each parameter of the best model against the value
    with which the twin is made, in bars: on a synthetic draw that value is the truth, on his traces it is not."""
    twin = model.truth(config)
    parameters = {}
    for name, cell in out["output"]["parameters"].items():
        bar, theory = cell["bar_of_the_fit"], cell["bar_of_theory"]["together"]
        whole = None if not bar else float(np.hypot(bar, theory))          # the bar of the fit and that of the held waist, added in squares
        of = twin.get(name)          # a term that the toy lacks has no value of the twin: on his traces with the noise of neighbouring samples the loop kept the pedestal offered last, 0.0170 +- 0.0042, and the count stopped on its name (0.1.0rc21)
        parameters[name] = {"value": cell["value"], "bar": bar, "bar_of_theory": theory, "of_the_twin": of,
                            "from_the_twin_in_bars": None if not bar or of is None else (cell["value"] - of) / bar,
                            "from_the_twin_in_whole_bars": None if not whole or of is None else (cell["value"] - of) / whole}
    return {"kept": [k["file"] for k in out["kept"]], "refused": {r["file"]: r["why"] for r in out["refused"]}, "best": out["best"], "parameters": parameters}


def freed(home: Path, config: dict, fit_of=None) -> dict:
    """The terms that the loop holds at theory, freed as a test once it has ended, his message 340 and row M480: the model
    that the folder's term files give, `loop.model_of`, fitted again in a copy with each held term free, and each freed
    value against the value of the twin in its own bar. None for a held term that the criterion leaves out of the model."""
    import shutil
    import tempfile
    twin = model.truth(config)
    fit_of = fit_of or loop.model_of
    with tempfile.TemporaryDirectory() as tmp:
        copy = Path(tmp) / "freed"
        shutil.copytree(home, copy, ignore=shutil.ignore_patterns("loop", "*.pdf", "*.png"))
        for path in sorted((copy / "terms").glob("[0-9][0-9].py")):
            path.write_text(path.read_text().replace("hold = True", "hold = False"))
        problem, result = fit_of(copy)[:2]
    names, out = list(problem.shared), {}
    for name in HELD:
        k = names.index(name) if name in names else None
        bar = float(result.sd[k]) if k is not None and result.sd is not None else None
        out[name] = None if not bar else {"value": float(result.shared[k]), "bar": bar, "of_the_twin": twin[name], "from_the_twin_in_bars": (float(result.shared[k]) - twin[name]) / bar}
    return out


def reported(home: Path) -> dict:
    """What the report and the proposals say on each model of a run with them: the first answer of each question beyond
    the twin, the terms proposed for the model, and the levers proposed for the experiment with their error over that of
    the records."""
    out = {}
    for path in sorted((Path(home) / "loop").glob("model_*.json")):
        held = json.loads(path.read_text())
        six = held.get("proposals") or {}
        out[path.stem.replace("_", " ")] = {"first": (held.get("where") or {}).get("first_of_each_question"),
                                            "for_the_model": [r["term"] for r in six.get("for_the_model", [])],
                                            "for_the_experiment": {r["lever"]: r["error_over_that_under_the_records"] for r in six.get("for_the_experiment", [])}}
    return out


CONDITIONS = {"one": ("one", False, False), "his": ("his", False, False), "draws": ("draws", False, False), "one law": ("one", True, False), "his law": ("his", True, False), "draws law": ("draws", True, False),
              "one report": ("one", True, True), "his report": ("his", True, True)}          # what the queue runs: a set of traces, a way of giving their noise, and whether the report and the proposals are read


def unit(condition: str, k: int, root: Path = ROOT) -> dict:
    """One run of the loop: on one draw, on his traces, or on draw `k` of the fifth detail, with the noise of neighbouring
    samples or, in a condition named with law, with the toy's noise law as a column of the description."""
    kind, law, report = CONDITIONS[condition]
    const, config = model.load(Path(root) / "toy", root)
    base = Path(root) / "runs" / ("chosen_report" if report else "chosen_law" if law else "chosen")
    if kind == "his" and not law:
        home = folder(base / "his", his())
    else:
        home = base / (kind if kind != "draws" else f"draw_{k}")
        home.mkdir(parents=True, exist_ok=True)
        rows = real.rows(const, config, root) if kind == "his" else fits.synthetic(5, const, config, 1 if kind == "one" else k)
        name = "his 8 recorded traces" if kind == "his" else f"the synthetic traces of the fifth detail, draw {1 if kind == 'one' else k}"
        folder(home, described(rows, home, name, const if law else None))
    description.read(home / "description.toml")          # a description that the package refuses is a fault of the count and ends the wave, and not a result of a run: on 2026-10-17 a hundred refusals were kept as rows
    try:
        got = read(loop.run(home, n=len(ORDER) + 1, sets=SETS, report=report), config)
        got["freed"] = freed(home, config)          # the held terms freed as a test, his message 340, row M480
        return dict(got, report=reported(home)) if report else got
    except (Exception, SystemExit) as why:
        return {"failed": f"{type(why).__name__}: {why}"[:300]}


def units_of(condition: str) -> int:
    """How many units a condition has: its draws, or one run."""
    kind, law, _ = CONDITIONS[condition]
    return (LAW_DRAWS if law else DRAWS) if kind == "draws" else 1


def unit_path(root: Path, condition: str, k: int) -> Path:
    return Path(root) / "runs" / "chosen_units" / f"{condition.replace(' ', '_')}_{k}.json"


def done(root: Path, condition: str, made_from: str) -> dict:
    """The units of a condition kept under the digest of the files as they are, by their number."""
    out = {}
    for k in range(1, units_of(condition) + 1):
        path = unit_path(root, condition, k)
        if path.exists():
            kept = json.loads(path.read_text())
            if kept.get("digest") == made_from:
                out[k] = kept["result"]
    return out


def units(root: Path = ROOT) -> list:
    """What the queue asks first: each condition with the units that are not yet kept."""
    made_from = digest(root)
    return [{"condition": c, "units": [k for k in range(1, units_of(c) + 1) if k not in done(root, c, made_from)], "after": [], "kind": kind + (" law" if law else "") + (" report" if report else "")}
            for c, (kind, law, report) in CONDITIONS.items()]


def slice_of(condition: str, first: int, last: int, limit: float = float("inf"), root: Path = ROOT) -> dict:
    """The units of a condition from `first` to `last` that are not yet kept, each kept as it ends, until `limit` seconds
    have passed: a wave of the queue, or a pre-wave."""
    import resource
    made_from, began, ids, seconds = digest(root), time.time(), [], []
    for k in range(first, last + 1):
        if k in done(root, condition, made_from):
            continue
        t = time.time()
        result = unit(condition, k, root)
        path = unit_path(root, condition, k)
        path.parent.mkdir(parents=True, exist_ok=True)
        aside = path.with_suffix(".part")
        aside.write_text(json.dumps({"digest": made_from, "result": result}, indent=1) + "\n")
        aside.replace(path)
        ids.append(k)
        seconds.append(time.time() - t)
        if time.time() - began > limit:
            break
    return {"done": ids, "seconds": seconds, "rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


def progress(root: Path = ROOT) -> str:
    """The preliminary result after a wave: the units kept, and on the draws so far the share of intervals of one bar that
    hold the value of the twin."""
    made_from, said = digest(root), []
    for c, (kind, _, _) in CONDITIONS.items():
        got = done(root, c, made_from)
        if kind == "draws" and got:
            s = summary({f"draw {k}": v for k, v in got.items()})
            said.append(f"{c}: {s['draws']} draws, within one bar " + ", ".join(f"{n} {v['share']:.2f}" for n, v in s["within_one_bar"].items() if v))
        else:
            said.append(f"{c}: {len(got)} kept")
    return "; ".join(said)


RUNS = ("toy/chosen.py",)          # the program whose imports say which files of the package a run runs: the loop and what it imports
MAKES = ("freed", "ORDER", "HELD", "WHY", "ABSENT", "ABSENT_TERM", "SETS", "design_of", "mechanism_of", "MODEL", "TERM", "folder", "law_sd", "described", "his", "read", "reported", "CONDITIONS", "unit")          # what makes a run


def digest(root: Path = ROOT) -> str:
    """What a run depends on: the kept fits of the last detail, the parts of this program that make a run, the
    description of his traces and the files of the package that the loop imports, as code with no docstring (rule 16,
    his message 175: every file whole made a word of the tutorial a day of counts)."""
    from toy import draws, power
    root = Path(root)
    h = hashlib.sha256(draws.digest(root).encode())
    h.update(power.sources((root / "toy" / "chosen.py").read_text(), MAKES))
    draws.hashed(h, root, ["toy/description.toml"] + draws.package_of(root, RUNS))
    return h.hexdigest()


def summary(draws_of: dict) -> dict:
    """Of the runs on many draws: the share of draws on which each term file was kept, and for each parameter of the best
    model the share on which the interval of one bar, and of two, holds the true value. A draw whose run failed counts in
    neither and is named."""
    ran = {k: v for k, v in draws_of.items() if "failed" not in v}
    out = {"draws": len(ran), "failed": sorted(k for k in draws_of if k not in ran), "kept": {}, "within_one_bar": {}, "within_two_bars": {}, "within_one_whole_bar": {}, "reached_the_fifth": None, "refused_the_absent": None}
    if not ran:
        return out
    absent = f"{len(ORDER) + 1:02d}.py"
    out["reached_the_fifth"] = sum(sum(f"{k:02d}.py" in v["kept"] for k in range(1, len(ORDER) + 1)) == len(ORDER) for v in ran.values()) / len(ran)
    out["refused_the_absent"] = sum(absent in v["refused"] for v in ran.values()) / len(ran)          # the term that the toy lacks
    whys = [v["refused"][absent] for v in ran.values() if absent in v["refused"]]          # two things: a fit of its size that ends at its bound, and the criterion (review 22, minor 13)
    out["refused_the_absent_by"] = {key: sum(w.startswith(opening) for w in whys) / len(ran) for key, opening in REASONS.items()}          # by the opening that the package writes, and not by a word: "bound" stands in a third reason too (review 24, minor 6)
    out["refused_the_absent_by"]["another_reason"] = sum(not any(w.startswith(opening) for opening in REASONS.values()) for w in whys) / len(ran)
    for k, name in enumerate(ORDER, start=1):
        out["kept"][name] = sum(f"{k:02d}.py" in v["kept"] for v in ran.values()) / len(ran)
    for name in ORDER:
        off = [v["parameters"][name]["from_the_twin_in_bars"] for v in ran.values() if name in v["parameters"] and v["parameters"][name]["from_the_twin_in_bars"] is not None]
        if name in HELD:          # a term held at theory is scored by the fit that frees it as a test, his message 340, row M480
            off = [v["freed"][name]["from_the_twin_in_bars"] for v in ran.values() if (v.get("freed") or {}).get(name)]
        out["within_one_bar"][name] = None if not off else {"share": sum(abs(o) < 1.0 for o in off) / len(off), "of": len(off)}
        out["within_two_bars"][name] = None if not off else {"share": sum(abs(o) < 2.0 for o in off) / len(off), "of": len(off)}
        whole = [v["parameters"][name].get("from_the_twin_in_whole_bars") for v in ran.values() if name in v["parameters"]]
        whole = [o for o in whole if o is not None]
        out["within_one_whole_bar"][name] = None if not whole else {"share": sum(abs(o) < 1.0 for o in whole) / len(whole), "of": len(whole)}
    return out


def table(root: Path = ROOT) -> dict:
    """toy/chosen.json, from the units kept in runs/chosen_units under the digest of the files as they are. Refused while
    a unit is not kept: the queue has not ended."""
    from toy import power
    root, made_from = Path(root), digest(root)
    got = {c: done(root, c, made_from) for c in CONDITIONS}
    lacking = {c: n for c, n in ((c, units_of(c) - len(v)) for c, v in got.items()) if n}
    if lacking:
        raise SystemExit("REFUSED: units are not kept under the files as they are: " + ", ".join(f"{c} {n}" for c, n in lacking.items()) + ". Add toy/chosen.py to the queue again")
    draws_of = lambda c: {f"draw {k}": v for k, v in sorted(got[c].items())}          # noqa: E731
    from twinfit.version import __version__
    out = {"digest": made_from, "package": __version__, "one": got["one"][1], "his": got["his"][1], "draws": draws_of("draws"), "summary": summary(draws_of("draws")),          # the version under the table (the red team of his message 331, severe 1). Until 2026-10-23 this comment stood before the his run, the draws and their summary, and left them out of the table
           "law": {"one": got["one law"][1], "his": got["his law"][1], "draws": draws_of("draws law"), "summary": summary(draws_of("draws law"))},
           "report": {"one": got["one report"][1], "his": got["his report"][1]}, "queue": power.queue_of(root, "toy/chosen.py", made_from)}          # with the waves that the queue spent on its runs (review 22, minor 18)
    (root / "toy" / "chosen.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def check(root: Path = ROOT) -> list:
    path = Path(root) / "toy" / "chosen.json"
    if not path.exists():
        return ["toy/chosen.json is not written: run `python toy/chosen.py one`, `his` and `draws`, then `table`"]
    held = json.loads(path.read_text())
    if held.get("digest") != digest(root):
        return ["toy/chosen.json was made from other files than those on disk now: run the loop on the toy again"]
    same = held["summary"] == summary(held["draws"]) and held["law"]["summary"] == summary(held["law"]["draws"])
    return [] if same else ["toy/chosen.json does not hold the summary of its own draws"]


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["units"]:
        print(json.dumps(units()))
    elif args[:1] == ["slice"] and len(args) >= 4:
        if not os.environ.get("TOY_QUEUE"):
            sys.exit("REFUSED: a slice is a wave of the queue: `python3 tools/autoqueue.py add toy/chosen.py` (his message 147)")
        print(json.dumps(slice_of(args[1], int(args[2]), int(args[3]), float(args[5]) if args[4:5] == ["--for"] else float("inf"))))
    elif args == ["progress"]:
        print(progress())
    elif args == ["digest"]:          # what the queue judges a wave by: a change that moves no unit makes none stale
        print(digest())
    elif args[:1] and args[0] in ("one", "his", "draws"):
        sys.exit("REFUSED: the loop on the toy runs through the queue: `python3 tools/autoqueue.py add toy/chosen.py` (his message 147)")
    elif args[:1] == ["table"]:
        print(json.dumps(table()["summary"], indent=1))
    elif args[:1] == ["--check"]:
        wrong = check()
        print("\n".join(wrong) or "the loop on the toy: made from the files as they are")
        sys.exit(1 if wrong else 0)
    else:
        print(__doc__)
        sys.exit(1)
