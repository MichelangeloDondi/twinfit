"""The rehearsal of the blind test (the owner's answer of 2026-10-13 to "May a rehearsal of the blind test run now, on
the toy of four or five terms?": "Yes, all three arms").

The aim of the package is that a fresh chat, with the package and the traces alone, names the terms of the model and
their sizes. Here that is tried in small, on the first tiers of the toy, where the truth is known, in three arms:

    A    an agent with the package and its protocol
    B    an agent with a plain fit and what the fit leaves, and no package
    C    a script with no agent: the terms of the toy fitted in turn and chosen by BIC. It knows the list of the terms

A case is one draw of the 8 synthetic traces of a detail. An agent gets the traces, what an experimenter knows of them,
a first model that is a line of one free Lorentzian width, and a list of kinds of term to name: more kinds than the toy
holds, and the true ones among them. So this is the level with hints, where the terms are named for the agent. It does
not get the detail of the case, the values, or any file of this repository.

    python toy/rehearsal.py make OUT       writes the folders of the cases into OUT, which is outside the repository, and
                                           the truth of each case into runs/rehearsal/truth.json, which no agent is given
    python toy/rehearsal.py script         arm C on every case: writes runs/rehearsal/script.json
    python toy/rehearsal.py score OUT      reads the answer of each folder of OUT and writes runs/rehearsal/scores.json
    python toy/rehearsal.py keep OUT       writes toy/rehearsal.json, the record that the repository keeps: for each case
                                           its truth, and for each arm its answer as it was handed back, its score, and
                                           for an agent the time and the tokens that it took
    python toy/rehearsal.py --check        ends without error when the record holds every case in its three arms, was
                                           scored by the scorer as it is, and gives the same scores when counted again

No number of the physics is written here.
"""
from __future__ import annotations

import hashlib
import json
import math
import shutil
import sys
import tomllib
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import twinfit  # noqa: E402
from twinfit import noise as N  # noqa: E402

from toy import fits, model, traces  # noqa: E402

FIRST = 9001          # the draws of the cases: none is a draw of the check on 100 draws, nor of the counts
KINDS = {          # the kinds of term that an answer may name: the four of the first tiers of the toy, and five that the toy does not hold
    "collisions": "a Lorentzian width that grows with the density of the vapour. Size: MHz for 1e15 atoms per cm^3",
    "constant": "a constant Lorentzian width beyond the natural width, the same for all traces. Size: MHz",
    "power_width": "a Lorentzian width that grows with the power. Size: MHz at the highest power",
    "gaussian": "a Gaussian spread of the frequency, the same for all traces. Size: its standard deviation in MHz",
    "gaussian_power": "a Gaussian spread that grows with the power. Size: its standard deviation in MHz at the highest power",
    "cusp": "a spread with exponential wings, exp(-|nu| / b), whose rate b follows the square root of the temperature in kelvin. Size: b in MHz at the highest temperature",
    "asymmetry": "a line that leans to one side. Size: the share of its height by which the two sides differ at half height",
    "second_line": "a second line beside the first. Size: its height as a share of the first",
    "baseline": "a baseline that is not a straight line. Size: its largest departure from a straight line, as a share of the height",
}
OF_TERM = {"collisions": "collisions", "gas_width": "constant", "laser_noise": "gaussian", "transit": "cusp"}          # the kind of each term of the toy
ALIKE = [{"gaussian", "cusp"}, {"gaussian_power"}, {"constant"}, {"collisions"}, {"power_width"}]          # the kinds that 8 traces cannot tell apart stand in one class: `power.map_of` gives the Gaussian spread and the spread with exponential wings a cosine of 0.999
AGENTS = ("A", "B")
ARMS = AGENTS + ("C",)          # and the script
RECORD = "toy/rehearsal.json"


def cases(config: dict) -> list:
    """The cases, each as (its name, the detail of its traces, its draw). A name says nothing of the detail: the letters
    follow the digest of the detail and the draw."""
    made = [(d, FIRST + k) for k, d in enumerate(config["folder"]["details"])]
    order = sorted(made, key=lambda c: hashlib.sha256(f"{c[0]} {c[1]}".encode()).hexdigest())
    return [(f"case {chr(ord('a') + k)}", d, draw) for k, (d, draw) in enumerate(order)]


def truth_of(detail: int, config: dict) -> dict:
    """The terms that the traces of a detail hold beyond the natural width, by their kind, with the size of each."""
    return {OF_TERM[name]: float(model.truth(config)[name]) for name, e in config["parameter"].items() if 1 < e["enters"] <= detail}


def known(const: dict, config: dict) -> dict:
    """What an experimenter knows of the traces and an agent is given: the conditions, the noise law of each trace, the
    natural width from the lifetime, and the density of the vapour at each temperature."""
    natural = config["parameter"]["natural_width"]
    rows = []
    for P, T in traces.conditions(const):
        law = next(r for r in const["noise"]["trace"] if (r["power_mW"], r["temperature_C"]) == (P, T))
        rows.append({"file": traces.name_of(P, T), "power_mW": P, "temperature_C": T, "density_per_cm3": float(model.density(T, const)),
                     "noise": {"a_V": law["a_V"], "b_V": law["b_V"], "c": law["c"], "d_MHz2": const["noise"]["d_MHz2"]}})
    return {"axis": "the frequency of the two photon transition, in MHz from the middle of the sweep",
            "noise_law": "variance = a_V^2 + b_V h + c h^2 + d_MHz2 s^2, with h the height of the line above its baseline in V and s its slope in V per MHz",
            "natural_width_MHz": float(natural["true"]), "natural_width_uncertainty_MHz": float(natural["uncertainty"]),
            "natural_width_source": "the lifetime of the upper state, from the literature: the full width at half maximum of a Lorentzian",
            "traces": rows, "kinds": KINDS}


def make(out: Path, root: Path = ROOT, texts: dict | None = None) -> dict:
    """The folders of the cases, one a case and an arm of an agent, in `out`, and the truth of each case, which goes into
    no folder. `texts`: the instructions of each arm. Refuses a folder inside the repository."""
    out, root = Path(out).resolve(), Path(root).resolve()
    if out == root or root in out.parents:
        raise SystemExit(f"REFUSED: {out} is inside the repository: an agent of the rehearsal reads no file of it")
    const, config = model.load(root / "toy", root)
    texts = texts or {}
    truth = {}
    for name, detail, draw in cases(config):
        truth[name] = {"detail": detail, "draw": draw, "terms": truth_of(detail, config)}
        for arm in AGENTS:
            home = out / f"{name}, arm {arm}"
            (home / "traces").mkdir(parents=True, exist_ok=True)
            for P, T in traces.conditions(const):
                t = traces.trace(detail, P, T, const, config, draw)
                (home / "traces" / traces.name_of(P, T)).write_text("frequency_MHz,signal_V\n" + "".join(f"{x:.6f},{y:.9e}\n" for x, y in zip(t["nu"], t["y"])))
            (home / "known.json").write_text(json.dumps(known(const, config), indent=1) + "\n")
            (home / "READ FIRST.md").write_text(texts.get(arm, ""))
            if arm == "A":
                if (home / "twinfit").exists():
                    shutil.rmtree(home / "twinfit")
                shutil.copytree(root / "twinfit", home / "twinfit", ignore=shutil.ignore_patterns("__pycache__"))
    return truth


def with_terms(const: dict, config: dict, free: list):
    """The model with the terms of `free` alone: every other term of the toy is absent, at zero."""
    last = max(int(k) for k in config["detail"])
    base = {name: 0.0 for name in config["parameter"]}

    def of(x, shared, own, design):
        values = dict(base, **{name: float(v) for name, v in zip(free, shared)})
        u = np.asarray(x, dtype=float) - own[1]
        top = float(model.shape(last, 0.0, design["T"], values, const))
        return own[0] * model.shape(last, u, design["T"], values, const) / top + own[2] + own[3] * u
    return of


def problem_of(rows: list, const: dict, config: dict, free: list) -> twinfit.Problem:
    """The joint problem of some traces under the model with the terms of `free` alone."""
    entries = [config["parameter"][name] for name in free]
    records = []
    for row in rows:
        y = np.asarray(row["y"], dtype=float)
        first = float(np.std(np.concatenate([y[:40], y[-40:]]), ddof=1))
        records.append(twinfit.Record(row["nu"], y, max(first, 1e-12), {"P": row["P"], "T": row["T"]}, name=traces.name_of(row["P"], row["T"])))
    return twinfit.Problem(records, with_terms(const, config, free), tuple(free), fits.OWN, [e["start"] for e in entries], fits.start_of,
                           centre=lambda s, o, d: o[1], name="with " + ", ".join(free), bounds_shared=([e["lower"] for e in entries], [e["upper"] for e in entries]))


def selection(rows: list, const: dict, config: dict, start: tuple = ("natural_width",), steps: int = 40) -> dict:
    """Arm C: forward selection by BIC over the terms of the toy. From the line of one free width, every term that the
    model lacks is fitted in turn at the weights of the model as it stands, the one with the largest fall of the chi
    square is added when that fall passes the logarithm of the number of samples, and the weights are settled again.
    Returns the terms in the order in which they were added, the fall of each, and the last fit."""
    free, added = list(start), []
    while True:
        settled, result, _ = N.settle(problem_of(rows, const, config, free), fits.variance_of(const))
        gains = {}
        for term, entry in config["parameter"].items():
            if term in free:
                continue
            entries = [config["parameter"][name] for name in free] + [entry]
            larger = settled.with_model(with_terms(const, config, free + [term]), tuple(free) + (term,), settled.own, list(result.shared) + [entry["start"]], list(result.own),
                                        bounds_shared=([e["lower"] for e in entries], [e["upper"] for e in entries]))
            gains[term] = float(result.chi2 - twinfit.fit(larger, list(result.shared) + [entry["start"]], result.own, max_iter=steps).chi2)
        best = max(gains, key=gains.get) if gains else None
        if best is None or gains[best] <= math.log(settled.samples):
            bars = np.sqrt(np.clip(np.diag(result.cov), 0.0, None))
            return {"added": added, "left_out": gains, "values": {n: float(v) for n, v in zip(free, result.shared)}, "bars": {n: float(b) for n, b in zip(free, bars)},
                    "chi2": float(result.chi2), "dof": int(result.dof)}
        added.append({"term": best, "gain": gains[best], "others": {k: v for k, v in gains.items() if k != best}})
        free.append(best)


def answer_of(selected: dict, config: dict) -> dict:
    """The answer of arm C in the form of an agent's: the kinds that it added, each with its fitted size. The free
    natural width, less the one of the lifetime, is its constant width, named when it is beyond two of its uncertainties."""
    natural = config["parameter"]["natural_width"]
    terms = [{"kind": OF_TERM[a["term"]], "size": selected["values"][a["term"]], "uncertainty": selected["bars"][a["term"]]} for a in selected["added"] if a["term"] in OF_TERM]
    more, bar = selected["values"]["natural_width"] - natural["true"], math.hypot(selected["bars"]["natural_width"], natural["uncertainty"])
    if more > 2.0 * bar and not any(t["kind"] == "constant" for t in terms):
        terms.append({"kind": "constant", "size": more, "uncertainty": bar})
    return {"terms": terms}


def score(answer: dict, truth: dict) -> dict:
    """An answer against the truth of its case. For each term of the truth: whether the answer names its kind, or a kind
    of its class, and how far its size lies from the true one, in its own uncertainty. And the terms that the answer
    names and that stand with no term of the truth."""
    named = {t["kind"]: t for t in answer.get("terms", []) if t.get("kind") in KINDS}
    of_class = lambda kind: next((c for c in ALIKE if kind in c), {kind})          # noqa: E731
    found = {}
    for kind, size in truth.items():
        got = named.get(kind)
        found[kind] = {"named": got is not None, "by_class": bool(of_class(kind) & set(named)), "true": size, "size": got.get("size") if got else None,
                       "off_in_uncertainties": float((got["size"] - size) / got["uncertainty"]) if got and got.get("uncertainty") else None}
    false = sorted(kind for kind in named if not any(of_class(kind) & of_class(true) for true in truth))
    return {"terms": found, "named": int(sum(f["named"] for f in found.values())), "by_class": int(sum(f["by_class"] for f in found.values())), "of": len(truth), "false": false,
            "right": bool(all(f["by_class"] for f in found.values()) and not false)}


def scorer() -> str:
    """The digest of what scores an answer: the scorer, the kinds that an answer may name and their classes. The record
    names it, so that an answer scored by another scorer than the one of the repository shows."""
    import inspect
    held = [inspect.getsource(score), sorted(sorted(c) for c in ALIKE), sorted(KINDS), sorted(OF_TERM.items())]
    return hashlib.sha256(json.dumps(held).encode()).hexdigest()


def gathered(out: Path, held: Path, config: dict) -> dict:
    """The record of the rehearsal, from the folders of the agents in `out` and from what `held` keeps: the truth of
    each case, the answers of the script, and what each agent took. For each case its truth, and for each arm that
    handed back its answer as it stands, with its score. An arm with no answer is left out, and the check says so."""
    truth = json.loads((Path(held) / "truth.json").read_text())
    script = json.loads((Path(held) / "script.json").read_text()) if (Path(held) / "script.json").exists() else {}
    took = json.loads((Path(held) / "took.json").read_text()) if (Path(held) / "took.json").exists() else {}
    record = {"scorer": scorer(), "kinds_named_to_an_agent": sorted(KINDS), "values": {name: float(v) for name, v in model.truth(config).items()}, "cases": {}}
    for name, entry in sorted(truth.items()):
        arms = {}
        for arm in AGENTS:
            home = Path(out) / f"{name}, arm {arm}"
            if (home / "answer.json").exists():
                answer = json.loads((home / "answer.json").read_text())
                arms[arm] = dict({"answer": answer, "score": score(answer, entry["terms"])}, **took.get(home.name, {}))
        if name in script:
            arms["C"] = {"answer": script[name]["answer"], "score": score(script[name]["answer"], entry["terms"])}
        record["cases"][name] = dict(entry, arms=arms)
    return record


def check(root: Path = ROOT) -> list:
    """What is wrong with the record of the rehearsal: an empty list when it holds every case in its three arms, when the
    truth of each case is the one that its detail gives at the values of the record, when the scorer is the one of this
    file, and when every answer gives its score again. A score that a hand wrote differs here."""
    path = Path(root) / RECORD
    if not path.exists():
        return [f"{RECORD} is not there"]
    record, wrong = json.loads(path.read_text()), []
    config = tomllib.loads((Path(root) / "toy" / "config.toml").read_text())          # the details and where each term enters: the check needs no trace
    if record["scorer"] != scorer():
        wrong.append("the answers were scored by another scorer than the one of toy/rehearsal.py")
    for name, detail, draw in cases(config):
        entry = record["cases"].get(name)
        if entry is None or (entry["detail"], entry["draw"]) != (detail, draw):
            wrong.append(f"{name} is not in the record with its detail and its draw")
            continue
        true = {OF_TERM[term]: float(record["values"][term]) for term, e in config["parameter"].items() if 1 < e["enters"] <= detail}
        if entry["terms"] != true:
            wrong.append(f"the truth of {name} is not the one that its detail gives at the values of the record")
        for arm in ARMS:
            if arm not in entry["arms"]:
                wrong.append(f"{name} has no answer of arm {arm}")
            elif entry["arms"][arm]["score"] != score(entry["arms"][arm]["answer"], entry["terms"]):
                wrong.append(f"the answer of {name}, arm {arm}, does not give its score when it is scored again")
    return wrong


COMMON = """# Eight traces of one spectral line

You are asked to analyse eight recorded curves of one spectral line. Work only in this folder. Read no other folder of
this machine and do not use the network. The Python to use is `{python}`: it has numpy, scipy and matplotlib.

## The measurement

One line of an atomic vapour, recorded while a frequency is swept, at eight conditions: four temperatures at one power,
and five powers at one temperature. `traces/` holds one file a trace, with the columns `frequency_MHz` and `signal_V`.
`known.json` holds what the experimenter knows: the power and the temperature of each trace, the density of the vapour
at each temperature, the noise law of each trace, and the natural width of the line, which is a Lorentzian width known
from the lifetime of the upper state.

## The question

The first model is a line of the natural width alone: a Lorentzian, with for each trace its own height, centre, offset
and slope of the baseline. Which physical terms does the line hold beyond that, and how large is each? `known.json`
lists under `kinds` the kinds of term that you may name, with the unit of the size of each. The traces may hold none of
them, one, or several. Name only what the traces support, and say what they cannot tell apart.

## How to work

{how}

## What to hand back, within 30 minutes

- `answer.json`: {{"terms": [{{"kind": "one of the kinds", "size": 0.0, "uncertainty": 0.0}}], "summary": "one paragraph"}}.
  An empty list of terms when the first model lacks nothing.
- `notes.md`: what you did, in its order, and what you could not tell apart.

Read the clock with a command. When 25 minutes have passed, write the answer that you have.
"""
HOW = {"A": "The package `twinfit` is in this folder. Its protocol is `twinfit/agent/AGENT.md`: read it and follow it. The researcher is not here: you open each step yourself, and you say in your notes why.",
       "B": "Use numpy and scipy alone, and no other package: fit the eight traces jointly, with the noise law as their weights, look at what the fit leaves, and decide."}


def texts(python: str) -> dict:
    """The instructions of each arm of an agent."""
    return {arm: COMMON.format(python=python, how=how) for arm, how in HOW.items()}


def main(argv: list) -> int:
    held = ROOT / "runs" / "rehearsal"
    const, config = model.load()
    if len(argv) == 2 and argv[0] == "make":
        truth = make(Path(argv[1]), ROOT, texts(sys.executable))
        held.mkdir(parents=True, exist_ok=True)
        (held / "truth.json").write_text(json.dumps(truth, indent=1) + "\n")
        print(f"{len(truth)} cases, {len(AGENTS)} arms of an agent each, in {Path(argv[1]).resolve()}")
        return 0
    if argv == ["script"]:
        out = {}
        for name, detail, draw in cases(config):
            selected = selection(fits.synthetic(detail, const, config, draw), const, config)
            out[name] = {"selected": selected, "answer": answer_of(selected, config), "score": score(answer_of(selected, config), truth_of(detail, config))}
            print(name, [a["term"] for a in selected["added"]], "right" if out[name]["score"]["right"] else "not right")
        held.mkdir(parents=True, exist_ok=True)
        (held / "script.json").write_text(json.dumps(out, indent=1) + "\n")
        return 0
    if len(argv) == 2 and argv[0] == "score":
        truth = json.loads((held / "truth.json").read_text())
        out = {}
        for home in sorted(Path(argv[1]).iterdir()):
            if (home / "answer.json").exists() and home.name.split(", arm ")[0] in truth:
                out[home.name] = score(json.loads((home / "answer.json").read_text()), truth[home.name.split(", arm ")[0]]["terms"])
                print(home.name, out[home.name]["named"], "of", out[home.name]["of"], "named,", out[home.name]["by_class"], "by class, false:", out[home.name]["false"])
        (held / "scores.json").write_text(json.dumps(out, indent=1) + "\n")
        return 0
    if len(argv) == 2 and argv[0] == "keep":
        record = gathered(Path(argv[1]), held, config)
        (ROOT / RECORD).write_text(json.dumps(record, indent=1) + "\n")
        print(f"{RECORD}: " + ", ".join(f"{name} {''.join(sorted(entry['arms']))}" for name, entry in record["cases"].items()))
        return 0
    if argv == ["--check"]:
        wrong = check(ROOT)
        for line in wrong:
            print("REFUSED: " + line)
        print("the rehearsal: " + ("every case has its three arms, and every answer gives its score again" if not wrong else f"{len(wrong)} thing(s) wrong"))
        return 1 if wrong else 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
