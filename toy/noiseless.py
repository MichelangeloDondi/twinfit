"""The package on noiseless traces, the first small step of his message 156: "start with noiseless traces and validate the
package on noiseless traces. Once you'll show me the results and after I'll tell you to do so, we'll go ahead also with
the noisy traces, then with more terms of the model".

The traces are the 8 of the fifth detail of the toy, the line of the model with no noise drawn. The package needs a
standard deviation for each point, to weight its fits and to make its twin: each point is given the one that the toy's
law gives it, and nothing of that noise is in the traces. So the report reads the traces against a twin with the noise
of the experiment, and a fit on them has no draw to be lucky or unlucky with.

First, the end of a fit: the first model fitted to the same traces with the standard deviation of each point at
`SCALES` times the toy's law. A search ends when the step that it still asks of each parameter is under a ten-thousandth
of that parameter's standard deviation. On the first run, with a thousandth of the law, the fit of the first model did
not end, and the loop stopped there. From the first model, the natural width alone,
the loop of the package runs with the four terms of the toy offered in their order and a fifth that the toy lacks, as
toy/chosen.py offers them, and reads at each model the report of where to search and the six proposals.

What a noiseless run must show, written before it ran:

    the report       under each model that lacks terms, something beyond the twin, and first the place and the place
                     that changes along a setting that `EXPECTED` admits for one of the terms that it lacks. Under the
                     fifth, nothing
    the loop         the four terms kept, the fifth refused
    the output       each parameter of the best model at the value that made the traces, within its own bar

    python toy/noiseless.py run       the end of the fit at four scales, then one run of the loop, kept in runs/noiseless/
    python toy/noiseless.py make OUT [--noise DRAW | --recorded] [--iterations N] [--minutes M]
                                      the folder of one agent, outside the repository: the noiseless traces, what the
                                      experimenter knows, the package, and the brief. No name of a term of the toy
    python toy/noiseless.py score OUT the answer of the agent against the truth, by `SCORE`, written before it ran
    python toy/noiseless.py keep OUT TOKENS SECONDS
                                      keeps the answer, the notes and the score of the agent in toy/noiseless.json, with the
                                      tokens and the seconds that its call took
    python toy/noiseless.py report RUN OUT
                                      writes into OUT a page for each step of the agent's run RUN, with its plots, the
                                      report of the package against the twin, the bias analysis and the agent's words
    python toy/noiseless.py run-report RUN OUT.html [--brief READ_FIRST]
    python toy/noiseless.py run-slides RUN OUT [--recorded RUN2] [--visual]       the same run as slides that show each decision, his messages 239 and 243, with --visual its figures alone, his message 297
                                      the whole run as one page for a reader who was not there: each turn, what the
                                      agent did and how it used the package, its step, then its answer beside the truth
    python toy/noiseless.py loop-report HOME OUT.html
                                      the loop of the package on the toy as one page, model by model, from the folder
                                      HOME that toy/chosen.py keeps, with the best model against the truth
    python toy/noiseless.py keep-blind OUT NAME
                                      keeps a blind run of the agent under its name, with what its run.json says, and
                                      refuses one whose notes do not open with `Context: none`
    python toy/noiseless.py table     writes toy/noiseless.json from the run kept
    python toy/noiseless.py --check   ends without error when toy/noiseless.json was made from the files as they are
"""
from __future__ import annotations

import hashlib
import itertools
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from twinfit import loop  # noqa: E402

from toy import chosen, fits, model  # noqa: E402

SCALE = 1.0          # the standard deviation given to each point, as a share of the toy's law: weights alone, since the traces hold no noise
SCALES = (1.0, 0.1, 0.01, 0.001)          # the scales at which the end of the fit of the first model is tried
ALONG_T = ("the core along T", "the flanks along T", "the wings along T")          # a width that grows with the temperature changes the core, the flanks or the wings along it
EXPECTED = {          # what the physics of each term says of where it shows and which place changes along which setting, written before the noiseless run, and again for the blind report before any count of it (the red team of his message 343, S1): a width of the toy leaves no slow trend and no periodic structure over the span, and no law of the noise, so the shape and the noise are answered by none of them
    "collisions": {"where": "the core", "changes along": ALONG_T, "physics": "a Lorentzian width proportional to the density of the vapour, which grows steeply with the temperature"},
    "gas_width": {"where": "the core", "changes along": None, "physics": "one Lorentzian width, the same in every trace"},
    "laser_noise": {"where": "the core", "changes along": None, "physics": "one Gaussian width, the same in every trace"},
    "transit": {"where": "the flanks", "changes along": ALONG_T, "physics": "a line exp(-|nu| / b) whose rate grows as the root of the temperature: its wings fall faster than the Lorentzian wings of the model, so a fit without it widens its other widths and leaves the difference within one width of the centre. toy/pointing.py wrote the wings, which his message 154 questioned"},
}


def rows_of(const: dict, config: dict, scale: float = SCALE) -> list:
    """The 8 traces of the fifth detail with no noise: each its line, and the standard deviation that weights it."""
    rows = fits.synthetic(5, const, config, 1)
    return [dict(r, y=r["signal"], sd=r["sd"] * scale) for r in rows]


NOISE_SAID = {False: ("with no noise", "These traces are synthetic and hold no noise at all: each is the line of a model."),
              True: ("with the noise of the experiment", "These traces are synthetic: each is the line of a model with the noise of the experiment, drawn by the law that `known.json` gives."),
              "recorded": ("recorded in the laboratory", "These traces were recorded in the laboratory: each is the line as it was measured, with the noise of the experiment, whose law `known.json` gives.")}          # the title and the sentence of the brief, his messages 211 and 228


def noisy_rows(const: dict, config: dict, draw: int) -> list:
    """The 8 traces of the fifth detail of draw `draw`, with the toy's law of the noise, his messages 156 and 211: the
    traces of the agentic version at the toy's noise."""
    return fits.synthetic(5, const, config, int(draw))


def ends(root: Path = ROOT) -> list:
    """The fit of the first model on the noiseless traces, the standard deviation of each point at each of `SCALES`
    times the toy's law: whether its search ended, in how many steps, its chi square times the square of the scale,
    which is the same at one minimum, and whether a search started at the point where the fit at the toy's law ended
    ends too."""
    from twinfit.fit import fit
    const, config = model.load(Path(root) / "toy", root)
    out, start = [], None
    for scale in SCALES:
        home = Path(root) / "runs" / "noiseless" / f"scale_{scale:g}"
        home.mkdir(parents=True, exist_ok=True)
        chosen.folder(home, chosen.described(rows_of(const, config, scale), home, "the synthetic traces of the fifth detail, with no noise", const))
        _, problem, _, _ = loop.initial(home, with_module=True)
        got = fit(problem)
        if start is None:
            start = (got.shared.copy(), got.own.copy())
        again = fit(problem, start[0], start[1])
        out.append({"scale": scale, "ended": bool(got.converged), "steps": int(got.iterations), "chi2_times_scale_squared": float(got.chi2 * scale ** 2),
                    "from_the_point_at_the_law": {"ended": bool(again.converged), "steps": int(again.iterations), "chi2_times_scale_squared": float(again.chi2 * scale ** 2)}})
    return out


def run(root: Path = ROOT) -> dict:
    """One run of the loop on the noiseless traces, with the report and the proposals at each model."""
    const, config = model.load(Path(root) / "toy", root)
    home = Path(root) / "runs" / "noiseless" / digest(root)[:12]
    home.mkdir(parents=True, exist_ok=True)
    chosen.folder(home, chosen.described(rows_of(const, config), home, "the synthetic traces of the fifth detail, with no noise", const))
    began = time.time()
    scales = ends(root)
    try:
        got = chosen.read(loop.run(home, n=len(chosen.ORDER) + 1, sets=chosen.SETS, report=True), config)
        got = dict(got, report=chosen.reported(home))
    except (Exception, SystemExit) as why:
        got = {"failed": f"{type(why).__name__}: {why}"[:300]}
    got["seconds"], got["ends"] = round(time.time() - began), scales
    (home / "run.json").write_text(json.dumps(got, indent=1) + "\n")
    return got


SPELLED = {"T": "temperature", "P": "power"}          # a covariate of EXPECTED, as toy/pointing.py writes it, by the name that description.toml gives it and that the report of the loop puts first


def admitted(term: str) -> dict:
    """The answers that the physics of a term of the toy admits to the questions of the blind report that it answers,
    the place and what changes along a setting, with the setting spelled as the report of the loop spells it. None: the
    term changes along no setting, and the report is right where it names none."""
    e = EXPECTED[term]
    spell = lambda a: a.rpartition(" along ")[0] + " along " + SPELLED.get(a.rpartition(" along ")[2], a.rpartition(" along ")[2])          # noqa: E731
    return {"where": (e["where"],), "changes along": None if e["changes along"] is None else tuple(spell(a) for a in e["changes along"])}


def against(got: dict) -> dict:
    """What the run shows against what a noiseless run must show: for each model, the terms of the toy that it lacks,
    and the first answers of the report beside those that `EXPECTED` writes for each of them: the terms whose place,
    kind and covariate the report gives first. Then the terms kept and refused, and how far each parameter of the best
    model lies from the value that made the traces, in its own bar. Model k lacks the terms after the k-th of the toy's
    order when the loop kept every term before it."""
    if "failed" in got:
        return {"failed": got["failed"]}
    out = {"models": {}, "kept": got["kept"], "refused": sorted(got["refused"]), "parameters": {}}
    for name, said in sorted(got["report"].items(), key=lambda item: int(item[0].split()[-1])):
        lacks = list(chosen.ORDER[int(name.split()[-1]):])          # model k, by the number in its name
        first = said.get("first") or {}
        agrees = {term: {q: first.get(q) is None if want is None else first.get(q) in want for q, want in admitted(term).items()} for term in lacks}          # "T" and "temperature" are one covariate: the run of 2026-10-15 put "temperature" first and was read as a miss
        out["models"][name] = {"lacks": lacks, "first": first, "agrees": agrees, "matches": [term for term in lacks if all(agrees[term].values())]}
    for name, cell in got["parameters"].items():
        out["parameters"][name] = {"value": cell["value"], "of_the_twin": cell["of_the_twin"], "bar": cell["bar"], "from_the_twin_in_bars": cell["from_the_twin_in_bars"]}
    return out


SCORE = {          # how a term of the toy is found in an answer, written before any agent ran: its shape, and what it changes with
    "collisions": {"shape": {"lorentzian"}, "changes_with": {"density", "temperature"}},
    "gas_width": {"shape": {"lorentzian"}, "changes_with": {"nothing"}},
    "laser_noise": {"shape": {"gaussian"}, "changes_with": {"nothing"}},
    "transit": {"shape": {"exponential"}, "changes_with": {"temperature"}},
}
ALIKE = ({"gaussian", "exponential"},)          # the shapes that 8 traces of this line hardly tell apart: power.map_of gives the two a cosine of 0.999
SHAPES = ("lorentzian", "gaussian", "exponential", "other")
CHANGES = ("nothing", "temperature", "density", "power", "other")
BRIEF = """# Eight traces of one spectral line, {noise_title}

You are asked to analyse eight curves of one spectral line. Work only in this folder. Read no other folder of this
machine and do not use the network. The Python to use is `{python}`: it has numpy, scipy and matplotlib.

## The measurement

One line of an atomic vapour, the two photon line of rubidium from 5S to 6S, recorded while a frequency is swept, at
eight conditions: four temperatures at one power, and five powers at one temperature. `traces/` holds one file a trace,
with the columns `frequency_MHz` and `signal_V`. {noise_said} `known.json` holds what the experimenter knows: the power and the temperature of each trace, the density of
the vapour at each temperature, the radius of the beam, the law of the noise that the experiment would have, to weight
the fits, and the natural width of the line, a Lorentzian width known from the lifetime of the upper state.

## The question

The first model is a line of the natural width alone: a Lorentzian of that width, with for each trace its own height,
centre, offset and slope of the baseline. Which physical terms does the line hold beyond it, and how large is each? No
list of terms is given: find them, from what the fits leave and from the physics of this vapour, and add each as a term
of physical origin. The traces may hold none, one or several. Name only what the traces support, and say what they
cannot tell apart.

## How to work

The package `twinfit` is in this folder. Its protocol is `twinfit/agent/AGENT.md`: read it and follow it, with at most
{iterations} iterations, each adding one term. The researcher is not here: you open each step yourself, and you say in your
notes why.

The point of the package is to go beyond the chi square: its report reads what each fit leaves through generalized
moments, their ratios, the moments of the derivative and of the running integral, windows along the frequency, bands
along the height, crossings and the chi of every sample, each against the twin of the fit. So:

- Every fit goes through the package: `python -m twinfit loop`, `python -m twinfit scan` and its functions. Write no
  fitting code of your own, no Voigt, no least squares and no chi square of your own. The noise of each sample comes
  from the law of `known.json` at the mean of the package's own fit, as `twinfit/noise.py` sets it.
- Each term is chosen by the answers of the report beyond the twin that it fills, the place, the kind, the setting and
  the family of statistics, and by the physics of this vapour. A chi square of your own is no reason to add a term or
  to stop. The loop's own rule of keeping a term stands, and two forms of one term are compared in copies of the folder
  through the loop, the report and the twin deciding between them.
- A shape that the package lacks is written as a Mechanism of your own in a term file: that is the model, not a fit. A red team may check a candidate with code of its own, as step 7 of the protocol
  allows, and the loop must bear out what it finds. Where this brief and `twinfit/agent/AGENT.md` differ, this brief
  holds.
- Stop when the twin admits the model and the report asks for no term, or puts one answer only just past its line, or
  after the last iteration.
- After each run of the loop, copy `DIR/loop/`, `DIR/model.py` and `DIR/terms/` into `steps/NN/`, NN the iteration
  from 00, before the next run, and write `steps/NN/step.md`. It holds, under the heading `## Proposals for the model`,
  the three proposals for the model that the report of the step gives in `loop/model_K.md`, under `## Proposals for
  the experiment` its three proposals for the experiment, and under `## Choice and why` which proposal you take and
  why. Then the answers of the report beyond the twin, family by family, the term that you choose, the answers that it fills and its physics, its size against your
  estimate, what the red team said, and what the twin says of the bias, the spread and the coverage of each parameter.
  The researcher red teams each step from this folder.

## What to hand back, within {minutes}

- `answer.json`: {{"terms": [{{"origin": "its physical origin, in words", "shape": "the shape of what it adds to the line, in your words", "changes_with":
  "one of {changes}", "width_at_130C_MHz": 0.0, "uncertainty_MHz": 0.0}}], "summary": "one paragraph"}}. The width
  is the full width at half maximum of the line that the term gives alone, at 130 C. An empty list of terms when the
  first model lacks nothing.
- `notes.md`: its first line `Context:` and the files of instructions or rules that your session carries, or `Context: none`, then what you did, in its order, and what you could not tell apart.

Read the clock with a command. When {last} minutes have passed, write the answer that you have.
"""


def truth(const: dict, config: dict) -> dict:
    """The terms of the fifth detail beyond the natural width, each with the full width at half maximum of the line
    that it gives alone at 130 C, from its own file in toy/terms."""
    import importlib
    values = model.truth(config)
    files = {d: name for name, d in const["details"].items()}          # the file of a term is named by the detail at which it enters
    return {name: float(importlib.import_module(f"toy.terms.{files[config['parameter'][name]['enters']]}").alone(130.0, values, const)) for name in chosen.ORDER}


HIDDEN = ("RUBIDIUM.md", "CONTINUING.md", "rubidium", "dictionary.py", "walkthrough.py", "survey.py", "two_scales.py", "against_the_twin.png", "turn_by_turn.png", "rubidium_twin")          # the page of his line in the package, and its figures, which name the five terms of the toy, the package's own list of terms with the two examples that read it (the review of candidate 95, moderate 2), and the two figures of the README, of a blind run on the line of the toy (lesson 35): left out of the folder of an agent
UNSAID = ((r"(?s)\s*With --dictionary,.*?(?:never in a blind run|nobody suspects)", ""), (r" \[--dictionary\]", ""), (r"(?m)^.*python -m twinfit walk .*\n", ""), (r"(?m)^.*examples/walkthrough\.py.*\n", ""),
          (r'dictionary = "--dictionary" in argv\[2:\][^\n]*', "dictionary = False"), (r'(?m)^[ ]*argv = argv\[:2\] \+ \[a for a in argv\[2:\] if a != "--dictionary"\]\n', ""),
          (r" And the package's own list, `twinfit/dictionary\.py`, read only on `--dictionary`, never in a blind run \(337 to 343\)\.", ""), (r"(?m)^ *dictionary\.py .*\n", ""),
          (r"(?m)^!\[[^\]]*\]\([^)]*(?:against_the_twin|turn_by_turn)\.png\)\n\n\*[^\n]*\n\n?", ""), (r"(?m)^.*python -m twinfit guide rubidium.*\n", ""), (r"(?s)## Continuing this work\n\n.*?\n\n", ""),
          (r"(?m)^ *python -m twinfit example rubidium DIR .*\n(?: {40,}.*\n){3}", ""),
          (r'(?s)    if argv\[:1\] == \["example"\] and argv\[1:2\] == \["rubidium"\].*?        return 0\n', ""))          # the twin of his traces, the owner's message 366, which the folder of a blind agent leaves out with the rest of `examples/rubidium`          # what the pages and the command line of the package say of that list and of the walk that reads it, taken out of the copy that an agent reads, where the flag is no more (his message 340)

ITERATIONS, MINUTES = 6, 120          # six iterations in two hours, the brief of the run of his message 310: his message 223 asked twelve iterations, and a run of the loop takes about 8 minutes
NUMBERS = ("no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve")


WORD_ON_HIS_TRACES = "runs/his_word_on_his_traces.txt"          # the number of the message of his that lets an agent run blind on his recorded traces: none yet, his message 229


def his_word(root: Path = ROOT) -> int | None:
    """The number of his message that lets an agent run blind on his recorded traces, from the file WORD_ON_HIS_TRACES,
    where OWNER_WORDS.md holds a message of that number, or None."""
    path = Path(root) / WORD_ON_HIS_TRACES
    said = path.read_text().strip() if path.exists() else ""
    words = Path(root) / "OWNER_WORDS.md"
    return int(said) if said.isdigit() and words.exists() and f"\n## {said}. " in words.read_text() else None


def make(out: Path, root: Path = ROOT, python: str = sys.executable, draw: int | None = None, iterations: int = ITERATIONS, minutes: int = MINUTES, recorded: bool = False) -> Path:
    """The folder of one agent: the noiseless traces, or with `draw` the traces of that draw with the toy's law of the
    noise, what the experimenter knows without the list of kinds of the rehearsal with hints, a copy of the package, and
    the brief, which says which, with its iterations and its minutes. The brief asks every fit through the package and
    every term from its report, and a folder of each step: two blind runs of his message 223 chose their terms by their
    own fits of line shapes and a chi square, outside the package (his message 224). With `recorded` the traces are his
    8 recorded traces of the peak of config.toml, read as `toy/real.py` reads them, his message 228: "At the end we would
    like to use the package to analyse the real traces!". Refuses a draw beside them, and a folder inside the repository. The copy of the package leaves out the page of
    his recorded traces and its figures: since 0.0.5 they name the five terms of the toy and the three that change with
    the power, the answer that the agent is asked for."""
    import shutil
    from toy import rehearsal, traces
    out, root = Path(out).resolve(), Path(root).resolve()
    if out == root or root in out.parents:
        raise SystemExit(f"REFUSED: {out} is inside the repository: the agent reads no file of it")
    if recorded and his_word(root) is None:          # his message 229: no blind run on his recorded traces before his word, which a file names (review 39, minor 19)
        raise SystemExit(f"REFUSED: no blind run on his recorded traces before his word, his message 229: {WORD_ON_HIS_TRACES} names none of his messages. "
                         "Write there the number of his message that allows such a run, once he sends one, and no other")          # the file holds a number, and the session alone writes it (review 40, minor 11)
    if recorded and draw is not None:
        raise SystemExit("REFUSED: the recorded traces are no draw: ask for one or the other")
    const, config = model.load(root / "toy", root)
    from toy import real
    (out / "traces").mkdir(parents=True, exist_ok=True)
    for r in (real.rows(const, config, root) if recorded else rows_of(const, config) if draw is None else noisy_rows(const, config, draw)):
        (out / "traces" / traces.name_of(r["P"], r["T"])).write_text("frequency_MHz,signal_V\n" + "".join(f"{x:.6f},{y:.9e}\n" for x, y in zip(r["nu"], r["y"])))
    known = {k: v for k, v in rehearsal.known(const, config).items() if k != "kinds"}
    known["beam_radius_um"], known["beam_radius_uncertainty_um"] = float(const["beam"]["waist_um"]), float(const["beam"]["uncertainty_um"])
    (out / "known.json").write_text(json.dumps(known, indent=1) + "\n")
    title, said = NOISE_SAID["recorded" if recorded else draw is not None]
    said_minutes = f"{minutes // 60} hours" if minutes >= 120 and minutes % 60 == 0 else f"{minutes} minutes"
    (out / "READ FIRST.md").write_text(BRIEF.format(python=python, changes=", ".join(CHANGES), noise_title=title, noise_said=said,
                                                    iterations=NUMBERS[iterations] if iterations < len(NUMBERS) else iterations, minutes=said_minutes, last=minutes - max(5, minutes // 8)))
    if (out / "twinfit").exists():
        shutil.rmtree(out / "twinfit")
    shutil.copytree(root / "twinfit", out / "twinfit", ignore=shutil.ignore_patterns("__pycache__", *HIDDEN))
    for blind in (out / "twinfit" / "TUTORIAL.md", out / "twinfit" / "tutorial.py"):          # the folder of a blind agent names no flag of the dictionary, row M471, where the tutorial for people and its writer name it since the repairs of his message 406
        if blind.exists():
            blind.write_text(blind.read_text().replace("`--dictionary`", "the flag of the dictionary").replace("--dictionary", "the flag of the dictionary"))
    import re
    for page in [p for p in (out / "twinfit").rglob("*") if p.suffix in (".md", ".py", ".toml")]:
        text = page.read_text()
        for pattern, instead in UNSAID:
            text = re.sub(pattern, instead, text)
        page.write_text(blind_history(text) if page.name == "version.py" and page.parent.name == "twinfit" else text)
    return out


def blind_history(text: str) -> str:
    """`twinfit/version.py` for the folder of a blind agent: each version keeps its digest, and its changes are one line
    that says they are left out. Eight of their lines named what the package's own list of terms held, the ladders of
    widths, the second component, the slope of the level, in a folder that holds no list (reviews 98 to 101)."""
    import ast
    tree = ast.parse(text)
    node = next(n for n in tree.body if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "VERSIONS" for t in n.targets))
    for value in node.value.values:
        for k, v in zip(value.keys, value.values):
            if isinstance(k, ast.Constant) and k.value == "changes":
                value.values[value.keys.index(k)] = ast.Tuple(elts=[ast.Constant("left out of the folder of a blind agent")], ctx=ast.Load())
    lines = text.splitlines(keepends=True)
    return "".join(lines[:node.lineno - 1]) + ast.unparse(node) + "\n" + "".join(lines[node.end_lineno:])


def shape_kind(said) -> str:
    """The kind of a shape that an agent wrote in its own words, as the scorer reads it: one of SHAPES, the one that it
    names first, since a description names its kernel before what its flanks look like. The blind run of 0.1.0rc16 wrote
    "a cusped, near Gaussian kernel with Lorentzian looking flanks", which a fixed order of the kinds read as Lorentzian."""
    text = str(said or "").lower()
    at = [(text.find(key), kind) for key, kind in (("lorentz", "lorentzian"), ("gauss", "gaussian"), ("expon", "exponential")) if key in text]
    return min(at)[1] if at else "other"


def score(answer: dict, true: dict) -> dict:
    """An answer against the truth by `SCORE`: for each term of the toy, the term of the answer of its shape, or of a
    shape of its class, that changes with what it changes with, and how far its width lies from the true one in its own
    uncertainty. And the terms of the answer that match no term of the toy."""
    terms = [dict(a, shape=shape_kind(a.get("shape"))) for a in answer.get("terms", []) if isinstance(a, dict)]          # the shape in the agent's own words, read for its kind: the brief offers no list of shapes (his message 340)
    of_class = lambda shape: next((c for c in ALIKE if shape in c), {shape})          # noqa: E731

    def fit(name: str, k: int) -> int:          # 2 a term of its shape, 1 of a shape of its class, 0 none, each changing with what the toy's term changes with
        a, rule = terms[k], SCORE[name]
        if a.get("changes_with") not in rule["changes_with"]:
            return 0
        return 2 if a.get("shape") in rule["shape"] else 1 if of_class(a.get("shape")) & rule["shape"] else 0

    def off(name: str, k: int | None) -> float:
        a = terms[k] if k is not None else {}
        u, w = a.get("uncertainty_MHz"), a.get("width_at_130C_MHz")
        return abs((w - true[name]) / u) if isinstance(u, (int, float)) and u and isinstance(w, (int, float)) else 0.0

    names, best = list(SCORE), None          # the terms of the toy matched to those of the answer all at once, the most named first, then the most by class, then the nearest: matched one after the other, the first term of the toy took the transit for the collisions (0.1.0rc16)
    for pick in itertools.permutations(list(range(len(terms))) + [None] * len(names), len(names)):
        if any(k is not None and fit(n, k) == 0 for n, k in zip(names, pick)):
            continue
        key = (sum(k is not None and fit(n, k) == 2 for n, k in zip(names, pick)), sum(k is not None for k in pick), -sum(off(n, k) for n, k in zip(names, pick)))
        if best is None or key > best[0]:
            best = (key, pick)
    chosen = dict(zip(names, best[1])) if best else dict.fromkeys(names)
    used = {k for k in chosen.values() if k is not None}
    found = {}
    for name in names:
        k = chosen[name]
        exact = [k] if k is not None and fit(name, k) == 2 else []
        a = terms[k] if k is not None else {}
        u = a.get("uncertainty_MHz")
        found[name] = {"named": bool(exact), "by_class": k is not None, "true_width_MHz": true[name], "width_MHz": a.get("width_at_130C_MHz"),
                       "off_in_uncertainties": float((a["width_at_130C_MHz"] - true[name]) / u) if k is not None and u else None}
    false = [terms[k].get("origin") for k in range(len(terms)) if k not in used]
    return {"terms": found, "named": sum(f["named"] for f in found.values()), "by_class": sum(f["by_class"] for f in found.values()), "of": len(SCORE),
            "within_two": sum(1 for f in found.values() if f["off_in_uncertainties"] is not None and abs(f["off_in_uncertainties"]) <= 2.0), "false": false}


MAKES = ("SCALE", "SCALES", "rows_of", "ends", "run")          # the parts of this program that make the run
RUNS = ("toy/noiseless.py", "toy/chosen.py")          # the programs whose imports say which files of the package the run runs
DEPENDS = ("toy/chosen.py", "toy/model.py", "toy/fits.py", "toy/traces.py", "toy/const.toml", "toy/config.toml")


def digest(root: Path = ROOT) -> str:
    """What the run depends on: the parts of this program that make it, the loop on the toy, the model, its values, the
    terms and the files of the package. Another part of this program, a reading or the folder of an agent, moves none."""
    from toy import power
    root = Path(root)
    h = hashlib.sha256(power.sources((root / "toy" / "noiseless.py").read_text(), MAKES))
    from toy import draws
    draws.hashed(h, root, list(DEPENDS) + sorted(str(p.relative_to(root)) for p in (root / "toy" / "terms").glob("*.py")) + draws.package_of(root, RUNS))          # the files of the package that the loop imports, as code (rule 16, his message 175)
    return h.hexdigest()


AGENTS = ("agent", "agent_at_the_noise")          # the keys of toy/noiseless.json where `keep` puts the answer of an agent
BLIND = "blind_runs"          # the key of toy/noiseless.json where `keep-blind` puts each blind run of the agent, by its name (his messages 219 and 223)


def table(root: Path = ROOT) -> Path:
    path = Path(root) / "runs" / "noiseless" / digest(root)[:12] / "run.json"
    if not path.exists():
        raise SystemExit("REFUSED: no run kept from the files as they are: `python toy/noiseless.py run` first")
    got = json.loads(path.read_text())
    out = Path(root) / "toy" / "noiseless.json"
    held = json.loads(out.read_text()) if out.exists() else {}
    table = {"digest": digest(root), "scale": SCALE, "seconds": got.get("seconds"), "ends": got.get("ends"), "run": got, "against": against(got)}
    table.update({key: held[key] for key in AGENTS + (BLIND,) if key in held})          # the answers of the agents, kept by `keep` and `keep-blind`, outlive a new table: review 33 found them dropped
    out.write_text(json.dumps(table, indent=1) + "\n")
    return out


def blind(notes: str) -> bool:
    """Whether an agent's run is blind: its notes open with `Context: none`, the line of the brief where it names the
    instructions its session carries. A subagent of his cloud chat carries the rules of this repository, which name the
    terms of the toy, and the agent at the toy's noise said so (review 32): a run that names any, or says nothing, is not.
    A clause after `none` is read only where it begins with no or nothing, so that it says none again."""
    import re
    head = [line.strip() for line in notes.splitlines() if line.strip()][:5]          # under a heading such as "# Notes", as both agents wrote (review 33)
    said = next((line for line in head if line.lower().startswith("context:")), "").lower().rstrip(".")
    return said == "context: none" or bool(re.match(r"context: none\s*[\u2014\u2013:,(.-]+\s*(no|nothing)\b", said))          # a clause that says none again, as the notes of a blind run of Opus wrote, and not one that names what was carried (review 38, minor 5)


def keep(out: Path, tokens: int, seconds: float, root: Path = ROOT, draw: int | None = None) -> dict:
    """The answer of the agent kept in toy/noiseless.json beside the run: its answer and its notes as it handed them
    back, its score by `SCORE`, and the tokens and the seconds that it took. Refused before the table is written."""
    path = Path(root) / "toy" / "noiseless.json"
    if not path.exists():
        raise SystemExit("REFUSED: toy/noiseless.json is not written: `python toy/noiseless.py table` first")
    held = json.loads(path.read_text())
    const, config = model.load(Path(root) / "toy", root)
    answer = json.loads((Path(out) / "answer.json").read_text())
    notes = (Path(out) / "notes.md").read_text() if (Path(out) / "notes.md").exists() else ""
    key = "agent" if draw is None else "agent_at_the_noise"          # the agent at the toy's noise beside the noiseless one, his message 211
    held[key] = {"answer": answer, "notes": notes, "truth": truth(const, config), "score": score(answer, truth(const, config)), "tokens": int(tokens), "seconds": round(float(seconds)),
                 "blind": blind(notes)}          # review 32: the second agent said that its session carried the rules of this repository, which name the terms
    if draw is not None:
        held[key]["draw"] = int(draw)
    path.write_text(json.dumps(held, indent=1) + "\n")
    return held[key]


RAN = ("nonce", "draw", "model", "effort", "tokens", "seconds", "iterations", "stopped", "own_fits")          # what the run.json of a blind run says of it: since review 38, moderate 2, how it ended and the files of its folder that fit outside the package


def keep_blind(out: Path, name: str, root: Path = ROOT) -> dict:
    """A blind run of the agent kept in toy/noiseless.json under its name: its answer and its notes as it handed them
    back, its score by `SCORE`, and what its run.json says of it, of `RAN`. Refused where its notes do not open with
    `Context: none`: a run whose session carried a rules file, or that says nothing of it, is not blind (his message 219:
    "Why you didn't run the blind run? What are you waiting for? It's the core point of the whole project"). Refused
    before the table is written, and for a run.json without its draw."""
    path = Path(root) / "toy" / "noiseless.json"
    if not path.exists():
        raise SystemExit("REFUSED: toy/noiseless.json is not written: `python toy/noiseless.py table` first")
    notes = (Path(out) / "notes.md").read_text() if (Path(out) / "notes.md").exists() else ""
    if not blind(notes):
        raise SystemExit(f"REFUSED: the notes in {out} do not open with `Context: none`: the run is not kept as blind")
    ran = json.loads((Path(out) / "run.json").read_text())
    if "draw" not in ran:
        raise SystemExit(f"REFUSED: the run.json in {out} names no draw: the traces of the run are not known")
    held = json.loads(path.read_text())
    const, config = model.load(Path(root) / "toy", root)
    answer = json.loads((Path(out) / "answer.json").read_text())
    true = truth(const, config)
    held.setdefault(BLIND, {})[name] = dict({"answer": answer, "notes": notes, "truth": true, "score": score(answer, true), "blind": True}, **{k: ran[k] for k in RAN if k in ran})
    path.write_text(json.dumps(held, indent=1) + "\n")
    return held[BLIND][name]


def steps_of(run: Path) -> list:
    """The steps of an agent's run, in order, each as (its name, its loop folder, its folder): the folders `steps/NN/`
    that the brief asks for, from the folder `steps` nearest to the top of the run that holds them, or where the run kept
    none, every loop folder that it left, a folder holding `model_0.json`."""
    run = Path(run)
    homes = sorted((p for p in run.rglob("steps") if p.is_dir() and any(c.is_dir() and c.name.isdigit() for c in p.iterdir())), key=lambda p: (len(p.parts), str(p)))          # where the run kept them: one run kept its steps a folder deeper than the brief asked (2026-10-18)
    kept = sorted(p for p in homes[0].iterdir() if p.is_dir() and p.name.isdigit()) if homes else []
    if kept:
        return [(p.name, loop_in(p), p) for p in kept]
    loops = sorted({p.parent for p in run.rglob("model_0.json")}, key=lambda p: (p.stat().st_mtime, str(p)))
    return [(str(p.relative_to(run)).replace("/", " "), p, p.parent) for p in loops]


def loop_in(step: Path) -> Path:
    """The loop folder of a step: its own `loop`, or where the step kept copies of the project each with its loop, as the
    run on the hidden project of message 235 did, the first of them by name, or the step itself."""
    step = Path(step)
    if (step / "loop" / "model_0.json").exists():
        return step / "loop"
    nested = sorted(p.parent for p in step.rglob("model_0.json"))
    return nested[0] if nested else step


def loops_of(loop: Path, folder: Path) -> list:
    """Every loop folder of a step: the one of `loop_in` alone, or where the step kept copies of the project each with
    its loop, all of them by name."""
    loop, folder = Path(loop), Path(folder)
    if loop in (folder / "loop", folder):
        return [loop]
    return sorted(p.parent for p in folder.rglob("model_0.json"))


def last_model(loop: Path) -> int | None:
    """The highest K of the `model_K.json` of a loop folder, the model that the step ended on, or None."""
    ks = [int(p.stem.split("_")[1]) for p in Path(loop).glob("model_*.json") if p.stem.split("_")[1].isdigit()]
    return max(ks) if ks else None


def _png(fig) -> str:
    import base64
    import io
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=90, bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return "<img src='data:image/png;base64," + base64.b64encode(buf.getvalue()).decode() + "' style='max-width:100%'>"


def _escape(text: str) -> str:
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def prose(text: str) -> str:
    """The notes of an agent, written in Markdown with hard breaks, as HTML for a reader: headings, lists, bold and code,
    and the lines of one paragraph or one item joined again. No package renders Markdown here, and none is installed
    without his word."""
    import re

    def inline(s):
        s = _escape(s)
        s = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", s)
        return re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    blocks, current, kind = [], [], None

    def close():
        nonlocal current, kind
        if current:
            body = " ".join(x.strip() for x in current)
            blocks.append(("li" if kind == "item" else "p", body))
        current, kind = [], None
    for line in str(text).splitlines():
        heading = re.match(r"^(#{1,4})\s+(.*)", line)
        if line.strip().startswith("|"):          # a row of a table of the loop's pages, as best.md holds them
            close()
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if not all(re.fullmatch(r":?-{2,}:?", c) for c in cells if c):
                blocks.append(("tr", cells))
            continue
        if heading:
            close()
            blocks.append((f"h{min(len(heading.group(1)) + 3, 6)}", heading.group(2)))
        elif not line.strip():
            close()
        elif re.match(r"^\s*([-*]|\d+\.)\s+", line) and not line.startswith("    "):
            close()
            current, kind = [re.sub(r"^\s*([-*]|\d+\.)\s+", "", line)], "item"
        else:
            current.append(line)
            kind = kind or "para"
    close()
    out, open_tag = [], None
    for tag, body in blocks:
        wanted = {"li": "ul", "tr": "table"}.get(tag)
        if open_tag and wanted != open_tag:
            out.append(f"</{open_tag}>")
            open_tag = None
        if wanted and not open_tag:
            out.append("<ul>" if wanted == "ul" else "<table border=1 cellpadding=3>")
            open_tag = wanted
        out.append("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in body) + "</tr>" if tag == "tr" else f"<{tag}>{inline(body)}</{tag}>")
    return "\n".join(out + ([f"</{open_tag}>"] if open_tag else []))


def in_scatters(cell: dict | None) -> str:
    """A bias of a known parameter on the twin, as `proposals` keeps it, in its own unit beside its scatter, said in
    scatters: the page printed the unit under a heading of scatters (review 38, moderate 3)."""
    cell = cell or {}
    bias, scatter = cell.get("bias"), cell.get("scatter")
    if bias is None or not scatter:
        return "not read"
    return f"{bias / scatter:+.2g}"


def shares(loop: Path, k: int) -> bool:
    """Whether model k of a loop folder shares a parameter, read from its fit in the notebook: a model that shares none
    has a twin that gives no value back, and its admission says nothing (review 38, minor 6)."""
    fits = sorted(Path(loop).glob(f"notebook/*_model_{k}_fit.json"))
    return bool(fits) and bool(json.loads(fits[-1].read_text()).get("shared"))


def step_page(name: str, loop: Path, folder: Path, k: int | None = None) -> str:
    """One step of an agent's run as a page for him to red team (his message 224): the eight fits of the model that it
    ended on and their residuals in units of the noise, every answer of the report of the package against the line of
    its twin, the bias, the spread and the coverage of each parameter on the twin, the proposals with the bias of each
    known parameter without and with the term, the model and the terms as written, and the agent's own words."""
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    loop, folder = Path(loop), Path(folder)
    one = k is not None          # one model of a loop folder, as the report of the loop on the toy reads each, his message 230
    k = last_model(loop) if k is None else k
    if k is None:
        return f"<h2>Step {_escape(name)}</h2><p>No model of the loop is kept in this step.</p>"
    table = json.loads((loop / f"model_{k}.json").read_text())
    parts = [f"<h2>Step {_escape(name)}: model {k}</h2>",
             f"<p>Chi square {table['chi_square']:.1f} on {table['degrees_of_freedom']} degrees of freedom, {table['chi_square'] / max(table['degrees_of_freedom'], 1):.3f} a degree. "
             + (f"The twin admits the fit: {'yes' if table.get('admitted_on_the_twin') else 'no'}. {_escape(table.get('measured_on_the_twin', ''))}</p>" if shares(loop, k)
                else "The model shares no parameter, so its twin has no value to give back and admits nothing: the sentence of the loop is left out.</p>")]
    words = folder / "step.md"
    if words.exists():
        parts.append("<h3>The agent's words for this step</h3><div style='border-left:3px solid #56B4E9;padding-left:12px'>" + prose(words.read_text()) + "</div>")
    rows = list(csv.DictReader(open(loop / f"residuals_{k}.csv"))) if (loop / f"residuals_{k}.csv").exists() else []
    curves = sorted({r["curve"] for r in rows})
    if curves:
        fig, axes = plt.subplots(2, len(curves), figsize=(2.4 * len(curves), 5.2), sharex=True, squeeze=False)
        for j, c in enumerate(curves):
            got = [r for r in rows if r["curve"] == c]
            x = np.array([float(r["abscissa"]) for r in got])
            axes[0][j].plot(x, [float(r["value"]) for r in got], ".", ms=1.5, color="#999999")
            axes[0][j].plot(x, [float(r["fitted"]) for r in got], "-", lw=1.2, color="#0072B2")
            axes[0][j].set_title(Path(c).stem.replace("trace_", ""), fontsize=8)
            axes[1][j].plot(x, [float(r["residual"]) / float(r["noise"]) for r in got], ".", ms=1.5, color="#D55E00")
            axes[1][j].axhline(0, color="black", lw=0.6)
            axes[1][j].set_xlabel("abscissa")
        axes[0][0].set_ylabel("value, and the fit")
        axes[1][0].set_ylabel("residual over the noise")
        parts.append("<h3>The fits of the model and what they leave</h3>" + _png(fig))
    where = table.get("where") or {}
    answers = where.get("answers") or []
    line = where.get("largest_excess_on_the_twin")
    if answers:
        questions = list(dict.fromkeys(a["question"] for a in answers))
        fig, axes = plt.subplots(len(questions), 1, figsize=(9, 1.2 + 0.32 * len(answers)), squeeze=False,
                                 gridspec_kw={"height_ratios": [sum(a["question"] == q for a in answers) for q in questions]})
        for ax, q in zip(axes[:, 0], questions):
            mine = sorted((a for a in answers if a["question"] == q), key=lambda a: a["excess"])
            ax.barh([a["answer"][:70] for a in mine], [max(a["excess"], 1e-3) for a in mine], color=["#0072B2" if a["beyond_the_twin"] else "#999999" for a in mine])
            ax.set_xscale("symlog", linthresh=1.0)
            if line is not None:
                ax.axvline(line, color="black", ls="--", lw=1)
            ax.set_title(q, fontsize=9, loc="left")
            ax.tick_params(labelsize=7)
        axes[-1, 0].set_xlabel("excess of each answer over the twin, with the line of the twin dashed")
        parts.append(f"<h3>The report of where to search, on {where.get('sets')} sets of the twin</h3><p>Blue: beyond the twin. First answer of each question: "
                     + _escape(json.dumps(where.get("first_of_each_question"))) + "</p>" + _png(fig))
    twin = sorted(loop.glob(f"notebook/*_model_{k}_twin.json"))
    admission = json.loads(twin[-1].read_text()).get("admission", {}) if twin else {}
    params = admission.get("parameters") or {}
    if params:
        names = list(params)
        bands = admission.get("bands", {})
        fig, axes = plt.subplots(1, 3, figsize=(10, 0.8 + 0.35 * len(names)), squeeze=False, sharey=True)          # the names once, at the left: they printed over the panels beside them (review 38, minor 6)
        for ax, (key, title) in zip(axes[0], (("bias_over_spread", "bias over spread"), ("spread_over_bound", "spread over the bar of the fit"), ("coverage", "coverage of the interval"))):
            ax.barh(names, [params[n].get(key) or 0.0 for n in names], color="#56B4E9")
            band = bands.get({"bias_over_spread": "bias_over_spread_under"}.get(key, key))
            if isinstance(band, list):
                ax.axvspan(band[0], band[1], color="#009E73", alpha=0.15)
            elif band is not None:
                ax.axvspan(-band, band, color="#009E73", alpha=0.15)
            ax.set_title(title, fontsize=9)
            ax.tick_params(labelsize=7)
        parts.append(f"<h3>The bias analysis on the twin of model {k}</h3><p>" + _escape(admission.get("rule", "")) + "</p>" + _png(fig))
    fits = sorted(loop.glob(f"notebook/*_model_{k}_fit.json"))
    if fits and json.loads(fits[-1].read_text()).get("shared"):
        shared = json.loads(fits[-1].read_text())
        parts.append("<h3>The shared parameters</h3><table border=1 cellpadding=3><tr><th>parameter</th><th>value</th><th>bar of the fit</th></tr>"
                     + "".join(f"<tr><td>{_escape(n)}</td><td>{v:.6g}</td><td>{shared.get('standard_deviations', {}).get(n, float('nan')):.2g}</td></tr>" for n, v in shared.get("shared", {}).items()) + "</table>")
    if (loop / "best.md").exists():          # the five parts of the uncertainty of each value and its verdict, which the loop writes and the page read nowhere (review 38, minor 6)
        parts.append("<h3>The best model of the step, as the loop writes it</h3><div style='font-size:90%'>" + prose((loop / "best.md").read_text()) + "</div>")
    proposals = (table.get("proposals") or {}).get("for_the_model") or []
    if proposals:
        cells = []
        for p in proposals:
            known = p.get("known_parameters") or {}
            moved = ", ".join(f"{_escape(n)}: {in_scatters(v.get('without_the_term'))} to {in_scatters(v.get('with_the_term'))}" for n, v in known.items())
            cells.append(f"<tr><td>{_escape(p.get('term'))}</td><td>{_escape(p.get('stands_for'))}</td><td>{p.get('excess_in_the_report', float('nan')):.3g}</td>"
                         f"<td>{p.get('size', float('nan')):.4g} +- {p.get('bar_of_the_size', float('nan')):.2g}</td><td>{p.get('gain_in_chi_square', float('nan')):.4g}</td><td>{moved}</td></tr>")
        parts.append("<h3>The proposals of the package for the model</h3><table border=1 cellpadding=3><tr><th>term</th><th>for the answer</th><th>excess</th><th>size</th><th>gain in chi square</th>"
                     "<th>bias of each known parameter on the twin, in its scatters, without and with the term</th></tr>" + "".join(cells) + "</table>")
    code = [p for p in [folder / "model.py"] + sorted((folder / "terms").glob("*.py")) if p.exists()] if (folder / "model.py").exists() else [p for p in [loop.parent / "model.py"] + sorted((loop.parent / "terms").glob("*.py")) if p.exists()]
    if one:          # the file that this model adds, and no other
        code = [p for p in code if (p.name == "model.py" if k == 0 else p.name == f"{k:02d}.py")]
    for p in code:
        parts.append(f"<h3>{_escape(p.parent.name + '/' + p.name if p.parent.name == 'terms' else p.name)}</h3><pre style='white-space:pre-wrap'>{_escape(p.read_text())}</pre>")
    return "\n".join(parts)


def steps_report(run: Path, out: Path) -> list:
    """The report of each step of an agent's run, his message 224: "make sure that I can review the log, creating a
    report (with all the plots, the terms of the model, the theoretical reasoning, the bias analysis, etc) after each
    step, so that I can red team it". One page a step by `step_page`, and an index with the notes of the run. Refuses a
    run with no step. Returns the pages written, the index first."""
    out = Path(out)
    steps = steps_of(run)
    if not steps:
        raise SystemExit(f"REFUSED: {run} holds no step: no folder steps/NN and no loop folder")
    out.mkdir(parents=True, exist_ok=True)
    pages = []
    for name, loop, folder in steps:
        page = out / f"step {name}.html"
        page.write_text(f"<!doctype html><meta charset='utf-8'><title>Step {_escape(name)}</title><body style='font-family:Arial;max-width:1100px;margin:auto'>" + step_page(name, loop, folder) + "</body>")
        pages.append(page)
    notes = Path(run) / "notes.md"
    index = out / "index.html"
    index.write_text("<!doctype html><meta charset='utf-8'><title>The steps of the run</title><body style='font-family:Arial;max-width:1100px;margin:auto'><h1>The steps of the run " + _escape(Path(run).name) + "</h1><ol>"
                     + "".join(f"<li><a href='{_escape(p.name)}'>{_escape(p.stem)}</a></li>" for p in pages) + "</ol>"
                     + ("<h2>The notes of the agent</h2><pre style='white-space:pre-wrap'>" + _escape(notes.read_text()) + "</pre>" if notes.exists() else "") + "</body>")
    return [index] + pages


STAMP_UTC = r"\w{3} \w{3} +\d+ (\d\d:\d\d:\d\d) UTC (\d{4})"          # the clock as `date` prints it in the container of his cloud chat


def executable(command: str) -> tuple:
    """A command of the shell split into what runs and the files that it writes through here documents: the lines of
    each here document are taken out of the first and kept, by the line that opens them, in the second."""
    import re
    run, bodies, closing, body = [], [], None, []
    for line in str(command).splitlines():
        if closing is not None:
            if line.strip() == closing:
                bodies.append("\n".join(body))
                closing, body = None, []
            else:
                body.append(line)
            continue
        run.append(line)
        opened = re.search(r"<<-?\s*['\"]?(\w+)['\"]?", line)
        closing = opened.group(1) if opened else None
    if closing is not None:
        bodies.append("\n".join(body))
    return "\n".join(run), bodies


OWN_FIT = r"least_squares|curve_fit|minimize\(|leastsq|lstsq|polyfit|\.fit\(\s*model"          # a fit made outside the package, which the brief refuses: the scan of his cloud chat found lstsq and polyfit in the red team of a run that this read as none (2026-10-18)


def kind_of(name: str, given: dict) -> tuple:
    """What one call of a tool did, as (its kind, its words): a command of the package from the shell or from Python,
    a fit of its own, a term or the model written, a red team, a reading, a writing, or another command."""
    import re
    raw = str(given.get("command") or "")
    run, bodies = executable(raw)
    first = " ".join(run.split())[:220]
    path = str(given.get("file_path") or "")
    written = str(given.get("content") or given.get("new_string") or "")
    code = "\n".join(bodies) + "\n" + written
    package = list(dict.fromkeys(re.findall(r"-m twinfit (\w+)", run)))
    imported = sorted(set(re.findall(r"from twinfit(?:\.(\w+))? import", code)) | set(re.findall(r"twinfit\.(\w+)\.", code)) - {""})
    if name == "Agent":
        return "red team", f"a subagent: {given.get('description', '')}"
    if re.search(OWN_FIT, code) and not re.search(r"terms/\d+\.py|model\.py", path):
        return "own fit", f"a fit outside the package: {first or path}"
    if package:
        return "package", "twinfit " + ", twinfit ".join(package) + f": {first}"
    if re.search(r"terms/\d+\.py$", path) or re.search(r"(>|cp\s.*)\s*\S*terms/\d+\.py", run):
        return "term", f"writes the term {path or first}"
    if path.endswith("model.py") or re.search(r">\s*\S*model\.py", run):
        return "model", f"writes the model {path or first}"
    if imported and re.search(r"python3? -|python3? \S+\.py", run + " " + path):
        return "package", f"the package from Python, twinfit {', '.join(m or 'itself' for m in imported)}: {first or path}"
    if name in ("Write", "Edit") or bodies:
        return "writes", f"writes {path or first}"
    if name in ("Read", "Grep", "Glob") or re.match(r"(cd \S+;\s*)?(cat|sed|head|tail|grep|ls|wc|find)\b", first):
        return "reads", f"reads {path or given.get('pattern') or first}"
    return "other", first or json.dumps(given)[:200]


def actions_of(stream: Path) -> list:
    """The actions of a run in their order, read from its stream: each call of a tool with its kind and its words by
    `kind_of`, its time in UTC, from the stamp of its message or else the clock that the run last showed, and the folders
    steps/NN that its command names. The words that the agent wrote between its calls are kept as actions of the kind
    "says". What a subagent of the run says and does, a red team, carries "subagent" before its kind: on 2026-10-18 the
    fits of the red team of a run read as fits of the agent."""
    import gzip
    import re
    opener = gzip.open if str(stream).endswith(".gz") else open
    out, clock = [], None
    with opener(stream, "rt", errors="replace") as fh:
        for line in fh:
            try:
                e = json.loads(line)
            except json.JSONDecodeError:
                continue
            message = e.get("message") if isinstance(e.get("message"), dict) else {}
            stamp = str(e.get("timestamp") or "")[11:19] or None
            by = "subagent " if e.get("parent_tool_use_id") else ""
            for b in message.get("content", []) if isinstance(message.get("content"), list) else []:
                if b.get("type") == "tool_result":
                    c = b.get("content")
                    text = c if isinstance(c, str) else " ".join(x.get("text", "") for x in c if isinstance(x, dict)) if isinstance(c, list) else ""
                    found = re.findall(STAMP_UTC, text)
                    clock = found[-1][0] if found else clock
                elif b.get("type") == "text" and e.get("type") == "assistant" and b.get("text", "").strip():
                    out.append({"kind": by + "says", "what": " ".join(b["text"].split())[:600], "at": stamp or clock, "steps": []})
                elif b.get("type") == "tool_use":
                    given = b.get("input") or {}
                    kind, what = kind_of(b.get("name"), given)
                    run, _ = executable(str(given.get("command") or ""))
                    out.append({"kind": by + kind, "what": what, "at": stamp or clock, "steps": [] if by else sorted({int(x) for x in re.findall(r"steps/(\d\d)\b", run + " " + str(given.get("file_path") or ""))})})
    return out


def turns_of(actions: list) -> dict:
    """The actions of each step NN: those after the last that names the folder of the step before and up to the last
    that names its own, where it is copied and its words are written. What comes after the last step is "after"."""
    last = {}
    for k, a in enumerate(actions):
        for n in a["steps"]:
            last[n] = k
    out, begin = {}, 0
    for n in sorted(last):
        out[f"{n:02d}"] = actions[begin:last[n] + 1]
        begin = last[n] + 1
    out["after"] = actions[begin:]
    return out


KEY_KINDS = ("says", "package", "own fit", "term", "model", "red team", "subagent says", "subagent own fit")          # the actions that tell a turn, in the body of the report: the rest stand in its appendix


GLOSSARY = (("The twin", "synthetic records drawn from the fitted model with the noise of the experiment, each fitted again as the records were. What the records show beyond what the twin shows is what the model lacks."),
            ("A family of statistics", "a way to read what a fit leaves: moments of powers, fractional powers, inverse and logarithmic kernels, of the derivative and of the running integral, ratios of moments, windows along the frequency, bands along the height, crossings of levels, and the chi of every sample."),
            ("The excess of an answer", "its statistic on the records less the median of the same statistic on the twin, over the span from that median to the 95th percentile of the twin. An excess of 1 is the 95th percentile of the twin."),
            ("The line of the twin", "the excess that the largest answer reaches on 95 of 100 synthetic sets of a model that lacks nothing. An answer past it is evidence that the model lacks something of that place, kind or setting."),
            ("Bias, spread and coverage", "the fit of each model is made again on its own twin, whose truth is known: the bias is the mean distance of the fitted value from that truth over its spread, the spread is over the bar that the fit gives, and the coverage is how often the interval of the fit holds the truth. Each has a band, and the twin admits the fit when all three stand in their bands."))


def probe_said(path: Path) -> str:
    """The answer of the probe of a run, the lines between `# output:` and the next comment of its file, in words. The
    report printed the last two words of the file, which were its exit status (2026-10-18)."""
    if not Path(path).exists():
        return "That question was not recorded."
    lines, inside = [], False
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            inside = line.strip().lower() == "# output:"
            continue
        if inside and line.strip():
            lines.append(line.strip())
    if lines[:2] == ["NONE", "NO"]:
        return "It answered NONE to the first and NO to the second."
    return "It answered: " + _escape(" / ".join(lines) or "nothing") + "."


def at_a_glance(steps: list) -> str:
    """The steps of a run in one table, the first thing a red team of a run of twelve steps asks (review 38, minor 6):
    for each, the model it ended on, its chi square a degree of freedom, whether the twin admits it, the first answers of
    the report to where and what, and the term files that it holds."""
    rows = []
    for step, shown, folder in steps:
        every = loops_of(shown, folder)          # a step that kept copies of the project has a row for each
        for loop in every:
            name = step if len(every) == 1 else f"{step} {Path(loop).relative_to(folder)}"
            k = last_model(loop)
            if k is None:
                rows.append(f"<tr><td>{_escape(name)}</td><td colspan=6>no model of the loop kept</td></tr>")
                continue
            table = json.loads((Path(loop) / f"model_{k}.json").read_text())
            first = (table.get("where") or {}).get("first_of_each_question") or {}
            admitted = ("yes" if table.get("admitted_on_the_twin") else "no") if shares(loop, k) else "nothing to admit"
            terms = ", ".join(p.name for p in sorted((Path(folder) / "terms").glob("*.py"))) or "none"
            rows.append(f"<tr><td>{_escape(name)}</td><td>{k}</td><td>{table['chi_square'] / max(table['degrees_of_freedom'], 1):.3f}</td><td>{admitted}</td>"
                        f"<td>{_escape(first.get('where') or 'nothing beyond the twin')}</td><td>{_escape(first.get('what') or 'nothing beyond the twin')}</td><td>{_escape(terms)}</td></tr>")
    return ("<h2>The steps at a glance</h2><table border=1 cellpadding=3><tr><th>step</th><th>model</th><th>chi square a degree</th><th>the twin admits it</th>"
            "<th>first place of the report</th><th>first kind</th><th>term files</th></tr>" + "".join(rows) + "</table>")


def hidden_truth(project: str, root: Path = ROOT) -> str:
    """The truth of a hidden project of the cold tests, the docstring of its file in tools/cold, which no folder of an
    agent holds."""
    import ast
    path = Path(root) / "tools" / "cold" / f"{project}.py"
    return " ".join((ast.get_docstring(ast.parse(path.read_text())) or "").split()) if path.exists() else f"tools/cold/{project}.py is not found"


def _many(n: int, one: str, more: str) -> str:
    return f"{n} {one if n == 1 else more}"


def run_report(run: Path, out: Path, brief: Path | None = None, root: Path = ROOT) -> Path:
    """A whole blind run as one page for a reader who was not there, his message 226: "I would like to send to my PI a
    full report (I mean log with plot and with explanations of what happened at each turn and how the package has been
    used) of a blind run". What the test is and how it was made blind, how to read the plots, then each turn: what the
    agent did in its order and at what time, the commands of the package that it ran, and the page of its step by
    `step_page`. Then its answer beside the truth that the toy hid, scored by `SCORE`, its notes, and the brief. The
    truth comes last, as the agent never saw it. Refuses a run with no step or no stream."""
    run, out = Path(run), Path(out)
    stream = next((p for p in (run / "stream.jsonl", run / "stream.jsonl.gz") if p.exists()), None)
    steps = steps_of(run)
    if stream is None or not steps:
        raise SystemExit(f"REFUSED: {run} holds {'no stream' if stream is None else 'no step'}: the report reads both")
    ran = json.loads((run / "run.json").read_text()) if (run / "run.json").exists() else {}
    acts_all = actions_of(stream)
    turns = turns_of(acts_all)
    times = [a["at"] for acts in turns.values() for a in acts if a["at"]]
    tally = {k: sum(a["kind"] == k for a in acts_all) for k in ("package", "own fit", "subagent own fit", "term", "red team")}
    scored, project = None, ran.get("project")          # a run on a hidden project of the cold tests carries its name: its truth is that file's, and no toy scores it
    if (run / "answer.json").exists() and not project:
        const, config = model.load(Path(root) / "toy", root)
        scored = score(json.loads((run / "answer.json").read_text()), truth(const, config))
    said = lambda key: _escape(str(ran.get(key, "not recorded")))          # noqa: E731
    title = run.parent.name if run.name == "live" else run.name
    what = ((f"The hidden project {_escape(project)} of the cold tests of the package: curves of another field, written by tools/cold/{_escape(project)}.py, whose truth no file of the agent's folder holds. "
             "A fresh session of Claude Code was given the curves, what the experimenter knows, the package twinfit and a brief, and no name of any term. "
             "It was asked what the curves hold beyond their first model and how large each part is, with every fit made by the package and every term chosen from its report. ") if project else
            ("Eight synthetic traces of the two photon line of rubidium from 5S to 6S, at four temperatures and five powers, made by the toy model of this repository with the noise of the experiment. "
             "The toy holds four terms beyond the natural width: collisions between rubidium atoms, a permeated gas, the noise of the laser and the transit of the atoms through the beam. "
             "A fresh session of Claude Code was given the traces, what the experimenter knows, the package twinfit and a brief, and no name of any term. "
             "It was asked which terms the line holds beyond the natural width and how large each is, with every fit made by the package and every term chosen from its report. "))
    parts = [f"<h1>A blind run of the agentic version of twinfit: {_escape(title)}</h1>",
             f"<p>Run {_escape(title)}: model {said('model')} at the effort {said('effort')}, " + (f"the hidden project {_escape(project)} of the cold tests" if project else f"draw {said('draw')} of the toy") + f", nonce {said('nonce')}. "
             f"The clock of the run read {times[0] if times else 'nothing'} to {times[-1] if times else 'nothing'} UTC. {said('tokens')} tokens, {round(float(ran['seconds']) / 60) if 'seconds' in ran else 'unknown'} minutes. "
             f"{len(steps)} steps kept.</p>"
             + (f"<p><b>How the run ended</b>: {_escape(ran['stopped'])}.</p>" if ran.get("stopped") else "")
             + (f"<p><b>Code of fitting outside the package</b>, found after the run in the files of its folder by a scan of the session that ran it, beside what this report reads in the stream: {_escape(', '.join(ran['own_fits']))}.</p>" if ran.get("own_fits") else ""),
             "<h2>In one paragraph</h2><p>" + (f"The agent named {scored['named']} of the {scored['of']} terms that the toy hid, {scored['by_class']} of their kind, {scored['within_two']} within two of the uncertainty that it gave, and {len(scored['false'])} that the toy does not hold. "
                                               if scored else "Its answer stands at the end beside the truth of the project, which the agent never saw. " if project and (run / "answer.json").exists() else "The run has not ended: it holds no answer yet. ")
             + f"In {len(steps)} steps it used the package {_many(tally['package'], 'time', 'times')}, wrote {_many(tally['term'], 'term file', 'term files')} and called {_many(tally['red team'], 'red team', 'red teams')}. "
             f"It made {_many(tally['own fit'], 'fit', 'fits')} outside the package, which the brief refuses. Its red teams made {tally['subagent own fit']}, which step 7 of the protocol of the package allows a red team, as a candidate that the loop must bear out. The detail of each turn follows, and every action stands in the appendix.</p>",
             at_a_glance(steps),
             "<h2>The package</h2><p>twinfit fits several curves together with one model, and reads what each fit leaves through many families of statistics: generalized moments and their ratios, moments of the derivative and of the running integral, windows along the frequency, bands along the height, crossings of levels and the chi of every sample. "
             "Each statistic is set against its own spread over a twin of the fit, synthetic records drawn from the fitted model with the noise of the experiment and fitted the same way. From what stands beyond the twin its report says where in the line the model misses, of which kind the missing term is and along which setting it changes, and proposes terms. "
             "Its loop fits the first model, then each term that a researcher, or here an agent, writes as a file, keeps or refuses each, and writes the report of each model.</p>",
             "<h2>What this test is</h2><p>" + what +
             "Before the run the same command asked the session which files of rules or instructions it carried, and whether any instruction named the line or its terms. "
             + probe_said(run / "probe.txt") + " "
             f"The numbers here come from synthetic {'curves' if project else 'traces'} whose truth is known, and from one run: they are a test of the method, not a measurement.</p>",
             "<h2>How to read the plots</h2><dl>" + "".join(f"<dt><b>{_escape(a)}</b></dt><dd>{_escape(b)}</dd>" for a, b in GLOSSARY) + "</dl>"]
    table = lambda acts: ("<table border=1 cellpadding=3 style='font-size:90%'><tr><th>clock, UTC</th><th>kind</th><th>what</th></tr>"          # noqa: E731
                          + "".join(f"<tr><td>{a['at'] or ''}</td><td>{a['kind']}</td><td>{_escape(a['what'])}</td></tr>" for a in acts) + "</table>")
    for name, loop, folder in steps:
        acts = turns.get(name, [])
        counts = {k: sum(a["kind"] == k for a in acts) for k in ("package", "own fit", "term", "model", "red team", "reads", "writes", "other", "subagent own fit")}
        counts["other"] += sum(a["kind"].startswith("subagent ") and a["kind"] not in ("subagent says", "subagent own fit") for a in acts)
        parts.append(f"<hr><h2>Turn {_escape(name)}</h2>")
        many = _many
        parts.append(f"<p>{many(counts['package'], 'use', 'uses')} of the package, {many(counts['own fit'], 'fit', 'fits')} outside it by the agent and {counts['subagent own fit']} by its subagents, {many(counts['term'], 'writing', 'writings')} of a term, {counts['model']} of the model, "
                     f"{many(counts['red team'], 'red team', 'red teams')}, and {many(counts['reads'] + counts['writes'] + counts['other'], 'other action', 'other actions')}, every one of them in the appendix.</p>")
        parts.append("<h3>The turn in the agent's words, with its uses of the package and its terms</h3>" + table([a for a in acts if a["kind"] in KEY_KINDS]))
        several = [str(q.relative_to(folder)) for q in loops_of(loop, folder)]
        if len(several) > 1:          # a step that kept copies of the project, each with its loop: the page shows the first
            parts.append(f"<p>This step kept {len(several)} loop folders, one a copy of the project. The page shows {_escape(several[0])}, and the others are {_escape(', '.join(several[1:]))}.</p>")
        parts.append(step_page(name, loop, folder).replace("<h2>", "<h3>").replace("</h2>", "</h3>"))
    appendix = ["<hr><h2>Appendix: every action of the run, turn by turn</h2>"] + [f"<h3>Turn {_escape(name)}</h3>" + table(acts) for name, acts in turns.items() if acts]
    if project and (run / "answer.json").exists():
        answer = json.loads((run / "answer.json").read_text())
        parts.append("<hr><h2>The answer of the agent</h2><table border=1 cellpadding=3><tr><th>origin</th><th>in the model</th><th>parameters</th></tr>"
                     + "".join(f"<tr><td>{_escape(a.get('origin'))}</td><td>{_escape(a.get('in the model'))}</td><td>{_escape(json.dumps(a.get('parameters')))}</td></tr>" for a in answer.get("terms", []))
                     + f"</table><p>{_escape(answer.get('summary', ''))}</p>")
        from toy.agent_example import beside_the_truth
        rows = beside_the_truth(project, answer)
        parts.append("<h2>Beside the truth of the project</h2><p>" + _escape(hidden_truth(project, root)) + "</p>"
                     + ("<table border=1 cellpadding=3><tr><th>parameter of the answer</th><th>value</th><th>its uncertainty</th><th>the truth</th><th>off, in its uncertainties</th></tr>"
                        + "".join(f"<tr><td>{_escape(n)}</td><td>{v:g}</td><td>{b:g}</td><td>{tr:g}</td><td>{o:+.1f}</td></tr>" for n, v, b, tr, o in rows) + "</table>" if rows else ""))
    if scored is not None:
        answer = json.loads((run / "answer.json").read_text())
        s = scored
        parts.append("<hr><h2>The answer of the agent</h2><table border=1 cellpadding=3><tr><th>origin</th><th>shape</th><th>changes with</th><th>full width at 130 C, MHz</th></tr>"
                     + "".join(f"<tr><td>{_escape(a.get('origin'))}</td><td>{_escape(a.get('shape'))}</td><td>{_escape(a.get('changes_with'))}</td><td>{a.get('width_at_130C_MHz')} +- {a.get('uncertainty_MHz')}</td></tr>" for a in answer.get("terms", []))
                     + f"</table><p>{_escape(answer.get('summary', ''))}</p>")
        parts.append("<h2>Beside the truth that the toy hid</h2><table border=1 cellpadding=3><tr><th>term of the toy</th><th>true width, MHz</th><th>named</th><th>of its kind</th><th>width given</th><th>off, in its uncertainties</th></tr>"
                     + "".join(f"<tr><td>{_escape(n)}</td><td>{f['true_width_MHz']:.3f}</td><td>{'yes' if f['named'] else 'no'}</td><td>{'yes' if f['by_class'] else 'no'}</td><td>{f['width_MHz'] if f['width_MHz'] is not None else ''}</td>"
                               f"<td>{'' if f['off_in_uncertainties'] is None else round(f['off_in_uncertainties'], 1)}</td></tr>" for n, f in s["terms"].items())
                     + f"</table><p>{s['named']} of {s['of']} terms named, {s['by_class']} of their kind, {s['within_two']} within two of the uncertainty given, {len(s['false'])} false.</p>")
    for title, path in (("The notes of the agent", run / "notes.md"), ("The brief that the agent was given", Path(brief) if brief else None)):
        if path is not None and path.exists():
            parts.append(f"<hr><h2>{title}</h2>" + prose(path.read_text()))
    parts += appendix
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("<!doctype html><meta charset='utf-8'><title>A blind run of twinfit</title><body style='font-family:Arial;max-width:1100px;margin:auto;line-height:1.35'>" + "\n".join(parts) + "</body>")
    return out


def chosen_loop(folder: Path) -> Path | None:
    """Of the loop folders of a step, the one whose last model leaves a chi square nearest one a degree of freedom: where
    a step kept several copies of the project, as the run on the hidden project did, the first by name was a form that
    the agent set aside. None where the step kept no model."""
    best = None
    for one in loops_of(loop_in(folder), folder):
        k = last_model(one)
        if k is None:
            continue
        table = json.loads((Path(one) / f"model_{k}.json").read_text())
        gap = abs(table["chi_square"] / max(table["degrees_of_freedom"], 1) - 1.0)
        if best is None or gap < best[0]:
            best = (gap, Path(one))
    return None if best is None else best[1]


def withheld(short: str, entry: dict, width: float) -> str:
    """A withheld term as the slide of the design of the test says it: its symbol with the value of its parameter, where
    that is not the width, and its full width at half maximum. "the collisions β_self: 0.124 MHz" put a full width under the
    symbol of a coefficient, and "σ_L: 0.869 MHz" one under a standard deviation of 0.369 (the red team of his message 328,
    moderate 7)."""
    symbol, true, unit = entry["symbol"], entry.get("true"), entry.get("unit", "MHz")
    if isinstance(true, (int, float)) and abs(float(true) - width) > 1e-3 * width:
        return f"{short}: {width:.3g} MHz, ${symbol}$ = {float(true):.3g}" + (" MHz" if unit == "MHz" else "")
    return f"{short} ${symbol}$: {width:.3g} MHz"


def short_name(name: str, words: int = 6) -> str:
    """A name of a term or a model as a slide shows it: its words before the first comma, at most `words` of them."""
    head = str(name).split(",")[0].split(" (")[0].strip()
    cut = head.split()
    return " ".join(cut[:words]) + (" ..." if len(cut) > words else "")


def chain_of(loop: Path, k: int) -> list:
    """The names of the first model and of each term that a loop kept up to its model k, from its best.json, or the
    model's own words where that file is not written."""
    best = Path(loop) / "best.json"
    if best.exists():
        steps = json.loads(best.read_text()).get("iterations") or []
        return [s.get("added") for s in steps[:k + 1] if s.get("added")]
    return [f"model {k}"]


def news(name: str, old: str) -> str:
    """What a name says that an older name of the same model did not: the words after those that both begin with, as a
    model whose name grew from one line to "one line, and a second line" gives "a second line"."""
    a, b = str(name).split(), str(old).split()
    k = 0
    while k < min(len(a), len(b)) and a[k].strip(",") == b[k].strip(","):
        k += 1
    rest = " ".join(a[k:]).lstrip(", ")
    return rest[4:] if rest.startswith("and ") else rest or str(name)


def turns_of_the_run(run: Path) -> list:
    """Each step of a run as the deck reads it: its name, the loop of `chosen_loop`, its last model, the chi square a
    degree of freedom of that model, the first answer of its report beyond the twin over the line of the twin where
    the loop wrote a report, or the answer nearest to that line where none stands beyond it, and the names that the step
    added to the chain of the step before."""
    out, before = [], []
    for name, _, folder in steps_of(run):
        chosen = chosen_loop(folder)
        if chosen is None:
            continue
        k = last_model(chosen)
        table = json.loads((chosen / f"model_{k}.json").read_text())
        where = table.get("where") or {}
        answers = [a for a in where.get("answers", []) if a.get("excess") is not None]
        lead = max(answers, key=lambda a: a["excess"]) if answers else None          # beyond the twin where it stands over its line, else the nearest to it
        chain = chain_of(chosen, k)
        added = [news(c, before[0]) if before and j == 0 else c for j, c in enumerate(chain) if c not in before]          # a model written again says what is new in its name, and a turn that added nothing adds no name
        out.append({"step": name, "loop": chosen, "model": k, "chi": table["chi_square"] / max(table["degrees_of_freedom"], 1), "added": added,
                    "lead": None if lead is None else (f"{lead['question']}: {lead['answer']}", lead["excess"] / max(where.get("largest_excess_on_the_twin", 1.0), 1e-12)),
                    "reported": bool(where)})
        before = chain
    return out


def curves_of(loop: Path, k: int) -> dict:
    """What a model of a loop leaves on each curve, from its residuals_k.csv: the abscissa less the place of the top of
    the fitted curve, the value and the fitted value over the height of the fitted curve, and the residual over the noise."""
    import csv
    import numpy as np
    rows = list(csv.DictReader(open(Path(loop) / f"residuals_{k}.csv")))
    out = {}
    for c in dict.fromkeys(r["curve"] for r in rows):
        got = [r for r in rows if r["curve"] == c]
        x, y, f = (np.array([float(r[key]) for r in got]) for key in ("abscissa", "value", "fitted"))
        z = np.array([float(r["residual"]) / float(r["noise"]) for r in got])
        low, span = float(f.min()), float(f.max() - f.min()) or 1.0
        out[Path(c).stem.replace("trace_", "")] = (x - x[int(np.argmax(f))], (y - low) / span, (f - low) / span, z)
    return out


def sizes_of_the_answer(run: Path, root: Path = ROOT) -> tuple:
    """The values of an answer beside the truth, as `meeting.of_the_sizes` draws them, and one line of what is missed:
    for a run on the toy each term of the toy that the answer holds by its kind, its width at 130 C and the true one,
    and for a run on a hidden project each value of TRUTHS that the answer names."""
    answer = json.loads((Path(run) / "answer.json").read_text())
    ran = json.loads((Path(run) / "run.json").read_text()) if (Path(run) / "run.json").exists() else {}
    if ran.get("project"):
        from toy.agent_example import beside_the_truth
        rows = beside_the_truth(ran["project"], answer)
        return [(short_name(n, 5), v, b, t, "", k % 10 + 1) for k, (n, v, b, t, _) in enumerate(rows)], f"{len(rows)} reported values compared with the true values of the project."
    import tomllib
    const, config = model.load(Path(root) / "toy", root)
    words = tomllib.loads((Path(root) / "toy" / "words.toml").read_text())
    s = score(answer, truth(const, config))
    sizes, missed = [], []
    for name, cell in s["terms"].items():
        k = config["parameter"][name]["enters"]
        match = next((a for a in answer.get("terms", []) if a.get("width_at_130C_MHz") == cell["width_MHz"]), None)
        if cell["by_class"] and match and match.get("uncertainty_MHz"):
            sizes.append((words["detail"][str(k)]["short"], float(cell["width_MHz"]), float(match["uncertainty_MHz"]), float(cell["true_width_MHz"]), "MHz", k))
        else:
            missed.append(words["detail"][str(k)]["short"])
    said = (f"{said_number(s['by_class']).capitalize()} of {said_number(s['of'])} terms identified, " + ("all" if s["within_two"] == s["by_class"] > 1 else said_number(s["within_two"]))
            + f" within two standard uncertainties, {'none' if not s['false'] else said_number(len(s['false']))} false.")          # false against the truth, where "spurious" stood over a term that the analyst called marginal, which its own slide now says (the red team of his message 328)          # where "none spurious" stood over a term that the analyst itself called marginal (the red team of his message 328)
    return sizes, said + (f" Missing: {', '.join(missed)}." if missed else "")


def biased_said(sizes: list, record: dict, said: str, mine: dict) -> list:
    """What the comparison with the true values shows, review 46 for his PI, his message 257: each reported value that
    stands farther than two of its standard uncertainties from its true value, by how many and on which side, then why
    the replicas of the final model, which covered every value, could not see it, and the terms that the run did not
    identify, from the words of `sizes_of_the_answer`, the last two joined by and. Empty where every value agrees within
    two: in BLIND4 the transit time stood 6.6 of them above its truth, on a slide that drew it and said nothing of it."""
    from twinfit.loop import bands
    off = [(name, v, b, true, unit) for name, v, b, true, unit, _ in sizes if b and abs(v - true) > 2.0 * b]
    if not off:
        return []

    def truly(v: float, b: float, true: float) -> str:          # the truth at the decimals of the value, so that the distance reads from the two (review 57, minor 2)
        shown = pages_pm(v, b).split(" ")[0]
        return f"{true:.{len(shown.split('.')[1]) if '.' in shown else 0}f}"
    out = [f"{name[:1].upper() + name[1:]}: {pages_pm(v, b)}{' ' + unit if unit else ''} against a true {truly(v, b, true)}{' ' + unit if unit else ''}, "
           f"{abs(v - true) / b:.1f} standard uncertainties {'above' if v > true else 'below'} it." for name, v, b, true, unit in off]
    cover = sorted({f"{float(p['coverage_on_the_twin']):.2f}" for p in (record or {}).get("parameters", {}).values() if isinstance(p, dict) and p.get("coverage_on_the_twin") is not None})
    if cover:
        sets = int((record or {}).get("sets_of_the_twin") or 0)
        lo, hi = bands(sets or 30)["coverage"]          # right only inside the band of the twin's rule, the package's own: a coverage of 0.20 was said right (review 57, moderate 1)
        right = all(lo <= float(c) <= hi for c in cover)
        out.append(mine["biased_twin" if right else "biased_wrong"].format(sets=sets or "its", covered=cover[0] if len(cover) == 1 else ", ".join(cover[:-1]) + " and " + cover[-1], low=f"{lo:.2f}", high=f"{hi:.2f}"))
    if "Missing: " in said:
        names = said.split("Missing: ", 1)[1].rstrip(".").split(", ")
        out.append(mine["biased_missing"].format(missing=names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]))
    out.append(mine["biased_meaning"])
    return out


def held_said(rows: list, unit: str, mine: dict) -> list:
    """The check of Mercury in sentences, row M314: each size over its theory with its distance in its own bars, and
    where that passes 2, each family kept beside the term held at its theory with its fitted value and the fall of the
    criterion, the families refused in one sentence, then each family kept beside the final model with every term free,
    and the rule. A width carries the unit of the abscissa."""
    out, what = [], {"a Gaussian width": "standard deviation ", "a Lorentzian width": "full width "}          # the number that each width is
    for r in rows:
        name = r["parameter"][:1].upper() + r["parameter"][1:]
        theory = pages_pm(r.get("theory", 1.0), r["theory_bar"]) if r.get("theory_bar") else f"{r.get('theory', 1.0):g}, given with no uncertainty"
        out.append(f"{name}: {pages_pm(r['value'], r['bar'])} against its theory, {theory}: {abs(r['distance']):.1f} of their uncertainties together {'above' if r['distance'] > 0 else 'below'} it. "
                   + (mine["held_near"] if r["held"] is None else mine["held_far"]))
        if r["held"] is None:
            continue
        for b in [b for b in r["beside_held"] if b["ended"] and b["criterion_less"] > 2.0]:
            out.append(f"Held, with {b['family']} beside it, {what.get(b['family'], '')}{pages_pm(b['value'], b['bar'])}{' ' + unit if unit and 'width' in b['family'] else ''}: the criterion lower by {b['criterion_less']:.1f}, kept.")
        refused = [b for b in r["beside_held"] if not (b["ended"] and b["criterion_less"] > 2.0)]
        if refused:
            read = [b for b in refused if b["ended"]]
            moves = ([("higher" if read[0]["criterion_less"] < 0 else "lower") + " by " + " and ".join(f"{abs(b['criterion_less']):.1f}" for b in read)] if len({b["criterion_less"] < 0 for b in read}) == 1 and read
                     else [f"{'lower' if b['criterion_less'] > 0 else 'higher'} by {abs(b['criterion_less']):.1f}" for b in read])
            out.append("Held, with " + " or ".join(b["family"] for b in refused) + " beside it: refused" + (", the criterion " + " and ".join(moves) if moves else "")
                       + (", a fit not ended" if len(read) < len(refused) else "") + ".")
        for b in [b for b in r["beside_free"] if b["ended"] and b["criterion_less"] > 2.0]:
            out.append(f"With every term free, {b['family']} beside them, {what.get(b['family'], '')}{pages_pm(b['value'], b['bar'])}{' ' + unit if unit and 'width' in b['family'] else ''}: the criterion lower by {b['criterion_less']:.1f}, kept.")
        if not any(b["ended"] and b["criterion_less"] > 2.0 for b in r["beside_held"] + r["beside_free"]):          # every family refused: the theory's value, or a shape outside the families (the red team of his message 261)
            out.append(mine["held_none"])
    return out + [mine["held_rule"]]


def curve_label(stem: str) -> tuple:
    """The name of a curve on a slide and the key that orders it: a trace of the toy, trace_P025mW_T130C, as 130 C, 25 mW
    in the order of his message 118, and any other curve by its own name."""
    import re
    found = re.search(r"P(\d+)mW_T(\d+)C", stem)
    if not found:
        return stem.replace("_", " "), (2, stem)
    power, temperature = int(found.group(1)), int(found.group(2))
    return f"{temperature} °C, {power} mW", (0, temperature) if power == 225 else (1, -power)


def unit_of(run: Path) -> str:
    """The unit of the abscissa of the curves of a run, from the first description that its folder holds, or none."""
    import tomllib
    for path in sorted(Path(run).rglob("description.toml")):
        try:
            return str(tomllib.loads(path.read_text()).get("problem", {}).get("abscissa_unit", ""))
        except (tomllib.TOMLDecodeError, OSError):
            continue
    return ""


ASKED = {"where": "location", "shape": "shape", "what": "noise", "changes along": "change", "seen by": "statistic", "along": "dependence on"}          # the questions of the report in the words of the slides, his message 246: in the blind report "what" asks the law of the noise, and the shape has a question of its own. "along" for a run kept before the blind report
BODY = 24.0          # the letters of the sentences of a slide of a run, his message 243: more words than on a slide of the group meeting, and never under its least size
QUOTE = 120          # the words of the agent's reasoning that a slide of a decision quotes: at 70 the physics of the term was cut off


def stop_in_notes(run: Path) -> str:
    """The sentences of the analyst's notes, `notes.md` of the run, that say why it stopped, where no page of a step has a
    section on it: the run of message 310 gave its reason there, and its slide said that it had recorded none (the blind
    red team of 0.1.0rc12, stumble 3)."""
    import re
    path = Path(run) / "notes.md"
    if not path.exists():
        return ""
    text = " ".join(re.sub(r"[*`#]", "", path.read_text()).split())
    return " ".join(s for s in re.split(r"(?<=[.!?])\s+", text) if re.search(r"\b(I|we)\s+(stopped|stop|ended|end)\b", s, re.I))          # a sentence in which the analyst stopped, and not one on the stopping rule of its brief (the red team of his message 320, minor 4)


def reason_of(path: Path) -> str:
    """The reasoning that a term file or a model file gives: the string `why` of a term, or the first paragraph of the
    docstring of a model, read from the file as Python reads it, so that a string written over several lines is whole."""
    import ast
    tree = ast.parse(Path(path).read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "why" for t in node.targets):
            try:
                return " ".join(str(ast.literal_eval(node.value)).split())
            except ValueError:
                return ""
    return " ".join((ast.get_docstring(tree) or "").split("\n\n")[0].split())


def term_file(run: Path, name: str) -> Path | None:
    """The term file of a run that names a term, the last of the run's steps first, or None."""
    found = [p for p in Path(run).rglob("terms/*.py") if name in p.read_text(errors="ignore")] or [p for p in Path(run).rglob("terms/*.py") if built_name(p) == name]          # a file that calls a term of the package by its function names no name, as `mechanisms.gaussian_width()` of the walk by hand of 2026-10-20
    return sorted(found, key=lambda p: (str(p.parent.parent).count("steps"), str(p)))[-1] if found else None


def built_name(path: Path) -> str | None:
    """The name of the term that a term file builds, read by running it as the loop does, or None where it builds none."""
    import runpy
    try:
        term = runpy.run_path(str(path)).get("term")
    except Exception:          # a file of a run that no longer runs names no term
        return None
    return getattr(term, "name", None)


def final_loop(run: Path) -> Path | None:
    """The loop folder of the last step of a run that keeps a model, read by `chosen_loop`: the longest chain of the run
    was, on the hidden project, a form that the agent set aside, with one term more."""
    steps = [folder for _, _, folder in steps_of(run) if chosen_loop(folder) is not None and (chosen_loop(folder) / "best.json").exists()]          # a loop with no output is no final loop (review 42, minor 10)
    return chosen_loop(steps[-1]) if steps else None


def section_of(page: Path, pattern: str, lines: bool = False) -> str:
    """The text of the first section of an agent's page of a step whose heading matches `pattern`, as plain words: its
    marks of Markdown taken out. With `lines`, its lines as the page writes them, so that its paragraphs and its tables
    can be read as such (his message 311). Nothing where the page or the section is not there."""
    import re
    if not Path(page).exists():
        return ""
    parts = re.split(r"(?m)^#{2,3} ", Path(page).read_text())
    for part in parts[1:]:
        head, _, body = part.partition("\n")
        if re.search(pattern, head, re.I):
            plain = re.sub(r"`|\*\*|\[([^\]]*)\]\([^)]*\)", lambda m: m.group(1) or "", body)
            return plain.strip() if lines else " ".join(plain.split())          # the marks of Markdown, and not the underscores of a name of a file
    return ""


CHOICE = r"chosen|choice|why this test|which proposal|\btaken?\b"          # the headings under which the brief asks the agent to say what it chose at a step, and why: since his message 310 it asks "which proposal you take and why" (the red team of message 310, S2)


def report_of(run: Path, k: int, chain: list) -> dict | None:
    """The report of where to search on model k of a chain: the latest loop of the run whose model k holds a report and
    whose chain up to k is that one. None where the agent ran that model without its report."""
    out = None
    for _, _, folder in steps_of(run):
        for one in loops_of(loop_in(folder), folder):
            path = Path(one) / f"model_{k}.json"
            if path.exists() and chain_of(one, k) == chain[:k + 1]:
                table = json.loads(path.read_text())
                if table.get("where"):
                    out = dict(table["where"], chi=table["chi_square"] / max(table["degrees_of_freedom"], 1))
    return out


def fit_of(loop: Path, k: int) -> dict:
    """The shared values of model k of a loop with their bars of the fit, from its notebook, or nothing."""
    found = sorted((Path(loop) / "notebook").glob(f"*_model_{k}_fit.json")) if (Path(loop) / "notebook").is_dir() else []
    if not found:
        return {}
    kept = json.loads(found[-1].read_text())
    return {name: (float(v), float(kept.get("standard_deviations", {}).get(name, float("nan")))) for name, v in kept.get("shared", {}).items()}


def decisions_of(run: Path) -> list:
    """Each decision of a run in its order: the model it was taken on, the term or the model that the agent wrote, its
    reasoning in the agent's words, whether twinfit kept it and why not, the chi square a degree of freedom before and
    after, the fall of the Bayesian criterion, and the values that the new term brought with their bars. A run whose loops
    hold terms is read from the loop of `final_loop`, and a run that wrote its model anew at each step from its steps."""
    import re
    final = final_loop(run)
    its = [it for it in json.loads((final / "best.json").read_text()).get("iterations", []) if re.fullmatch(r"model \d+", str(it.get("model", "")))] if final else []
    out = []
    if len(its) > 1:
        chain = [it["added"] for it in its]
        tables = [json.loads((final / f"model_{k}.json").read_text()) for k in range(len(its))]
        for k, it in enumerate(its):
            for name in it.get("tried", []):
                refused = next((r for r in it.get("refused", []) if r.get("term") == name), None)
                kept = k + 1 < len(its) and its[k + 1]["added"] == name
                path = term_file(run, name)
                before, after = fit_of(final, k), fit_of(final, k + 1) if kept else {}
                out.append({"on": k, "chain": chain, "wrote": name, "reason": reason_of(path) if path else "", "kept": kept, "refused": refused.get("why") if refused else None,
                            "chi": (tables[k]["chi_square"] / max(tables[k]["degrees_of_freedom"], 1), tables[k + 1]["chi_square"] / max(tables[k + 1]["degrees_of_freedom"], 1) if kept else None),
                            "fall": its[k]["criterion"] - its[k + 1]["criterion"] if kept else None, "new": {n: v for n, v in after.items() if n not in before},
                            "moved": {n: (before[n], after[n]) for n in before if n in after}})
        return out
    turns = turns_of_the_run(run)
    folders = {name: folder for name, _, folder in steps_of(run)}
    for j in range(1, len(turns)):          # a run that wrote its model or its data anew at each step: each step is a decision, in the agent's words of the step before
        page = Path(folders[turns[j - 1]["step"]]) / "step.md"
        heading = next((line[2:].split(":", 1)[-1].strip() for line in (Path(folders[turns[j]["step"]]) / "step.md").read_text().splitlines() if line.startswith("# ")), "") if (Path(folders[turns[j]["step"]]) / "step.md").exists() else ""
        model_file = Path(turns[j]["loop"]).parent / "model.py"
        reason = section_of(page, CHOICE) or (reason_of(model_file) if model_file.exists() else "")
        out.append({"on": j - 1, "chain": None, "step": turns[j]["step"], "wrote": turns[j]["added"][0] if turns[j]["added"] else heading or "no new term", "reason": reason,
                    "kept": True, "refused": None, "chi": (turns[j - 1]["chi"], turns[j]["chi"]), "fall": None, "new": {}, "moved": {}, "loop_before": turns[j - 1]["loop"]})
    return out


def _frame(title: str, foot: str, phase: int | None = None):
    """A slide of a run with its title alone, in the look of the slides of the group meeting, and with `phase` the marks of
    the phases of a turn at its foot."""
    from toy import meeting
    fig = meeting.slide(title, "", foot)
    if phase is not None:
        phase_strip(fig, phase)
    return fig


WORDS_SLIDE = 60          # the words of the sentences of one slide of a run at most, as on a page of the folder, his messages 149 and 311: "try to reduce the number of words per slide ... and/or split them on multiple slides"


def split_points(points: list, cap: int = WORDS_SLIDE) -> list:
    """The points of a slide in parts of at most `cap` words each, in their order, a point never cut: a point longer
    than `cap` stands alone. One part where they all hold."""
    parts, words = [[]], 0
    for point in points:
        n = len(str(point).split())
        if parts[-1] and words + n > cap:
            parts.append([])
            words = 0
        parts[-1].append(point)
        words += n
    return parts


def _points(fig, points: list, top: float, left: float = 0.06, width: float = 0.88, size: float = BODY, gap: float = 0.03) -> float:
    """Sentences as points, each wrapped to `width` of the slide, from `top` down. Returns where they end."""
    from toy import meeting, pages
    y = top
    for k, point in enumerate(points):
        body = pages.wrapped(point, size, width * meeting.SIZE[0] / pages.SIZE[0])
        fig.text(left, y - 0.004, "●", fontsize=size, color=pages.GREY, va="top")          # one colour: the bullets took the colours that the slide of the phases gives the package and the analyst (the red team of his message 328, minor)
        fig.text(left + 0.03, y, body, fontsize=size, color=pages.INK, va="top", linespacing=1.25)
        y -= 0.058 * size / 24.0 * (body.count("\n") + 1) + gap
    return y


def _boxes(fig, boxes: list, y: float, colours: list | None = None) -> None:
    """A chain of boxes joined by arrows across the slide at height `y`, each with its words."""
    import textwrap
    from matplotlib.patches import FancyBboxPatch
    from toy import meeting, pages
    ax = fig.add_axes([0.03, y - 0.17, 0.94, 0.30])
    ax.axis("off")
    ax.set_xlim(0.0, len(boxes))
    ax.set_ylim(0.0, 1.0)
    for k, words in enumerate(boxes):
        ax.add_patch(FancyBboxPatch(((k + 0.07) / len(boxes), 0.05), 0.86 / len(boxes), 0.9, boxstyle="round,pad=0.01", fc=(colours or meeting.PALETTE)[k % len(colours or meeting.PALETTE)], ec="none", transform=ax.transAxes))
        ax.text(k + 0.5, 0.5, pages.fill(words, max(8, 75 // len(boxes))), ha="center", va="center", fontsize=meeting.LEAST, color="white", fontweight="bold", linespacing=1.1)
        if k:
            ax.annotate("", xy=(k + 0.07, 0.5), xytext=(k - 0.07, 0.5), arrowprops={"arrowstyle": "-|>", "lw": 2.5, "color": pages.INK})


def _chart(fig, nodes: list, arrows: list, bottom: float = 0.10) -> object:
    """A figure of boxes and arrows in the room of a slide below its title, his message 311: "use visual explanations to
    decrease the number of words per slide". Each node is (x, y, width, height, words, colour), its centre and size as
    shares of the room, its words on lines of its own, and each arrow (from, to, words) between two points of the room,
    its words beside its middle. Returns the axes."""
    from matplotlib.patches import FancyBboxPatch
    from toy import meeting, pages
    ax = fig.add_axes([0.03, bottom, 0.94, fig.below - bottom - 0.02])
    ax.axis("off")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    for x, y, w, h, words, colour in nodes:
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.008", fc=colour, ec="none", transform=ax.transAxes))          # a box of a figure and no bar of a plot, as `_boxes` draws them
        ax.text(x, y, words, ha="center", va="center", fontsize=meeting.LEAST, color="white", fontweight="bold", linespacing=1.1)
    for start, end, words, *bent in arrows:          # an arrow with a bend, as a turn that goes back to the loop, curves by it
        ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "-|>", "lw": 2.5, "color": pages.INK, "mutation_scale": 22, "connectionstyle": f"arc3,rad={bent[0] if bent else 0.0}"})
        level = abs(end[1] - start[1]) < abs(end[0] - start[0])          # the words of an arrow across stand above it, those of an arrow down or up beside it
        if words:
            ax.text((start[0] + end[0]) / 2 + (0.0 if level else 0.012), (start[1] + end[1]) / 2 + (0.03 if level else 0.0), words, ha="center" if level else "left", va="bottom" if level else "center", fontsize=meeting.LEAST, color=pages.INK)
    return ax


def _algorithm_slide(title: str, foot: str, words: dict):
    """The procedure of a turn as a figure, his messages 248 and 311: the fit, the twin and the report, the stop where no
    deviation stands beyond the twin, then the three proposals for the model, the analyst's term with its reason, the
    decision with the three proposals for the experiment, and the next turn. Blue boxes are twinfit's, orange the
    analyst's, green the end."""
    from toy import meeting
    blue, orange, green = meeting.PALETTE[0], meeting.PALETTE[1], meeting.PALETTE[2]
    box = words["algorithm_boxes"]
    fig = _frame(title, foot)
    w, h, top, low = 0.21, 0.30, 0.80, 0.27
    cx = [0.13, 0.38, 0.63, 0.88]
    nodes = [(cx[0], top, w, h, box[0], blue), (cx[1], top, w, h, box[1], blue), (cx[2], top, w, h, box[2], blue), (cx[3], top, w, h, box[3], green),
             (cx[2], low, w, h, box[4], blue), (cx[1], low, w, h, box[5], orange), (cx[0], low, w, h, box[6], blue)]
    arrows = [((cx[0] + w / 2, top), (cx[1] - w / 2, top), ""), ((cx[1] + w / 2, top), (cx[2] - w / 2, top), ""), ((cx[2] + w / 2, top), (cx[3] - w / 2, top), ""),
              ((cx[2], top - h / 2), (cx[2], low + h / 2), words["algorithm_yes"]), ((cx[2] - w / 2, low), (cx[1] + w / 2, low), ""), ((cx[1] - w / 2, low), (cx[0] + w / 2, low), ""),
              ((cx[0], low + h / 2), (cx[0], top - h / 2), words["algorithm_next"])]
    _chart(fig, nodes, arrows)
    meeting.inside(fig, title)
    return fig


def _modes_slide(title: str, foot: str, words: dict, by_hand: bool = False):
    """The two modes of use as two lanes of boxes, his messages 247, 310 and 311: without an agent, the commands and what
    the researcher does, and with one, what the agent and its red team do, each turn going back to the loop. The run
    shown is named under its own lane."""
    from toy import meeting, pages
    blue, orange = meeting.PALETTE[0], meeting.PALETTE[1]
    fig = _frame(title, foot)
    w, h, ys = 0.40, 0.17, [0.84, 0.62, 0.40, 0.18]
    nodes, arrows = [], []
    for x, head, boxes, last in ((0.25, words["modes_without"], words["modes_without_boxes"], blue), (0.75, words["modes_with"], words["modes_with_boxes"], orange)):
        fig.text(x * 0.94 + 0.03, fig.below - 0.005, head, ha="center", va="top", fontsize=meeting.LEAST + 2, fontweight="bold", color=last)
        for k, (y, said) in enumerate(zip(ys, boxes)):
            nodes.append((x, y - 0.03, w, h, said, blue if k < 3 else orange))
            if k:
                arrows.append(((x, ys[k - 1] - 0.03 - h / 2), (x, y - 0.03 + h / 2), ""))
        arrows.append(((x + w / 2, ys[3] - 0.03), (x + w / 2, ys[1] - 0.03), "", 0.6))          # the next turn, back to the loop
    _chart(fig, nodes, arrows, bottom=0.21)
    _say(fig, 0.19, words["modes_foot_by_hand" if by_hand else "modes_foot"], colour=pages.GREY)
    meeting.inside(fig, title)
    return fig


def _traces(ax, loop: Path, k: int, unit: str, one: bool = False) -> None:
    """The curves of model k of a loop, every sample faint and the fit as a line, each over its fitted height, and with `one`
    the hottest and strongest curve alone."""
    from toy import meeting, pages
    curves = sorted(curves_of(loop, k).items(), key=lambda item: curve_label(item[0])[1])
    if one:
        curves = [curves[len(curves) // 2 - 1] if curve_label(curves[0][0])[1][0] == 2 else next(c for c in curves if curve_label(c[0])[0] == "130 °C, 225 mW") if any(curve_label(c[0])[0] == "130 °C, 225 mW" for c in curves) else curves[-1]]
    for j, (name, (x, y, f, z)) in enumerate(curves):
        colour = pages.OF_CONDITION[j % len(pages.OF_CONDITION)] if not one else pages.INK
        ax.plot(x, y, ".", ms=3 if one else 2.5, alpha=0.35 if one else 0.2, color=pages.GREY if one else colour, rasterized=True)
        ax.plot(x, f, "-", lw=2.4 if one else 2.0, color=meeting.WRONG if one else colour, label=curve_label(name)[0])
    ax.set_xlabel(f"distance from the\nline centre{' (' + unit + ')' if unit else ''}", fontsize=meeting.LEAST)          # two lines: the plot takes a third of the slide
    ax.set_ylabel("signal / peak height", fontsize=meeting.LEAST)
    meeting.plain(ax)


def said_number(n: int) -> str:
    """A count as a word, as formal prose writes the counts under ten, and a larger one in figures."""
    return NUMBERS[n] if 0 <= n < min(10, len(NUMBERS)) else str(n)


def pages_pm(value: float, bar: float) -> str:
    """A value with its bar at two significant digits, as every page writes it."""
    from toy import pages
    return pages.pm(value, bar)


def gave_of(run: Path, decided: list, sizes: list) -> list:
    """What twinfit gave the researcher in one run, as sentences with its numbers: the chi square of the first model and
    where the report put the miss, the terms that the report asked for and the chi square they reached, the terms
    refused, how the values of the last model came back from synthetic traces that hold them, and how many of the
    values given lie within two bars of the truth."""
    import re
    out, first = [], decided[0]
    if first["chain"]:
        report = report_of(run, 0, first["chain"])
    else:          # a run that wrote its model anew: the report of the loop of the step before the first change
        before = Path(first["loop_before"]) / f"model_{last_model(first['loop_before'])}.json"
        report = json.loads(before.read_text()).get("where") if before.exists() else None
    lead = [f"{ASKED.get(a['question'], a['question'])}: {a['answer']}" for a in sorted((report or {}).get("answers", []), key=lambda a: -a["excess"]) if a.get("beyond_the_twin")][:1]
    out.append(f"The initial model gave a reduced chi-squared of {first['chi'][0]:.3g}" + (f", and the diagnostic report placed the deviation far above the detection threshold ({' and '.join(lead)})." if lead else "."))
    kept = [d for d in decided if d["kept"]]
    refused = [d for d in decided if not d["kept"]]
    if kept and first["chain"]:
        out.append(f"{said_number(len(kept)).capitalize()} {'terms' if len(kept) > 1 else 'term'}, each written by the analyst after the report of the model before it, reduced it to {kept[-1]['chi'][1]:.3f}.")          # written after the report, and not motivated by it: the term of turn 1 of the run of message 310 came from the analyst's physics (the blind red team of 0.1.0rc12, stumble 2)
    elif kept:
        out.append(f"{said_number(len(kept)).capitalize()} {'changes' if len(kept) > 1 else 'change'} of the model or of the treatment of the data, each written by the analyst after the report of the model before it, reduced it to {kept[-1]['chi'][1]:.3f}.")
    final = final_loop(run)
    introduced = next((d for d in kept if d["new"]), None)
    if introduced and final:          # a value that looked settled when its term was written, moved by the terms written after it
        name, (v1, b1) = next(iter(introduced["new"].items()))
        k_last = last_model(final)
        v2, b2 = fit_of(final, k_last).get(name, (v1, b1))
        if b1 and b2 and abs(v2 - v1) > 2.0 * b2:          # counted in the bar that the value was quoted with then: what a reader of that fit would have been told
            out.append(f"{name.split(',')[0][:1].upper() + name.split(',')[0][1:]} was estimated as {pages_pm(v1, b1)} when its term was introduced, at a reduced chi-squared of {introduced['chi'][1]:.3f}, and as {pages_pm(v2, b2)} in the final model: a shift of {abs(v2 - v1) / b1:.0f} times the uncertainty quoted at that stage.")
    if refused:
        out.append(f"{said_number(len(refused)).capitalize()} proposed {'terms' if len(refused) > 1 else 'term'} not supported by the data {'were' if len(refused) > 1 else 'was'} rejected by the Bayesian information criterion.")
    best = json.loads((final / "best.json").read_text()).get("parameters", {}) if final and (final / "best.json").exists() else {}
    cover = [float(v["coverage_on_the_twin"]) for v in best.values() if isinstance(v, dict) and v.get("coverage_on_the_twin") is not None]
    bias = [abs(float(v["bias_on_the_twin_in_bars"])) for v in best.values() if isinstance(v, dict) and v.get("bias_on_the_twin_in_bars") is not None]
    if cover and bias:
        k_last = last_model(final)
        sets = int(json.loads((final / "best.json").read_text()).get("sets_of_the_twin") or phase_file(final, k_last, "twin").get("admission", {}).get("sets") or 0)          # the replicas of the output, as best.json counts them (review 41, minor)
        values = sorted({f"{c:.2f}" for c in cover})
        covered = values[0] if len(values) == 1 else ", ".join(values[:-1]) + " and " + values[-1]
        out.append(f"On {sets} replicas drawn for the output, apart from those of its checks, every parameter came back with a coverage of {covered}" + (f", one replica being {1.0 / sets:.2f}," if sets else "")
                   + f" and a bias of at most {max(bias):.2f} standard uncertainties.")          # a number keeps its condition (review 39, moderate 9)
    if sizes:
        within = sum(1 for _, v, b, true, _, _ in sizes if abs(v - true) <= 2.0 * b)
        if len(sizes) == 1:
            out.append(f"The single reported value {'agrees' if within else 'does not agree'} with the true value within two standard uncertainties.")
        else:
            out.append(("All " + said_number(len(sizes)) + " reported values agree" if within == len(sizes) else f"{said_number(within).capitalize()} of the {said_number(len(sizes))} reported values {'agrees' if within == 1 else 'agree'}") + " with the true values within two standard uncertainties.")
    return [re.sub(r"\s+", " ", s) for s in out]


def model_parts(path: Path) -> dict:
    """SHARED, OWN and HELD of a model file, as it writes them, read without running it: what is joint, what is free in
    each curve and what is held from theory. Each is empty where the file does not write it."""
    import ast
    out = {"SHARED": (), "OWN": (), "HELD": {}}
    if not Path(path).exists():
        return out
    for node in ast.parse(Path(path).read_text()).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if getattr(target, "id", None) in out:
                    try:
                        out[target.id] = ast.literal_eval(node.value)
                    except ValueError:
                        pass
    return out


def sweeps_of(names: list) -> list:
    """The curves of a run in the sweeps of his message 245, each as (its words, its names in order): for the traces of
    the line, the temperature sweep at the highest power and the power sweep at the highest temperature, in the order of
    his message 118, and for other curves every curve in one."""
    import re
    found = {n: re.search(r"P(\d+)mW_T(\d+)C", n) for n in names}
    if names and all(found.values()):
        power = {n: int(m.group(1)) for n, m in found.items()}
        heat = {n: int(m.group(2)) for n, m in found.items()}
        top, hot = max(power.values()), max(heat.values())
        return [(f"temperature series at {top} mW", sorted((n for n in names if power[n] == top), key=lambda n: heat[n])),
                (f"power series at {hot} °C", sorted((n for n in names if heat[n] == hot), key=lambda n: -power[n]))]
    return [("all spectra", sorted(names))]


def joint_said(parts: dict, fitted: dict, chi: float, count: int, unit: str = "") -> str:
    """The words under the plots of a joint fit: its chi square, what is one value for all the curves with its value and
    its bar, what is held from theory, and what is free in each curve."""
    joint = "; ".join(f"{name.split(',')[0]} = {pages_pm(v, b)}" for name, (v, b) in fitted.items()) or "none"
    held = ", ".join(f"{name.replace('_', ' ')} = {value:.5g}" + (f" {unit}" if unit and name.endswith("width") else "") for name, value in parts["HELD"].items()) or "none"
    free = ", ".join(parts["OWN"]) or "none"
    return (f"Joint fit of all {count} spectra, reduced chi-squared {chi:.3f}. Global parameters: {joint}. "
            f"Fixed at theoretical values: {held}. Local parameters, one per spectrum: {free}.").replace("; ", ", ")


def rebuilt(run: Path, final: Path | None = None) -> list:
    """Each model of the chain of a run fitted again by the package, as (problem, fit), from the project whose noise is
    the one that the loop of `final_loop` read and from the model file and the term files of its step: what the windows
    of the deck are read on. Nothing where that project is not in the run or the chain was not made of term files."""
    from twinfit.refine import develop_one
    first = chain_terms(run, final)
    if first is None:
        return []
    problem, result, terms = first
    out = [(problem, result)]
    for _, m in terms:
        grown = develop_one(*out[-1], m)
        if not grown["ended"]:
            return []
        out.append((grown["_problem"], grown["_fit"]))
    return out


def chain_terms(run: Path, final: Path | None = None):
    """The first model of the chain of a run fitted by the package, from the project whose noise is the one that the loop
    of `final_loop` read, with the term of each later model as (its name, its Mechanism), read from the term files of its
    step: (problem, fit, terms). None where that project is not in the run or a term of the chain has no file."""
    import csv
    import shutil
    import tempfile
    from twinfit import loop
    final = Path(final) if final is not None else final_loop(run)
    if final is None or not (final / "residuals_0.csv").exists() or not (final.parent / "model.py").exists():
        return None
    first = next(csv.DictReader(open(final / "residuals_0.csv")))

    def noise(description: Path) -> float | None:
        """The noise of the first sample of the first curve, as the project gives it: the third column of its file, or the
        one value of its description."""
        import tomllib
        row = list(csv.reader(open(description.parent / first["curve"])))[1]
        if len(row) >= 3:
            return float(row[2])
        value = tomllib.loads(description.read_text()).get("noise", {}).get("value")
        return float(value) if value is not None else None
    samples = sum(1 for row in csv.DictReader(open(final / "residuals_0.csv")) if row["curve"] == first["curve"])

    def same(description: Path) -> bool:
        """Whether a project holds the data that the loop read: its first curve with as many samples as the loop's, the
        same first value and the same noise, as the two files print them, to six digits. The noise alone took the
        project of 3001 samples a curve for its copies thinned to 376, which print the same first sample (review 44,
        moderate 1)."""
        path = description.parent / first["curve"]
        if not path.exists():
            return False
        rows = [row for row in list(csv.reader(open(path)))[1:] if row]
        value, sd = noise(description), float(first["noise"])
        return (len(rows) == samples and abs(float(rows[0][1]) - float(first["value"])) <= 1e-6 * max(abs(float(first["value"])), 1e-12)
                and value is not None and abs(value - sd) < 1e-5 * abs(sd))
    project = next((d.parent for d in sorted(Path(run).rglob("description.toml")) if same(d)), None)          # and no project where none holds them, never the first that prints the noise
    chain = chain_of(final, last_model(final))
    if project is None:
        return None
    home = Path(tempfile.mkdtemp(prefix="deck_"))
    try:
        shutil.copy(project / "description.toml", home)
        shutil.copytree(project / "traces", home / "traces")
        shutil.copy(final.parent / "model.py", home)
        (home / "terms").mkdir()
        _, problem, _, source = loop.initial(home, with_module=True)
        built = {}
        for path in sorted((final.parent / "terms").glob("*.py")):          # each term file by the name of the term it builds: a name that the package makes from a mechanism stands in no file, and the windows of every turn of the run of his message 310 were left out
            try:
                term = loop._module(path, f"deck_{path.stem}", getattr(source, "HELD", None), source).term
            except Exception:          # a file that does not build: no term of the chain is read from it
                continue
            built.setdefault(getattr(term, "name", None), term)
        terms = []
        for name in chain[1:]:
            if name not in built:
                return None
            terms.append((name, built[name]))
        return problem, loop.searched(problem), terms
    finally:
        shutil.rmtree(home, ignore_errors=True)


OVER_THEORY = "over its theory"          # a parameter that is a size over its value of theory, one there, as the blind runs wrote the transit
FAMILIES = ("a Gaussian width", "a Lorentzian width", "a pedestal")          # what is offered beside a term held at its theory: the widths of another shape, and the generic wings


def width_beside(shape: str, along=None, name: str | None = None, parameter: str | None = None):
    """A width of one value for every curve beside the line, a Gaussian or a Lorentzian: the families of the package
    since version 0.0.9, `mechanisms.gaussian_width` and `mechanisms.lorentzian_width`, which lived here while the
    count of row M250 ran on 0.0.8."""
    from twinfit import mechanisms
    return (mechanisms.gaussian_width if shape == "gaussian" else mechanisms.lorentzian_width)(along=along, name=name, parameter=parameter)


def held_at_theory(run: Path, final: Path | None = None, limit: float = 2.0) -> list:
    """The check of Mercury, his cloud chat's red team of message 259, row M314: each parameter of the final model of a
    run that is a size over its value of theory, with its distance from one in its own bars, and where that passes
    `limit` the chain fitted again with its term held at one, alone and with each of `FAMILIES` beside it, and with each
    family beside the final model, all free. Each row gives the change of the loop's criterion against the final model,
    so that a family is kept where it falls by more than 2, the loop's own rule. In BLIND4 the transit stood at 1.150 of
    its theory and the term that it took up was the laser noise, a Gaussian width, which the run never offered: a
    distance alone is no detector, since a fifth of the draws with no term missing raise it, and the criterion is the
    test. Nothing where the chain cannot be fitted again."""
    first = chain_terms(run, final)
    return [] if first is None else checked(*first, limit, beam_bar_of(run))


def beam_bar_of(run: Path) -> float | None:
    """The bar that the radius of the beam puts on a size over its theory, from the first description of a run that
    gives the beam with its uncertainty: twice the relative uncertainty of the radius, since the size of the transit
    goes as one over its square. None where no description gives it. Without it the distance of BLIND4's transit from
    its theory was said in the fit's bar alone, 7.6 where the beam of 43 +- 1 um makes it 3.0 (the red team of his
    message 260)."""
    import tomllib
    for path in sorted(Path(run).rglob("description.toml")):
        try:
            beam = tomllib.loads(path.read_text()).get("setup", {}).get("beam", {})
        except (tomllib.TOMLDecodeError, OSError):
            continue
        if beam.get("value") and beam.get("uncertainty") is not None:
            return 2.0 * float(beam["uncertainty"]) / float(beam["value"])
    return None


def checked(problem, result, terms: list, limit: float = 2.0, theory_bar: float | None = None) -> list:
    """`held_at_theory` on a chain given as its first model, its fit and its terms, (name, Mechanism): what toy/mercury.py
    counts on draws. A term that gives its theory, `theory = (centre, bar)`, is read against it, and else against one
    with `theory_bar`: the distance counts the fit's bar and the theory's together, since a theory known to 4.7 per
    cent read in the fit's bar alone raised the check on a quarter of the draws that lack no term (the red team of his
    message 260). Nothing where a term of the chain does not end."""
    from twinfit import fit, loop, mechanisms
    from twinfit.refine import develop_one
    chain = [(problem, result)]
    for _, m in terms:
        grown = develop_one(*chain[-1], m)
        if not grown["ended"]:
            return []
        chain.append((grown["_problem"], grown["_fit"]))
    last = loop._criterion(*chain[-1])
    offered = {"a Gaussian width": lambda: width_beside("gaussian"), "a Lorentzian width": lambda: width_beside("lorentzian"), "a pedestal": lambda: mechanisms.pedestal()}

    def beside(on, family):
        got = develop_one(*on, offered[family]())
        if not got["ended"]:
            return {"family": family, "ended": False}
        k = len(on[0].shared)
        return {"family": family, "ended": True, "value": float(got["_fit"].shared[k]), "bar": float(got["_fit"].sd[k]), "criterion_less": float(last - loop._criterion(got["_problem"], got["_fit"]))}
    import math
    out = []
    for name, m in terms:
        if OVER_THEORY not in m.parameter:
            continue
        k = list(chain[-1][0].shared).index(m.parameter)
        value, bar = float(chain[-1][1].shared[k]), float(chain[-1][1].sd[k])
        centre, tb = (float(v) for v in getattr(m, "theory", None) or (1.0, theory_bar or 0.0))
        row = {"term": name, "parameter": m.parameter, "value": value, "bar": bar, "theory": centre, "theory_bar": tb, "distance": (value - centre) / math.hypot(bar, tb),
               "held": None, "beside_held": [], "beside_free": []}
        if abs(row["distance"]) > limit:
            on = (problem, result)
            for other, mo in terms:          # the chain again, this term held at its theory and every other one free
                if other == name:
                    held = mo.hold(on[0], centre)
                    on = (held, fit(held, on[1].shared, on[1].own))
                else:
                    grown = develop_one(*on, mo)
                    on = (grown["_problem"], grown["_fit"])
            row["held"] = float(last - loop._criterion(*on))
            row["beside_held"] = [beside(on, family) for family in FAMILIES]
            row["beside_free"] = [beside(chain[-1], family) for family in FAMILIES]
        out.append(row)
    return out


def _sweep_slide(title: str, curves: dict, names: list, said: str, foot: str, unit: str, phase: int | None = None):
    """One sweep of a model: each curve in a row, every sample faint and the joint fit as a line, what the fit leaves
    under each in units of its noise with its local mean as a dense curve, and below them the words of the joint fit."""
    import numpy as np
    from toy import meeting, pages
    fig = _frame(title, foot, phase)
    n = len(names)
    left, right, gap = 0.08, 0.97, 0.012
    width = (right - left - gap * (n - 1)) / n
    body = pages.wrapped(said, meeting.LEAST, 0.9 * meeting.SIZE[0] / pages.SIZE[0])
    words_top = 0.085 + 0.047 * (body.count("\n") + 1)          # the words above the foot, as many lines as they take, and the plots in the room left
    low_base = words_top + 0.125
    top_base = low_base + 0.20
    for j, name in enumerate(names):
        x, y, f, z = curves[name]
        colour = pages.OF_CONDITION[j % len(pages.OF_CONDITION)]
        top = fig.add_axes([left + j * (width + gap), top_base, width, fig.below - 0.08 - top_base])
        top.plot(x, y, ".", ms=2, alpha=0.25, color=pages.GREY, rasterized=True)
        top.plot(x, f, "-", lw=2.2, color=colour)
        top.set_title(curve_label(name)[0], fontsize=meeting.LEAST, color=pages.INK, pad=6)
        top.set_xticks([])
        top.set_ylim(-0.12, 1.15)
        low = fig.add_axes([left + j * (width + gap), low_base, width, 0.18])
        low.plot(x, z, ".", ms=1.5, alpha=0.15, color=pages.GREY, rasterized=True)
        low.plot(*smooth_mean(x, z), "-", lw=2.2, color=colour)          # 400 points, his message 247: 30 bins were "so unpleasant"
        low.axhline(0.0, color=pages.INK, lw=0.8)
        reach = max(4.0, float(np.percentile(np.abs(z), 99.5)))
        low.set_ylim(-reach, reach)
        span = float(max(abs(x.min()), abs(x.max())))
        tick = 10.0 * max(1, int(span // 20))
        low.set_xticks([-tick, 0.0, tick], labels=[f"{-tick:g}", "0", f"{tick:g}"])
        low.set_xlabel(unit or "abscissa", fontsize=meeting.LEAST, labelpad=2)          # each axis of numbers named, as every page asks: the distance from the line centre
        for ax in (top, low):
            meeting.plain(ax)
            ax.tick_params(labelsize=meeting.LEAST)
            if j:
                ax.set_yticks([])
        if j == 0:
            top.set_ylabel("signal", fontsize=meeting.LEAST)
            low.set_ylabel("residual", fontsize=meeting.LEAST)
    fig.text(0.06, words_top, body, fontsize=meeting.LEAST, color=pages.INK, va="top", linespacing=1.15)
    meeting.inside(fig, title)
    return fig


BANDS_SAID = ("0 to 0.25", "0.25 to 0.5", "0.5 to 1", "1 to 2", "2 to 4", "4 to 8")          # the pairs of windows of `report.BANDS`, in full widths of the line from its centre


def _windows_slide(title: str, problem, result, foot: str, text: str, phase: int | None = None):
    """The windows of what one model leaves, as the report of the package reads them: one row a curve, one column a pair
    of windows at equal distances from the centre, the number of each the residuals summed over it over their spread
    after the fit, both sides added on the left and one side less the other on the right. A complete model leaves noise
    about zero."""
    import numpy as np
    from twinfit.report import by_position_and_record
    from toy import meeting, pages
    rows = by_position_and_record(problem, result)
    names = [Path(rec.name).stem for rec in problem.records]
    order = sorted(range(len(names)), key=lambda k: curve_label(names[k])[1])
    fig = _frame(title, foot, phase)
    body = pages.wrapped(text, meeting.LEAST, 0.9 * meeting.SIZE[0] / pages.SIZE[0])
    fig.text(0.06, fig.below - 0.01, body, fontsize=meeting.LEAST, color=pages.INK, va="top", linespacing=1.2)
    base, height = 0.20, fig.below - 0.32 - 0.05 * (body.count("\n") + 1)          # the name of the axis above the foot, and the heads of the maps below the words
    image = None
    for c, (side, head) in enumerate((("both_sides", "symmetric sum"), ("one_side_less_the_other", "antisymmetric difference"))):
        z = np.array([[row[side][k] for row in rows] for k in order], dtype=float)
        ax = fig.add_axes([0.23 + c * 0.34, base, 0.30, height])
        image = ax.imshow(np.clip(z, -5.0, 5.0), cmap="RdBu_r", vmin=-5.0, vmax=5.0, aspect="auto")
        for i in range(z.shape[0]):
            for j in range(z.shape[1]):
                if np.isfinite(z[i, j]) and abs(z[i, j]) >= 3.0:          # the numbers beyond 3 written, the rest left to the colour
                    ax.text(j, i, f"{z[i, j]:.0f}", ha="center", va="center", fontsize=meeting.LEAST, color="white", fontweight="bold")
        ax.set_xticks(range(len(rows)), labels=[f"{row['from']:g}" for row in rows], fontsize=meeting.LEAST)
        ax.set_yticks(range(len(order)), labels=[curve_label(names[k])[0] for k in order] if c == 0 else [], fontsize=meeting.LEAST)
        ax.set_title(head, fontsize=meeting.LEAST, color=pages.INK, pad=6)
        ax.set_xlabel("inner edge (line widths)", fontsize=meeting.LEAST)
        ax.tick_params(length=0)
    bar = fig.colorbar(image, cax=fig.add_axes([0.885, base, 0.012, height]), ticks=[-5, 0, 5])
    bar.ax.tick_params(labelsize=meeting.LEAST)
    bar.set_label("standard deviations", fontsize=meeting.LEAST)
    meeting.inside(fig, title)
    return fig


def _how_windows_slide(problem, result, foot: str, words: dict):
    """How a window is read, on one curve of the first model: what the fit leaves against the distance from the centre in
    widths of the line, the pairs of windows shaded, and the number of each written in it."""
    import numpy as np
    from twinfit.linear import Linear
    from twinfit.report import BANDS, by_position_and_record
    from toy import meeting, pages
    lin = Linear(problem, result)
    names = [Path(rec.name).stem for rec in problem.records]
    k = max(range(len(names)), key=lambda j: (curve_label(names[j])[0] == "130 °C, 225 mW", -j))
    rec = problem.records[k]
    u = (np.asarray(rec.x) - lin.centre_of(k)) / lin.width_of(k)
    z = (np.asarray(rec.y) - np.asarray(problem.predict(k, result.shared, result.own[k]))) / np.asarray(rec.sd)
    rows = by_position_and_record(problem, result)
    fig = _frame(words["windows_title"], foot)
    largest = max((abs(rows[j]["both_sides"][k]) for j in range(len(rows)) if np.isfinite(rows[j]["both_sides"][k])), default=0.0)
    body = pages.wrapped(words["windows"].format(largest=f"{largest:.0f}"), meeting.LEAST, 0.9 * meeting.SIZE[0] / pages.SIZE[0])          # the largest window of this spectrum, from the run
    fig.text(0.06, fig.below - 0.01, body, fontsize=meeting.LEAST, color=pages.INK, va="top", linespacing=1.2)
    ax = fig.add_axes([0.10, 0.20, 0.86, fig.below - 0.30 - 0.05 * (body.count("\n") + 1)])          # the plot under the words, as many lines as they take
    for j, (lo, hi) in enumerate(zip(BANDS, BANDS[1:])):
        for sign in (1, -1):
            ax.axvspan(sign * lo, sign * hi, color=meeting.PALETTE[j % len(meeting.PALETTE)], alpha=0.15, lw=0)
        value = rows[j]["both_sides"][k] if j < len(rows) else float("nan")
        if np.isfinite(value) and hi <= 4.0:
            ax.annotate(f"{value:.0f}", (0.5 * (lo + hi), 1.0), xycoords=("data", "axes fraction"), xytext=(0, -4 - 30 * (j % 2)), textcoords="offset points",
                        ha="center", va="top", fontsize=meeting.LEAST, color=pages.INK, fontweight="bold")          # the narrow inner windows at two heights, 30 points apart, so that their numbers stand apart
    ax.plot(u, z, ".", ms=2.5, alpha=0.5, color=pages.INK, rasterized=True)
    ax.axhline(0.0, color=pages.INK, lw=0.8)
    ax.set_xlim(-4.0, 4.0)
    shown = z[np.abs(u) <= 4.0]
    ax.set_ylim(float(shown.min()) * 1.05, float(shown.max()) + 1.2 * float(shown.max() - shown.min()))          # the top half free of points, for the numbers of the windows at their two heights
    ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])          # no number at the corner, where it met those of the residual
    ax.set_xlabel(f"distance from the line centre (line widths), spectrum at {curve_label(names[k])[0]}", fontsize=meeting.LEAST)
    ax.set_ylabel("normalized residual", fontsize=meeting.LEAST)
    meeting.plain(ax)
    ax.tick_params(labelsize=meeting.LEAST)
    meeting.inside(fig, words["windows_title"])
    return fig


def smooth_mean(x, z, points: int = 400, share: float = 0.012) -> tuple:
    """The local mean of the normalized residuals as a smooth curve on `points` abscissae, his message 247: "the line of
    the residual average has too little points". Each point is the mean of the residuals weighted by a Gaussian of a
    standard deviation `share` of the span of the abscissa, so the curve is dense and follows a structure of the width of
    the line core."""
    import numpy as np
    x, z = np.asarray(x, dtype=float), np.asarray(z, dtype=float)
    grid = np.linspace(float(x.min()), float(x.max()), points)
    sigma = share * float(x.max() - x.min()) or 1.0
    weights = np.exp(-0.5 * ((grid[:, None] - x[None, :]) / sigma) ** 2)
    return grid, (weights @ z) / weights.sum(axis=1)


def phase_file(loop: Path, k: int, phase: str) -> dict:
    """What the notebook of a loop kept of one phase of its model k, `intake`, `identify`, `twin`, `fit`, `read`,
    `propose`, `report` or `candidate_on_the_twin`, the latest of them, or nothing where the loop kept no such file."""
    book = Path(loop) / "notebook"
    found = sorted(book.glob(f"*_model_{k}_{phase}.json")) if book.is_dir() else []
    return json.loads(found[-1].read_text()) if found else {}


def proposed(run: Path, loop: Path, k: int) -> tuple:
    """The report of where to search on model k of a loop and the generic candidates that twinfit fitted for it, from the
    proposal phase of its notebook: of the loop itself, or where it ran that model for a count with no report, of the
    latest loop of the run whose chain up to k is the same. A loop with no best.json is read alone, since its chain is no
    more than its number. Returns (report, candidates), each empty where none was made."""
    chain = chain_of(loop, k)
    found = [Path(loop)]
    if (Path(loop) / "best.json").exists():
        found += [Path(one) for _, _, folder in reversed(steps_of(run)) for one in loops_of(loop_in(folder), folder) if Path(one) != Path(loop)]
    for one in found:
        kept = phase_file(one, k, "propose")
        if kept.get("where") and chain_of(one, k) == chain:
            return kept["where"], kept.get("candidates", [])
    return {}, []


def threshold_of(report: dict) -> float:
    """The detection threshold of a report: the 95th percentile over the replicas of the largest excess of any
    statistic, one where the report does not give it."""
    return float(report.get("largest_excess_on_the_twin") or 1.0) or 1.0


SAID = {"where": "located in", "shape": "shaped as", "what": "with", "changes along": "a change of", "seen by": "seen by", "along": "dependent on"}          # the questions of the report as the opening of a phrase of a sentence


def detected(report: dict, most: int | None = 3, questions: tuple | None = None) -> list:
    """The deviations that a report detects, as phrases with their excess over the detection threshold: the answers beyond
    the threshold, the largest first, at most `most` of them, all where `most` is None, and of the `questions` alone where
    they are named."""
    line = threshold_of(report)
    beyond = sorted((a for a in report.get("answers", []) if a.get("beyond_the_twin") and (questions is None or a["question"] in questions)), key=lambda a: -a["excess"])[:most]
    return [f"{SAID.get(a['question'], a['question'])} {a['answer']}, {times_threshold(a['excess'] / line)} times the threshold" for a in beyond]


def times_threshold(ratio: float) -> str:
    """An excess over the detection threshold as words: two significant digits, and three near the threshold itself, where
    "1.0 times the threshold" was called past it (review 47)."""
    return f"{ratio:.2f}" if 0.95 <= ratio < 1.1 else plain_number(ratio)


def phase_strip(fig, current: int, count: int = 8) -> None:
    """The phases of a turn as numbered marks at the foot of a slide, the current one filled, so that a reader knows where
    in the turn the slide stands."""
    from toy import meeting, pages
    for j in range(count):
        on = j + 1 == current
        fig.text(0.665 + 0.042 * j, 0.03, str(j + 1), fontsize=meeting.LEAST, ha="center", va="bottom", fontweight="bold", color="white" if on else pages.GREY,
                 bbox={"boxstyle": "circle,pad=0.22", "fc": meeting.PALETTE[0] if on else "white", "ec": meeting.PALETTE[0] if on else pages.GREY, "lw": 1.4})


def model_said(loop: Path, k: int, words: dict, unit: str = "", why: bool = False) -> list:
    """The model of one turn as sentences, his message 247: "the model used are not clearly states in each turn..., it's
    not clear what is fixed, how the joint fits are performed, why". Its components with the turn that added each, what is
    held at a value of theory, its global and its local parameters, the samples and the degrees of freedom, and how and
    why the fit is joint."""
    parts = model_parts(Path(loop).parent / "model.py")
    chain = chain_of(loop, k)
    fitted = fit_of(loop, k)
    intake, read = phase_file(loop, k, "intake"), phase_file(loop, k, "read")
    records = int(intake.get("records") or len(curves_of(loop, k)))
    out = [f"Component {j + 1}, " + ("the initial model" if j == 0 else f"added in turn {j}") + f": {name if j == 0 else short_name(name, 14)}." for j, name in enumerate(chain)]          # a term by its name before its first comma, whose whole name stands on the slide of the turn that wrote it, and the model whole, which a run that wrote its model anew names in full
    held = [f"{name.replace('_', ' ')} = {value:.5g}" + (f" {unit}" if unit and name.endswith("width") else "") for name, value in parts["HELD"].items()]
    out.append("Fixed at theoretical values: " + (", ".join(held) if held else "none") + ".")
    out.append("Global parameters, one value for all spectra: " + (", ".join(fitted) if fitted else "none, so the model has no free physical parameter") + ".")
    own = list(parts["OWN"])
    out.append("Local parameters, one value per spectrum: " + (", ".join(own) + f", that is {len(own)} for each of the {records} spectra, {len(own) * records} in all" if own else "none") + ".")
    if intake.get("samples") and read.get("degrees_of_freedom"):
        out.append(f"Data: {intake['samples']} samples in {records} spectra, leaving {read['degrees_of_freedom']} degrees of freedom.")
    return out + ([words["joint_how"], words["joint_why"]] if why else [])          # said at the first turn alone, and on the slide of its definitions (his message 311: fewer words a slide)


def _say(fig, y: float, text: str, colour: str | None = None, weight: str = "normal", size: float | None = None, gap: float = 0.015) -> float:
    """One sentence across the slide from `y` down, wrapped to its width. Returns where the next one starts."""
    from toy import meeting, pages
    size = size or meeting.LEAST
    body = pages.wrapped(text, size, 0.90 * meeting.SIZE[0] / pages.SIZE[0], weight)
    fig.text(0.05, y, body, fontsize=size, color=colour or pages.INK, va="top", fontweight=weight, linespacing=1.15)
    return y - 0.047 * size / 21.5 * (body.count("\n") + 1) - gap


def plain_number(value: float) -> str:
    """A number at two significant digits, written without an exponent: 5800, 1.5, 0.27, and 1.0 where the second digit
    is a zero, which "1 times the threshold" hid."""
    said = f"{float(f'{value:.2g}'):g}" if abs(value) < 1e6 else f"{value:.2g}"
    return f"{float(said):.1f}" if said.lstrip("-").isdigit() and abs(value) < 10 else said


def _model_slide(title: str, points: list, foot: str) -> list:
    """Phase 1 of a turn: the model stated, as `model_said` writes it, on as many slides as its words ask at sixty each
    (his message 311)."""
    from toy import meeting
    parts = split_points(points)
    figs = []
    for j, part in enumerate(parts, 1):
        head = title if len(parts) == 1 else f"{title}, {j} of {len(parts)}"
        fig = _frame(head, foot, 1)
        _points(fig, part, fig.below - 0.03, width=0.90, size=meeting.LEAST, gap=0.012)
        meeting.inside(fig, head)
        figs.append(fig)
    return figs


def _table(fig, top: float, heads: list, rows: list, columns: list, first_width: float, colours: list | None = None) -> float:
    """Words in columns from `top` down: a row of heads, then each row, the first column wrapped to `first_width` of the
    slide and the others written as they are. Returns where the table ends."""
    from toy import meeting, pages
    size = meeting.LEAST
    step = 0.046
    lines = max(h.count("\n") + 1 for h in heads)
    for x, h in zip(columns, heads):
        fig.text(x, top, h, fontsize=size, color=pages.GREY, va="top", linespacing=1.05, fontweight="bold")
    y = top - step * lines - 0.012
    for r, row in enumerate(rows):
        first = pages.wrapped(row[0], size, first_width * meeting.SIZE[0] / pages.SIZE[0])
        height = step * (first.count("\n") + 1)
        for c, (x, cell) in enumerate(zip(columns, row)):
            fig.text(x, y, first if c == 0 else cell, fontsize=size, va="top", linespacing=1.05,
                     color=(colours[r] if colours and colours[r] and (c > 0 or colours[r] == pages.GREY) else pages.INK))
        y -= height + 0.008
    return y


def _checks_slide(title: str, loop: Path, k: int, foot: str, words: dict, first: bool = True):
    """Phase 3 of a turn: what the identification and the replicas say of the global parameters of the model, each with
    the band that the rule of admission sets, and whether the model is admitted. Its definition at the `first` turn alone
    (his message 311: fewer words a slide)."""
    from toy import meeting, pages
    ident = phase_file(loop, k, "identify")
    twin = phase_file(loop, k, "twin").get("admission", {})
    params = ident.get("parameters", {})
    bands = twin.get("bands", {})
    fig = _frame(title, foot, 3)
    y = fig.below - 0.02
    if first and len(params) <= 4:          # the definition again where the table leaves room for it: a model of seven global parameters has none, and the slide of the phases defines them
        y = _say(fig, y, words["checks"].format(sets=twin.get("sets", 30)), gap=0.02)
    if params:
        spread = bands.get("spread_over_bound", [float("nan")] * 2)
        cover = bands.get("coverage", [float("nan")] * 2)
        rows = [["required", "", f"< {bands.get('bias_over_spread_under', float('nan')):.2f}", f"{spread[0]:.2f} to {spread[1]:.2f}", f"{cover[0]:.2f} to {cover[1]:.2f}"]]
        colours = [pages.GREY]
        for name, cell in params.items():
            back = twin.get("parameters", {}).get(name, {})
            rows.append([name, pages_pm(float(cell["value"]), float(cell["standard_deviation"])), f"{abs(back.get('bias_over_spread', float('nan'))):.2f}",
                         f"{back.get('spread_over_bound', float('nan')):.2f}", f"{back.get('coverage', float('nan')):.2f}"])
            colours.append(None)
        y = _table(fig, y, ["global parameter", "fitted value", "|bias| /\nspread", "spread /\nuncertainty", "coverage"], rows, [0.05, 0.32, 0.57, 0.68, 0.835], 0.25, colours)
    else:
        fig.text(0.05, y, pages.wrapped(words["checks_none"], meeting.LEAST, 0.90 * meeting.SIZE[0] / pages.SIZE[0]), fontsize=meeting.LEAST, color=pages.INK, va="top", linespacing=1.2)
        y -= 0.12
    admitted = bool(twin.get("admitted"))
    verdict = (f"{twin.get('sets', 30)} replicas, fits converged on {100.0 * float(twin.get('ended', 1.0)):.0f} per cent. "
               + ("No value to return: the model has no global parameter." if not params and admitted else "Admitted: the fit returns the values put in." if admitted else "Not admitted: " + ", ".join(twin.get("failed", [])) + "."))
    _say(fig, max(0.16, min(y, 0.22)), verdict, colour=pages.INK if not params and admitted else meeting.RIGHT if admitted else meeting.WRONG, weight="bold")          # a check with nothing to check is not said to pass (review 39, minor 12)
    meeting.inside(fig, title)
    return fig


def _diagnosis_slide(title: str, report: dict, chi: float, foot: str, words: dict):
    """Phase 5 of a turn: the largest excesses of the diagnostic statistics over the detection threshold, beyond it in red,
    or the words that say the report was not made for this model."""
    import textwrap
    from toy import meeting, pages
    fig = _frame(title, foot, 5)
    if report and report.get("answers"):
        line = threshold_of(report)
        answers = sorted(report["answers"], key=lambda a: -a["excess"])[:5]
        y = _say(fig, fig.below - 0.02, words["diagnosis"].format(chi=f"{chi:.3f}" if chi < 2.0 else f"{chi:.3g}"))
        ax = fig.add_axes([0.50, 0.20, 0.40, y - 0.24])          # room on its left for the whole name of each bar
        for j, a in enumerate(answers):
            value = max(a["excess"] / line, 0.02)
            ax.barh(-j, value, height=0.62, color=meeting.WRONG if a.get("beyond_the_twin") else meeting.RIGHT, lw=0)
            ax.text(max(value, 1.0) * 1.3, -j, times_threshold(value), va="center", fontsize=meeting.LEAST, color=pages.INK)          # the bars as the sentences after them say them (review 48, minor 2)
        ax.axvline(1.0, color=pages.INK, lw=1.6, ls="--")
        ax.set_xscale("log")
        ax.set_xlim(0.05, 40.0 * max(1.0, max(a["excess"] / line for a in answers)))
        ax.set_yticks([-j for j in range(len(answers))])
        ax.set_yticklabels([pages.fill(short_name(f"{ASKED.get(a['question'], a['question'])}: {a['answer']}", 12), 44) for a in answers], fontsize=meeting.LEAST, linespacing=1.0)          # whole: cut at eight words, five bars read "a Gaussian width of 0.3 widths beside ..." (the red team of his message 328, moderate 6)
        ax.set_xlabel("excess / threshold", fontsize=meeting.LEAST)
        ax.tick_params(axis="y", length=0)
        meeting.plain(ax)
        ax.tick_params(labelsize=meeting.LEAST)
    else:
        _points(fig, [f"Reduced chi-squared {chi:.3f}.", words["diagnosis_none"]], fig.below - 0.04)
    meeting.inside(fig, title)
    return fig


def _families_slide(title: str, report: dict, foot: str, words: dict) -> list:
    """Phase 5 of a turn, the evidence of every family of statistics (the owner's messages 262 to 265): the excess over the
    threshold of each family's reading, the generalised moments of every kernel, the windows, the crossings, the
    derivatives, the ratios and the transformations of the value, beyond the threshold in red, those that read the most
    as many as the slide holds with their names apart, and on a slide of their own the candidates that the readings
    together point to, the report's matched answers, with their excess. Nothing where the report holds no reading of a
    family. Sixteen families and three shapes on one slide put their names on one another (the walk by hand of
    2026-10-20)."""
    from toy import meeting, pages
    fig = _frame(title, foot, 5)
    line = threshold_of(report)
    every = sorted([a for a in report["answers"] if a["question"] == "seen by"], key=lambda a: -a["excess"])
    matched = sorted([a for a in report["answers"] if a["question"] == "matched" and a.get("beyond_the_twin")], key=lambda a: -a["statistic"])[:3]
    y = fig.below - 0.02
    room = y - 0.11 - meeting.BASE          # under the sentence and its second one, of how many are shown
    seen = every[:max(3, int(room / FAMILY_ROW))]
    y = _say(fig, y, words["families"] + (f" The {len(seen)} that read the most, of {len(every)}." if len(seen) < len(every) else ""))
    ax = fig.add_axes([0.50, meeting.BASE, 0.42, y - meeting.BASE - 0.03])
    for j, a in enumerate(seen):
        value = max(a["excess"] / line, 0.02)
        ax.barh(-j, value, height=0.62, color=meeting.WRONG if a.get("beyond_the_twin") else meeting.RIGHT, lw=0)
    ax.axvline(1.0, color=pages.INK, lw=1.6, ls="--")
    ax.set_xscale("log")
    ax.set_xlim(0.05, 40.0 * max(1.0, max(a["excess"] / line for a in seen)))
    ax.set_ylim(-len(seen) + 0.4, 0.6)
    ax.set_yticks([-j for j in range(len(seen))])
    ax.set_yticklabels([pages.fill(short_name(a["answer"].split(":")[0] + ":" + a["answer"].split(":")[-1], 6), 30) for a in seen], fontsize=meeting.LEAST, linespacing=0.95)
    ax.set_xlabel("excess / threshold", fontsize=meeting.LEAST)
    ax.tick_params(axis="y", length=0)
    meeting.plain(ax)
    ax.tick_params(labelsize=meeting.LEAST)
    meeting.inside(fig, title)
    out = [fig]
    if matched:
        head = title + ", the shapes they point to"
        fig = _frame(head, foot, 5)
        _points(fig, [words["families_lead"]] + [f"{short_name(a['answer'], 16)}: {times_threshold(a['excess'] / line)} times the threshold." for a in matched], fig.below - 0.04, size=meeting.LEAST, gap=0.012)          # whole: two of them read alike when cut at nine words (the red team of his message 328, moderate 6)
        meeting.inside(fig, head)
        out.append(fig)
    return out


FAMILY_ROW = 0.10          # the height of the row of one family on its slide, as a share of the sheet: a name of two lines at the least letters and the gap to the next


def candidate_said(j: int, c: dict, line: float, short: bool = False) -> str:
    """One generic term that twinfit proposes for the model, as a sentence with its numbers, his message 250: "you didn't
    explain the proposal of 3 candidate terms for the model before choosing". What it stands for and how far beyond the
    threshold its deviation lies, its size fitted to the spectra with what it takes off the chi square, and the largest
    bias of a global parameter on the replicas without it and with it, but with `short`, where the table of the next slide
    holds the biases (his message 311: fewer words a slide)."""
    package = c.get("the_term_of_the_package") or {}
    said = f"{j}. {c['term']}, of the kind {c.get('stands_for') or 'that the residuals give'}, at {times_threshold(c.get('excess_in_the_report', float('nan')) / line)} times the threshold" + (", beyond the replicas." if c.get("beyond_the_twin") else ", within the replicas.")
    if package.get("stands") and package.get("bar"):
        said += f" Its one free coefficient fitted {abs(package['size']) / package['bar']:.1f} uncertainties from zero, chi-squared lower by {plain_number(package['gain_in_chi_square'])}."
    else:
        said += " No term of one free size in the package fits it."
    biases = [(abs(v["without_the_term"]["bias"]) / v["without_the_term"]["scatter"], abs(v["with_the_term"]["bias"]) / v["with_the_term"]["scatter"], name)
              for name, v in c.get("known_parameters", {}).items() if v.get("without_the_term", {}).get("scatter") and v.get("with_the_term", {}).get("scatter")]
    if biases and not short:
        without, with_it, name = max(biases)
        said += f" Bias of {name} on the replicas: {plain_number(without)} scatters without it, {plain_number(with_it)} with it."
    return said


def _proposals_slide(title: str, report: dict, candidates: list, foot: str, words: dict, part: int = 0):
    """Phase 6 of a turn, before the table that compares them: how twinfit makes its three proposals for the model, and
    each of them in a sentence by `candidate_said`, his message 250. With `part` 1 the first two, with 2 the third, each
    on a slide of its own (his message 311)."""
    from toy import meeting, pages
    if part:
        title = f"{title}, {part} of 2"
    fig = _frame(title, foot, 6)
    said = [candidate_said(j, c, threshold_of(report), short=True) for j, c in enumerate(candidates[:3], 1)] if candidates else [words["candidates_none"] if report else words["diagnosis_none"]]
    _points(fig, {0: said, 1: said[:2], 2: said[2:]}[part], fig.below - 0.03, width=0.90, size=meeting.LEAST, gap=0.012)
    meeting.inside(fig, title)
    return fig


def experiment_said(proposal: dict) -> str:
    """One proposal for the experiment in a sentence: its kind and lever, what it changes, and by how much the twin says
    that it moves the error of the known parameters, where the package measured it."""
    factor = proposal.get("error_over_that_under_the_records")
    moved = f" On the twin, the error of the known parameters times {factor:.2f}." if isinstance(factor, (int, float)) else ""
    return f"{str(proposal.get('his_kind', 'a change')).capitalize()}, {proposal.get('lever', '')}: {proposal.get('words', '')}.{moved}"


def _experiment_slide(title: str, table: dict, foot: str, words: dict):
    """The end of a turn, his message 310: the three proposals for the experiment that the report of the model gives, a
    better setup, a better characterisation of it or a better data collection, each measured on the twin."""
    from toy import meeting, pages
    fig = _frame(title, foot, 8)
    chosen = [p for p in (((table or {}).get("proposals") or {}).get("for_the_experiment") or [])[:3] if isinstance(p.get("error_over_that_under_the_records"), (int, float))]
    if not chosen:
        _points(fig, [words["experiment_none"]], fig.below - 0.03, width=0.90)
        meeting.inside(fig, title)
        return fig
    y = _say(fig, fig.below - 0.01, words["experiment_read"])          # a figure in the place of three sentences, his message 311
    ax = fig.add_axes([0.50, meeting.BASE, 0.42, y - meeting.BASE - 0.02])          # half the sheet for the labels, a lever of 33 letters at 21.5 points left a sheet at 0.42
    ratios = [float(p["error_over_that_under_the_records"]) for p in chosen]
    ax.barh(range(len(chosen)), ratios, color=[meeting.PALETTE[j] for j in range(len(chosen))])
    ax.axvline(1.0, color=pages.INK, lw=2)
    ax.set_yticks(range(len(chosen)))
    ax.set_yticklabels([f"{j}. {p.get('lever', '')}\n" + str(p.get("his_kind", "")).replace("a better ", "") for j, p in enumerate(chosen, 1)], fontsize=meeting.LEAST)          # two lines a proposal, its lever whole: three lines a proposal lay on one another under a sentence of two lines (the report of the blind run of 0.1.0rc16)          # the change on two lines at most, and its kind, setup, characterization or data collection, under it
    ax.invert_yaxis()
    for j, r in enumerate(ratios):          # each value right of its bar and of the line of no change, on its row
        ax.text(max(r, 1.0) + 0.04, j, f"× {r:.2f}", va="center", fontsize=meeting.LEAST, color=pages.INK)
    ax.set_xlim(0.0, max(1.45, max(ratios) + 0.45))
    ax.set_xlabel(words["experiment_axis"], fontsize=meeting.LEAST)
    meeting.plain(ax)
    meeting.inside(fig, title)
    return fig


def package_slides(where: Path, k: int) -> list:
    """The slides of the package's own visual log for model `k` of a step, its equation and what each of its terms does
    to the curves, his message 309, read from the loop folder `where` that the step kept: none where that loop wrote no
    equations.json, as a run of a version before 0.1.0rc8 did."""
    from twinfit import log as package_log
    import matplotlib.pyplot as plt
    if not (Path(where) / "equations.json").is_file():
        return []
    home = Path(where).parent          # steps/NN of a run keeps no description: the project's, from the run's own folder (the dry run of the report of message 310)
    run = package_log.read_run(home, next(iter(sorted(home.parent.parent.glob("*/description.toml"))), None))
    with plt.rc_context(package_log._rc()):
        figs = package_log.slides_equation(run, k) + [package_log.slide_effects(run, k)]          # every term of the equation, on as many slides as it asks
    for fig in figs:          # the package numbers its own eight phases from the records, and a turn of the analyst its eight from the model: the corner says whose (the red team of message 310)
        for x in fig.texts:
            if " phase " in x.get_text() and x.get_text().startswith("model "):
                x.set_text("twinfit's own visual log, " + x.get_text().replace("  ·  phase", ", its phase"))
    return figs


def _candidates_slide(title: str, report: dict, candidates: list, foot: str, words: dict, first: bool = True):
    """Phase 6 of a turn, his message 247: "how the new terms are chosen to be attempted". Each generic term that twinfit
    built for a detected deviation, with its excess over the threshold, its fitted size over its uncertainty, its gain in
    chi square, and the largest bias of a global parameter on the replicas without the term and with it, in scatters. Its
    explanation at the `first` turn, and at the others its key alone (his message 311: fewer words a slide)."""
    from toy import meeting, pages
    fig = _frame(title, foot, 6)
    intro = pages.wrapped(words["candidates"] if first else words["candidates_key"], meeting.LEAST, 0.90 * meeting.SIZE[0] / pages.SIZE[0])
    fig.text(0.05, fig.below - 0.02, intro, fontsize=meeting.LEAST, color=pages.INK, va="top", linespacing=1.2)
    y = fig.below - 0.04 - 0.05 * (intro.count("\n") + 1)
    if candidates:
        line = threshold_of(report)
        rows, colours = [], []
        for c in candidates:
            package = c.get("the_term_of_the_package") or {}
            biases = [(abs(v["without_the_term"]["bias"]) / v["without_the_term"]["scatter"], abs(v["with_the_term"]["bias"]) / v["with_the_term"]["scatter"])
                      for v in c.get("known_parameters", {}).values() if v.get("without_the_term", {}).get("scatter") and v.get("with_the_term", {}).get("scatter")]
            worst = max(biases) if biases else None
            rows.append([c["term"], times_threshold(c.get("excess_in_the_report", float("nan")) / line),
                         plain_number(package["gain_in_chi_square"]) if package.get("stands") and package.get("gain_in_chi_square") is not None else "none", "none" if worst is None else f"{plain_number(worst[0])} to {plain_number(worst[1])}"])          # the gain of the term that is fitted, or none where no term of one size is, as the proposal says (review 40, moderate 6)
            colours.append(meeting.WRONG if c.get("beyond_the_twin") else None)
        y = _table(fig, y, ["generic term", "excess /\nthreshold", "$\\chi^2$\ngain", "bias of a\nparameter"], rows, [0.05, 0.50, 0.65, 0.81], 0.43, colours)          # chi-squared written as LaTeX writes it (his message 314)
        tail = (words["candidates_foot"] if any(row[3] != "none" for row in rows) else words["candidates_foot_unread"]) if first else ""          # "none, none, none" said where the model has no global parameter (review 44, minor 6)
    else:
        tail = words["candidates_none"] if report else words["diagnosis_none"]
    if tail:
        _say(fig, max(0.18, min(y, 0.24)), tail, colour=pages.GREY)
    meeting.inside(fig, title)
    return fig


def fitted_quote(quote: str, room: float, size: float) -> str:
    """A quotation cut at a word to the lines that `room`, a height of the slide, holds at `size`, with three dots where it
    was cut."""
    from toy import meeting, pages
    lines = max(1, int(room / (0.047 * size / 21.5)))
    cut = quote.split()
    for n in range(min(len(cut), QUOTE), 0, -1):
        text = "\u201c" + " ".join(cut[:n]) + (" ..." if n < len(cut) else "") + "\u201d"
        if pages.wrapped(text, size, 0.90 * meeting.SIZE[0] / pages.SIZE[0]).count("\n") + 1 <= lines:
            return text
    return "\u201c...\u201d"


CLOSE = 0.8          # two settings whose excesses lie within this share of each other are ranked almost alike by the data


def close_call(report: dict, problem=None) -> str:
    """Where the data rank two settings of the design almost alike, the sentence that says so, with how closely the two
    settings go together across the records: the choice between them rests on the physics of the term, his message 257.
    Empty where one setting leads clearly, or where the report names one alone."""
    import numpy as np
    keys = [] if problem is None else list(problem.records[0].design)
    named = settings_named(report, keys) if problem is not None else [(a["answer"], float(a["excess"])) for a in sorted(report.get("answers", []), key=lambda a: -a["excess"]) if a["question"] == "along" and a.get("beyond_the_twin")]
    if len(named) < 2 or named[1][1] < CLOSE * named[0][1] or named[0][1] <= 0:          # both past the threshold (review 42, minor 8)
        return ""
    (a, ea), (b, eb) = named[:2]
    line = threshold_of(report)
    said = f"The data rank {a} and {b} almost alike, {times_threshold(ea / line)} and {times_threshold(eb / line)} times the threshold"
    if problem is not None and all(a in rec.design and b in rec.design for rec in problem.records):
        r = float(np.corrcoef([float(rec.design[a]) for rec in problem.records], [float(rec.design[b]) for rec in problem.records])[0, 1])
        said += f", and across the spectra the two go together with a correlation of {r:.2f}"
    return said + ": the choice between them rests on the physics of the term."


def report_asked(report: dict, candidates: list) -> bool | None:
    """Whether the report of a model asked for a term: one of its proposals stands beyond the replicas. None where the
    report was not computed. Turn 3 of BLIND4 showed its deviations and its term as the algorithm's choice where every
    proposal stood within the replicas, and turn 3 of BLIND3 where no report was computed (review 44, moderate 2)."""
    if not report:
        return None
    return any(c.get("beyond_the_twin") for c in candidates)


def page_of(run: Path, wrote: str, on: int | None = None, lines: bool = False, pattern: str = CHOICE) -> str:
    """The analyst's own words on the step that wrote a term, the section of its page under which the brief asks why: the
    page whose first line names the term. The reason in a term's file was written while a report was read, and on BLIND4
    turn 3 the page corrects it, "The report of model 2 asks for no term" (review 45, minor 3). Else, with `on`, the page
    of the step of that model, where the brief asks the choice of the next term: the run of his message 310 named its
    terms by their files. Nothing where no page names the term."""
    steps = steps_of(run)
    for _, _, folder in steps:
        page = Path(folder) / "step.md"
        if page.exists() and wrote and wrote.lower() in page.read_text().split("\n", 1)[0].lower():
            return section_of(page, pattern, lines)
    own = next((Path(folder) / "step.md" for name, _, folder in steps if on is not None and name == f"{on:02d}"), None)
    return section_of(own, pattern, lines) if own is not None and own.exists() else ""


def run_said(run: Path) -> list:
    """The version of the package that the run used, the model and effort of its agent, and its tokens, from its own copy
    of the package, `run.json` and `metrics.json`: nothing where the run kept none."""
    import re
    run = Path(run)
    version = re.search(r'^__version__ = "([^"]+)"', (run / "twinfit" / "version.py").read_text(), re.M) if (run / "twinfit" / "version.py").exists() else None
    ran = json.loads((run / "run.json").read_text()) if (run / "run.json").exists() else {}
    m = json.loads((run / "metrics.json").read_text()) if (run / "metrics.json").exists() else {}
    tokens = sum(int(m.get(k) or 0) for k in ("tokens_in", "tokens_out", "tokens_cache_read", "tokens_cache_made"))
    said = []
    if version or ran.get("model"):
        said.append("Run with twinfit " + (version.group(1) if version else "of an unknown version") + (f", by {ran['model']} at {ran.get('effort', 'its default')} effort" if ran.get("model") else "") + ".")
    if tokens:
        said.append(f"Its tokens: {tokens / 1e6:.2f} million, of which {100 * int(m.get('tokens_cache_read') or 0) / max(tokens, 1):.0f} per cent read again from the cache, and {int(m.get('tokens_out') or 0) / 1e3:.0f} thousand written.")
    return said


def experiment_words(table: dict) -> list:
    """The three proposals for the experiment of a model's report, each in the words of the package's log: the change,
    what it is, and what the twin measured of it. They are for the experimenter, and a blind run follows none."""
    from twinfit.log import experiment_line
    chosen = (((table or {}).get("proposals") or {}).get("for_the_experiment") or [])[:3]
    return [f"{j}. " + experiment_line(p).replace("_", " ") + "." for j, p in enumerate(chosen, 1) if isinstance(p, dict) and p.get("lever")]          # a name of the code in words


def chosen_among(named: list, wrote: str) -> int | None:
    """Which of the proposals the term written is, counted from 1, by its name, or None where it is none of them."""
    w = " ".join(str(wrote or "").lower().split())
    if len(w.split()) < 2:          # an empty or one-word name is no term: it read as proposal 1 (the red team of his message 333)
        return None
    for i, n in enumerate(named, 1):
        m = " ".join(str(n).lower().split())
        if m == w or (len(m.split()) >= 3 and (m in w or w in m)):
            return i
    return None


def answer_in_its_words(run: Path) -> list:
    """The summary of the analyst's answer, `answer.json` of the run, in its sentences: its terms, their widths and what
    it says that the traces could not tell apart. The closing slides gave the procedure's own reading and left the
    analyst's caveats out (the red team of his message 328)."""
    import re
    path = Path(run) / "answer.json"
    if not path.exists():
        return []
    summary = json.loads(path.read_text()).get("summary") or ""
    return reasoning_points(summary)


def red_team_of(run: Path, n: int, lines: bool = True) -> str:
    """What the red team of the analyst said of term file `n`, the section of the page of a step whose heading names it,
    "Red team of term 01": the slides of the run cited it and did not show it (the red team of his message 328)."""
    for _, _, folder in steps_of(run):
        said = section_of(Path(folder) / "step.md", rf"\bred team of term {n:02d}\b", lines)
        if said:
            return said
    return ""


def clauses(sentence: str, most: int = 60) -> list:
    """A sentence of the analyst in parts of at most `most` words, cut after a comma or a semicolon where it runs longer:
    one sentence of 86 words stood alone on a slide of the run of message 310."""
    words = sentence.split()
    if len(words) <= most:
        return [sentence]
    out, part = [], []
    for w in words:
        part.append(w)
        if len(part) >= most // 2 and w.endswith((",", ";")) or len(part) >= most:
            out.append(" ".join(part))
            part = []
    if part:
        out.append(" ".join(part))
    return out


REASONING = None          # the words of the analyst's reasoning that its slides carry: all of them, on as many slides as they take. At 180, three slides, a derivation was cut in the middle of its sentence, and his message 319 asks the report "really full, in particular with the clear reasoning"


def reasoning_points(text: str, cap: int | None = REASONING) -> list:
    """The analyst's reasoning of a step as points of slides, his message 311: "the clear reasoning and clear decision
    making process". Each paragraph a point, each row of a table a point of its cells, the marks of the page taken off,
    every word of it, or up to `cap` words with the last point cut at a word with three dots where the page goes on."""
    import re
    points = []
    runs = []
    for block in re.split(r"\n\s*\n", text or ""):          # a block in runs of table lines and of text lines: a table written right under a paragraph is a table all the same (the run of message 310)
        for x in (x.strip() for x in block.splitlines() if x.strip()):
            if runs and runs[-1][0] == x.startswith("|") and runs[-1][2] is block:
                runs[-1][1].append(x)
            else:
                runs.append((x.startswith("|"), [x], block))
    for table, lines, _ in runs:
        if table:
            rows = [[c.strip() for c in x.strip("|").split("|")] for x in lines if not re.fullmatch(r"\|[\s|:-]+\|?", x)]
            if rows and len(rows[0]) > 2:          # each value with the name of its column, where a row of six numbers read as numbers alone (the red team of his message 328)
                points += [", ".join(f"{h} {c}" if h else c for h, c in zip(rows[0], row)) for row in rows[1:]]
            else:
                points += [": ".join(row) for row in rows[1:]]          # the head of a table names its columns, which each row then says
        elif lines:
            group = []
            for sentence in (part for whole in re.split(r"(?<=[.!?])\s+", " ".join(lines)) for part in (clauses(whole) if len(whole.split()) > 90 else [whole])):          # a long paragraph in points of whole sentences, of 60 words at most where its sentences allow, a longer sentence whole on its own up to 90 words, more than a slide holds only then cut at its commas: cut at 60, sentences ran across two slides (the red team of his message 333, severe 2)
                if group and len(" ".join(group + [sentence]).split()) > 60:
                    points.append(" ".join(group))
                    group = []
                group.append(sentence)
            points.append(" ".join(group))
    joined = []
    for point in points:          # a point that opens a list, on a colon, keeps its items with it on one slide: the list began on the next slide (the red team of his message 333)
        if joined and joined[-1].rstrip().endswith(":") or joined and joined[-1].count("\n") and not joined[-1].rstrip().endswith((".", "!", "?")) and len(point.split()) <= 15:
            joined[-1] = joined[-1] + "\n" + point
        else:
            joined.append(point)
    points = joined
    out, left = [], float("inf") if cap is None else cap
    for point in (re.sub(r"[*`#]", "", x).strip() for x in points):
        words = point.split()
        if not words or left <= 0:
            break
        out.append(point if len(words) <= left else " ".join(words[:left]) + " ...")
        left -= len(words)
    return out


def quote_room(y: float, gloss_lines: int) -> float:
    """The height that a quote of the choice may take below `y`: down to 0.10 of the slide, and higher by the lines of a
    gloss below it with their gap. The gloss of BLIND4 turn 3 ran onto the foot under a quote that took the room."""
    return y - 0.055 - 0.10 - (0.047 * gloss_lines + 0.01 if gloss_lines else 0.0)


def _choice_slide(title: str, d: dict, report: dict, candidates: list, foot: str, words: dict, problem=None, page: str = "", part: int = 0):
    """Phase 7 of a turn: what the report detected and what twinfit proposed, where the data leave the choice to the
    physics, the term that the analyst wrote, and its physical rationale in the analyst's own words. With `part` 1 what
    the report detected alone, with 2 the term and its rationale, each on a slide of its own (his message 311)."""
    from toy import meeting, pages
    if part:
        title = f"{title}, {part} of 2"
    fig = _frame(title, foot, 7)
    from twinfit.loop import ACTIONABLE
    found = detected(report, 3, ACTIONABLE)          # the deviations that name a term, as the record names them (review 45, moderate 1)
    lines = [("Deviations detected by the report: " + ", and ".join(found) + ".") if found else (words["choice_unread"] if not report else "No statistic exceeds the detection threshold.")]
    asked = report_asked(report, candidates)
    if asked is not True:          # a term that the algorithm does not ask for is said to be the analyst's own, first (review 44, moderate 2)
        lines.insert(0, words["choice_not_computed"] if asked is None else words["choice_not_asked"])
    if close_call(report, problem):
        lines.append(close_call(report, problem))          # the generic terms stand on the two slides of phase 6 before it, and the room goes to the physics of the choice (his message 257)
    y = _points(fig, lines, fig.below - 0.03, width=0.90, size=meeting.LEAST, gap=0.012) if part != 2 else fig.below - 0.03
    if part == 1:          # what the report detected alone: the term and its rationale on the next slide
        meeting.inside(fig, title)
        return fig
    wrote = pages.wrapped(f"Term written by the analyst: {d['wrote']}.", meeting.LEAST, 0.90 * meeting.SIZE[0] / pages.SIZE[0], "bold")
    fig.text(0.05, y - 0.005, wrote, fontsize=meeting.LEAST, fontweight="bold", color=pages.INK, va="top", linespacing=1.15)
    y -= 0.05 * (wrote.count("\n") + 1) + 0.03
    fig.text(0.05, y, words["choice_page"] if page else "Physical rationale, verbatim:", fontsize=meeting.LEAST, color=pages.GREY, va="top")
    whole = page or d["reason"] or "No rationale was recorded."
    gloss = pages.wrapped(words["choice_gloss"], meeting.LEAST, 0.90 * meeting.SIZE[0] / pages.SIZE[0]) if "asks for no term" in whole and found else ""
    quote = fitted_quote(whole, quote_room(y, gloss.count("\n") + 1 if gloss else 0), meeting.LEAST)          # where the report asked for no term, the analyst's page of the step, which says why it went on, with room for the gloss below it
    said = pages.wrapped(quote, meeting.LEAST, 0.90 * meeting.SIZE[0] / pages.SIZE[0])
    fig.text(0.05, y - 0.055, said, fontsize=meeting.LEAST, color=pages.INK, va="top", linespacing=1.15, style="italic")
    if gloss:          # the analyst's reading set beside the rule of the package, which reads a term asked by the answers past the threshold above (review 48, minor 6)
        _say(fig, y - 0.055 - 0.047 * (said.count("\n") + 1) - 0.01, words["choice_gloss"], colour=pages.GREY)
    meeting.inside(fig, title)
    return fig


def bic_said(d: dict) -> tuple:
    """What the Bayesian information criterion said of a term, as a sentence, and whether it met the rule: a fall of more
    than 2. None where the run did not compute it."""
    import re
    if d["fall"] is not None:
        return f"Bayesian information criterion: decreases by {d['fall']:.1f}, where a decrease of more than 2 is required.", d["fall"] > 2.0
    found = re.search(r"Bayesian criterion (rises|falls) by (?:only )?(-?[\d.]+)", d.get("refused") or "")
    if found:
        return f"Bayesian information criterion: {'increases' if found.group(1) == 'rises' else 'decreases'} by {float(found.group(2)):.1f}, where a decrease of more than 2 is required.", False
    return None, None


def facts_of(d: dict, loop: Path, k: int, after: tuple | None, asked: bool | None = True) -> dict:
    """What decided a term, read once for the slide of its decision and for the record of its turn: the validation of the
    term on replicas that contain it, the criterion and whether it met its rule, whether the replicas admit the enlarged
    model, and the verdict with its reason and its colour. `after` is the loop and the model of the enlarged model."""
    from toy import meeting, pages
    book = Path(loop) / "notebook"
    cand = next((c for c in reversed([json.loads(p.read_text()) for p in sorted(book.glob(f"*_model_{k}_candidate_on_the_twin.json"))]) if c.get("the_planted_candidate") == d["wrote"]), None) if book.is_dir() else None
    bic, met = bic_said(d)
    if bic is None and d["chain"] is None and after is not None and d.get("loop_before"):
        before = (Path(d["loop_before"]), last_model(d["loop_before"]))
        same = phase_file(*before, "intake").get("samples") == phase_file(*after, "intake").get("samples")
        a, b = (json.loads((Path(p) / f"model_{m}.json").read_text()).get("criterion") for p, m in (before, after))
        if same and a is not None and b is not None:
            bic, met = f"Bayesian information criterion: {'decreases' if a > b else 'increases'} by {abs(a - b):.1f}, on the same samples.", a - b > 2.0
        elif not same:
            bic = "Bayesian information criterion: not comparable, since the data of the two models differ."
    adm = phase_file(*after, "twin").get("admission", {}) if after is not None else {}
    if d["chain"] is None and asked is False and not met:          # EXAMPLE1 adopted three turns on evidence that shows nothing, with a criterion that did not fall (review 44, moderate 3)
        verdict, why, colour = "Not a term", "a test of the model. No proposal of the report stood beyond the replicas at this turn, and the criterion does not fall by more than 2.", pages.GREY
    elif d["chain"] is None and asked is False:
        verdict, why, colour = "Adopted", "by the analyst, by its own protocol. No proposal of the report stood beyond the replicas at this turn, and the criterion falls by more than 2.", meeting.RIGHT
    elif d["chain"] is None and asked is None:
        verdict, why, colour = "Adopted", "by the analyst, with no report computed on the model before it.", meeting.RIGHT
    elif d["chain"] is None:
        verdict, why, colour = "Adopted", "by the analyst, on the evidence above.", meeting.RIGHT
    elif d["kept"]:
        verdict, why, colour = "Accepted", "the criterion and the replicas both support the new term.", meeting.RIGHT
    else:
        verdict, why, colour = "Rejected", "the criterion does not support the new term." if met is False else str(d["refused"] or "") + ".", meeting.WRONG
    return {"cand": cand, "bic": bic, "met": met, "admitted": bool(adm.get("admitted")) if adm else None, "verdict": verdict, "why": why, "colour": colour}


def settings_named(report: dict, keys: list) -> list:
    """The settings of the design that the report names past its threshold, through the questions "along" and "changes
    along", each with its largest excess, the largest first."""
    best = {}
    for a in report.get("answers", []):
        key = a["answer"].split(" along ")[-1] if a["question"] == "changes along" else a["answer"] if a["question"] == "along" else None
        if key in keys and a.get("beyond_the_twin"):
            best[key] = max(best.get(key, 0.0), float(a["excess"]))
    return sorted(best.items(), key=lambda kv: -kv[1])


def leading_covariate(report: dict, problem, result) -> str | None:
    """The setting of the design along which the deviation is read, his message 257 on the steps of a decision and the
    analysis that chooses the models: the one that the report names past its threshold with the largest excess, through
    the questions "changes along" and "along". None where it names none, and no evidence along a setting is then drawn."""
    from twinfit.candidates import covariates_of
    named = settings_named(report, covariates_of(problem))
    return named[0][0] if named else None          # no setting that the report does not name past its threshold (review 42, moderate 2)


def trends_of(problem, result, key: str) -> tuple:
    """The windowed sums of each record in each band, against the covariate `key`, place by place, as the report reads a
    change along a setting: for each band the sum over the records of the covariate, centred and of unit length, times
    the reading, its trend, a standard normal number where nothing changes along it. Returns the values of the
    covariate and, for each place, its bands as (from, to, readings, trend)."""
    import numpy as np
    from twinfit.report import by_position_and_record
    from twinfit.where import PLACES
    z = np.array([float(rec.design[key]) for rec in problem.records])
    unit = (z - z.mean()) / (np.sqrt((z - z.mean()) @ (z - z.mean())) or 1.0)
    rows = by_position_and_record(problem, result)
    places = {place: [(r["from"], r["to"], np.asarray(r["both_sides"], dtype=float), float(np.nansum(unit * np.asarray(r["both_sides"], dtype=float))))
                      for r in rows if pick(r)] for place, pick in PLACES}
    return z, places


def largest_trend(places: dict) -> tuple:
    """The band whose trend is the largest in size, as (its place, its inner edge, its outer edge, its trend)."""
    return max(((place, lo, hi, t) for place, bands in places.items() for lo, hi, _, t in bands), key=lambda b: abs(b[3]))


def trend_at(places: dict, place: str, lo: float) -> float | None:
    """The trend of the band of `place` whose inner edge is `lo`, or None where the places hold no such band."""
    return next((t for b_lo, _, _, t in places.get(place, []) if b_lo == lo), None)


def straight_fit(z, reading, at=None):
    """The least squares straight line of the readings that are numbers against the setting, at each value of `at`, the
    values of the setting where it is None, or None where fewer than two readings are numbers: a reading is no number
    where the fit takes a record's share of a band whole (review 43, minor 5)."""
    import numpy as np
    z, reading = np.asarray(z, dtype=float), np.asarray(reading, dtype=float)
    at = z if at is None else np.asarray(at, dtype=float)
    ok = np.isfinite(reading)
    if ok.sum() < 2:
        return None
    zc = z[ok] - z[ok].mean()
    slope = float(np.sum(zc * (reading[ok] - reading[ok].mean())) / (np.sum(zc * zc) or 1.0))
    return reading[ok].mean() + slope * (at - z[ok].mean())


def value_scale(ax, reach: float, size: float = 17.0) -> None:
    """The scale of the windowed sums of a panel of evidence, from -reach to reach: linear, or where the sums reach past
    20, linear within the band of 2 and logarithmic beyond it, so the band is seen beside sums of hundreds (review 43,
    minor 4). The marks of the value stand clear of the corner, where the marks of the setting are."""
    import numpy as np
    ax.set_ylim(-reach, reach)
    if reach > 20.0:
        ax.set_yscale("symlog", linthresh=2.0, linscale=0.6)          # the band of 2 seen, and a decade tall enough for its marks
        marks = [10.0 ** e for e in range(1, int(np.log10(0.75 * reach)) + 1)]          # a mark at each power of ten, so that a sum of 30 or 60 has marks to be read by (review 44, moderate 4)
        decade = ax.get_position().height * ax.figure.get_size_inches()[1] * 72.0 / 2.0 / (0.6 / 0.9 + np.log10(reach / 2.0))          # the height of a decade in points, on this panel
        if decade < 1.25 * size:          # a decade too short for the words of its mark keeps the largest mark alone: a caption that grew took two marks onto one another
            marks = marks[-1:]
        ax.set_yticks([-m for m in reversed(marks)] + [0.0] + marks)
        ax.yaxis.set_major_formatter(__import__("matplotlib").ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    else:
        mark = next(m * 10.0 ** e for e in range(int(np.floor(np.log10(0.75 * reach))), -2, -1) for m in (5.0, 2.0, 1.0) if m * 10.0 ** e <= 0.75 * reach)
        ax.set_yticks([-mark, 0.0, mark])


def edges_said(places: dict) -> tuple:
    """The edges of the inner and of the outer band of each place, in widths, as words: "0 to 0.25, 0.5 to 1 and 2 to 4"
    and "0.25 to 0.5, 1 to 2 and 4 to 8". A title gave the span of a place alone (review 44, minor 2)."""
    def said(j: int) -> str:
        parts = [f"{bands[j][0]:g} to {bands[j][1]:g}" for bands in places.values() if len(bands) > j]
        return ", ".join(parts[:-1]) + (" and " if len(parts) > 1 else "") + parts[-1] if parts else ""
    return said(0), said(-1)


def ranking_said(report: dict, keys: list) -> str:
    """The settings that the report names past its threshold, the largest first, each with its excess over the
    threshold: on turn 1 of BLIND4 a slide drew the deviation along the temperature with no word of why that setting
    (review 44, minor 4). Nothing where it names none."""
    named = settings_named(report, keys)
    if not named:
        return ""
    line = threshold_of(report)
    parts = [f"{key} {times_threshold(excess / line)}" for key, excess in named]
    return "Past the threshold: " + (", ".join(parts[:-1]) + " and " if len(parts) > 1 else "") + parts[-1] + " times it."


def scale_before(scales: dict, seen: tuple, places: dict) -> float:
    """The half height of the slide of evidence drawn for a model, where one was drawn, at its own phase 5 or as the slide
    after the term of the turn before, or else that of its own sums: the slide after the next term is drawn on it. EXAMPLE1
    slide 37 took the reach of a slide of phase 5 that was never drawn (review 47, moderate 1)."""
    return scales.get(seen, reach_of(places))


def reach_of(places: dict) -> float:
    """The half height of the panels of a slide of evidence: a little beyond the largest windowed sum, and at least 4."""
    import numpy as np
    return max(4.0, 1.15 * max(float(np.nanmax(np.abs(r))) for bands in places.values() for _, _, r, _ in bands))


def _evidence_slide(title: str, problem, result, key: str, foot: str, text: str, phase: int, scale: str = "", ranked: str = "", at_least: float = 0.0):
    """The evidence of a deviation along a setting of the design, his message 257: for the core, the flanks and the wings,
    the windowed sum of each spectrum in each band against the setting, the grey band of plus or minus 2 that a complete
    model leaves them in, and the trend of each band as a line and as a number above its panel."""
    import numpy as np
    from toy import meeting, pages
    z, places = trends_of(problem, result, key)
    reach = max(reach_of(places), at_least)          # the slide after a term on the scale of the slide before it, so that the two are read together (review 44, moderate 4)
    inner, outer = edges_said(places)
    text = text.format(key=key, inner=inner, outer=outer) + (" " + scale if reach > 20.0 and scale else "") + (" " + ranked if ranked else "")          # the scale said where it is not plain (review 44, moderate 4)
    fig = _frame(title, foot, phase)
    top = _say(fig, fig.below - 0.02, text) - 0.07          # room for the one line of each panel's title, and no more
    logx = bool(np.all(z > 0) and z.max() / z.min() > 10.0)
    along = np.geomspace(z.min(), z.max(), 40) if logx else np.linspace(z.min(), z.max(), 40)          # each line drawn at forty points, so that it stays straight on any scale (review 44, moderate 4)
    for j, (place, bands) in enumerate(places.items()):
        ax = fig.add_axes([0.09 + j * 0.31, 0.19, 0.25, top - 0.19])          # panels as tall as the foot allows: a decade must hold the words of its mark (review 46, moderate 3)
        ax.axhspan(-2.0, 2.0, color=pages.GREY, alpha=0.18, lw=0)
        ax.axhline(0.0, color=pages.INK, lw=0.8)
        for b, (_, _, reading, t) in enumerate(bands):
            colour = meeting.PALETTE[b % len(meeting.PALETTE)]
            ax.plot(z, reading, "os"[b % 2], ms=10, color=colour, ls="none")
            line = straight_fit(z, reading, along)
            if line is not None:          # the straight fit of the band to the setting (review 42, moderate 3)
                ax.plot(along, line, "-", lw=2.4, color=colour)
        if logx:
            ax.set_xscale("log")
        value_scale(ax, reach, meeting.LEAST)
        ax.set_title(f"{place}: " + " and ".join(f"{t:+.1f}" for *_, t in bands), fontsize=meeting.LEAST, color=pages.INK, loc="left", pad=8)          # one line, the place and its trends, inner first: the edges stand in the caption, and a second line took the room of the marks at 10 (review 49, moderates 1 and 2)
        ax.set_xlabel(key, fontsize=meeting.LEAST)
        if j == 0:
            ax.set_ylabel("windowed sum", fontsize=meeting.LEAST)
        meeting.plain(ax)
        ax.tick_params(labelsize=meeting.LEAST)
        if logx:
            ax.xaxis.set_minor_formatter(__import__("matplotlib").ticker.NullFormatter())
    meeting.inside(fig, title)
    return fig


def deviation_after(report: dict, later: dict, k: int, names: tuple | None = None) -> str:
    """What a term did to the deviation that leads the report it was written for, where that deviation is a place or a
    shape and no setting: its excess over the threshold under the model and under the model that adds the term, each
    from its own report. The record gave the trend of one band along one setting where the deviation it named was the
    pedestal, at 6000 times the threshold and 0.73 after (review 44, minor 5). Nothing where the leading deviation is
    along a setting, which the trend says, or where the later report does not read it."""
    from twinfit.loop import ACTIONABLE
    lead = max((a for a in report.get("answers", []) if a.get("beyond_the_twin") and a["question"] in ACTIONABLE), key=lambda a: a["excess"], default=None)
    if lead is None or lead["question"] in ("along", "changes along") or not later:
        return ""
    same = next((a for a in later.get("answers", []) if (a["question"], a["answer"]) == (lead["question"], lead["answer"])), None)
    if same is None:
        return ""
    before_name, after_name = names or (f"model {k}", f"model {k + 1}")
    ratio = same["excess"] / threshold_of(later)
    return (f"The leading deviation, {SAID.get(lead['question'], lead['question'])} {lead['answer']}, stands at {times_threshold(lead['excess'] / threshold_of(report))} times the threshold "
            f"under {before_name} and " + (f"{times_threshold(ratio)} under {after_name}." if ratio > 0 else f"shows no excess over the replicas under {after_name}."))          # an excess below none said as such, not as a negative number of thresholds


def trend_said(key: str, before: tuple, after: float | None, k: int, largest_after: tuple | None = None, names: tuple | None = None) -> str:
    """The trend of the band that leads along a setting, under a model and under the model that adds the term, as a
    sentence: what the term did to the evidence that it was written for. `before` is a band of `largest_trend`, and
    `largest_after` the band of `largest_trend` after the term, said where it is larger than the band quoted (review 42,
    minor 7). `names` names the two models where they are no models of one chain (minor 9)."""
    place, lo, hi, t = before
    first, second = names or (f"model {k}", f"model {k + 1}")
    said = f"Along {key}, the band from {lo:g} to {hi:g} widths has a trend of {t:+.1f} under {first}"
    if after is None:
        return said + "."
    said += f", {after:+.1f} under {second}"
    if largest_after is not None and abs(largest_after[3]) > abs(after) + 0.05 and (largest_after[0], largest_after[1]) != (place, lo):
        said += f", where the largest trend is {largest_after[3]:+.1f}, from {largest_after[1]:g} to {largest_after[2]:g} widths in {largest_after[0]}"
    return said + "."


FOOT_CLEAR = 0.13          # the height of a slide, from its bottom, that its words leave to the foot and the line above it
FOOT_TOP = 0.085          # the height of a slide, from its bottom, that its foot and the phases of a turn take


def record_gap(heights: list, top: float) -> float:
    """The room between the rows of a record of a decision: 0.014 of the slide's height, or less where the rows would
    then come within `FOOT_CLEAR` of the bottom, never under none. A record of six lines of test ended a line's height
    above the foot of the slide (review 43, minor 6)."""
    room = top - FOOT_CLEAR - sum(heights)
    return max(0.0, min(0.014, room / max(len(heights) - 1, 1)))


def quoted_within(words: list, room: float, size: float, share: float) -> str:
    """The words of a quote between quotation marks, wrapped to `share` of the slide's width at `size` points, at most
    QUOTE of them and fewer where their lines would take more than `room` of the slide's height, with ... where they are
    cut. The record of the run added to the slide of the termination took the reason of BLIND4 onto its foot."""
    from toy import pages
    most = max(int(room / (0.058 * size / 24.0 * 1.2 / 1.25)), 1)
    n = min(len(words), QUOTE)
    while True:
        text = pages.wrapped("\u201c" + " ".join(words[:n]) + (" ..." if n < len(words) else "") + "\u201d", size, share)
        if text.count("\n") + 1 <= most or n <= 1:
            return text
        n -= 1


def deviations_said(report: dict, words: dict) -> str:
    """The deviations of a report as the record of a decision says them: the two largest past the threshold, and how many
    more the slide of the diagnosis shows, so that two are not read as the whole list. On BLIND4 the record of turn 2
    named the density and the temperature, and its diagnosis showed the wings past the threshold as well."""
    from twinfit.loop import ACTIONABLE
    found = detected(report, 2, ACTIONABLE)          # the answers that name a term, as the package counts them: a family of statistics seen past its line names none
    more = len(detected(report, None, ACTIONABLE)) - len(found)
    if found:
        return "The report detects a deviation " + ", and ".join(found) + (f", and {said_number(more)} more {'answers' if more > 1 else 'answer'} of the report past the threshold" if more else "") + "."          # the diagnosis draws five statistics, and not every answer (review 45, moderate 1)
    return words["choice_unread"] if not report else "The report detects no deviation above the threshold."


SPLIT_ABOVE = 60          # the words of sentences above which a slide of a turn goes on two halves, each on a slide of its own (his message 311), 60 as his message 149 asks of a page, where 80 left 31 of 140 slides above it (review 87, moderate 2)


def in_halves(make, *args, **kwargs) -> list:
    """A slide of a turn as `make` draws it, whole where its sentences hold `SPLIT_ABOVE` words at most, else its two
    halves, `make(..., part=1)` and `make(..., part=2)`, each on a slide of its own: every sentence kept, fewer on each
    slide (his message 311: "split them on multiple slides")."""
    import matplotlib.pyplot as plt
    whole = make(*args, **kwargs)
    if slide_words([whole])["sentences"][0] <= SPLIT_ABOVE:
        return [whole]
    plt.close(whole)
    return [make(*args, part=1, **kwargs), make(*args, part=2, **kwargs)]


def _record_slide(title: str, d: dict, k: int, turn: int, report: dict, candidates: list, chi: float, facts: dict, foot: str, words: dict, trend: str = "", name: str | None = None, part: int = 0):
    """The decision of a turn in one view, his message 248: "The log has to tell a (clear, and data-driven and physics-driven
    decision making!!) story". The evidence that the data give, the hypothesis that the physics gives, the test that the
    statistics give, and the verdict with the model that the next turn starts from. With `part` 1 the evidence and the
    hypothesis alone, with 2 the test and the verdict, each on a slide of its own where the whole holds more words than
    one slide should (his message 311)."""
    from toy import meeting, pages
    if part:
        title = f"{title}, {part} of 2"
    fig = _frame(title, foot)
    evidence = [f"Reduced chi-squared of {name or f'model {k}'}: {chi:.3f}.", deviations_said(report, words)]          # a step's model by its step: every loop of EXAMPLE1 starts at a model 0 (review 44)
    lead = next((c for c in candidates if c.get("beyond_the_twin")), None)          # a proposal within the replicas is no deviation of the data, and its bias is not quoted (review 44, moderate 2)
    if report and lead is None:
        evidence.append(words["record_not_asked"])
    biases = [(abs(v["without_the_term"]["bias"]) / v["without_the_term"]["scatter"], abs(v["with_the_term"]["bias"]) / v["with_the_term"]["scatter"], name)
              for name, v in (lead or {}).get("known_parameters", {}).items() if v.get("without_the_term", {}).get("scatter") and v.get("with_the_term", {}).get("scatter")]
    if biases:
        without, with_it, name = max(biases)
        evidence.append(f"On the replicas, a deviation of this kind left unmodelled biases {name} by {plain_number(without)} scatters, and by {plain_number(with_it)} with a generic term of its kind.")
    test = []
    back = (facts["cand"] or {}).get("what_the_whole_procedure_returns", {})
    if back:
        test.append(f"Every value comes back within {max(abs(v['bias_in_scatters']) for v in back.values()):.2f} scatters on {facts['cand'].get('read', 30)} replicas that contain the term.")
    if facts["bic"]:
        test.append(facts["bic"] if "not comparable" in facts["bic"] else facts["bic"].replace("Bayesian information criterion: ", "The Bayesian information criterion ").replace(", where a decrease of more than 2 is required.", ", where a fall of more than 2 is required."))          # a rise reads as a rise against the rule (review 42, moderate 5)
    if facts["admitted"] is not None:
        test.append("The replicas admit the enlarged model." if facts["admitted"] else "The replicas do not admit the enlarged model.")
    if trend:          # what the term did to the evidence that it was written for, his message 257
        test.append(trend)
    if d["chain"] is None:
        nxt = "The model stays as it was." if facts["verdict"] == "Not a term" else "The next step starts from the model that the analyst adopted."          # a test of the model adopts none (review 45, minor 4)
    else:
        nxt = f"Turn {turn + 1} starts from model {k + 1}." if d["kept"] else f"The model stays model {k}."
    def laid(test: list) -> tuple:
        rows = [(words["record_heads"][0], " ".join(evidence), pages.INK), (words["record_heads"][1], f"{d['wrote'][:1].upper() + d['wrote'][1:]}.".replace("`", ""), pages.INK),
                (words["record_heads"][2], " ".join(test) or "No statistical test was recorded.", pages.INK), (words["record_heads"][3], f"{facts['verdict']}: {facts['why']} {nxt}", facts["colour"])]
        bodies = [pages.wrapped(text, meeting.LEAST, 0.68 * meeting.SIZE[0] / pages.SIZE[0], "bold" if j == 3 else "normal") for j, (_, text, _) in enumerate(rows)]
        heads = [pages.wrapped(head, meeting.LEAST, 0.22 * meeting.SIZE[0] / pages.SIZE[0], "bold") for head, _, _ in rows]
        return rows, bodies, heads, [0.044 * max(body.count("\n") + 1, head.count("\n") + 1) for body, head in zip(bodies, heads)]          # each row as tall as its body or its wrapped head, at the height of a line as it is drawn, measured on BLIND4
    if part == 1:          # the first half leaves the biases to the table of the proposals, which shows them (his message 311)
        evidence = [x for x in evidence if not x.startswith("On the replicas, a deviation of this kind left unmodelled biases")]
    rows, bodies, heads, heights = laid(test if not (part and trend) else test[:-1] + [trend.split(". ")[0].rstrip(".") + "."])          # a half keeps the first sentence of the trend, which the slide of evidence before it holds whole
    for shorter in ([trend.split(". ")[0].rstrip(".") + "."], []) if trend and not part else ():          # a record too long for its slide keeps the first sentence of its trend, then none: the whole stands on the slide of evidence before it
        if sum(heights) <= fig.below - 0.03 - FOOT_CLEAR:
            break
        rows, bodies, heads, heights = laid(test[:-1] + shorter)
    shown = {0: range(4), 1: range(2), 2: range(2, 4)}[part]          # the rows of this slide: all, or one half of them
    y = fig.below - 0.03
    gap = record_gap([heights[j] for j in shown], y)
    for j, ((_, _, colour), head, body, height) in ((j, (rows[j], heads[j], bodies[j], heights[j])) for j in shown):
        fig.text(0.05, y, head, fontsize=meeting.LEAST, fontweight="bold", color=meeting.PALETTE[j], va="top", linespacing=1.1)
        fig.text(0.29, y, body, fontsize=meeting.LEAST, color=colour, va="top", linespacing=1.15, fontweight="bold" if j == 3 else "normal")
        y -= height + gap
    meeting.inside(fig, title)
    return fig


def _decision_slide(title: str, d: dict, loop: Path, k: int, after: tuple | None, foot: str, words: dict, asked: bool | None = True, part: int = 0):
    """Phase 8 of a turn: the new term on replicas that contain it, analysed by the whole procedure, the gain in chi square
    on the data beside the replicas, the Bayesian information criterion against its rule, whether the replicas admit the
    enlarged model, and the verdict. `after` is the loop and the model of the enlarged model, where it was fitted. With
    `part` 1 the replicas that contain the term alone, with 2 the criterion, the admission and the verdict (his message
    311)."""
    from toy import meeting, pages
    if part:
        title = f"{title}, {part} of 2"
    fig = _frame(title, foot, 8)
    y = fig.below - 0.02
    facts = facts_of(d, loop, k, after, asked)
    cand, bic, verdict, why, colour = facts["cand"], facts["bic"], facts["verdict"], facts["why"], facts["colour"]
    if part == 2:          # the replicas of the term stand on the half before
        pass
    elif cand and cand.get("what_the_whole_procedure_returns"):
        y = _say(fig, y, words["decision"].format(sets=cand.get("read", cand.get("sets", 30))))
        rows = [[name, f"{v['put_in']:.4g}", f"{v['mean']:.4g}", f"{round(v['bias_in_scatters'], 2) + 0.0:.2f}", f"{v['inside_one_bar']:.2f}"] for name, v in cand["what_the_whole_procedure_returns"].items()]
        y = _table(fig, y, ["parameter", "value\nput in", "mean\nrecovered", "bias in\nscatters", "within one\nuncertainty"], rows, [0.05, 0.40, 0.53, 0.67, 0.80], 0.33)
        gain = cand.get("gain", {})
        if gain.get("on_the_records") is not None:
            y = _say(fig, y, f"Gain in chi-squared from the new term: {gain['on_the_records']:.1f} on the data, {pages_pm(gain['mean'], gain['scatter'])} on the replicas.")
    elif d["chain"] is None:
        y = _say(fig, y, f"Model of step {d.get('step')} fitted by twinfit: reduced chi-squared {d['chi'][0]:.3f} before, {d['chi'][1]:.3f} after.")
    if part == 1:          # the criterion, the admission and the verdict on the half after
        meeting.inside(fig, title)
        return fig
    if bic:
        y = _say(fig, y, bic)
    if facts["admitted"] is not None:
        y = _say(fig, y, "Replicas of the enlarged model admit its fit: " + ("yes." if facts["admitted"] else "no."))
    fig.text(0.955, 0.945, verdict, ha="right", va="top", fontsize=meeting.TITLE, fontweight="bold", color=colour)          # the verdict where the eye lands first, beside the title
    _say(fig, y, f"{verdict}: {why}", colour=colour, weight="bold")
    meeting.inside(fig, title)
    return fig


def recorded_model(loop: Path, where: Path, k: int) -> int | None:
    """The model of a loop on his recorded traces whose chain of terms is that of model k of `where`, the loop of the
    synthetic traces, or None where the loop on his traces kept no such model: his message 288, each model of a run shown
    on his traces too."""
    chain = chain_of(where, k)
    found = [j for j in range(len(chain) + 2) if (Path(loop) / f"model_{j}.json").exists() and chain_of(loop, j) == chain]
    return found[0] if found else None


def _normalized_slide(title: str, stages: list, foot: str, unit: str):
    """The normalized residuals of two models side by side, every sample over its noise with its dense local mean: each
    stage is (its loop, its model, its chi square a degree of freedom, its label)."""
    import numpy as np
    from toy import meeting, pages
    fig = meeting.slide(title, "Residual of every sample over its noise standard deviation, with its local mean.", foot)
    for c, (where, k, chi, label) in enumerate(stages):
        ax = fig.add_axes([0.09 + 0.46 * c, meeting.BASE, 0.40, fig.below - meeting.BASE - 0.17])
        x = np.concatenate([v[0] for v in curves_of(where, k).values()])
        z = np.concatenate([v[3] for v in curves_of(where, k).values()])
        ax.plot(x, z, ".", ms=2, alpha=0.08, color=pages.GREY, rasterized=True)
        inner = (x >= np.percentile(x, 1)) & (x <= np.percentile(x, 99))
        ax.plot(*smooth_mean(x[inner], z[inner]), "-", lw=2.5, color=meeting.PALETTE[0])          # the dense local mean, his message 247
        ax.axhline(0.0, color=pages.INK, lw=1.0)
        reach = max(5.0, float(np.percentile(np.abs(z), 99.5)))
        ax.set_ylim(-reach, reach)
        ax.set_title(f"{label}, reduced\nchi-squared {chi:.3f}" if chi < 2.0 else f"{label}, reduced\nchi-squared {chi:.3g}", fontsize=meeting.LEAST, loc="left", color=pages.INK, pad=10)
        ax.set_xlabel(f"distance from the line centre{' (' + unit + ')' if unit else ''}", fontsize=meeting.LABEL)
        if c == 0:
            ax.set_ylabel("normalized residual", fontsize=meeting.LABEL)
        meeting.plain(ax)
    meeting.inside(fig, title)
    return fig


def said_of(artist) -> str:
    """What an artist of a slide says: its text, or the text of the text it holds, as a cell of a table answers get_text
    with a text and not a string (the slides of the package in a blind report, 2026-10-22)."""
    x = artist.get_text()
    x = x.get_text() if hasattr(x, "get_text") else x
    return x if isinstance(x, str) else ""


def slide_words(figs: list) -> dict:
    """The words of each slide, every text drawn on it numbers and names of axes included, as a reader counts them, and
    the foot of each: what row M370 reads of a deck."""
    import re
    words, feet, titles, sentences = [], [], [], []
    for fig in figs:
        drawn = [said_of(a) for a in fig.findobj(lambda a: hasattr(a, "get_text") and hasattr(a, "get_visible") and a.get_visible())]
        words.append(sum(len(x.split()) for x in drawn))
        feet.append(next((t.get_text() for t in fig.texts if t.get_text().startswith("twinfit")), ""))
        titles.append(max(fig.texts, key=lambda t: t.get_fontsize()).get_text().replace("\n", " ") if fig.texts else "")
        own = [x.get_text() for x in fig.texts if x.get_visible() and x.get_text() not in (titles[-1].replace(" ", "\n"), feet[-1]) and x.get_text().replace("\n", " ") != titles[-1]
               and 0.08 <= x.get_position()[1] <= 0.92]          # the frame of a slide, its corner, the strip of its phases and its foot, in the bands at its top and its bottom, is no sentence
        own += [x.get_text() for ax in fig.axes for x in ax.texts if x.get_visible()] + [said_of(c) for ax in fig.axes for table in ax.tables for c in table.get_celld().values()]
        sentences.append(sum(1 for x in own for w in x.split() if re.search(r"[A-Za-z]", w)))          # a word holds a letter: the signs of a formula and the numbers are no words          # what a reader reads as sentences and tables: the title, the foot, the marks of the phases, the numbers along the axes and the keys left out (his message 311)
    return {"titles": titles, "words": words, "feet": feet, "sentences": sentences}


def run_slides(run: Path, out: Path, root: Path = ROOT, recorded: Path | None = None, visual: bool = False) -> list:
    """A blind run as slides for a reader who knows neither the line nor the package, his messages 239, 243 and 247: "the
    whole point of the report is to make clear the decision making process... don't give anything for granted". In the
    look of the slides of the group meeting, by `toy/meeting.py`: the question, the physical problem, what the agent was
    told and what was hidden, the two modes of use, the eight phases of a turn, the first model and why, and the terms
    the phases use, each defined. Then each turn phase by phase: its model with what is fixed, global and local, its joint
    fit by sweep with what it leaves, its checks, its windows, its diagnosis, the generic candidates of twinfit, and each
    term the analyst wrote with its reason and its verdict. Then why the run ended, what the first and the last models
    leave, what the procedure established, the answer beside the truth, and the cost. The whole log stays in
    `run_report`. Returns the paths of the slides and of their PDF, which holds the slides alone."""
    import tomllib
    from twinfit.candidates import covariates_of
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages
    from toy import meeting, pages
    run, out = Path(run), Path(out)
    ran = json.loads((run / "run.json").read_text()) if (run / "run.json").exists() else {}
    turns = turns_of_the_run(run)
    if not turns:
        raise SystemExit(f"REFUSED: {run} keeps no model in any step: no slide has anything to show")
    words = tomllib.loads((Path(root) / "toy" / "words.toml").read_text())
    mine = dict(words["run_slides"], **(words["run_slides_by_hand"] if ran.get("by_hand") else {}))
    mine_recorded = final_loop(Path(recorded)) if recorded is not None else None          # the loop of the same run on his recorded traces, his message 288          # a walk by hand says what it was, and is no blind test (the walk of 2026-10-20)
    title = run.parent.name if run.name == "live" else run.name
    foot = "twinfit, a blind test on a hidden project." if ran.get("project") else mine["foot"] if ran.get("by_hand") else f"twinfit, a blind test on synthetic traces, draw {ran.get('draw')}."
    his_foot = mine.get("foot_his", "twinfit, on his 8 recorded traces.")          # a slide of his recorded traces says so at its foot (review 71, moderate 5)
    texts = []          # the slides made of sentences alone, which the visual form leaves out (his message 297)
    project = ran.get("project")
    unit = unit_of(run)
    acts = actions_of(next(p for p in (run / "stream.jsonl", run / "stream.jsonl.gz") if p.exists()))
    tally = {key: sum(a["kind"] == key for a in acts) for key in ("package", "own fit", "subagent own fit", "term", "red team")}
    if ran.get("by_hand"):          # a walk by hand keeps no stream of commands: its record says how often it ran the loop (the walk of 2026-10-20)
        tally["package"] = int(ran.get("loops", 0))
    minutes = round(float(ran["seconds"]) / 60) if "seconds" in ran else None
    ended = "ended by the analyst" if not ran.get("stopped") else "stopped before the analyst ended it, as the slide of its termination says"          # the record of the run says why (review 42, moderate 6)
    what = "a hidden project of another field" if project else f"8 synthetic traces of the 5S to 6S line, draw {ran.get('draw')}"
    described = {}
    for path in sorted(run.rglob("description.toml")):
        try:
            described = tomllib.loads(path.read_text())
            break
        except (tomllib.TOMLDecodeError, OSError):
            continue
    slides = []
    fig = _frame(mine["title"], foot)
    opening = mine["opening"] if not project else [f"The curves: {str(described.get('problem', {}).get('what', 'curves of another field')).split('. ')[0].strip()}, whose physics is known to us and hidden from the analysis.",
                                                   "The test: find what the curves hold beyond their first model, and how large each part is, from the curves alone.", mine["opening"][2]]
    _points(fig, opening, fig.below - 0.04)
    meeting.inside(fig, mine["title"])
    slides.append(fig)
    texts.append(fig)
    decided = decisions_of(run)
    final = final_loop(run)
    chained = bool(decided) and decided[0]["chain"] is not None
    stages = [(final, k) for k in range(last_model(final) + 1)] if chained and final else [(t["loop"], t["model"]) for t in turns]

    def stage(where: Path, k: int) -> dict:          # a model as the slides before and after the turns read it
        table = json.loads((Path(where) / f"model_{k}.json").read_text())
        return {"loop": where, "model": k, "chi": table["chi_square"] / max(table["degrees_of_freedom"], 1)}
    first, last = stage(*stages[0]), stage(*stages[-1])
    fig = _frame(mine["physics_title"], foot)
    physics = mine["physics"] if not project else [s.strip() + "." for s in str(described.get("problem", {}).get("what", "")).split(". ") if s.strip()][:4]
    _points(fig, physics, fig.below - 0.03, width=0.44, size=BODY - 2.0)
    ax = fig.add_axes([0.64, meeting.BASE, 0.33, fig.below - meeting.BASE - 0.04])
    _traces(ax, first["loop"], first["model"], unit)
    ax.plot([], [])
    meeting.inside(fig, mine["physics_title"])
    slides.append(fig)
    if not project:          # why a physicist cares about the widths, before what the analysis was told
        fig = _frame(mine["matters_title"], foot)
        _points(fig, mine["matters"], fig.below - 0.04)
        meeting.inside(fig, mine["matters_title"])
        slides.append(fig)
        texts.append(fig)
    fig = _frame(mine["known_title"], foot)
    if project:
        told = [str(v.get("words", k)) for k, v in described.get("setup", {}).items() if isinstance(v, dict)][:4] or ["the curves and the first model alone"]
        hidden = ["what the curves hold beyond the first model, and its size: it stands after the answer, at the end"]          # the truth after the answer alone (the fresh red team of 0.1.0rc16, mislead 5)
    else:
        told = mine["known"]
        hidden = mine.get("hidden") or [withheld(words["detail"][str(config["parameter"][n]["enters"])]["short"], config["parameter"][n], w) for const, config in [model.load(Path(root) / "toy", root)] for n, w in truth(const, config).items()]          # what was withheld, and the truth after the answer alone (the fresh red team of 0.1.0rc16, mislead 5). A walk by hand, whose researcher knew the truth, names each term by its symbol (his message 317)
    fig.text(0.06, fig.below - 0.02, mine["known_head"], fontsize=BODY, fontweight="bold", color=pages.INK, va="top")
    left_end = _points(fig, told, fig.below - 0.09, width=0.40, size=BODY - 2.0)
    fig.text(0.53, fig.below - 0.02, mine["hidden_head"], fontsize=BODY, fontweight="bold", color=pages.INK, va="top")
    end = _points(fig, hidden, fig.below - 0.09, left=0.53, width=0.40, size=BODY - 2.0)
    fig.text(0.06, max(0.21, min(end, left_end) - 0.01), pages.wrapped(mine["hidden_foot"] if not project else "The curves are synthetic, drawn from a model that holds these effects with the noise of the experiment. The analysis saw no list and no name.", BODY - 2.0, 0.88 * meeting.SIZE[0] / pages.SIZE[0]),
             fontsize=BODY - 2.0, color=pages.GREY, va="top", linespacing=1.25)
    meeting.inside(fig, mine["known_title"])
    slides.append(fig)
    texts.append(fig)
    slides.append(_modes_slide(mine["use_title"], foot, mine, by_hand=bool(ran.get("by_hand"))))          # with an agent and without one, as two lanes, his messages 247 and 311
    fig = _frame(mine["phases_title"], foot)          # the phases of a turn, at the start, his message 247
    who = [meeting.PALETTE[1] if p[1] == "analyst" else meeting.PALETTE[0] for p in mine["phases"]]
    _boxes(fig, [f"{j + 1}. {p[0]}" for j, p in enumerate(mine["phases"][:4])], fig.below - 0.14, colours=who[:4])
    _boxes(fig, [f"{j + 5}. {p[0]}" for j, p in enumerate(mine["phases"][4:])], fig.below - 0.43, colours=who[4:])
    _say(fig, 0.22, mine["phases_foot"])
    meeting.inside(fig, mine["phases_title"])
    slides.append(fig)
    fig = _frame(mine["first_title"], foot)
    first_words = mine["first"] if not project else [f"{short_name(chain_of(first['loop'], 0)[0], 12)}.", reason_of(Path(first['loop']).parent / 'model.py') or "The first model that the description gives."]
    _points(fig, first_words, fig.below - 0.03, width=0.44, size=BODY - 2.0)
    ax = fig.add_axes([0.64, meeting.BASE, 0.33, fig.below - meeting.BASE - 0.10])
    _traces(ax, first["loop"], first["model"], unit, one=True)
    ax.set_title(f"one spectrum, reduced\nchi-squared {first['chi']:.3g}", fontsize=meeting.LEAST, loc="left", color=pages.INK, pad=8)
    meeting.inside(fig, mine["first_title"])
    slides.append(fig)
    def points_slide(head: str, points: list, size: float = BODY, at_foot: str | None = None, kept: bool = False) -> None:          # a slide of sentences: a definition, or one step of the argument, on as many slides as its words ask (his message 311)
        parts = split_points(points)
        for j, part in enumerate(parts, 1):
            title = head if len(parts) == 1 else f"{head}, {j} of {len(parts)}"
            fig = _frame(title, at_foot or foot)
            _points(fig, part, fig.below - 0.04, size=size)
            meeting.inside(fig, title)
            slides.append(fig)
            if not kept:          # a slide of sentences that the visual form keeps all the same
                texts.append(fig)
    points_slide(mine["notation_title"], mine["notation"], BODY - 2.0)
    points_slide(mine["for_title"], mine["for"])
    points_slide(mine["replicas_title"], mine["replicas"])
    if rebuilt_first := initial_of(run, chained, stages):          # how a window is read, on the initial model of the run, his message 245: "show how the windowing is performed and used to detect deviation from a complete model", and not on the first model of the last step's loop (review 45, moderate 2)
        slides.append(_how_windows_slide(*rebuilt_first[0], foot, mine))
    points_slide(mine["line_title"], mine["line"])
    points_slide(mine["report_title"], mine["report"])
    for half in (0, 4):          # each phase defined, with who performs it, now that its terms are
        points_slide(mine["phases_title"] + f", {half + 1} to {half + 4}", [f"{half + j + 1}. {p[0][:1].upper() + p[0][1:]}, by {'the analyst' if p[1] == 'analyst' else 'twinfit'}: {p[2]}" for j, p in enumerate(mine["phases"][half:half + 4])], meeting.LEAST)
    slides.append(_algorithm_slide(mine["algorithm_title"], foot, mine))          # the procedure as a figure, his messages 248 and 311
    drawn = set()          # the figures of evidence drawn already
    scales = {}          # the half height that each figure of evidence was drawn at, so that the slide after a term reads beside the slide drawn before it
    loops = {}          # the fits of the models of each loop, made again once: the windows, the evidence and the record read them. Each turn phase by phase, his message 247: "logging turn by turn, phase by phase, focusing in particular to the core issues of the package: how the analysis is performed and how the new terms are chosen to be attempted"

    def fitted(where: Path, k: int):
        if Path(where) not in loops:
            loops[Path(where)] = rebuilt(run) if chained else rebuilt(run, where)
        return next(iter(loops[Path(where)][k:k + 1]), None)
    names = [p[0] for p in mine["phases"]]
    each = []          # whether the report asked for a term, turn by turn: the turn from which the first clause held is said at the termination (review 45, minor 4)
    for at, (where, k) in enumerate(stages):
        def head(p: int, what: str = "") -> str:
            return f"Turn {at + 1}, phase {p}: {names[p - 1]}" + (f", {what}" if what else "")
        curves = curves_of(where, k)
        parts = model_parts(Path(where).parent / "model.py")
        table = json.loads((Path(where) / f"model_{k}.json").read_text())
        chi = table["chi_square"] / max(table["degrees_of_freedom"], 1)
        slides.extend(_model_slide(f"Turn {at + 1}, phase 1: " + (f"model {k}" if chained else f"the model of step {turns[at]['step']}"), model_said(where, k, mine, unit, why=at == 0), foot))
        slides.extend(package_slides(where, k))          # its equation and what each term does, from the package's own log, his message 309          # turns from 1 and models from 0 (review 39, moderate 7)
        said = joint_said(parts, fit_of(where, k), chi, len(curves), unit)
        for words_of, order in sweeps_of(list(curves)):
            slides.append(_sweep_slide(head(2, words_of), curves, order, said, foot, unit, 2))
        if mine_recorded is not None:          # the same model on his recorded traces, its terms decided on the synthetic ones, his message 288
            j = recorded_model(mine_recorded, where, k)
            if j is None:
                points_slide(head(2, "his recorded traces"), [mine["recorded_none"]], at_foot=his_foot)
            else:
                got = curves_of(mine_recorded, j)
                held = json.loads((Path(mine_recorded) / f"model_{j}.json").read_text())
                theirs = mine["recorded_said"] + " " + joint_said(model_parts(Path(mine_recorded).parent / "model.py"), fit_of(mine_recorded, j), held["chi_square"] / max(held["degrees_of_freedom"], 1), len(got), unit)
                for words_of, order in sweeps_of(list(got)):
                    slides.append(_sweep_slide(head(2, words_of + ", his traces"), got, order, theirs, his_foot, unit, 2))
                his_report = proposed(recorded, mine_recorded, j)[0]          # and what the report reads in his traces against their own twin, which the statistics of the package are for (his message 289)
                if his_report:
                    slides.append(_diagnosis_slide(head(5, "his traces"), his_report, held["chi_square"] / max(held["degrees_of_freedom"], 1), his_foot, mine))
        slides.append(_checks_slide(head(3), where, k, foot, mine, first=at == 0))
        windows = fitted(where, k)
        if windows:
            slides.append(_windows_slide(head(4), *windows, foot, mine["model_windows"].format(model=f"model {k}" if chained else f"the model of step {turns[at]['step']}"), 4))          # every loop of a step starts at a model 0 (review 45, minor 6)
        read_at = turns[at]["loop"] if at < len(turns) and turns[at]["model"] == k and chain_of(turns[at]["loop"], k) == chain_of(where, k) else where
        report, candidates = proposed(run, read_at, k)          # the report and the proposals that the analyst read at this turn, in the loop of its own step: the last run held other proposals for model 0, since a term file then stood beside it (the red team of his message 333, severe 1)
        if not report and read_at != where:          # a step that kept no report of the model: the last run's, as before
            report, candidates = proposed(run, where, k)
        slides.append(_diagnosis_slide(head(5), report, chi, foot, mine))
        if report and any(a["question"] == "seen by" for a in report.get("answers", [])):          # the evidence of every family of statistics, his messages 262 to 265
            slides.extend(_families_slide(head(5, "every family of statistics"), report, foot, mine))
        key = leading_covariate(report, *windows) if windows else None
        seen = (key, tuple(round(t, 1) for bands in trends_of(*windows, key)[1].values() for *_, t in bands)) if key else None
        if key and seen not in drawn:          # the deviation along the setting that the report names, spectrum by spectrum, his message 257, and a model unchanged since an earlier turn drawn once (review 42, minor 12)
            drawn.add(seen)
            scales[seen] = reach_of(trends_of(*windows, key)[1])
            slides.append(_evidence_slide(f"Turn {at + 1}, phase 5: along {key}", *windows, key, foot, mine["evidence"], 5, mine["evidence_scale"], ranking_said(report, covariates_of(windows[0]))))
        slides.extend(in_halves(_proposals_slide, f"Turn {at + 1}, phase 6: three proposals", report, candidates, foot, mine))          # each proposal explained before the choice, his message 250, on two slides where their words ask it (his message 311)
        slides.append(_candidates_slide(f"Turn {at + 1}, phase 6: the proposals compared", report, candidates, foot, mine, first=at == 0))
        for d in [d for d in decided if d["on"] == at]:
            asked = report_asked(report, candidates)
            each.append((at + 1, clause_held(report, phase_file(where, k, "twin").get("admission", {}).get("admitted"))))
            page = page_of(run, d["wrote"], d["on"])
            slides.extend(in_halves(_choice_slide, head(7), d, report, candidates, foot, mine, windows[0] if windows else None, page if asked is not True else ""))
            if reasoned := reasoning_points(page_of(run, d["wrote"], d["on"], lines=True)):          # and the analyst's reasoning in its own words, which choice it took and why (his messages 310 and 311)
                points_slide(head(7, "the analyst's reasoning"), reasoned, BODY - 2.0, kept=True)
            named = [c["term"] for c in candidates[:3] if isinstance(c, dict) and c.get("term")]
            if named:          # the term written against the three proposals, where the slides never said that it was none of them (the red team of his message 328)
                which = chosen_among(named, d["wrote"])
                points_slide(head(7, "the term chosen against the three proposals"),
                             ["The three proposals of the report: " + ", ".join(f"{i + 1}. {n}" for i, n in enumerate(named)) + ".",
                              f"The term written, {d['wrote']}, is proposal {which}." if which else f"The term written, {d['wrote']}, is none of the three: the analyst chose it on its own reading, as its reasoning says."],
                             BODY - 2.0, kept=True)
            if read := reasoning_points(page_of(run, d["wrote"], d["on"], lines=True, pattern=r"\bparameters of the best model\b")):          # the table that the reasoning cites as below (the red team of his message 328)
                points_slide(head(7, "the values that the analyst read"), read, BODY - 2.0, kept=True)
            if team := reasoning_points(red_team_of(run, at + 1)):          # and what its red team said of the term, before its fit, which the slides cited and did not show
                points_slide(head(7, "the red team of the term, in its words"), team, BODY - 2.0, kept=True)
            else:          # a turn with no red team says so, where the deck went on in silence (the red team of his message 333, P6)
                points_slide(head(7, "the red team of the term"), ["The analyst wrote no red team of this term before its fit: its page of the step holds none."], BODY - 2.0, kept=True)
            after = ((final, k + 1) if d["kept"] else None) if d["chain"] else ((turns[at + 1]["loop"], turns[at + 1]["model"]) if at + 1 < len(turns) else None)
            slides.extend(in_halves(_decision_slide, head(8), d, where, k, after, foot, mine, asked))
            enlarged = fitted(*after) if after is not None else None
            shape = deviation_after(report, proposed(run, *after)[0] if after is not None else {}, k, None if d["chain"] else (f"step {turns[at]['step']}", f"step {d.get('step')}"))
            trend = ""
            if key and windows:
                before = largest_trend(trends_of(*windows, key)[1])
                places_after = trends_of(*enlarged, key)[1] if enlarged and key in enlarged[0].records[0].design else None
                later = trend_at(places_after, before[0], before[1]) if places_after else None
                trend = trend_said(key, before, later, k, largest_trend(places_after) if places_after else None,
                                   None if d["chain"] else (f"step {turns[at]['step']}", f"step {d.get('step')}"))
                if later is not None:          # the same evidence on the enlarged model: the term is borne out where the trend that it was written for has gone
                    at_least = scale_before(scales, seen, trends_of(*windows, key)[1])
                    drawn.add((key, tuple(round(t, 1) for bands in places_after.values() for *_, t in bands)))          # and its figure is not drawn again at the next turn (review 43, minor 3)
                    scales[(key, tuple(round(t, 1) for bands in places_after.values() for *_, t in bands))] = max(reach_of(places_after), at_least)
                    slides.append(_evidence_slide(f"Turn {at + 1}, phase 8: along {key}, " + (f"model {k + 1}" if d["chain"] else f"step {d.get('step')}"), *enlarged, key, foot, mine["evidence_after"], 8, mine["evidence_scale"], "", at_least))
            trend = " ".join(s for s in (shape, trend) if s)          # what the term did to the deviation that it was written for, first, then to the trend along the setting (review 44, minor 5)
            slides.extend(in_halves(_record_slide, f"Turn {at + 1}: the decision in one view", d, k, at + 1, report, candidates, chi, facts_of(d, where, k, after, asked), foot, mine, trend, None if d["chain"] else f"the model of step {turns[at]['step']}"))          # his message 248: the story of the decision, from the data, the physics and the test
        before = len(slides)
        slides.append(_experiment_slide(f"Turn {at + 1}, phase 8: three proposals for the experiment", table, foot, mine))          # at the end of each turn, decided or not, his message 310
        if said_of_them := experiment_words(table):          # each proposal in words, where the bars gave a label alone (the red team of his message 328)
            points_slide(f"Turn {at + 1}, phase 8: the three proposals for the experiment, in words", said_of_them, BODY - 2.0, kept=True)
        experiment_slides = len(slides) - before
    stop = next((section_of(Path(folder) / "step.md", r"\bstop\b") for _, _, folder in reversed(steps_of(run)) if section_of(Path(folder) / "step.md", r"\bstop\b")), "") or stop_in_notes(run)
    final_turn = len(stages)
    closing = [slides.pop() for _ in range(experiment_slides)][::-1]          # the three proposals for the experiment end the last turn too, after its choice to stop and its decision: the last turn had no slide of reasoning nor of decision of its own, and its stop stood in the closing slides (his message 328)
    reason_title, record_title = f"Turn {final_turn}, phase 7: choice, the analyst's reason to stop", f"Turn {final_turn}, phase 8: decision, the procedure ends"
    cut = stop.split()
    if cut:          # the analyst's reason on a slide of its own: under the record of the run it took the slide to 102 words (the red team of his message 320)
        fig = _frame(reason_title, foot, 7)
        fig.text(0.05, fig.below - 0.03, "Reason recorded by the analyst, verbatim:", fontsize=meeting.LEAST, color=pages.GREY, va="top")
        fig.text(0.05, fig.below - 0.075, quoted_within(cut, fig.below - 0.075 - FOOT_TOP, BODY - 2.0, 0.90 * meeting.SIZE[0] / pages.SIZE[0]),
                 fontsize=BODY - 2.0, color=pages.INK, va="top", linespacing=1.25, style="italic")
        meeting.inside(fig, reason_title)
        slides.append(fig)
    fig = _frame(record_title, foot, 8)          # why the analyst ended the run, in its own words where its page says so
    final_report = proposed(run, *stages[-1])[0]
    record = json.loads((Path(stages[-1][0]) / "best.json").read_text()) if (Path(stages[-1][0]) / "best.json").exists() else {}
    end = _points(fig, [mine["stopping"], ended_by(final_report, stop, ran, minutes, record, phase_file(*stages[-1], "twin").get("admission", {}).get("admitted"))] + ([run_ended_said(ran["ended"])] if isinstance(ran.get("ended"), str) else [])
                  + ([mine["held_from"].format(turn=held_from(each)[0], more=said_number(held_from(each)[1]), turns="turns" if held_from(each)[1] > 1 else "turn")] if held_from(each) else []), fig.below - 0.03, size=meeting.LEAST)          # how the run ended, beside the analyst's words of the time (review 43, minor 8)
    if not cut and not ran.get("stopped"):          # a stopped run has said why, from its record
        _say(fig, end - 0.02, "The analyst ended the run without recording a reason.")
    meeting.inside(fig, record_title)
    slides.append(fig)
    slides.extend(closing)
    decided = decisions_of(run)
    if decided:          # every decision in one view, before what it left: the story that the slides before told one by one
        fig = _frame("Summary of the iterations", foot)
        _boxes(fig, ["initial model"] + [short_name(d["wrote"], 9) + ("" if d["chain"] is None else ", accepted" if d["kept"] else ", rejected") for d in decided][:5], fig.below - 0.20,          # a name that names its term (review 40, minor 7)
               colours=[meeting.PALETTE[0]] + [meeting.RIGHT if d["kept"] else meeting.WRONG for d in decided][:5])
        _points(fig, [f"Reduced chi-squared: {decided[0]['chi'][0]:.3g} for the initial model, then " + ", ".join(f"{d['chi'][1]:.3f}" for d in decided if d["kept"] and d["chi"][1] is not None) + ".",
                      "Green: term accepted. Orange: term rejected, with the reason stated on its slide."], fig.below - 0.45)
        meeting.inside(fig, "Summary of the iterations")
        slides.append(fig)
    slides.append(_normalized_slide("Normalized residuals, initial and final models", [(first["loop"], first["model"], first["chi"], "initial model"), (last["loop"], last["model"], last["chi"], "final model")], foot, unit))
    ends = [(stage, recorded_model(mine_recorded, stage["loop"], stage["model"])) for stage in (first, last)] if mine_recorded is not None else []
    if ends and all(j is not None for _, j in ends):          # and the same two models on his recorded traces, his message 288
        chis = [json.loads((Path(mine_recorded) / f"model_{j}.json").read_text()) for _, j in ends]
        slides.append(_normalized_slide("Normalized residuals of his recorded traces, initial and final models", [(mine_recorded, j, c["chi_square"] / max(c["degrees_of_freedom"], 1), label) for (_, j), c, label in zip(ends, chis, ("initial model", "final model"))], his_foot, unit))
    sizes, said = sizes_of_the_answer(run, root) if (run / "answer.json").exists() else ([], "")
    if decided:          # what the package gave, in the numbers of this run: the slide that a reader who skips the others should read, his message 244
        if said_by_the_analyst := answer_in_its_words(run):          # the analyst's own answer, with what it says that it could not tell apart, before what the procedure established (the red team of his message 328, severe 2)
            points_slide("The analyst's answer, in its words", said_by_the_analyst, BODY - 2.5, kept=True)
        points_slide(mine["gave_title"], gave_of(run, decided, sizes), BODY - 2.5, kept=True)          # on as many slides as its words ask: a run that kept five values ran past the page (the run of his message 310). A slide that a reader who skips the others reads, kept in the visual form
    if sizes:
        slides.append(meeting.of_the_sizes("Comparison with the true values", said, sizes, foot, "reported / true, bars of the answer"))          # a bar of 0.040 drawn where the fit gave 0.021, and the slide did not say whose it was (the blind red team of 0.1.0rc16, mislead 8)
        if read := biased_said(sizes, record, said, mine):          # what the comparison shows, said and not only drawn (review 46)
            points_slide(mine["biased_title"], read, BODY - 2.0)
    if checked := held_at_theory(run):          # a size far from its theory, held there with a family beside it: the check of Mercury, which reads no truth (row M314)
        points_slide(mine["held_title"], held_said(checked, unit, mine), meeting.LEAST)
    fig = _frame("Computational cost and limitations", foot)
    _points(fig, cost_said(minutes, ended, tally, decided) + run_said(run), fig.below - 0.04)          # the package, the model and the tokens of the run, named nowhere before (the red team of his message 328)
    meeting.inside(fig, "Computational cost and limitations")
    slides.append(fig)
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    if visual:          # his message 297: mainly figures, the slides of sentences alone left out, the recorded slide that has no model kept
        slides = [f for f in slides if not any(f is t for t in texts) or any(his_foot == x.get_text() for x in f.texts)]
    (out / "slides.json").write_text(json.dumps(slide_words(slides), indent=1) + "\n")          # the words and the foot of each slide, which row M370 reads
    with PdfPages(out / f"{title}.pdf", metadata={"Creator": None, "Producer": None, "CreationDate": None}) as pdf:
        for j, fig in enumerate(slides, 1):
            path = out / f"slide {j:02d}.png"
            fig.savefig(path, dpi=110, metadata={"Software": None})
            pdf.savefig(fig)
            plt.close(fig)
            paths.append(path)
    return paths + [out / f"{title}.pdf"]


def ended_by(final_report: dict, stop: str, ran: dict, minutes, record: dict | None = None, admitted: bool | None = None) -> str:
    """Which clause of the stopping rule of the brief ended a run (review 40, moderate 5), read from the record of the
    run and from the package's own record of its final loop, `best.json`, which says whether one answer only just
    passes the threshold (review 42, moderate 1), and never from the analyst's words alone: a run stopped before the
    analyst ended it, no deviation above the threshold, one answer only just past it, or neither."""
    from twinfit.loop import ACTIONABLE
    found = detected(final_report, None, ACTIONABLE)          # every answer of the questions that name a term, as the package counts them (review 43, minor 2)
    record = record or {}
    if ran.get("stopped"):
        return f"The run was stopped before the analyst ended it, after {minutes} minutes, so no clause of the rule ended it. Its record says why: \u201c{ran['stopped']}\u201d" if isinstance(ran["stopped"], str) else f"The run was stopped before the analyst ended it, after {minutes} minutes, so no clause of the rule ended it."
    if final_report and not found and admitted is not False:          # the first clause holds with the admission of the replicas, as `clause_held` reads it (review 48, minor 3)
        return "The run ended on the first clause: the replicas admit the final model and its report detects no deviation above the threshold." if admitted else "The run ended on the first clause: the report on the final model detects no deviation above the threshold."
    if final_report and not found:
        return "The report on the final model detects no deviation above the threshold, but the replicas do not admit the model, so the first clause does not hold: the analyst ended the run on its own reading, below."
    if found and record.get("one_answer_just_past"):
        return "The report on the final model detects one answer only just past the threshold, " + found[0] + ", as the package's record of the loop says: the run ended on the second clause."
    if found:
        whose = (f"The package's record does not count {'them' if len(found) > 1 else 'it'} as one answer only just past" if record
                 else "No record of the final loop was read, so the second clause is not read from one")          # a record said read where none was (review 44, minor 1)
        return (f"The report on the final model detects {said_number(len(found))} {'answers' if len(found) > 1 else 'answer'} past the threshold: " + ", and ".join(found[:3]) + (", and more" if len(found) > 3 else "")
                + f". {whose}, so neither of the first two clauses holds, and the analyst ended the run on its own reading, below.")
    return "The report was not computed for the final model, so the rule was not read on it."


def initial_of(run: Path, chained: bool, stages: list) -> list:
    """The initial model of a run, refitted, as a list of one (problem, fit) or none: the first model of the chain of a
    chained run, and for a run that wrote its model anew at each step, the model of its first stage, from that stage's
    own loop. The slide that teaches a window drew the first model of the last step's loop of EXAMPLE1, which leaves
    nothing, as the initial model (review 45, moderate 2)."""
    if not stages:
        return []
    where, k = stages[0]
    return (rebuilt(run) if chained else rebuilt(run, where))[k:k + 1]


def clause_held(report: dict, admitted: bool | None = None) -> bool | None:
    """Whether the first clause of the stopping rule held on a model: the replicas admit it, where their admission is
    read, and its report finds no deviation above the threshold, no answer of the questions that name a term past it
    (review 47). None where the report was not computed. The termination of BLIND4 read the clause from the proposals,
    none beyond the replicas, where its report had two answers past the threshold (review 46, moderate 1)."""
    from twinfit.loop import ACTIONABLE
    if not report:
        return None
    return admitted is not False and not detected(report, None, ACTIONABLE)


def held_from(held: list) -> tuple | None:
    """From which turn the first clause of the stopping rule held at every decision to the last, and at how many
    decisions the analyst went on there, from (turn, held) pairs in their order: each decision from it was the analyst's
    own protocol. None where the clause did not hold at the last decision, or its report was not computed (review 45,
    minor 4)."""
    run = []
    for turn, one in held:
        run = run + [turn] if one is True else []
    return (run[0], len(run)) if run else None


def run_ended_said(ended: str) -> str:
    """How a run ended, from its record, in plain words: a record that ended by itself in 69 minutes of the 105 is "The
    analyst ended the run itself, after 69 of its 105 minutes." The words of a log, its exit code and its zone of time,
    stood on a slide for his PI (review 44, minor 8). A record of another form is quoted as it stands."""
    import re
    found = re.search(r"in (\d+) minutes of the (\d+)", ended)
    if ended.startswith("by itself") and found:
        return f"The analyst ended the run itself, after {found.group(1)} of its {found.group(2)} minutes."
    return f"The run ended at {ended}, by the record of the session." if re.fullmatch(r"\d{4}-\d\d-\d\d \d\d:\d\d(:\d\d)?", ended.strip()) else f"The record of the run: \u201c{ended}\u201d."          # a time alone, with no word of what it is, stood on the last slide of a turn (the red team of his message 328, minor)


def cost_said(minutes, ended: str, tally: dict, decided: list) -> list:
    """What a run cost and what it shows, in sentences: its minutes and how it ended, the runs of twinfit, the terms that
    the loop tried as the summary slide counts them, accepted and rejected, the fits made outside twinfit, and what one run
    on one data set can show. A command that writes a term file is no count of terms (review 39, moderate 8)."""
    chained = bool(decided) and decided[0]["chain"] is not None
    kept = sum(1 for d in decided if d["kept"])
    tried = (f"The loop tried {said_number(len(decided))} {'terms' if len(decided) != 1 else 'term'}: {said_number(kept)} accepted and {said_number(len(decided) - kept)} rejected."
             if chained else f"The analyst changed the model or the data {said_number(len(decided))} {'times' if len(decided) != 1 else 'time'}.")
    return [f"Duration {minutes} minutes, {ended}. twinfit was run {tally['package']} times. {tried}",
            f"Fits outside twinfit: {tally['own fit']} by the analyst, {tally['subagent own fit']} in independent checks of a candidate before twinfit evaluated it.",
            "A single test on one synthetic data set of known content. It illustrates the method, and is not a measurement of how often the method succeeds."]


def term_why(path: Path) -> str:
    """The `why` of a term file, the physics in the words of whoever wrote it, or nothing."""
    import re
    found = re.search(r"^why\s*=\s*(['\"])(.+)\1\s*$", Path(path).read_text(), re.M) if Path(path).exists() else None
    return found.group(2) if found else ""


def toy_said(const: dict, config: dict, root: Path = ROOT) -> str:
    """How the line of the toy is built and each of its terms, for a reader who red teams it (his message 230): the
    docstring of toy/terms, which writes the line as the Fourier transform of one function of time, then for each term
    its file, the value with which the traces are made, its uncertainty, where that value comes from, and whether the
    fit of the last detail holds it or frees it."""
    import ast
    files = {d: name for name, d in const["details"].items()}
    last = max(const["details"].values())
    free = set(config["detail"][str(last)].get("free", []))
    init = ast.get_docstring(ast.parse((Path(root) / "toy" / "terms" / "__init__.py").read_text())) or ""
    parts = ["<h2>The toy model: how its line is built</h2><pre style='white-space:pre-wrap'>" + _escape(init) + "</pre>",
             "<h2>Its terms, one file a term</h2>"]
    for name, cell in sorted(config["parameter"].items(), key=lambda item: item[1]["enters"]):
        path = Path(root) / "toy" / "terms" / f"{files[cell['enters']]}.py"
        status = "free in the fit of the last detail" if name in free else "held at its value in the fit of the last detail, its uncertainty a bar of its own"
        if cell["kind"] == "theory":          # the loop of toy/chosen.py gives the term its theory, centre and bar, which the records may move: the folder holds it, the loop does not
            status += ". In the loop of the toy, as toy/chosen.py writes the term, its theory enters the fit as a value with its bar, which the records may move"
        parts.append(f"<h3>Detail {cell['enters']}: {_escape(name)}</h3><p>Value with which the traces are made {cell['true']:.4g}, uncertainty {cell['uncertainty']:.2g}, of the kind {_escape(cell['kind'])}"
                     + (f": {_escape(cell['source'])}" if cell.get("source") else "") + f". It is {status}.</p><pre style='white-space:pre-wrap;font-size:85%'>" + _escape(path.read_text()) + "</pre>")
    return "\n".join(parts)


def loop_report(home: Path, out: Path, parameters: dict | None = None, title: str = "", root: Path = ROOT) -> Path:
    """The loop of the package on the toy as one page, his message 230: "I would like to see a full report (with the same
    level of details as the one to send to my PI) on the toy model before analysing the real traces". `home` is the
    folder of a run of the loop with the toy's terms written in its order, as `toy/chosen.py` keeps it. For each model:
    the term that it adds and its physics, the first answers of the report beside what the physics of each term that the
    model still lacks would leave, written in `EXPECTED` before any run, and the page of the model by `step_page`. It opens
    with what the toy holds and the models at a glance, and ends with each parameter of the best model against the value
    that made the traces, from `parameters`, as the counted file of the loop holds them. Refuses a folder with no loop."""
    home, out = Path(home), Path(out)
    loop = home / "loop"
    last = last_model(loop) if loop.is_dir() else None
    if last is None:
        raise SystemExit(f"REFUSED: {home} holds no loop folder with a model")
    const, config = model.load(Path(root) / "toy", root)
    true = truth(const, config)
    order = list(chosen.ORDER)
    rows = []
    for k in range(last + 1):
        table = json.loads((loop / f"model_{k}.json").read_text())
        first = (table.get("where") or {}).get("first_of_each_question") or {}
        rows.append(f"<tr><td>{k}</td><td>{_escape(order[k - 1] if 0 < k <= len(order) else ('none' if k == 0 else 'beyond the toy'))}</td>"
                    f"<td>{table['chi_square'] / max(table['degrees_of_freedom'], 1):.3f}</td><td>{(('yes' if table.get('admitted_on_the_twin') else 'no') if shares(loop, k) else 'nothing to admit')}</td>"
                    f"<td>{_escape(first.get('where') or 'nothing beyond the twin')}</td><td>{_escape(first.get('what') or 'nothing beyond the twin')}</td><td>{_escape(first.get('along') or 'none')}</td></tr>")
    parts = [f"<h1>The package on the toy, model by model{': ' + _escape(title) if title else ''}</h1>",
             "<p>The loop of twinfit run on the eight synthetic traces of the toy, at its noise, with the terms of the toy written one at a time in its order, "
             "and a fifth term that the toy does not hold, which the loop must refuse. Nothing here is blind: the terms are named, and what is measured is what the package says at each model, "
             "where the model lacks the terms that come after it. The numbers come from synthetic traces whose truth is known: they test the method, they measure nothing.</p>",
             "<h2>What the toy holds</h2><table border=1 cellpadding=3><tr><th>term</th><th>full width at half maximum at 130 C, MHz</th></tr>"
             + "".join(f"<tr><td>{_escape(n)}</td><td>{true[n]:.4f}</td></tr>" for n in order) + "</table>",
             toy_said(const, config, root),
             "<h2>The models at a glance</h2><table border=1 cellpadding=3><tr><th>model</th><th>the term it adds</th><th>chi square a degree</th><th>the twin admits it</th>"
             "<th>first place of the report</th><th>first kind</th><th>first setting</th></tr>" + "".join(rows) + "</table>",
             "<h2>How to read the plots</h2><dl>" + "".join(f"<dt><b>{_escape(a)}</b></dt><dd>{_escape(b)}</dd>" for a, b in GLOSSARY) + "</dl>"]
    for k in range(last + 1):
        table = json.loads((loop / f"model_{k}.json").read_text())
        first = (table.get("where") or {}).get("first_of_each_question") or {}
        lacks = order[k:]
        parts.append(f"<hr><h2>Model {k}{': ' + _escape(order[k - 1]) + ' added' if 0 < k <= len(order) else ''}</h2>")
        why = term_why(home / "terms" / f"{k:02d}.py") if k else ""
        if why:
            parts.append(f"<p><b>Its physics</b>, as the term file says it: {_escape(why)}.</p>")
        if lacks:
            cells = []
            for term in lacks:
                said = admitted(term) if term in EXPECTED else {"where": None, "changes along": None}
                cells.append(f"<tr><td>{_escape(term)}</td>" + "".join(f"<td>{_escape(' or '.join(said[q]) if said[q] else 'none')}: {'yes' if (first.get(q) is None if said[q] is None else first.get(q) in said[q]) else 'no'}</td>" for q in ("where", "changes along")) + "</tr>")
            parts.append("<h3>The first answers of the report, beside what each term that the model lacks would leave</h3>"
                         f"<p>The report put first: {_escape(json.dumps({q: first.get(q) for q in ('where', 'changes along')}))}. Each cell gives what the physics of the term would leave, written before any run, and whether the report put it first.</p>"
                         "<table border=1 cellpadding=3><tr><th>term still lacking</th><th>place</th><th>what changes along a setting</th></tr>" + "".join(cells) + "</table>")
        else:
            parts.append("<p>The model holds every term of the toy: the report should put nothing beyond its twin.</p>")
        parts.append(step_page(f"model {k}", loop, home, k=k).replace("<h2>", "<h3>").replace("</h2>", "</h3>"))
    if parameters:
        parts.append("<hr><h2>The best model against the values that made the traces</h2><table border=1 cellpadding=3><tr><th>parameter</th><th>fitted</th><th>bar of the fit</th><th>bar of theory</th><th>value that made the traces</th><th>off, in bars</th></tr>"
                     + "".join(f"<tr><td>{_escape(n)}</td><td>{c.get('value', float('nan')):.4g}</td><td>{c.get('bar', float('nan')):.2g}</td><td>{'none' if c.get('bar_of_theory') is None or c.get('bar_of_theory') != c.get('bar_of_theory') else format(c['bar_of_theory'], '.2g')}</td>"
                               f"<td>{c.get('of_the_twin', float('nan')):.4g}</td><td>{'' if c.get('from_the_twin_in_bars') is None else round(c['from_the_twin_in_bars'], 2)}</td></tr>" for n, c in parameters.items()) + "</table>")
    derivation = Path(root) / "toy" / "DERIVATION.md"
    if derivation.exists():          # the first three details derived from the physics alone, by a reader who did not see the code
        parts.append("<hr><h2>Appendix: the first three details from the physics alone</h2>" + prose(derivation.read_text()))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("<!doctype html><meta charset='utf-8'><title>The package on the toy</title><body style='font-family:Arial;max-width:1100px;margin:auto;line-height:1.35'>" + "\n".join(parts) + "</body>")
    return out


def unmet(got: dict) -> list:
    """What a noiseless run must show and its run does not, as written at the head of this program before it ran: the
    run ended, the four terms kept and the fifth refused, under each model that lacks terms a first answer of the report
    that one of them explains, under the right model nothing, and each parameter of the best model within two of its own
    bars of the value that made the traces. A run that falls short fails its proof (the red team of his message 160)."""
    if "failed" in got:
        return [f"the run did not end: {got['failed']}"]
    out = []
    said = against(got)
    if sorted(said["kept"]) != [f"{k:02d}.py" for k in range(1, len(chosen.ORDER) + 1)]:
        out.append(f"the terms kept are {said['kept']}, not the {len(chosen.ORDER)} of the toy")
    if f"{len(chosen.ORDER) + 1:02d}.py" not in said["refused"]:
        out.append("the term that the toy lacks is not refused")
    for name, m in said["models"].items():
        if m["lacks"] and not m["matches"]:
            out.append(f"under {name} the report puts first {m['first']}, which no term that it lacks explains")
        if not m["lacks"] and any(m["first"].get(q) for q in ("where", "shape", "what", "changes along")):
            out.append(f"under {name}, the right model, the report says {m['first']}")
    for name, cell in said["parameters"].items():
        off = cell["from_the_twin_in_bars"]
        if off is None or abs(off) > 2.0:
            out.append(f"{name} is {'without a bar' if off is None else f'{off:.2g} of its bars'} from the value that made the traces")
    return out


def check(root: Path = ROOT) -> list:
    path = Path(root) / "toy" / "noiseless.json"
    if not path.exists():
        return ["toy/noiseless.json is not written: `python toy/noiseless.py run`, then `table`"]
    held = json.loads(path.read_text())
    if held.get("digest") != digest(root):
        return ["toy/noiseless.json was made from other files than those on disk now: run it again"]
    wrong = [] if held["against"] == against(held["run"]) else ["toy/noiseless.json does not hold what its own run gives"]
    wrong += [f"the noiseless run falls short: {line}" for line in unmet(held["run"])]
    for key in ("agent", "agent_at_the_noise"):
        agent = held.get(key)
        if agent is not None and agent["score"] != score(agent["answer"], agent["truth"]):
            wrong.append(f"the score of the {key.replace('_', ' ')} in toy/noiseless.json is not what the scorer gives its answer")
    for name, kept in held.get(BLIND, {}).items():          # a blind run whose score was written by hand, or whose notes name a rules file
        if kept["score"] != score(kept["answer"], kept["truth"]) or not blind(kept["notes"]):
            wrong.append(f"the blind run {name} in toy/noiseless.json is not what its answer and its notes give")
    return wrong


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["run"]:          # a computation of about an hour, which goes through the queue (his messages 147 and 178): it ran bare on 2026-10-16
        print("REFUSED: the noiseless run is a computation of about an hour: `python3 tools/autoqueue.py add toy/noiseless.py`, the queue runs it as one unit of its own", file=sys.stderr)
        sys.exit(1)
    elif args == ["units"]:          # one unit, the whole loop: it keeps nothing until it ends, so no wave can cut it
        made = (ROOT / "runs" / "noiseless" / digest(ROOT)[:12] / "run.json").exists()
        print(json.dumps([{"condition": "loop", "units": [] if made else [1], "after": [], "kind": "loop"}]))
    elif args[:1] == ["slice"] and len(args) >= 4:
        if not os.environ.get("TOY_QUEUE"):
            print("REFUSED: a slice is a wave of the queue: `python3 tools/autoqueue.py add toy/noiseless.py` (his message 147)", file=sys.stderr)
            sys.exit(1)
        import resource
        began, done = time.time(), []
        if not (ROOT / "runs" / "noiseless" / digest(ROOT)[:12] / "run.json").exists():
            run()
            done = [1]
        print(json.dumps({"done": done, "seconds": [time.time() - began] if done else [], "rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))
    elif args == ["progress"]:
        path = ROOT / "runs" / "noiseless" / digest(ROOT)[:12] / "run.json"
        print(("the noiseless run is kept: " + json.dumps(against(json.loads(path.read_text())))[:250]) if path.exists() else "the noiseless run is not made yet")
    elif args == ["digest"]:
        print(digest(ROOT))
    elif len(args) >= 2 and args[0] == "make" and len([a for a in args[2:] if a != "--recorded"]) % 2 == 0 and set([a for a in args[2:] if a != "--recorded"][::2]) <= {"--noise", "--iterations", "--minutes"}:
        rest = [a for a in args[2:] if a != "--recorded"]
        flags = dict(zip(rest[::2], (int(v) for v in rest[1::2])))
        print(f"written: {make(Path(args[1]), draw=flags.get('--noise'), iterations=flags.get('--iterations', ITERATIONS), minutes=flags.get('--minutes', MINUTES), recorded='--recorded' in args)}")
    elif len(args) in (4, 6) and args[0] == "keep" and args[4:5] in ([], ["--noise"]):
        print(json.dumps(keep(Path(args[1]), int(args[2]), float(args[3]), draw=int(args[5]) if len(args) == 6 else None)["score"], indent=1))
    elif len(args) in (3, 5) and args[0] == "run-report" and args[3:4] in ([], ["--brief"]):
        print(run_report(Path(args[1]), Path(args[2]), brief=Path(args[4]) if len(args) == 5 else None))
    elif len(args) == 3 and args[0] == "loop-report":          # the loop of the package on the toy, draw 1 at its noise: the parameters of the run of the folder, else those that toy/chosen.py counted
        ran = Path(args[1]) / "run.json"
        held = json.loads(ran.read_text())["parameters"] if ran.exists() else json.loads((ROOT / "toy" / "chosen.json").read_text())["report"]["one"].get("parameters")
        print(loop_report(Path(args[1]), Path(args[2]), parameters=held, title="draw 1 at the toy's noise"))
    elif args[:1] == ["run-slides"] and len([a for a in args if a != "--visual"]) in (3, 5) and (len([a for a in args if a != "--visual"]) == 3 or [a for a in args if a != "--visual"][3] == "--recorded"):          # a blind run as slides, his message 239, with --recorded the same run on his recorded traces beside it, his message 288, and with --visual its figures alone, his message 297
        plain = [a for a in args if a != "--visual"]
        print("\n".join(str(p) for p in run_slides(Path(plain[1]), Path(plain[2]), recorded=Path(plain[4]) if len(plain) == 5 else None, visual="--visual" in args)))
    elif len(args) == 3 and args[0] == "report":
        print("\n".join(str(p) for p in steps_report(Path(args[1]), Path(args[2]))))
    elif len(args) == 3 and args[0] == "keep-blind":
        print(json.dumps(keep_blind(Path(args[1]), args[2])["score"], indent=1))
    elif len(args) == 2 and args[0] == "score":
        const, config = model.load()
        print(json.dumps(score(json.loads((Path(args[1]) / "answer.json").read_text()), truth(const, config)), indent=1))
    elif args == ["table"]:
        print(f"written: {table().relative_to(ROOT)}")
    elif args == ["--check"]:
        wrong = check()
        print("\n".join(wrong) or "the noiseless run: made from the files as they are")
        sys.exit(1 if wrong else 0)
    else:
        print(__doc__)
        sys.exit(1)
