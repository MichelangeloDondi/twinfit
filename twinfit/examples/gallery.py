"""The examples of the package, the owner's message 188: "a few synthetic examples from different fields [...] shown in
the package. they should be useful either for the researchers but also for the agent of a fresh chat to get inspired
about what to do, what not to do, to better understand what is it the supposed outcome".

Projects of other fields, each a folder that `python -m twinfit example NAME DIR` writes. Its truth is known and the
folder does not hold it. Six are the projects that a fresh agent worked with no hint in a round of the cold tests of the
package, their curves as that agent saw them, with the setup of their data collection added, and the others were written
for the gallery, so that no agent has worked them: `source` of each says which. For each, the gallery runs the
loop as a researcher would: the first model alone and what its report says first; then what to do, the move that the
report asks for; then, in a copy, what not to do, the move that a tester made and the report does not ask for. What each
run gave is kept beside the truth in gallery.json, and the page twinfit/EXAMPLES.md is drawn from it, so that no number
of the page is written by hand.

    python -m twinfit.examples.gallery --write    runs every example and writes gallery.json and EXAMPLES.md
    python -m twinfit.examples.gallery --check    says whether EXAMPLES.md is the page that gallery.json gives
    python -m twinfit.examples.gallery NAME DIR   runs one example in DIR and prints what it gave
"""
from __future__ import annotations

import json
import math
import shutil
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECTS = HERE / "projects"
KEPT = HERE / "gallery.json"
PAGE = HERE.parent / "EXAMPLES.md"
SETS = 30          # the synthetic sets of each twin of a run of the gallery
LEVELS = (1.0, 3.0, 0.0)          # his message 194: each project at its own noise, at three times it, and with no noise drawn
RHO = 0.8          # the correlation of neighbouring samples of a noise of the kind "correlated", as the filter of a detector gives it
KINDS = ("white", "correlated", "counts")          # the kinds of noise that a project draws: independent and normal, normal and correlated over neighbouring samples, or that of a counter
WORKERS = 6          # the examples run side by side by --write

SATURATION_MODEL = '''"""A model of one's own: the width of the line grows with the power by the law of saturation, its two numbers shared."""
import numpy as np

SHARED = ("width at no power", "saturation power")      # the full width at half maximum at no power, in MHz, and the power of saturation, in mW
OWN = ("height", "place", "background")                 # of each record
START_SHARED = [1.5, 30.0]
BOUNDS_SHARED = ([0.05, 1.0], [20.0, 1000.0])
CENTRE = "place"
NAME = "one Lorentzian line whose width saturates with the power"


def model(x, shared, own, design):
    width = shared[0] * np.sqrt(1.0 + design["power"] / shared[1])
    return own[0] / (1.0 + (2.0 * (x - own[1]) / width) ** 2) + own[2]


def start_own(record):
    base = float(np.median(record.y))
    return [float(record.y.max() - base), float(record.x[np.argmax(record.y)]), base]


def centre(shared, own, design):
    return own[1]
'''

FANO_MODEL = '''"""A model of one's own: a Fano profile, the line of a discrete state that interferes with a continuum, its half width
and its asymmetry q shared."""
import numpy as np

SHARED = ("half_width", "q")                 # the half width of the resonance, in meV, and the asymmetry of the profile
OWN = ("height", "place", "background")      # of each record
START_SHARED = [1.5, 2.0]
BOUNDS_SHARED = ([0.05, 0.1], [20.0, 50.0])
CENTRE = "place"
NAME = "one Fano profile"


def model(x, shared, own, design):
    e = (x - own[1]) / shared[0]
    q = shared[1]
    return own[0] * ((q + e) ** 2 / (1.0 + e ** 2)) / (q * q + 1.0) + own[2]


def start_own(record):
    base = float(np.median(record.y))
    top = int(np.argmax(np.abs(record.y - base)))
    return [float(record.y[top] - base), float(record.x[top]), base]


def centre(shared, own, design):
    return own[1]
'''

CLOCK_MODEL = '''"""A model of one's own: the centre of the line moves with the power of the probe, the light shift, from the frequency of
the clock at no power, which the model then measures."""
import numpy as np

SHARED = ("width", "clock frequency", "light shift")      # the full width, in Hz, the centre at no power less the nominal frequency, in Hz, and its shift, in Hz for each microwatt
OWN = ("height", "background")                           # of each record
START_SHARED = [1.0, 0.0, 0.0]
BOUNDS_SHARED = ([0.05, -3.0, -1.0], [10.0, 3.0, 1.0])
HELD = {"power scale": 1.0}                              # the powers read times this, calibrated within 2 per cent, which description.toml holds
NAME = "one Lorentzian line whose centre moves with the power"


def model(x, shared, own, design):
    centre = shared[1] + shared[2] * design["power"] * HELD["power scale"]
    return own[0] / (1.0 + (2.0 * (x - centre) / shared[0]) ** 2) + own[1]


def start_own(record):
    base = float(np.median(record.y))
    return [float(record.y.max() - base), base]

'''

OWN_CENTRE_MODEL = '''"""A centre of its own for each curve: each line fitted where it stands."""
import numpy as np

SHARED = ("width",)                         # the full width at half maximum, in Hz
OWN = ("height", "centre", "background")    # of each record
START_SHARED = [1.0]
BOUNDS_SHARED = ([0.05], [10.0])
HELD = {"power scale": 1.0}                 # held by description.toml, and read by no line of this model
CENTRE = "centre"
NAME = "one Lorentzian line of one width, its centre its own in each curve"


def model(x, shared, own, design):
    return own[0] / (1.0 + (2.0 * (x - own[1]) / shared[0]) ** 2) + own[2]


def start_own(record):
    base = float(np.median(record.y))
    return [float(record.y.max() - base), float(record.x[np.argmax(record.y)]), base]


def centre(shared, own, design):
    return own[1]
'''

TRANSIT_MODEL = '''"""A model of one's own: the dip of a star darker at its limb, by a coefficient that falls in proportion to the wavelength."""
import numpy as np

SHARED = ("half_duration", "limb darkening at 450 nm", "its fall to 800 nm")      # hours, and the coefficient u of the limb and its fall from 450 to 800 nm
OWN = ("depth", "centre", "baseline")                                            # of each record
START_SHARED = [1.0, 0.4, 0.1]
BOUNDS_SHARED = ([0.2, -1.0, -1.0], [2.5, 1.5, 1.5])
CENTRE = "centre"
NAME = "the dip of a star darker at its limb"
EDGE = 0.05                                 # the time over which the planet enters the disc of the star, in hours


def model(x, shared, own, design):
    u = shared[1] - shared[2] * (design["wavelength"] - 450.0) / 350.0
    t = x - own[1]
    inside = 1.0 / (1.0 + np.exp(-(shared[0] - np.abs(t)) / EDGE))
    rim = 1.0 - (t / shared[0]) ** 2
    mu = np.sqrt(0.5 * (rim + np.sqrt(rim ** 2 + 1e-4)))          # the cosine of the angle at the limb, smooth at the edge of the star, where a root clipped at zero has a corner that a search cannot cross
    return own[2] - own[0] * inside * (1.0 - u * (1.0 - mu)) / (1.0 - u / 3.0)


def start_own(record):
    base = float(np.median(record.y))
    return [float(base - record.y.min()), float(record.x[np.argmin(record.y)]), base]


def centre(shared, own, design):
    return own[1]
'''

EXAMPLES = {
    "pedestal": {
        "source": "cold",
        "field": "absorption lines of a gas cell at five pressures",
        "truth": {"width": 4.0, "pedestal 6": 0.04},
        "do": {"words": "scan the stretch of the shape that the report names first, a wider pedestal under the line, and write the term at the top of the scan", "scan": ("pedestal(x)", [3.0, 4.0, 6.0, 8.0, 12.0])},
        "dont": {"words": "make the width change along the pressure, a setting that the report does not name", "term": 'term = mechanisms.along("width", "pressure")'}},
    "satellite": {
        "source": "cold",
        "field": "absorption lines of a gas cell with a weak line beside each",
        "truth": {"width": 4.0, "second component at 12": 0.08},
        "do": {"words": "scan the distance of a second component, `python -m twinfit scan DIR \"second(x)\" ...`, and write the term at the top of the scan", "scan": ("second(x)", [6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0])},
        "dont": {"words": "widen the line with a kernel, as if the satellite were wings of the line", "term": "term = mechanisms.widened(lambda design: 1.0)"}},
    "saturation": {
        "source": "cold",
        "field": "fluorescence of a resonance driven at five powers of a laser",
        "truth": {"width at no power": 1.0, "saturation power": 50.0},
        "do": {"words": "write the law that the physics gives, the width saturating with the power, as shared parameters that read the power, in model.py", "model": SATURATION_MODEL},
        "dont": {"words": "a width that grows in proportion to the power, the generic term where the physics gives a law", "term": 'term = mechanisms.along("width", "power")'}},
    "filter": {
        "source": "cold",
        "kinds": ("correlated", "white"),
        "white truth": "truth: the model is whole, width 1.2 pm, and the noise is white: nothing is missing",
        "field": "absorption dips recorded behind the low-pass filter of a detector",
        "truth": {"width": 1.2},
        "do": {"words": "read the correlated noise as noise: keep one sample in every few with `python -m twinfit thin DIR N`, at the span that the report gives, and run the loop on the copy", "thin": 15},
        "dont": {"words": "hold the correlation with a term of the curve, a pedestal", "term": "term = mechanisms.pedestal(6.0)"}},
    "fano": {
        "source": "cold",
        "field": "photoionization across an autoionizing resonance at four pressures",
        "truth": {"half_width": 2.0, "q": 3.0},
        "do": {"words": "the asymmetry is the line's own shape, which no term of the package holds: write the profile that the physics gives, a Fano profile, as the model", "model": FANO_MODEL},
        "dont": {"words": "a tail on one side, the generic term of an asymmetry", "term": "term = mechanisms.tail(2.0, 'upper')"}},
    "counts": {
        "source": "cold",
        "kinds": ("counts", "white"),
        "white truth": "truth: the model is whole, half width 1.5 GHz, and the noise white, 6.3 at every sample as stated: nothing is missing",
        "field": "photon counts of an emission line at four powers of the pump",
        "truth": {"half_width": 1.5},
        "do": {"words": "give each sample the noise of a counter, the root of the fitted counts, with `python -m twinfit counts DIR`, and run the loop on the copy", "counts": True},
        "dont": {"words": "hold the misfit at the top of the lines with a term of the curve, a pedestal", "term": "term = mechanisms.pedestal(6.0)"}},
    "clock": {
        "source": "gallery",
        "field": "the resonance of an atomic clock read at five powers of the probe",
        "truth": {"width": 1.5, "clock frequency": 0.12, "light shift": 0.025},
        "do": {"words": "write the light shift that the physics gives, the centre moving in proportion to the power read through its calibration, which the description holds within 2 per cent, so that the model measures the frequency at no power and the light shift carries the bar of that calibration", "model": CLOCK_MODEL},
        "dont": {"words": "free the centre of each curve: the fit is right, and the frequency of the clock is then in no parameter, the shift hidden in five centres", "model": OWN_CENTRE_MODEL}},
    "gel": {
        "source": "gallery",
        "field": "a band of electrophoresis in a gel at four voltages",
        "truth": {"width": 1.0, "tail lower 3": 0.15},
        "do": {"words": "scan the stretch of the tail on the side of slower migration that the report points to, and write the term at the top of the scan", "scan": ("tail(x, 'lower')", [2.0, 3.0, 4.0, 5.0, 6.0])},
        "dont": {"words": "a pedestal under the band, wider on both sides alike", "term": "term = mechanisms.pedestal(3.0)"}},
    "transit": {
        "source": "gallery",
        "field": "the transit of an exoplanet seen through four filters",
        "truth": {"half_duration": 1.2, "limb darkening at 450 nm": 0.6, "its fall to 800 nm": 0.3},
        "do": {"words": "write the darkening of the limb of the star that the physics gives, its coefficient falling with the wavelength, as the model", "model": TRANSIT_MODEL},
        "dont": {"words": "a curvature of the baseline under every curve, the generic term of a round bottom", "term": "term = mechanisms.level(2)"}},
}


def kinds_of(name: str) -> tuple:
    """The two kinds of noise of a project: its own, then a second."""
    return EXAMPLES[name].get("kinds", ("white", "correlated"))


def drawing(level: float = 1.0, kind: str | None = None) -> dict:
    """What a project draws its noise with, his message 194: "for the synthetic projects for each one keep at least two
    noise levels and two different kinds of noise (maybe noiseless too?)". `level` multiplies the project's own noise, 0
    draws none, and `kind`, where given, takes the place of the project's own kind. At level 1 and its own kind a
    project draws the very numbers that it drew before: the normal numbers are drawn whatever the level, and scaled. The
    description states the noise at its level, and with no noise drawn the noise of level 1, since a noise is stated
    above zero."""
    if not math.isfinite(level) or level < 0 or (kind is not None and kind not in KINDS):          # nan and inf wrote traces of nan and a noise of inf (review 29, minor 9)
        raise ValueError(f"--noise takes a multiple of the example's own noise, 0 or more, and --kind one of {', '.join(KINDS)}, and they were given {level:g} and {kind}")          # "1.0, counter" read as one value (round 20 of the tests of the tutorial)

    def noise(rng, sd, size, own="white"):
        made = kind or own
        if made == "white":
            return level * sd * rng.standard_normal(size)
        if made == "correlated":          # each sample keeps RHO of the last, its spread that of one sample
            e = np.empty(size)
            e[0] = rng.normal(0.0, sd)
            for i in range(1, size):
                e[i] = RHO * e[i - 1] + rng.normal(0.0, sd * (1.0 - RHO ** 2) ** 0.5)
            return level * e
        raise ValueError(f"a noise of the kind {made!r} is drawn by counted, on the counts themselves")

    def counted(rng, mean, sd):          # the counts of a counter at the level: a variance of level squared times the counts
        if kind in (None, "counts"):
            return np.array(mean, dtype=float) if level == 0 else level ** 2 * rng.poisson(np.asarray(mean) / level ** 2).astype(float)
        return np.asarray(mean, dtype=float) + noise(rng, sd, np.size(mean), "white")

    def stated(sd):
        return level * sd if level > 0 else sd
    return {"noise": noise, "counted": counted, "stated": stated}


def told_of(name: str, level: float, kind: str | None, truth: str) -> str:
    """The truth of a project at a level and a kind of noise, in words."""
    own = kinds_of(name)[0]
    if kind is not None and kind != own:
        truth = EXAMPLES[name].get(f"{kind} truth") or f"{truth}, and the noise is correlated, {RHO} from one sample to the next, which its description does not say"
    return truth + ("" if level == 1 else ", with no noise drawn" if level == 0 else f", at {level:g} times its noise")


def write(name: str, folder, level: float = 1.0, kind: str | None = None) -> str:
    """Writes the project `name` into `folder`, which must not hold one, at a level and a kind of noise, and returns its
    truth in words."""
    if name not in EXAMPLES:
        raise ValueError(f"no example {name!r}: the examples are {', '.join(EXAMPLES)}")
    folder = Path(folder)
    if (folder / "description.toml").exists():
        raise ValueError(f"{folder} holds a project already")
    space = dict(drawing(level, kind), home=folder, __name__=f"twinfit.examples.projects.{name}")
    exec(compile((PROJECTS / f"{name}.py").read_text(), str(PROJECTS / f"{name}.py"), "exec"), space)
    return told_of(name, level, kind, space["TRUTH"])


def _act(home: Path, act: dict) -> Path:
    """Applies a move to the project in `home` and returns the folder that the loop then runs on: a term file, a model of
    one's own, a scan whose top is written as a term, or the copy that thin or counts writes."""
    from .. import loop
    from ..__main__ import counts, thin
    if "term" in act:
        act["typed"] = act["term"]
        (home / "terms" / "01.py").write_text(f"from twinfit import mechanisms\n{act['term']}\nwhy = {json.dumps(act.get('words', 'the move of the gallery for this example'))}\n")
        return home
    if "model" in act:
        act["typed"] = "the model.py of the gallery"
        (home / "model.py").write_text(act["model"])
        return home
    if "scan" in act:          # as the command of the scan says: its top through three neighbours, and a finer scan about it where it asks for one. Until round 19 the best value of the grid was written, which no researcher is told to write
        expression, values = act["scan"]
        got = loop.scan(home, expression, values)
        said = [f"`python -m twinfit scan DIR \"{expression}\" {' '.join(format(v, 'g') for v in values)}` puts its top at {got['best']['x']:.4g}"]
        lo, hi = got["within_one"]
        step = hi - lo if hi > lo else 0.5 * min(abs(b - a) for a, b in zip(values, values[1:]))          # the spacing of the range that the command gives, kept: a spacing that shrank at each scan walked down an edge and never reached the top
        for _ in range(4):          # finer about the top, until the top stands inside the values scanned, as the command asks
            if not got.get("from_a_parabola") and not got.get("at_the_edge"):
                break
            values = [round(got["best"]["x"] + k * step, 6) for k in (-3, -2, -1, 0, 1, 2, 3)]
            got = loop.scan(home, expression, values)
            said.append(f"a finer scan about it, at {' '.join(format(v, 'g') for v in values)}, at {got['best']['x']:.4g}")
            if not got.get("at_the_edge"):
                break
        act["done"] = ", and ".join(said) + f": the term file holds `{got['best']['line']}`"
        act["typed"], act["range"] = got["best"]["line"], tuple(float(v) for v in got["within_one"])          # the term of each row and the range of its scan, which the page did not carry (round 20 of the tests of the tutorial)
        (home / "terms" / "01.py").write_text(f"from twinfit import mechanisms\n{got['best']['line']}\nwhy = {json.dumps(act.get('words', 'the top of the scan'))}\n")
        return home
    if "thin" in act:
        copies = thin(home, act["thin"])
        act["typed"] = f"python -m twinfit thin DIR {act['thin']}"
        act["done"] = f"`python -m twinfit thin DIR {act['thin']}` writes {len(copies)} copies, and the loop runs on the first, {Path(copies[0]).name}"          # which N and which copy: a newcomer could not tell (round 19)
        return copies[0]
    if "counts" in act:
        act["typed"] = "python -m twinfit counts DIR"
        act["done"] = "`python -m twinfit counts DIR` writes the copy, and the loop runs on it"
        return counts(home)
    raise ValueError(f"no move {act}")


def _read(folder: Path, truth: dict) -> dict:
    """What a run of the loop in `folder` gave: the best model, its terms, the chi square for each degree of freedom, the
    first answers of the last report that stand beyond its twin, and each parameter of the truth with its bar."""
    best = json.loads((folder / "loop" / "best.json").read_text())
    reports = sorted((folder / "loop").glob("model_*.json"), key=lambda p: int(p.stem.split("_")[1]))
    last = json.loads(reports[-1].read_text())["where"] if reports else {"answers": []}
    beyond = [f"{a['question']}: {a['answer']}" + (f", excess {a['excess']:.3g}" if a.get("excess") is not None else "")
              for a in last.get("answers", []) if a.get("beyond_the_twin") and a["question"] in ("where", "what", "along", "changes along")][:3]          # names alone, where a newcomer asked how far each stood past its line (round 20 of the tests of the tutorial)
    if beyond and last.get("largest_excess_on_the_twin") is not None:
        beyond.append(f"its line {last['largest_excess_on_the_twin']:.3g}")
    found = {}
    params = best.get("parameters") or {}
    for name in truth:
        key = name if name in params else next((k for k in params if k.split(" ")[:-1] == name.split(" ")[:-1] and name.split(" ")[:-1]), None)          # a term written at the top of a scan, "pedestal 6.7" for the truth's "pedestal 6"
        c = params.get(key) if key else None
        found[name] = None if c is None else {"value": c["value"], "bar": c["bar_of_the_fit"], "coverage": c.get("coverage_on_the_twin"), "as": key, "quotable": (best.get("quotable") or {}).get(key),
                                              "theory": ((c.get("bar_of_theory") or {}).get("together") or None) if isinstance(c.get("bar_of_theory"), dict) else None}          # the bar that a value held by an ansatz puts on it, as the clock's calibration of the power (round 21 of the tests of the tutorial)          # the output's own verdict beside the value: with a correlated noise it refuses the bars that the page showed bare (his message 194)
    first = (json.loads(reports[0].read_text()).get("proposals") or {}) if reports else {}
    model = [f"{m['term']}" + (f", `{m['as_a_term_file']}`" if m.get("as_a_term_file") else "") for m in (first.get("for_the_model") or [])[:3]]
    experiment = [f"{e['lever']}, {e['his_kind']}" + (f": the error of the known parameters times {e['error_over_that_under_the_records']:.2f}" if e.get("error_over_that_under_the_records") is not None else "") for e in (first.get("for_the_experiment") or [])[:3]]
    return {"best": best.get("best"), "terms": best.get("terms", []), "chi_square_per_degree": best.get("chi_square_per_degree"), "beyond": beyond,
            "parameters": found, "refused": [r.get("term") for r in best.get("tried_and_refused", [])], "by_a_hair": any(r.get("by_a_hair") for r in best.get("tried_and_refused", [])),
            "proposals": {"model": model, "experiment": experiment}}


def run(name: str, folder, sets: int = SETS) -> dict:
    """One example as a researcher runs it, in `folder`: the first model, then what to do, then in a copy what not to do."""
    from .. import loop
    began = time.time()
    folder = Path(folder)
    spec = EXAMPLES[name]
    truth = write(name, folder / "first")
    loop.run(folder / "first", n=0, sets=sets)
    out = {"field": spec["field"], "truth_in_words": truth, "truth": spec["truth"], "first": _read(folder / "first", spec["truth"])}
    for move in ("do", "dont"):
        act = dict(spec[move])
        if spec[move].get("counts"):          # counts reads the fitted values of the last run, so it writes its copy beside the folder that the first model ran in
            there = _act(folder / "first", act)
        else:
            home = folder / move
            shutil.copytree(folder / "first", home, ignore=shutil.ignore_patterns("loop", "loop.before", ".twinfit.lock"))
            there = _act(home, act)
        loop.run(there, n=1, sets=sets)
        out[move] = dict(_read(there, spec["truth"]), words=spec[move]["words"], done=act.get("done"), sets=sets)
    out["variants"] = [variant(name, folder, level, kind, sets) for level, kind in variants_of(name)]
    out["minutes"] = round((time.time() - began) / 60.0, 1)
    return out


def variants_of(name: str) -> list:
    """The other noises of a project, his message 194: three times its own, none drawn, and its second kind at its level
    and at three times it."""
    second = kinds_of(name)[1]
    return [(3.0, None), (0.0, None), (1.0, second), (3.0, second)]


def variant(name: str, folder: Path, level: float, kind: str | None, sets: int = SETS) -> dict:
    """The project at another noise, in a folder of its own: the first model alone, then what to do, as the main run
    makes them, so that a researcher reads how the first report and the outcome move with the noise."""
    from .. import loop
    spec = EXAMPLES[name]
    here = Path(folder) / f"noise_{level:g}_{kind or kinds_of(name)[0]}"
    truth = write(name, here / "first", level, kind)
    loop.run(here / "first", n=0, sets=sets)
    act = dict(spec["do"])
    if act.get("counts"):
        there = _act(here / "first", act)
    else:
        shutil.copytree(here / "first", here / "do", ignore=shutil.ignore_patterns("loop", "loop.before", ".twinfit.lock"))
        there = _act(here / "do", act)
    loop.run(there, n=1, sets=sets)
    return {"level": level, "kind": kind or kinds_of(name)[0], "truth_in_words": truth, "first": _read(here / "first", spec["truth"]), "do": dict(_read(there, spec["truth"]), done=act.get("done"), typed=act.get("typed"), range=act.get("range"))}


KIND_WORDS = {"white": "white", "correlated": "correlated over neighbouring samples", "counts": "that of a counter"}


def said_noise(level: float, kind: str) -> str:
    """A level and a kind of noise in words."""
    return "no noise drawn" if level == 0 else f"{KIND_WORDS[kind]}, " + ("as large as its own noise" if level == 1 else f"{level:g} times as large as its own noise")          # "its level" read as the baseline, which every report calls the level of a curve (round 20 of the tests of the tutorial)


def _value(cell, name: str = "") -> str:
    if not cell:
        return "not in the model"
    named = (f", as {cell['as']}" if cell.get("as") and cell["as"] != name else "") + (", not yet to be quoted" if cell.get("quotable") is False else "")          # the parameter of a term written at a scanned value, beside the truth's name, and the output's verdict
    if not cell.get("bar"):
        return f"{cell['value']:.4g}{named}"
    from ..loop import _pm
    theory = f" ± {_pm(cell['value'], cell['theory'])[1]} of theory" if (cell.get("theory") or 0.0) >= 0.01 * cell["bar"] else ""          # the bar of a held value beside that of the fit, as best.md quotes them, where it is not lost in the bar of the fit: the clock frequency at no power read 7e-9 from the calibration of the power
    return " ± ".join(_pm(cell["value"], cell["bar"])) + theory + named


def told(cold: list, written: list) -> str:
    """Which projects a fresh agent worked first, blind, and which were written for the gallery and
    worked by no agent: until the gallery held nine, its head said that six were, each worked by an agent."""
    names = lambda some: ", ".join(some[:-1]) + " and " + some[-1] if len(some) > 1 else "".join(some)          # noqa: E731
    said = ([f"{names(cold)} were each first worked by a fresh agent that was not told the truth, before they joined the gallery"] if cold else []) \
        + ([f"{names(written)} were written for the gallery and worked by no agent yet"] if written else [])
    return ("The projects " + ", and ".join(said) + ".") if said else ""


def page(kept: dict) -> str:
    """The page twinfit/EXAMPLES.md, drawn from what the runs gave."""
    lines = ["# Examples of the package", "",
             f"{len(kept)} projects of other fields, each written by `python -m twinfit example NAME DIR`, each drawn with a truth that its folder does not hold. "
             + told([n for n in kept if EXAMPLES.get(n, {}).get("source") == "cold"], [n for n in kept if EXAMPLES.get(n, {}).get("source") != "cold"]) + " "
             "For each the loop ran three times: on the first model alone, after what to do, and after what not to do. Then at four other noises, on the first model and after what to do: three times its own, none drawn, and a second kind of noise as large as its own and three times as large. "
             "`python -m twinfit example NAME DIR --noise L --kind K` writes each: L is a multiple of the example's own noise, 0 for none drawn, and K one of `white`, `correlated` and `counts`, the example's own kind where K is not given. "
             "Every number below is read from those runs, which `python -m twinfit.examples.gallery --write` makes again for those who build the package.", "",
             "Find the example whose first report leads as yours does, and read what to do there. Each answer of a report is given with its excess, 0 at the median of its statistic over the synthetic sets of the twin and 1 at their 95th percentile, "
             "and each report with its line, the excess that the largest answer of a synthetic set passes on 5 sets of 100: an answer past it stands beyond the twin. "
             "The excess says how clearly a term shows against the noise, and not how large the term is.", "",
             "| example | the field | what the first report puts first beyond its twin |", "|---|---|---|"]
    for name, e in kept.items():
        lines.append(f"| {name} | {e['field']} | {' · '.join(e['first']['beyond']) or 'nothing'} |")
    lines += ["", "## How to read the tables at other noises", "",
              "Each row is the project written with `--noise L --kind K`, its first model run, and then what to do, whose term the third column gives. For a term written at the top of a scan, it gives the range of the scan too: the values whose fall is within 1 of its top.", "",
              "At L times its noise the description states L times the noise, so the bars grow about L times. With none drawn it states the example's own noise: the bars stay those of the own row, the chi square falls near 0, "
              "and the check of the noise says that the noise is stated too large, which is why those values are not yet to be quoted.", "",
              "The rows of white and of correlated noise draw one stream of random numbers, scaled: a row at three times its noise holds the draw of the row at its own noise, three times larger. "
              "Their chi squares are then equal, unless the project also draws something of a fixed size, as the glitches of the saturation, and the two rows are not independent experiments.", "",          # the saturation read 1.1558 against 1.0062 (review 29, moderate 8)
              "A counter's noise at three times is drawn as nine times fewer counts, scaled up nine times, so the root of the counts is a third of the noise. What to do for a counter, `python -m twinfit counts DIR`, states that root: "
              "there it states the noise three times too small, and the chi square comes near 9.", "",
              "A parameter of a term written at the top of a scan is named by that term, as \"as pedestal 8\": its size is read at that stretch, distance or period, and a size at another one is another quantity. "
              "So the sizes of two rows are compared at one stretch only. Where the range of the scan is wide, the records do not fix the stretch, and the bar of the size, read at one stretch, is too small."]
    for name, e in kept.items():
        lines += ["", f"## {name}: {e['field']}", "", f"`python -m twinfit example {name} DIR`. The truth, which the folder does not hold: {e['truth_in_words'].removeprefix('truth: ')}.", ""]
        lines += ["| | the first model | after what to do | after what not to do |", "|---|---|---|---|",
                  f"| the best model | {e['first']['best']} | {e['do']['best']} | {e['dont']['best']} |",
                  f"| its terms | {', '.join(e['first']['terms']) or 'none'} | {', '.join(e['do']['terms']) or 'none'} | {', '.join(e['dont']['terms']) or 'none'} |",
                  f"| chi square for each degree of freedom | {e['first']['chi_square_per_degree']:.3f} | {e['do']['chi_square_per_degree']:.3f} | {e['dont']['chi_square_per_degree']:.3f} |",
                  f"| the last report, first beyond its twin | {' · '.join(e['first']['beyond']) or 'nothing'} | {' · '.join(e['do']['beyond']) or 'nothing'} | {' · '.join(e['dont']['beyond']) or 'nothing'} |"]
        for p, true in e["truth"].items():
            lines.append(f"| {p}, true {true:g} | {_value(e['first']['parameters'].get(p), p)} | {_value(e['do']['parameters'].get(p), p)} | {_value(e['dont']['parameters'].get(p), p)} |")
        props = e["first"].get("proposals") or {}
        if props.get("model") or props.get("experiment"):          # the end of an iteration, the owner's messages 142 and 191: three proposals for the model and three for the experiment, each measured on the twin
            lines += ["", "The first report proposes for the model: " + " · ".join(props.get("model") or ["nothing"]) + ". And for the experiment: " + " · ".join(props.get("experiment") or ["nothing"]) + "."]
        lines += ["", f"**What to do**: {e['do']['words']}." + (f" Here {e['do']['done']}." if e["do"].get("done") else "")]
        if EXAMPLES.get(name, {}).get("do", {}).get("model"):          # the model that what to do writes, which the page named and did not show (round 21 of the tests of the tutorial)
            lines += ["", "The `model.py` that it writes in the place of the first:", "", "```python", EXAMPLES[name]["do"]["model"].rstrip("\n"), "```"]
        lines += ["", f"**What not to do**: {e['dont']['words']}."]
        if e["dont"]["refused"]:          # on how many sets, and by how much: a refusal by a hair at 30 sets was kept at 100, and the page said only that the loop refused it (round 19). Under the table at other noises it read as a row of that table (round 20)
            lines[-1] += (f" The loop refused it on {e['dont'].get('sets', SETS)} synthetic sets: {', '.join(e['dont']['refused'])}"
                          + (", by a hair, which the loop itself asks to run again at 500 sets before the term is set aside: at more sets it may be kept, and that a wrong term is kept says that the noise check and the chi square, not the twin alone, judge it." if e["dont"].get("by_a_hair") else "."))
        if e.get("variants"):          # the same project at other noises, his message 194
            names = list(e["truth"])
            own, second = kinds_of(name) if name in EXAMPLES else ("white", "correlated")
            lines += ["", f"**At other noises**: the first model, then what to do. Its own noise is {KIND_WORDS[own]}, and its second kind {KIND_WORDS[second]}. "
                      f"The rows are written by `--noise 3`, `--noise 0`, `--kind {second}` and `--noise 3 --kind {second}`.", "",
                      "| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | " + " | ".join(f"{p}, true {e['truth'][p]:g}" for p in names) + " |",
                      "|---|---|---|---|" + "---|" * len(names)]
            for v in e["variants"]:
                typed = (f", `{v['do']['typed']}`" if v["do"].get("typed") else "") + (f", its scan within 1 of its top from {v['do']['range'][0]:.3g} to {v['do']['range'][1]:.3g}" if v["do"].get("range") else "")
                lines.append(f"| {said_noise(v['level'], v['kind'])} | {' · '.join(v['first']['beyond']) or 'nothing'} | {v['do']['best']}, {v['do']['chi_square_per_degree']:.3f}{typed} | {' · '.join(v['do']['beyond']) or 'nothing'} | "
                             + " | ".join(_value(v["do"]["parameters"].get(p), p) for p in names) + " |")
    return "\n".join(lines) + "\n"


def main(argv: list) -> int:
    if argv[1:] == ["--check"]:
        same = KEPT.exists() and PAGE.exists() and PAGE.read_text() == page(json.loads(KEPT.read_text()))
        print("twinfit/EXAMPLES.md is the page that gallery.json gives" if same else "REFUSED: twinfit/EXAMPLES.md is not the page that gallery.json gives: python -m twinfit.examples.gallery --write")
        return 0 if same else 1
    if argv[1:] == ["--write"]:
        import tempfile
        from concurrent.futures import ProcessPoolExecutor
        with tempfile.TemporaryDirectory() as tmp, ProcessPoolExecutor(min(WORKERS, len(EXAMPLES))) as pool:          # the examples side by side, one process each
            got = dict(zip(EXAMPLES, pool.map(run, list(EXAMPLES), [Path(tmp) / name for name in EXAMPLES])))
        kept = {name: got[name] for name in EXAMPLES}
        for name, e in kept.items():
            print(f"{name}: {e['minutes']} minutes")
        KEPT.write_text(json.dumps(kept, indent=1) + "\n")
        PAGE.write_text(page(kept))
        print(f"written: {KEPT} and {PAGE}")
        return 0
    if len(argv) == 3 and argv[1] in EXAMPLES:
        print(json.dumps(run(argv[1], Path(argv[2])), indent=1))
        return 0
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
