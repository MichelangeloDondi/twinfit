"""The scoreboard of the core feature, rows M150 and M151, his messages 142 and 165: "the core feature of the package
should be to use the machinery of the generalized moments ... to tell to a senior researcher where to search new terms
of the model (in the core? in the wings? in the asymmetry? adding a polynomial? adding a periodic function? changing the
noise model? adding glitches? etc, with a temperature behavior? with a power behavior? with a waist behavior? etc)".

Each question is a set of noiseless traces with something that its model lacks, and the answers that the report of the
package should put first, written in this file before any reading: where, the core, the wings or the asymmetry, of
which kind, and along what. Four groups:

    bench      on the 8 traces of the fifth detail, a shape that the toy does not hold, planted at the size at which the
               fit at the toy's noise leaves a chi square of `SEEN`: a polynomial, a ripple, wings, glitches, a tail, a
               width along the power, a width along the temperature, and nothing. Held out, written before a count
               that read them: ripples of two other periods, a cubic, a noise correlated from sample to sample, a
               noise that grows with the signal, and a dispersive part of the line
    toy        each term of details 2 to 5 hidden alone from the model that lacks it, at its size and at half, and each
               right model, as toy/pointing.py makes the cases
    field      on the curves of the package's own example of another field, twinfit/examples/two_scales.py, each term
               of its list, inside the package's dictionary and outside it, and nothing
    waist      the 8 conditions at three radii of the beam, the transit time hidden

Every question is asked at two weights, `WEIGHTS`. At the toy's noise the answer measures whether the report sees the
term. At a thousandth of it every planted shape stands far beyond the twin, so the answer measures whether the report
names it rightly: that is the core feature itself. The noise model is asked when noise is drawn. A report that puts
nothing beyond the twin is right on no question, but where nothing is missing.

    python3 tools/autoqueue.py add toy/scoreboard.py   the twin of each question at each weight in `PARTS` units, then its report
    python toy/scoreboard.py table                     writes toy/scoreboard.json from the readings kept
    python toy/scoreboard.py --check                   ends without error when toy/scoreboard.json was made from the files as they are

Each question is asked again on `DRAWS` draws of noise on its records, row M259: the report on one set of records is
one draw, and a share of how often it puts the right answer first needs several. A draw is read against the twin of
the noiseless records, made once, so the parts of every twin are kept.
"""
from __future__ import annotations

import copy
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from toy import draws, fits, model, noiseless, pointing, power  # noqa: E402

WEIGHTS = {"law": 1.0, "clean": 1e-3}          # the noise that weights each fit and draws its twin, as a share of the toy's law, or of the stated noise of another field
SAID = {"law": "at the toy's noise", "clean": "at a thousandth of the noise"}
SEEN = 100.0          # the chi square that the fit at the toy's noise leaves of a planted shape, over its traces
MOST = 100.0          # the largest factor by which a size is raised to leave `SEEN`
SETS = 100
PARTS = 10
WAISTS = (35.0, 43.0, 51.0)
QUESTIONS = ("where", "shape", "what", "changes along")          # the questions of the blind report that a key answers: the place, the shape of the residuals, the law of the noise, and which place changes along which setting


def blind_key(where=None, shape=None, what=None, along=None) -> dict:
    """A key of the blind report: the place, the shape, the noise, and the place that changes along a setting."""
    return {"where": where, "shape": shape, "what": what, "changes along": None if along is None else f"{where} along {along}"}


TREND, PERIODIC = "a slow trend over the span of each curve", "a periodic structure along the abscissa"
BENCH = {          # a shape planted on the traces of the fifth detail, and what the blind report should put first, written again from the physics of each shape before any count of the blind report (the red team of his message 343, S1)
    "nothing": blind_key(),
    "polynomial": blind_key("the wings", TREND),
    "ripple": blind_key(None, PERIODIC),
    "wings": blind_key("the wings"),
    "glitches": blind_key(what="glitches"),
    "tail": blind_key("the asymmetry"),
    "power": blind_key("the core", along="P"),
    "temperature": blind_key("the core", along="T"),
    "ripple of 0.7 widths": blind_key(None, PERIODIC),
    "ripple of 3.1 widths": blind_key(None, PERIODIC),
    "cubic": blind_key("the asymmetry", TREND),
    "correlated noise": blind_key(what="a noise correlated from sample to sample"),
    "noise that grows": blind_key(what="a noise that grows with the signal"),
    "dispersive part": blind_key("the asymmetry"),
}
HELD_OUT = ("ripple of 0.7 widths", "ripple of 3.1 widths", "cubic",          # written 2026-10-15 at 20:45, after the repair of version 0.0.4.dev8 and before any reading of it: shapes that the repair never saw (the red team of his message 167)
            "correlated noise", "noise that grows", "dispersive part")          # written 2026-10-23 at 01:44, before the count at 0.1.0rc14: what his message 327 asks and the bench lacked, a change of the noise model and an asymmetry of the line itself (row M439)
NOISE_OF = {"correlated noise": ("correlated", 0.5), "noise that grows": ("grows", 1.0)}          # a question of the noise: the lines hold nothing, and each draw is drawn with neighbouring samples correlated at 0.5, or with a deviation that grows to twice its own at the top of the line
FIELD = {          # a term of the example of another field, inside the package's dictionary or outside it, and what the blind report should put first, from the physics of each term
    "none": dict(blind_key(), inside=True),
    "second component": dict(blind_key("the asymmetry"), inside=True),
    "pedestal": dict(blind_key("the wings"), inside=True),
    "periodic": dict(blind_key(None, PERIODIC), inside=True),
    "scale grows with the square of the covariate": dict(blind_key("the core", along="z"), inside=True),
    "second component off the grid": dict(blind_key("the asymmetry"), inside=False),
    "step of the level": dict(blind_key("the asymmetry"), inside=False),
    "narrow dip": dict(blind_key("the core"), inside=False),
}
WAIST = blind_key("the core", along="waist")          # the transit time is a width in the core, and its rate goes as one over the waist
DRAWS = 5          # the draws of noise on the records of each question, row M259: the slides said of one draw what a share must say
DRAWS_IN = {"field": 20, "held out": 20, "waist": 20}          # the groups asked at 20 draws a question, row M412 and the red team of his message 308, R8: 5 left a share of a group of one question a fifth wide
DRAWN = ("DRAWS", "DRAWS_IN", "NOISE_OF", "draws_of", "drawn", "agreed", "draw_report")          # what makes a draw: its rows hang on these and on the digest of the twins, which do not read them


def questions(config: dict) -> list:
    """Every question of the scoreboard: its name, its group, and the answers that the report should put first, all
    None where nothing is missing."""
    out = [{"name": f"bench: {n}", "group": "held out" if n in HELD_OUT else "bench", "expected": e} for n, e in BENCH.items()]
    for name, made, d, free, sizes, term in pointing.cases(config):
        e = pointing.EXPECTED[term] if term else {"where": None, "setting": None}
        out.append({"name": f"toy: {name}", "group": "toy", "expected": blind_key(e["where"], along=e["setting"]) if term else blind_key()})          # a width of the toy: its place, and that place along its setting
    out += [{"name": f"field: {n}", "group": "field", "expected": {q: e[q] for q in QUESTIONS}, "inside": e["inside"]} for n, e in FIELD.items()]
    out.append({"name": "waist", "group": "waist", "expected": dict(WAIST)})
    return out


def missing(expected: dict) -> bool:
    return any(v is not None for v in expected.values())


def scaled(const: dict, share: float) -> dict:
    """The fixed values of the toy with the noise of each point `share` times its own: the constant term of the law by
    it, the terms of the variance by its square."""
    if share == 1.0:
        return const
    quiet = copy.deepcopy(const)
    for law in quiet["noise"]["trace"]:
        law["a_V"] *= share
        law["b_V"] *= share ** 2
        law["c"] *= share ** 2
    quiet["noise"]["d_MHz2"] *= share ** 2
    return quiet


def bench_lines(const: dict, config: dict, extra: float = 0.0, along: str | None = None) -> list:
    """The 8 noiseless traces of the fifth detail, with the natural width widened by `extra` MHz times the power over
    the highest power, or times the temperature above 70 C over 60 C."""
    values = dict(model.truth(config))
    rows = fits.synthetic(5, const, config, 1)
    if extra:
        out, top = [], max(r["P"] for r in rows)
        for r in rows:
            grow = r["P"] / top if along == "P" else (r["T"] - 70.0) / 60.0
            wide = dict(values, natural_width=values["natural_width"] + extra * grow)
            u = np.asarray(r["nu"]) - r["nu"][int(np.argmax(r["signal"]))]
            base, base_top = model.shape(5, u, r["T"], values, const), model.shape(5, 0.0, r["T"], values, const)
            broad, broad_top = model.shape(5, u, r["T"], wide, const), model.shape(5, 0.0, r["T"], wide, const)
            out.append(dict(r, signal=r["signal"] + float(np.max(r["signal"]) - np.min(r["signal"])) * (broad / broad_top - base / base_top)))
        rows = out
    return [dict(r, y=r["signal"]) for r in rows]


def shape_of(name: str, row: dict) -> np.ndarray:
    """The planted shape of a question of the bench on one trace, in volts for a line of height one."""
    nu, s = np.asarray(row["nu"], dtype=float), np.asarray(row["signal"], dtype=float)
    u = nu - nu[int(np.argmax(s))]
    above = u[s >= s.min() + 0.5 * (s.max() - s.min())]
    width, height = float(above.max() - above.min()), float(s.max() - s.min())
    lorentz = lambda w: 1.0 / (1.0 + (2.0 * u / w) ** 2)          # noqa: E731
    if name == "polynomial":
        return height * (u / np.max(np.abs(u))) ** 2
    if name.startswith("ripple"):
        widths = 1.3 if name == "ripple" else float(name.split()[2])
        return height * np.sin(2.0 * np.pi * u / (widths * width))
    if name == "cubic":
        return height * (u / np.max(np.abs(u))) ** 3
    if name == "wings":
        return height * lorentz(6.0 * width)
    if name == "glitches":
        spikes = np.zeros_like(u)
        spikes[np.random.default_rng([2026, 10, 15, int(row["P"]), int(row["T"])]).choice(u.size, 5, replace=False)] = height
        return spikes
    if name == "tail":
        return height * lorentz(2.0 * width) * 0.5 * (1.0 - np.tanh(u / (0.3 * width)))
    if name == "dispersive part":          # the dispersion profile of the line's own width, odd about its top: a Fano or a dispersive admixture in a spectrum
        return height * (2.0 * u / width) * lorentz(width)
    raise ValueError(f"no shape for {name!r}")


def bench_rows(name: str, size: float, const: dict, config: dict) -> list:
    """The noiseless traces of a question of the bench at a size: for a width, in MHz at the highest power or
    temperature, for another shape the factor that gives it a chi square of one against the toy's law."""
    if name == "nothing" or name in NOISE_OF:
        return bench_lines(const, config)
    if name in ("power", "temperature"):
        return bench_lines(const, config, size, "P" if name == "power" else "T")
    rows = bench_lines(const, config)
    shapes = [shape_of(name, r) for r in rows]
    unit = 1.0 / float(np.sqrt(sum(float(np.sum((s / r["sd"]) ** 2)) for s, r in zip(shapes, rows))))
    return [dict(r, y=r["y"] + size * unit * s) for r, s in zip(rows, shapes)]


def field_records(name: str, size: float, share: float) -> list:
    """The noiseless curves of the example of another field with one of its terms at a size, each weighted at `share`
    times its stated noise."""
    from twinfit import Record
    from twinfit.examples import two_scales as ex
    return [Record(r.x, r.y, ex.NOISE * share, r.design, name=r.name) for r in ex.records(name, size, noise=0.0)]


def waist_rows(const: dict, config: dict) -> list:
    """The 8 noiseless traces of the fifth detail at each radius of `WAISTS`, the rate of the transit time at each in
    proportion to one over the radius."""
    at, true = float(const["beam"]["waist_um"]), float(model.truth(config)["transit"])
    return [dict(r, y=r["signal"], waist=w) for w in WAISTS for r in fits.synthetic(5, const, power.drawn_with(config, {"transit": true * at / w}), 1)]


def at_weight(question: str, size: float, share: float, const: dict, config: dict) -> tuple:
    """The fit of a question at a size, weighted at `share` times the noise: its problem and its fit."""
    from twinfit import noise as N
    from toy import analysis
    group, _, name = question.partition(": ")
    if group == "bench":
        problem, result, _ = analysis.at(5, bench_rows(name, size, const, config), scaled(const, share), config, list(config["detail"]["5"]["free"]))
        return problem, result
    if group == "toy":
        case = next(c for c in pointing.cases(config) if c[0] == name)
        rows = [dict(r, y=r["signal"]) for r in fits.synthetic(case[1], const, power.drawn_with(config, case[4]), 1)]
        problem, result, _ = analysis.at(case[2], rows, scaled(const, share), config, case[3])
        return problem, result
    if group == "field":
        from twinfit import fit
        from twinfit.examples import two_scales as ex
        problem = ex.problem(field_records(name, size, share))
        return problem, fit(problem)
    rows = waist_rows(const, config)
    problem = fits.problem_of(4, rows, const, config, list(config["detail"]["4"]["free"]))
    for record, row in zip(problem.records, rows):
        record.design["waist"] = row["waist"]
        record.name = f"{record.name}, {row['waist']:g} um"
    settled, result, _ = N.settle(problem, fits.variance_of(scaled(const, share)))
    return settled, result


def size_of(question: str, const: dict, config: dict) -> dict:
    """The size at which the fit at the toy's noise leaves `SEEN`, for a planted shape of the bench or a term of another
    field: tried at a first size, then raised by the square root of what is still missing, at most `MOST` times, since
    what is left grows as the square of the size. The other questions take the size of their own case."""
    group, _, name = question.partition(": ")
    if group not in ("bench", "field") or name in ("nothing", "none") or name in NOISE_OF:
        return {"size": 0.0, "stopped": False}
    first = 0.01 if (group == "field" or name in ("power", "temperature")) else float(np.sqrt(SEEN))
    left = float(at_weight(question, first, 1.0, const, config)[1].chi2)
    factor = float(np.sqrt(SEEN / max(left, 1e-12)))
    size = first * min(factor, MOST)
    return {"size": size, "left_at_the_first": left, "left": float(at_weight(question, size, 1.0, const, config)[1].chi2), "stopped": factor > MOST}


def fitted(question: str, weight: str, const: dict, config: dict) -> tuple:
    return at_weight(question, size_of(question, const, config)["size"], WEIGHTS[weight], const, config)


MAKES = ("WEIGHTS", "SEEN", "MOST", "SETS", "PARTS", "WAISTS", "BENCH", "HELD_OUT", "FIELD", "WAIST", "questions", "scaled", "bench_lines", "shape_of", "bench_rows",
         "field_records", "waist_rows", "at_weight", "size_of", "fitted", "twin_part", "report_of")


RUNS = ("toy/scoreboard.py", "toy/analysis.py", "toy/power.py")          # the programs whose imports say which files of the package a reading runs


def digest(root: Path = ROOT) -> str:
    """What a reading depends on: the parts of this program that make it, the cases of toy/pointing.py and the places
    that it expects of each, which it reads, the analysis, the model, its values, the terms, and the package with its
    example of another field."""
    root = Path(root)
    h = hashlib.sha256(power.sources((root / "toy" / "scoreboard.py").read_text(), MAKES))
    h.update(power.sources((root / "toy" / "pointing.py").read_text(), ("cases", "EXPECTED")))
    toy = ["toy/analysis.py", "toy/model.py", "toy/fits.py", "toy/traces.py", "toy/power.py", "toy/const.toml", "toy/config.toml"] + sorted(str(p.relative_to(root)) for p in (root / "toy" / "terms").glob("*.py"))
    draws.hashed(h, root, toy + draws.package_of(root, RUNS))          # the files of the package that a reading imports, as code (rule 16, his message 175)
    return h.hexdigest()


def path_of(root: Path, question: str, weight: str, part: int | None = None) -> Path:
    stem = question.replace(": ", "_").replace(" ", "_").replace(",", "")
    return Path(root) / "runs" / "scoreboard" / digest(root)[:12] / (f"{stem}_{weight}_" + ("report.json" if part is None else f"twin_part{part}.json"))


def twin_part(question: str, weight: str, part: int, const: dict, config: dict, root: Path = ROOT) -> list:
    """One unit of the twin of a question at a weight: the statistics of the report on its sets from SETS / PARTS * (part
    - 1), each fitted again, kept as it ends."""
    from twinfit import twin, where
    from toy import analysis
    problem, result = fitted(question, weight, const, config)
    sets = twin.replicas(problem, result, SETS, seed=23)
    size = SETS // PARTS
    centre = None if question.startswith("field") else analysis.CENTRE
    got = [where._one_set((sets[i], result, 80, centre)) for i in range(size * (part - 1), size * part)]
    path = path_of(root, question, weight, part)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps([None if s is None else {"|".join(k): v for k, v in s.items()} for s in got]))
    return got


def report_of(question: str, weight: str, const: dict, config: dict, root: Path = ROOT) -> dict:
    """The report of a question at a weight against its twin, and its first answers against those written before. A
    report that puts nothing beyond the twin is right on no question, but where nothing is missing."""
    from twinfit import where
    from toy import analysis
    null = [{tuple(k.split("|", 1)): v for k, v in s.items()} for k in range(1, PARTS + 1) for s in json.loads(path_of(root, question, weight, k).read_text()) if s is not None]
    problem, result = fitted(question, weight, const, config)
    report = where.where(problem, result, null=null, centre=None if question.startswith("field") else analysis.CENTRE)
    first = report["first_of_each_question"]
    expected = next(q["expected"] for q in questions(config) if q["name"] == question)
    beyond = [[a["question"], a["answer"], power.said(a["excess"])] for a in report["answers"] if a["beyond_the_twin"]]
    if question.partition(": ")[2] in NOISE_OF:          # noiseless records hold no noise: only the draws ask a question of the noise
        agrees = {}
    elif missing(expected):
        agrees = {q: bool(beyond) and first.get(q) == e for q, e in expected.items() if e is not None}
        if expected["changes along"] is None:          # a setting named where the planted shape changes along none is a false answer (the red team of his message 167)
            agrees["unmoved"] = first.get("changes along") is None
    else:
        agrees = {"silent": not beyond}
    out = {"first": first, "beyond": beyond, "agrees": agrees, "sets_of_the_twin": len(null), "leading": [[a["question"], a["answer"], power.said(a["excess"])] for a in report["answers"][:5]]}
    path = path_of(root, question, weight)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out))
    return out


def draw_digest(root: Path = ROOT) -> str:
    """What the reading of a draw depends on: the digest of the twins, and the parts of this program that make a draw,
    which `digest` does not read, so that a change of a draw makes no twin again."""
    h = hashlib.sha256(digest(root).encode())
    h.update(power.sources((Path(root) / "toy" / "scoreboard.py").read_text(), DRAWN))
    return h.hexdigest()


def draw_path(root: Path, question: str, weight: str, draw: int, made: tuple | None = None) -> Path:
    """Where a draw is kept. `made` is the pair of the digests of the twins and of the draws, read once by a caller that
    asks for many: read for each, the 320 draws took minutes of the add of the queue."""
    stem = question.replace(": ", "_").replace(" ", "_").replace(",", "")
    twins, drawn = made or (digest(root), draw_digest(root))
    return Path(root) / "runs" / "scoreboard" / twins[:12] / f"draws_{drawn[:12]}" / f"{stem}_{weight}_draw{draw}.json"


def drawn(problem, draw: int, noise: tuple | None = None):
    """The records of a question with noise drawn on them at their own standard deviations, from the seed of the draw: the
    records of one experiment, which hold the planted shape, where the twin draws its noise on the fitted model. `noise`,
    from `NOISE_OF`, draws it otherwise than the twin does: ("correlated", rho), each sample correlated at rho with its
    neighbour along the abscissa, of the same deviation, or ("grows", g), the deviation 1 + g times its own at the top of
    the line. The records keep their stated deviations, which are what the fit and its twin assume."""
    from twinfit import Record
    rng = np.random.default_rng([2026, 10, 19, int(draw)])
    out = []
    for r in problem.records:
        y = np.asarray(r.y, dtype=float)
        sd = np.asarray(r.sd, dtype=float)
        z = rng.standard_normal(np.asarray(r.x).size)
        if noise and noise[0] == "correlated":
            order, rho = np.argsort(np.asarray(r.x, dtype=float), kind="stable"), float(noise[1])
            ar = z[order]
            for k in range(1, ar.size):
                ar[k] = rho * ar[k - 1] + np.sqrt(1.0 - rho ** 2) * ar[k]
            z[order] = ar
        elif noise and noise[0] == "grows":
            sd = sd * (1.0 + float(noise[1]) * (y - y.min()) / max(float(y.max() - y.min()), 1e-300))
        out.append(Record(r.x, y + sd * z, r.sd, r.design, r.name))
    return problem.with_records(out)


def agreed(first: dict, beyond: list, expected: dict) -> dict:
    """Whether a report puts first the answers written before, as `report_of` reads them: a report that puts nothing
    beyond the twin is right on no question but where nothing is missing, and a covariate named where the shape changes
    along none is a false answer. Written beside `report_of` and not called by it, since the twins hang on its text."""
    if not missing(expected):
        return {"silent": not beyond}
    out = {name: bool(beyond) and first.get(name) == answer for name, answer in expected.items() if answer is not None}
    if expected["changes along"] is None:
        out["unmoved"] = first.get("changes along") is None
    return out


def draw_report(question: str, weight: str, draw: int, const: dict, config: dict, root: Path = ROOT) -> dict:
    """The report of the package on one draw of a question at a weight, kept: the drawn records fitted from the fit of
    the noiseless ones, and read against the twin of that fit. A draw's own twin would centre its sets on its own fit,
    and cost ten parts a draw. A fit that does not end is kept as one, and counts as no answer."""
    from twinfit import fit, where
    from toy import analysis
    null = [{tuple(k.split("|", 1)): v for k, v in s.items()} for k in range(1, PARTS + 1) for s in json.loads(path_of(root, question, weight, k).read_text()) if s is not None]
    problem, result = fitted(question, weight, const, config)
    records = drawn(problem, draw, NOISE_OF.get(question.partition(": ")[2]))
    got = fit(records, result.shared, result.own, max_iter=80)
    out = {"ended": bool(got.converged)}
    if got.converged:
        report = where.where(records, got, null=null, centre=None if question.startswith("field") else analysis.CENTRE)
        first = report["first_of_each_question"]
        beyond = [[a["question"], a["answer"], power.said(a["excess"])] for a in report["answers"] if a["beyond_the_twin"]]
        out.update(first=first, beyond=beyond, agrees=agreed(first, beyond, next(q["expected"] for q in questions(config) if q["name"] == question)))
    path = draw_path(root, question, weight, draw)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out))
    return out


def counted_draws(draws: dict, config: dict) -> dict:
    """`counted` over the draws: for each weight and group, on how many draws of its questions the report puts first the
    answer written before. A draw whose fit did not end is asked and right on nothing."""
    out = {}
    for weight in WEIGHTS:
        for group in ("bench", "held out", "toy", "field", "waist"):
            cell = {q: [0, 0] for q in QUESTIONS + ("silent", "unmoved")}
            for q in [q for q in questions(config) if q["group"] == group and f"{q['name']}|{weight}" in draws]:
                asked = agreed({}, [], q["expected"])
                for d in draws[f"{q['name']}|{weight}"]:
                    for key in asked:
                        cell[key][0] += int(bool(d.get("agrees", {}).get(key)))
                        cell[key][1] += 1
            out.setdefault(weight, {})[group] = cell
    return out


def draws_of(question: dict) -> int:
    """The draws of noise on the records of a question: `DRAWS_IN` of its group, or `DRAWS`."""
    return DRAWS_IN.get(question["group"], DRAWS)


def kept_at(root: Path, twins: str, question: str, weight: str, part: int | None = None) -> Path:
    """`path_of` with the digest of the twins given, read once by a caller that asks for many paths. Read for each, the
    queue's calls of `units` and `progress` took up to 208 seconds each on 2026-10-19, after each wave, and fifteen of
    them held the machine at a load of 37. `path_of` itself is a part of what makes a twin, so it is not changed."""
    stem = question.replace(": ", "_").replace(" ", "_").replace(",", "")
    return Path(root) / "runs" / "scoreboard" / twins[:12] / (f"{stem}_{weight}_" + ("report.json" if part is None else f"twin_part{part}.json"))


def units(root: Path = ROOT) -> list:
    """What the queue asks: for each question at each weight the parts of its twin not yet made, then its report."""
    _, config = model.load(Path(root) / "toy", root)
    out, made = [], (digest(root), draw_digest(root))
    for q in questions(config):
        for weight in WEIGHTS:
            parts = [k for k in range(1, PARTS + 1) if not kept_at(root, made[0], q["name"], weight, k).exists()]
            kind = f"{q['group']} {weight}"
            out.append({"condition": f"twin|{weight}|{q['name']}", "units": parts, "after": [], "kind": f"twin {kind}"})
            out.append({"condition": f"report|{weight}|{q['name']}", "units": [] if kept_at(root, made[0], q["name"], weight).exists() else [1], "after": [f"twin|{weight}|{q['name']}"] if parts else [], "kind": f"report {kind}"})
            out.append({"condition": f"draw|{weight}|{q['name']}", "units": [k for k in range(1, draws_of(q) + 1) if not draw_path(root, q["name"], weight, k, made).exists()],
                        "after": [f"twin|{weight}|{q['name']}"] if parts else [], "kind": f"draw {kind}"})
    return out


def slice_of(condition: str, first: int, last: int, limit: float = float("inf"), root: Path = ROOT) -> dict:
    """The units of a condition from `first` to `last` not yet made, each kept as it ends, until `limit` seconds."""
    import resource
    const, config = model.load(Path(root) / "toy", root)
    kind, weight, question = condition.split("|", 2)
    began, ids, seconds = time.time(), [], []
    for k in range(first, last + 1):
        t = time.time()
        if kind == "twin" and not path_of(root, question, weight, k).exists():
            twin_part(question, weight, k, const, config, root)
        elif kind == "report" and not path_of(root, question, weight).exists():
            report_of(question, weight, const, config, root)
        elif kind == "draw" and not draw_path(root, question, weight, k).exists():
            draw_report(question, weight, k, const, config, root)
        else:
            continue
        ids.append(k)
        seconds.append(time.time() - t)
        if time.time() - began > limit:
            break
    return {"done": ids, "seconds": seconds, "rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


def counted(reports: dict, config: dict) -> dict:
    """For each weight and each group, on how many questions the report puts first the answer written before, place,
    kind and covariate apart, on how many where nothing is missing it says nothing, and on how many where the planted
    shape changes along nothing it names no covariate."""
    out = {}
    for weight in WEIGHTS:
        for group in ("bench", "held out", "toy", "field", "waist"):
            mine = [q for q in questions(config) if q["group"] == group and f"{q['name']}|{weight}" in reports]
            cell = {q: [0, 0] for q in QUESTIONS + ("silent", "unmoved")}
            for q in mine:
                for key, right in reports[f"{q['name']}|{weight}"]["agrees"].items():
                    cell.setdefault(key, [0, 0])[0] += int(right)          # a report kept before the blind report holds the question along
                    cell[key][1] += 1
            out.setdefault(weight, {})[group] = cell
    return out


def seen_by(reports: dict) -> dict:
    """For each weight, the family of statistics that the report names first as the one that sees most, and on how many
    questions: which part of the machinery reads the missing term."""
    out = {}
    for key, r in reports.items():
        family = (r.get("first") or {}).get("seen by")
        weight = key.rsplit("|", 1)[1]
        out.setdefault(weight, {})
        out[weight][family or "none"] = out[weight].get(family or "none", 0) + 1
    return {w: dict(sorted(c.items(), key=lambda kv: -kv[1])) for w, c in out.items()}


def progress(root: Path = ROOT) -> str:
    _, config = model.load(Path(root) / "toy", root)
    twins = digest(root)
    reports = {f"{q['name']}|{w}": json.loads(kept_at(root, twins, q["name"], w).read_text()) for q in questions(config) for w in WEIGHTS if kept_at(root, twins, q["name"], w).exists()}
    c = counted(reports, config)
    return f"{len(reports)} reports; " + "; ".join(f"{SAID[w]}: " + ", ".join(f"{q} {sum(c[w][g][q][0] for g in c[w])} of {sum(c[w][g][q][1] for g in c[w])}" for q in QUESTIONS + ("silent",)) for w in WEIGHTS)


def table(root: Path = ROOT) -> Path:
    _, config = model.load(Path(root) / "toy", root)
    twins = digest(root)
    reports = {f"{q['name']}|{w}": json.loads(kept_at(root, twins, q["name"], w).read_text()) for q in questions(config) for w in WEIGHTS if kept_at(root, twins, q["name"], w).exists()}
    if len(reports) != len(questions(config)) * len(WEIGHTS):
        raise SystemExit(f"REFUSED: {len(reports)} reports of {len(questions(config)) * len(WEIGHTS)} are kept: the queue has not ended")
    made = (digest(root), draw_digest(root))
    of = {q["name"]: draws_of(q) for q in questions(config)}
    kept = {key: [draw_path(root, key.rsplit("|", 1)[0], key.rsplit("|", 1)[1], k, made) for k in range(1, of[key.rsplit("|", 1)[0]] + 1)] for key in reports}
    held, asked = sum(p.exists() for paths in kept.values() for p in paths), sum(len(paths) for paths in kept.values())
    if held and held != asked:
        raise SystemExit(f"REFUSED: {held} draws of {asked} are kept: the queue has not ended")
    from twinfit.version import __version__
    board = {"digest": digest(root), "package": __version__, "weights": WEIGHTS, "seen": SEEN, "questions": questions(config), "reports": reports, "counted": counted(reports, config), "seen_by": seen_by(reports)}
    if held:          # each report with its draws, and the count over them (row M259)
        for key, paths in kept.items():
            reports[key]["draws"] = [json.loads(p.read_text()) for p in paths]
        board.update(draw_digest=draw_digest(root), counted_draws=counted_draws({key: r["draws"] for key, r in reports.items()}, config))
    out = Path(root) / "toy" / "scoreboard.json"
    out.write_text(json.dumps(board, indent=1) + "\n")
    return out


def check(root: Path = ROOT) -> list:
    path = Path(root) / "toy" / "scoreboard.json"
    if not path.exists():
        return ["toy/scoreboard.json is not written: `python3 tools/autoqueue.py add toy/scoreboard.py`, then `python toy/scoreboard.py table`"]
    held = json.loads(path.read_text())
    if held.get("digest") != digest(root):
        return ["toy/scoreboard.json was made from other files than those on disk now: count it again"]
    _, config = model.load(Path(root) / "toy", root)
    if held["counted"] != counted(held["reports"], config):
        return ["toy/scoreboard.json does not hold what its own reports add up to"]
    if "counted_draws" in held and held.get("draw_digest") != draw_digest(root):
        return ["the draws of toy/scoreboard.json were made by another reading of a draw than the one on disk now: count them again"]
    if "counted_draws" in held and held["counted_draws"] != counted_draws({key: r.get("draws", []) for key, r in held["reports"].items()}, config):
        return ["toy/scoreboard.json does not hold what its own draws add up to"]
    return []


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["units"]:
        print(json.dumps(units()))
    elif args[:1] == ["slice"] and len(args) >= 4:
        if not os.environ.get("TOY_QUEUE"):
            sys.exit("REFUSED: a slice is a wave of the queue: `python3 tools/autoqueue.py add toy/scoreboard.py` (his message 147)")
        print(json.dumps(slice_of(args[1], int(args[2]), int(args[3]), float(args[5]) if args[4:5] == ["--for"] else float("inf"))))
    elif args == ["progress"]:
        print(progress())
    elif args == ["digest"]:          # the digest of the twins, which the gate and the builder compare with the table's: a draw of another reading is kept under another folder, which `units` lists as missing, so the queue needs no second digest (row M259)
        print(digest())
    elif args == ["table"]:
        print(f"written: {table().relative_to(ROOT)}")
    elif args == ["--check"]:
        wrong = check()
        print("\n".join(wrong) or "the scoreboard: made from the files as they are, and its reports give the same counts")
        sys.exit(1 if wrong else 0)
    else:
        print(__doc__)
        sys.exit(1)
