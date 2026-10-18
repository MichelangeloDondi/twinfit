"""The six proposals that end an iteration, the owner's message 142: three ways to improve the model and three ways to
improve the design of the experiment, each measured on the twin by the bias and the scatter of each known parameter.

    out = proposals(problem, result, report, levers, sets=50)

`report` is the report of where to search (where.py) on the same fit, and `levers` the levers of the description
(description.py), each a dict with its name, its kind and its words.

For the model. The classes and the covariates that the report puts beyond the twin, then the others, each by its
excess, classes and covariates together: for each, the candidate of the dictionary that the residuals take most
of, made into a term of one free size with the shape that the residuals give it (`term_of`). The term is fitted, and the
twin is drawn from the model with the term at its fitted size: on each set the working model and the model with the term
are fitted. Each proposal gives, for each known parameter, its bias and its scatter on those sets without the term and
with it. The term is a shape and not a cause: which term of physical origin has that shape is the researcher's choice.

For the experiment. The levers are those of the description, and beside them the package adds
those that every design has, `DEFAULTS`, a quieter setup, more samples, more repeats and a wider span, of the kinds that
the description does not hold, and a better measured value for each value that it holds. Each lever is one of the three
ways of the owner's message 142, a better setup, a better characterization of the setup or a better data collection,
and the output says which. Each enters the twin in one of six ways:

    noise      the noise of every sample times `factor`: a better stabilisation, a better instrument, a better gain
    samples    `factor` times as many samples over the same span of each record: a slower scan, a faster recorder
    span       the span of the abscissa of each record times `factor`, about its middle, with the same step
    repeats    each record taken `factor` times, each with its own noise: more traces at the same conditions
    characterization  the uncertainty of the value `of` that the description holds, times `factor`, under one
    covariate  records added at the new `values` of one covariate, each with the own parameters of the record nearest
               in that covariate, among those at the conditions `at` when the lever names them: a new temperature, a
               new power, a sweep of a quantity that the records hold at one value

The twin of a lever is drawn under its design from the working model at its fit, with the term of the first model
proposal at its fitted size, and fitted by the working model, as the twin of the records is. On a new record the term
keeps the shape that it has on the nearest record, an assumption that the output names. The levers are ranked by the
error of the known parameters under each, the root of the bias squared, the scatter squared and the bar of theory of
the values that the description holds squared, over that under the design of the records, and the three of smallest
ratio are proposed. Where the working model shares no parameter, as a first model may, a lever is measured by the
size of the term that the report proposes first, fitted free on each set: how well the design measures what the model
lacks. A lever of kind characterization leaves the records as they are and narrows the bar of one held
value, so a better characterization of the setup is set against a better instrument on one scale. Where a covariate is held at one value, a lever
along it is the only way to read a behaviour along it.

The package reads no physics in a lever: its words are the researcher's, and so is the choice.
"""
from __future__ import annotations

import math
import re

import numpy as np

from .best import theory_moves
from .candidates import ALONG, ALONG_SQUARE, Candidate, dictionary
from .fit import fit
from .linear import Linear
from .mechanisms import Mechanism
from .problem import Problem, Record
from .propose import _one
from .refine import develop_one, on_the_twin
from .where import CLASSES

KINDS = {"noise": "the noise of every sample is multiplied by the factor",
         "samples": "each record has the factor times as many samples over the same span",
         "span": "the abscissa of each record spans the factor times its span about its middle, with the same step",
         "repeats": "each record is taken the factor times, each with its own noise",
         "covariate": "records are added at the new values of one covariate, each with the own parameters of the record nearest in that covariate",
         "characterization": "the uncertainty of a value held by the description is multiplied by the factor, and the records stay as they are",
         "setup": "a change of an item of the setup: where the model holds a value of that name, its uncertainty is multiplied by the factor, a half where none is given, and otherwise the change is named and not measured"}
SHAPES = dict(CLASSES)
HIS_KINDS = {"noise": "a better setup", "characterization": "a better characterization of the setup", "samples": "a better data collection",
             "span": "a better data collection", "repeats": "a better data collection", "covariate": "a better data collection", "setup": "a better setup"}          # the three ways of his message 142
DEFAULTS = ({"name": "a quieter setup", "kind": "noise", "factor": 0.5, "words": "the noise of every sample halved: a better stabilisation, a better instrument or a better gain, or, for counts, four times the time of counting"},
            {"name": "more samples", "kind": "samples", "factor": 2.0, "words": "twice the samples over the same span of each record: a slower scan or a faster recorder, and not, for counts, twice the channels of a histogram in the same time, which halves the counts of each"},
            {"name": "more repeats", "kind": "repeats", "factor": 2, "words": "each record taken twice, each with its own noise"},
            {"name": "a wider span", "kind": "span", "factor": 1.5, "words": "each record spans one and a half times its span about its middle, with the same step"})
NOISE_OF_THE_SAMPLES = ("glitches", "a noise that grows with the signal", "a noise correlated from sample to sample")          # the kinds that the noise of each sample holds, and no term
TERM = "size of the term"          # what a lever is measured by where the working model shares no parameter


def by_gain(model: list) -> list:
    """The proposals for the model by the gain in chi square that each brings, those beyond the twin first, then those
    whose term of the package does not end at its bound when it is fitted, and those with no gain last: in the order of
    the excess the first was an impostor that the loop kept (the cold tests of row M169, round 2), and a first proposal
    whose term ended at its bound, the residuals asking its shape with a sign that it does not take, was refused by the
    loop (round 3). A term whose search does not end at the distance or the stretch that the dictionary guessed keeps its
    place: the class is right and the number is rough, and a scan finds it."""
    def gain(r):          # the gain of the term of the package once fitted where it stands, else that of the shape that the residuals give: a periodic term of 78599 came first beside a second component of 2683953 (round 10 of the cold tests)
        pk = r.get("the_term_of_the_package") or {}
        return pk["gain_in_chi_square"] if pk.get("stands") and pk.get("gain_in_chi_square") is not None else r.get("gain_in_chi_square", -math.inf)
    kinds = [r for r in model if r["beyond_the_twin"] and "excess_in_the_report" in r]
    lead = max(kinds, key=lambda r: r.get("excess_in_the_report", -math.inf)) if kinds else None          # the proposal of the answer that leads the report comes first, of what or of along: a tail was put before the second component that the row of what named (round 11), and the width along the temperature of the worked spectra leads its report
    return sorted(model, key=lambda r: (not r["beyond_the_twin"], bool((r.get("the_term_of_the_package") or {}).get("at_its_bound")), r is not lead, -gain(r)))


def package_fit(problem, result, line: str, max_iter: int = 80) -> dict | None:
    """The term that the line of a proposal writes, fitted beside the working model as the loop fits a term file: its
    size and its bar, or why it does not stand. The proposal itself is the shape that the residuals give a candidate, of
    either sign, and a term of the package may bound its size: a tail on the lower side, proposed with a gain in chi
    square of 9230, ended at its bound of zero in the loop, which refused it, and the report had not said so (round 3 of
    the cold tests). None where the line names a parameter that the model holds already."""
    from . import mechanisms
    space = {"mechanisms": mechanisms}
    try:
        exec(line, space)          # a line that this file wrote, of a term of mechanisms.py
        term = space["term"]
        if term.parameter in problem.shared:
            return None
        got = develop_one(problem, result, term, max_iter)
    except (ValueError, KeyError, TypeError) as err:
        return {"stands": False, "why_not": str(err)}
    if got["ended"] and not got["at_its_bound"]:
        return {"stands": True, "parameter": got["parameter"]["name"], "size": got["parameter"]["value"], "bar": got["parameter"]["bar"],
                "gain_in_chi_square": got["criteria"]["gain_in_chi_square"]}
    if got["ended"]:
        return {"stands": False, "at_its_bound": True, "why_not": "its fit ends at the bound of its size, where the term is absent or as large as it may be: the residuals ask for this shape "
                                                                  "with a sign that the term does not take, and the loop refuses it"}
    chosen = re.match(r"term = mechanisms\.(second|pedestal|tail|periodic)\(", line)          # a term with a number that the dictionary guessed
    return {"stands": False, "at_its_bound": False, "why_not": got.get("why_not", "its fit did not end")
            + (f". Its {'distance' if chosen[1] == 'second' else 'period' if chosen[1] == 'periodic' else 'stretch'} is the rough guess of the dictionary: "
               f"`python -m twinfit scan DIR \"{chosen[1]}(x{', ' + repr(line.split(',', 1)[1].strip(' )').strip(chr(34) + chr(39))) if chosen[1] == 'tail' else ''})\" ...` tries others" if chosen else "")}


def as_a_term_file(name: str, unit: float, along: str | None = None) -> tuple:
    """The line of a term file that writes a proposal for the model with a term of the package, and None with the reason
    where no term of the package holds it. The dictionary counts a distance in widths, the full width of the fitted
    curves at half their height, `unit`, and a term file takes it in the unit of the abscissa: "0.35 widths" was read as
    0.35 standard deviations, and the scan took seconds (the cold tests of row M169). A background of its own size in
    each record has no term of one parameter, and where the report names a setting `along` which the misfit changes,
    the two term files that hold a background growing along it are said."""
    import re
    number = r"([-+0-9.e]+)"
    if m := re.fullmatch(f"a second component {number} widths (above|below) the centre", name):
        return f"term = mechanisms.second({float(m[1]) * unit * (1.0 if m[2] == 'above' else -1.0):.3g})", None
    if m := re.fullmatch(f"a second component {number} widths away on both sides", name):
        return None, f"two components, one each side: two term files, second({float(m[1]) * unit:.3g}) and second({-float(m[1]) * unit:.3g})"
    if m := re.fullmatch(f"a pedestal {number} times wider than the curve", name):
        return f"term = mechanisms.pedestal({float(m[1]):g})", None
    if m := re.fullmatch(f"a tail on the (upper|lower) side, {number} times wider", name):
        return f'term = mechanisms.tail({float(m[2]):g}, "{m[1]}")', None
    if m := re.fullmatch(f"a periodic term of period {number}", name):
        return f'term = mechanisms.periodic({float(m[1]):.3g}, "cosine")', None
    if name in ("a slope of the level", "a curvature of the level", "a cubic term of the level"):
        return f"term = mechanisms.level({1 if 'slope' in name else 2 if 'curvature' in name else 3})", None
    if name.endswith(", its own in each record"):          # the term of one size for all is offered too: a curvature of one size in every record was read as one of its own in each, and no line was given (round 6 of the cold tests)
        power = 1 if "slope" in name else 2 if "curvature" in name else 3
        grows = (f". Where it grows along {along}, two term files hold it: term = mechanisms.level({power}), then term = mechanisms.along(\"level power {power}\", \"{along}\")"
                 if along else "")
        return f"term = mechanisms.level({power})", ("the term of this line holds it at one size in every record, where the residuals give each record its own. Where each record carries its own size, no term of one parameter holds it: "
                                                     "a parameter of the model, own to each record, or a Mechanism of your own" + grows)
    for joint, square in ((ALONG_SQUARE, True), (ALONG, False)):
        if joint in name:
            parameter, key = name.split(joint, 1)
            return f'term = mechanisms.along("{parameter}", "{key}"' + (", power=2)" if square else ")"), None
    return None, "no term of the package holds it: a Mechanism of your own, as section 11 of TUTORIAL.md shows"


def term_of(problem, cand, amplitude) -> Mechanism:
    """A candidate of the dictionary as a term of the model: the shape that the residuals give it, the sum of its columns
    at their fitted amplitudes, with one free size that is one at that shape. The shape is known on the records that it
    was read from, by the design of each: between their samples it is interpolated, and outside their span it keeps its
    end values. On a record of another design it has no value, and is refused."""
    shapes = shapes_of(problem, cand, amplitude)

    def wrap(model, larger, n):
        def new(x, shared, own, design):
            return np.asarray(model(x, shared[:n], own, design), dtype=float) + float(shared[n]) * _shape(shapes, x, design)
        return new
    return Mechanism(cand.name, cand.family, f"size of {cand.name}", 1.0, wrap)


def shapes_of(problem, cand, amplitude) -> list:
    """The shape of a candidate on each record, as (abscissa, design, shape): the sum of its columns at `amplitude`."""
    out = []
    for k, rec in enumerate(problem.records):
        parts = [a * np.asarray(col[k], dtype=float) for col, a in zip(cand.columns, amplitude) if col[k] is not None]
        out.append((np.asarray(rec.x, dtype=float), rec.design, sum(parts) if parts else np.zeros(rec.x.size)))
    return out


def _shape(shapes, x, design):
    for xs, d, s in shapes:
        if d is design:
            return np.interp(np.asarray(x, dtype=float), xs, s)
    for xs, d, s in shapes:
        if d == design:
            return np.interp(np.asarray(x, dtype=float), xs, s)
    raise ValueError("the shape of a term is known on the records that it was read from: a record of another design has none")


def period_of(problem, result, periods: int = 600) -> dict | None:
    """The period of a ripple in what the fit leaves, read from a periodogram of the whitened residuals summed over the
    records, on periods from three steps of the abscissa to half its span, and in which frame the ripple stands: the
    coherence of the records at that period, the share of the summed power that stays when the records are added with
    their phases, about the abscissa itself and about the centre of each curve. A ripple of the instrument, fixed to the
    abscissa, agrees from curve to curve at equal abscissa, which the term of the package, placed about each centre, does
    not hold. The proposal gave the width of the curves, 2.82, as the period of a ripple of 0.35 (round 14 of the cold
    tests). None where the records are too short to read one."""
    series = []
    for k, rec in enumerate(problem.records):
        x, r = np.asarray(rec.x, dtype=float), np.asarray(problem.whitened(k, result.shared, result.own[k]), dtype=float)
        centre = problem.centre_of(k, result.shared, result.own[k]) if getattr(problem, "centre", None) is not None else None
        if x.size >= 8:
            series.append((x, r - r.mean(), centre))
    if not series:
        return None
    step = float(np.median([np.median(np.abs(np.diff(np.sort(x)))) for x, _, _ in series]))
    span = float(max(np.ptp(x) for x, _, _ in series))
    if not (step > 0 and span > 6.0 * step):
        return None
    grid = np.exp(np.linspace(np.log(3.0 * step), np.log(0.5 * span), periods))
    power, phasors = np.zeros(grid.size), []
    for x, r, _ in series:
        z = np.exp(-2j * np.pi * np.outer(1.0 / grid, x)) @ r
        power += np.abs(z) ** 2
        phasors.append(z)
    best = int(np.argmax(power))
    period, at = float(grid[best]), np.array([z[best] for z in phasors])
    out = {"period": period, "coherence_about_the_abscissa": float(abs(at.sum()) ** 2 / (len(at) * float(np.sum(np.abs(at) ** 2))))}
    if all(c is not None for _, _, c in series):
        turned = at * np.exp(2j * np.pi * np.array([c for _, _, c in series]) / period)
        out["coherence_about_the_centres"] = float(abs(turned.sum()) ** 2 / (len(turned) * float(np.sum(np.abs(turned) ** 2))))
    return out


def held_near(parameter: str | None, shared, unit: float) -> bool:
    """Whether the model holds a term of the kind of `parameter` at nearly its number: the same name but for its last
    number, within a tenth of it, or for a distance within a quarter of the width of the curves. A second component at
    -203 was proposed beside one kept at -204.6 (round 13 of the cold tests). Two powers of a level are two terms."""
    import re
    form = r"(.*?)\s*(-?\d+(?:\.\d+)?)"
    m = re.fullmatch(form, parameter or "")
    if not m:
        return False
    for name in shared:
        n = re.fullmatch(form, str(name))
        if n and n.group(1).strip() == m.group(1).strip():
            a, b = float(m.group(2)), float(n.group(2))
            if 0.0 < abs(a - b) and (abs(a - b) <= 0.1 * max(abs(a), abs(b)) or (m.group(1).strip().endswith(" at") and abs(a - b) <= 0.25 * unit)):          # an exact copy is left out by its name, before
                return True
    return False


def line_parameter(line: str | None) -> str | None:
    """The parameter that the line of a term file adds, or None where the line names no term of the package."""
    if not line:
        return None
    from . import mechanisms
    space = {"mechanisms": mechanisms}
    try:
        exec(line, space)          # a line that this file wrote, of a term of mechanisms.py
        return space["term"].parameter
    except (ValueError, KeyError, TypeError, NameError, SyntaxError):
        return None


def refused_before(line: str | None, refused) -> str | None:
    """The file and the reason of the last refusal of the term that a line writes, by the parameter that it adds, or
    None: level(1) was proposed again after its refusal, with no word (round 8 of the cold tests)."""
    added = line_parameter(line)
    before = [r for r in refused if added is not None and r.get("parameter") == added]
    return f"{before[-1]['file']}: {before[-1]['why']}" if before else None


def follows(problem, result, amplitude, bars) -> dict | None:
    """Whether the size of a shape that each record carries at its own size, one a record, is a share of one own
    parameter of the model, the same share in every record, within the bars of the sizes: the share with its bar and
    the chi square of that law on its degrees of freedom, for the own parameter that holds it best, where its chi square
    is within three of its standard deviations of its degrees of freedom. None elsewhere, and None where one value in every record holds the sizes as well, since such sizes follow any own value that does not change. A step under a peak, of a size
    in proportion to its height, was proposed as a level along the angle that the height followed, and the loop kept it
    and moved the width 30 bars (round 8 of the cold tests): the sizes read against each own value say what they
    follow, which the settings say only through it."""
    a, s = np.asarray(amplitude, dtype=float), np.asarray(bars, dtype=float)
    n = a.size
    if n < 3 or a.size != len(problem.records) or not np.all(s > 0):
        return None
    w = 1.0 / s ** 2
    if float(np.sum(w * (a - np.sum(w * a) / np.sum(w)) ** 2)) <= (n - 1) + 3.0 * math.sqrt(2.0 * (n - 1)):
        return None          # sizes of one value in every record follow any own value that does not change: a share of a baseline of 0.1 in every curve read as a lead (round 12 of the cold tests)
    best = None
    for j, name in enumerate(problem.own):
        v = np.array([float(result.own[k][j]) for k in range(n)])
        if np.any(v == 0.0) or np.ptp(v) <= 1e-9 * np.max(np.abs(v)) or (v.min() < 0.0 < v.max()):
            continue          # an own value that does not change from record to record tells nothing apart from a constant, and one that changes sign is a place and not a size: the worked NMR read a slope as a share of its places, the trace of a misfit at each place (round 9 of the cold tests)
        w = 1.0 / s ** 2
        share = float(np.sum(w * a * v) / np.sum(w * v * v))
        chi = float(np.sum(w * (a - share * v) ** 2))
        if chi <= (n - 1) + 3.0 * math.sqrt(2.0 * (n - 1)) and (best is None or chi < best["chi_square"]):
            best = {"own": name, "index": j, "share": share, "bar": float(1.0 / math.sqrt(np.sum(w * v * v))), "chi_square": chi, "degrees_of_freedom": n - 1}
    return best


def model_terms(problem, result, report, most: int = 3, exclude=()) -> list:
    """The candidates of the model proposals, in the order of the report, at most `most`: each with the class that it
    stands for and why the report names it. A candidate named as a term that the model holds, in `exclude`, is left out:
    the report proposed again a tail that the loop had kept (the cold tests of row M169, round 2)."""
    lin = Linear(problem, result)
    best = {}
    unit = float(np.median([lin.width_of(k) for k in range(len(problem.records))]))
    along = (report.get("first_of_each_question") or {}).get("along")
    for cand in dictionary(lin, covariates=[a["answer"] for a in report["answers"] if a["question"] == "along"]):
        named = line_parameter(as_a_term_file(cand.name, unit, along)[0])
        if cand.name in exclude or named in problem.shared or held_near(named, problem.shared, unit):          # a term that the model holds, by the parameter of its line: kept terms are named by their term and not by the candidate, and level(2) was proposed beside a kept level(2) (round 8)
            continue
        entry = _one(lin, cand)
        if entry["status"] == "ranked" and (cand.family not in best or entry["score"] > best[cand.family][1]["score"]):
            best[cand.family] = (cand, entry)
        if entry["status"] == "ranked" and cand.family == "grows with a covariate":
            key = cand.name.split(ALONG_SQUARE)[-1] if ALONG_SQUARE in cand.name else cand.name.split(ALONG)[-1]
            if ("along", key) not in best or entry["score"] > best[("along", key)][1]["score"]:
                best[("along", key)] = (cand, entry)
    family_of = {said: kind for kind, said in CLASSES}
    order = sorted([a for a in report["answers"] if (a["question"] == "what" and a["answer"] in family_of) or a["question"] == "along"],
                   key=lambda a: (not a["beyond_the_twin"], -a["excess"]))
    out = []
    for a in order:
        key = family_of[a["answer"]] if a["question"] == "what" else ("along", a["answer"])
        if key in best and all(o["candidate"] is not best[key][0] for o in out):
            cand, entry = best[key]
            out.append({"candidate": cand, "entry": entry, "stands_for": a["answer"] if a["question"] == "what" else f"a change along {a['answer']}",
                        "beyond_the_twin": a["beyond_the_twin"], "excess": a["excess"]})
        if len(out) == most:
            break
    return out


def _designed(problem, result, lever) -> list:
    """The records of the design of a lever, as (abscissa, standard deviation, design, own parameters, the record it is
    made from)."""
    kind = lever["kind"]
    out = []
    for k, rec in enumerate(problem.records):
        x, sd, own = np.asarray(rec.x, dtype=float), np.asarray(rec.sd, dtype=float), np.asarray(result.own[k], dtype=float)
        if kind == "noise":
            out.append((x, sd * float(lever["factor"]), rec.design, own, k))
        elif kind == "samples":
            many = max(int(round(x.size * float(lever["factor"]))), 2)
            new = np.linspace(x.min(), x.max(), many)
            out.append((new, np.interp(new, x, sd), rec.design, own, k))
        elif kind == "span":
            step, middle = (x.max() - x.min()) / max(x.size - 1, 1), 0.5 * (x.max() + x.min())
            half = 0.5 * (x.max() - x.min()) * float(lever["factor"])
            new = np.arange(middle - half, middle + half + 0.5 * step, step)
            order = np.argsort(x)
            out.append((new, np.interp(new, x[order], sd[order]), rec.design, own, k))          # the noise of the nearest sample, the edges carried past the ends: the median gave the peak of a curve of counts the noise of its tail (round 14 of the cold tests)
        elif kind == "repeats":
            out += [(x, sd, rec.design, own, k)] * int(lever["factor"])
        elif kind in ("covariate", "characterization"):
            out.append((x, sd, rec.design, own, k))
        else:
            raise ValueError(f"a lever of kind {kind!r}: the kinds are {', '.join(KINDS)}")
    if kind == "covariate":
        key = lever["covariate"]
        among = [k for k, rec in enumerate(problem.records) if all(rec.design.get(c) == v for c, v in (lever.get("at") or {}).items())]
        if not among:
            raise ValueError(f"the lever {lever.get('name', '')!r} starts from the records at {lever['at']}: none is at those conditions")
        for value in lever["values"]:
            k = among[int(np.argmin([abs(float(problem.records[k].design[key]) - float(value)) for k in among]))]
            rec = problem.records[k]
            out.append((np.asarray(rec.x, dtype=float), np.asarray(rec.sd, dtype=float), dict(rec.design, **{key: value}), np.asarray(result.own[k], dtype=float), k))
    return out


def _held_term(problem, designed, shapes, size) -> object:
    """The working model with the term of fixed shape and free size, on the records of a design: each record takes the
    shape of the record that it was made from."""
    like = {id(design): k for _, _, design, _, k in designed}
    n = len(problem.shared)

    def model(x, shared, own, design):
        return np.asarray(problem.model(x, shared[:n], own, design), dtype=float) + float(shared[n]) * _shape(shapes, x, problem.records[like[id(design)]].design)
    lo, hi = problem.bounds_shared
    centre = None if problem.centre is None else (lambda shared, own, design, c=problem.centre: c(shared[:n], own, design))
    return lambda records: Problem(records, model, problem.shared + ("size of the term",), problem.own, list(problem.start_shared) + [size],
                                   [own for _, _, _, own, _ in designed], centre, f"{problem.name} with the term", (list(lo) + [-np.inf], list(hi) + [np.inf]))


def _on_design(problem, result, designed, term, sets: int, seed: int, max_iter: int) -> dict | None:
    """The bias and the scatter of each known parameter on `sets` twins of one design, drawn from the working model at its
    fit with `term`, (its shapes on the records, its size), and fitted from the fitted values by the working model, and
    with a term by the working model with that term, of free size."""
    larger = None if term is None else _held_term(problem, designed, *term)
    est, held = [], []
    owns = [own for _, _, _, own, _ in designed]
    for s in range(sets):
        rng = np.random.default_rng([seed, s])
        recs = []
        for x, sd, design, own, k in designed:
            mean = np.asarray(problem.model(x, result.shared, own, design), dtype=float)
            if term is not None:
                shapes, size = term
                mean = mean + size * _shape(shapes, x, problem.records[k].design)
            recs.append(Record(x, mean + sd * rng.standard_normal(x.size), sd, design, ""))
        for into, make, start in ((est, problem.with_records, result.shared), (held, larger, None if term is None else np.append(result.shared, term[1]))):
            if make is None:
                continue
            try:
                got = fit(make(recs), start, owns, max_iter=max_iter)
            except (np.linalg.LinAlgError, FloatingPointError, ValueError):
                continue
            if got.converged:
                into.append(got.shared[:len(problem.shared) + (into is held)])          # with the term, its size too: the one parameter that a model sharing none can be measured by
    if len(est) < max(10, sets // 2) or (term is not None and len(held) < max(10, sets // 2)):
        return None
    out = {}
    for label, values, truth, names in (("with_the_working_model", est, list(result.shared), list(problem.shared)),
                                        ("with_the_term", held, list(result.shared) + ([] if term is None else [term[1]]), list(problem.shared) + [TERM])):
        if not values or not names:
            continue
        values = np.array(values)
        out[label] = {name: {"bias": float(values[:, j].mean() - truth[j]), "scatter": float(values[:, j].std(ddof=1)),
                             "error": float(np.hypot(values[:, j].mean() - truth[j], values[:, j].std(ddof=1)))} for j, name in enumerate(names)}
    return {"sets": len(est), **out}


def theory_of(moved: dict, j: int, narrowed: str | None = None, factor: float = 1.0) -> float:
    """The bar of theory of the shared parameter `j`: the bars that the held values put on it, added in squares, with
    the value `narrowed` at `factor` times its own."""
    return float(np.sqrt(sum((v[j] * (factor if name == narrowed else 1.0)) ** 2 for name, v in moved.items())))


def proposals(problem, result, report, levers=(), sets: int = 50, seed: int = 41, max_iter: int = 80, held: dict | None = None, exclude=(), refused=()) -> dict:
    """Three proposals for the model and three for the experiment, each measured on the twin. Fewer are given where the
    report or the description holds fewer, and the output says so. `held`: the values that the description holds,
    {name: (centre, bar, at)}. Their bar of theory is added in squares to the error of each known parameter under each
    lever, so that a quieter detector is not proposed where the theory sets the bar, and a lever of kind
    characterization, which narrows one of them, is measured by it. The bar of theory is that of the fit of the records,
    for every lever: a lever that adds records at other values of a held covariate would move it, and that is not
    measured here (review 22, minor 15). `refused`: the terms that the loop refused, each with its file, its parameter and
    why: a proposal of a term refused before says so, since level(1) was proposed again after its refusal with no word
    (round 8 of the cold tests)."""
    model = []
    lin = Linear(problem, result)
    unit = float(np.median([lin.width_of(k) for k in range(len(problem.records))]))          # the width in which the dictionary counts its distances
    for c in model_terms(problem, result, report, exclude=exclude):
        along = (report.get("first_of_each_question") or {}).get("along")
        line, no_line = as_a_term_file(c["candidate"].name, unit, along)
        ripple = period_of(problem, result) if line and line.startswith("term = mechanisms.periodic(") else None
        if ripple:          # the period of the residuals, where the dictionary counts it in widths of the curves: 2.82 for a ripple of 0.35 (round 14 of the cold tests)
            period = float(f"{ripple['period']:.4g}")
            centre = [lin.centre_of(k) for k in range(len(problem.records))]
            found = Candidate(f"a periodic term of period {period:g}", "periodic",
                              [[np.cos(2.0 * np.pi * (rec.x - centre[k]) / period) for k, rec in enumerate(problem.records)],
                               [np.sin(2.0 * np.pi * (rec.x - centre[k]) / period) for k, rec in enumerate(problem.records)]], "any")
            entry = _one(lin, found)
            if entry.get("amplitude") is not None:
                c = dict(c, candidate=found, entry=entry)
            line = f'term = mechanisms.periodic({period:g}, "cosine")'
        mechanism = term_of(problem, c["candidate"], c["entry"]["amplitude"])
        developed = develop_one(problem, result, mechanism, max_iter)
        row = {"term": c["candidate"].name, "stands_for": c["stands_for"], "beyond_the_twin": c["beyond_the_twin"], "excess_in_the_report": c["excess"],
               "as_a_term_file": line, "no_term_of_the_package": no_line}
        if ripple:
            row["ripple"] = ripple
        said = refused_before(line, refused)
        if said:
            row["refused_before"] = said
        if c["candidate"].name.endswith(", its own in each record"):
            got = follows(problem, result, c["entry"]["amplitude"], c["entry"]["amplitude_standard_deviation"])
            if got:
                row["follows"] = got
        if line:
            row["the_term_of_the_package"] = package_fit(problem, result, line, max_iter)
        if not developed["ended"]:
            row["why_not_measured"] = developed.get("why_not", "its fit did not end")
            model.append(row)
            continue
        twin = on_the_twin(problem, result, developed, sets=sets, seed=seed, max_iter=max_iter)
        row.update({"size": developed["parameter"]["value"], "bar_of_the_size": developed["parameter"]["bar"],
                    "gain_in_chi_square": developed["criteria"]["gain_in_chi_square"], "sets_of_the_twin": twin["read"]})
        if "what_the_whole_procedure_returns" in twin:
            back, shows = twin["what_the_whole_procedure_returns"], twin["what_the_working_model_shows"]
            row["known_parameters"] = {name: {"without_the_term": {"bias": shows[name]["mean_on_these_sets"] - back[name]["put_in"], "scatter": shows[name]["scatter"]},
                                              "with_the_term": {"bias": back[name]["mean"] - back[name]["put_in"], "scatter": back[name]["scatter"]},
                                              "the_records_in_scatters": shows[name]["the_records_in_scatters_from_that_mean"]}          # where the working model reads the records far from what it reads on this twin, the twin of the term does not make the records, and the cost is not the term's (round 20 of the tests of the tutorial)
                                       for name in problem.shared}
        row["_term"] = (shapes_of(problem, c["candidate"], c["entry"]["amplitude"]), developed["parameter"]["value"])
        model.append(row)
    model = by_gain(model)
    leading = sorted([a for a in report["answers"] if a["question"] == "what" and a["beyond_the_twin"]], key=lambda a: -a["excess"])
    noisy = bool(leading) and leading[0]["answer"] in NOISE_OF_THE_SAMPLES          # a second component measured in the twin of the levers, under a report that led with a correlated noise (round 15 of the cold tests)
    chosen = None if noisy else next((r for r in model if "_term" in r and r["beyond_the_twin"]), None)
    term = chosen["_term"] if chosen else None          # a proposal within the twin is not put into the twin of the levers: the report does not ask for it
    here = _on_design(problem, result, _designed(problem, result, {"kind": "noise", "factor": 1.0}), term, sets, seed, max_iter)
    moved = theory_moves(result, held, max_iter)
    levers = list(levers)
    kinds = {lever["kind"] for lever in levers}          # the levers of the design that every problem has, of each kind that the description does not name: three named levers left none of them, against the guide (round 16 of the cold tests)
    levers += [dict(d, from_the_package=True) for d in DEFAULTS if d["kind"] not in kinds]
    named = {lever.get("of") for lever in levers if lever["kind"] in ("characterization", "setup")}
    levers += [{"name": f"{name} measured better", "kind": "characterization", "of": name, "factor": 0.5, "words": "its uncertainty halved", "from_the_package": True} for name in moved if name not in named]
    names = list(problem.shared) or [TERM]          # where the model shares no parameter, a lever is measured by the size of the term that the report proposes first
    experiment = []
    for lever in levers:
        row = {"lever": lever["name"], "kind": lever["kind"], "his_kind": HIS_KINDS[lever["kind"]], "words": lever.get("words", ""), "how_it_enters_the_twin": KINDS[lever["kind"]],
               "from_the_package": bool(lever.get("from_the_package")), "measured_by": "the known parameters" if problem.shared else f"the {TERM} that the report proposes first, since the model shares no parameter"}
        row.update({k: lever[k] for k in ("factor", "covariate", "values", "of") if k in lever})          # which value a characterization narrows (review 22, minor 15)
        if lever["kind"] == "covariate":
            row["assumes"] = "the term of the first model proposal keeps, on a new record, the shape that it has on the record nearest in that covariate"
        narrowed = lever.get("of") if lever["kind"] in ("characterization", "setup") else None
        if narrowed is not None and narrowed not in moved:
            row["why_not_measured"] = (f"no value of the model reads {narrowed!r} of the setup: a term that reads it, with its value held in HELD of the model and in an ansatz, lets the twin measure this change"
                                       if lever["kind"] == "setup" else f"the model holds no value {narrowed!r} from the description: a better characterization of it narrows nothing here")          # a lever of the setup that no term reads, named and not measured (his message 184, row M216)
            experiment.append(row)
            continue
        got = here if narrowed is not None else _on_design(problem, result, _designed(problem, result, lever), term, sets, seed, max_iter)          # a characterization leaves the records as they are
        if got is None or here is None:
            row["why_not_measured"] = "fewer than half of the fits of its twin ended"
        else:
            factor = float(lever.get("factor", 0.5 if lever["kind"] == "setup" else 1.0)) if narrowed is not None else 1.0
            for label in ("with_the_term", "with_the_working_model"):
                if label in got and all(n in got[label] for n in names):
                    bar = lambda j, *a: theory_of(moved, j, *a) if j < len(problem.shared) else 0.0          # noqa: E731
                    row[label] = {"known_parameters": {name: dict(got[label][name], under_the_records=here[label][name], bar_of_theory=bar(j, narrowed, factor),
                                                                  bar_of_theory_under_the_records=bar(j)) for j, name in enumerate(names)}}
                    totals = [(math.hypot(got[label][n]["error"], bar(j, narrowed, factor)), math.hypot(here[label][n]["error"], bar(j))) for j, n in enumerate(names)]
                    ratios = [a / b for a, b in totals if b > 0]
                    row[label]["error_over_that_under_the_records"] = float(np.exp(np.mean(np.log(ratios)))) if ratios else None
                    row[label]["of_each"] = {n: a / b for (a, b), n in zip(totals, names) if b > 0}          # which known parameter the change helps (round 3 of the cold tests)
            first = "with_the_term" if "with_the_term" in row else "with_the_working_model"
            row["error_over_that_under_the_records"] = row[first]["error_over_that_under_the_records"] if first in row else None
            row["of_each"] = row[first]["of_each"] if first in row else {}
            if row["error_over_that_under_the_records"] is not None and lever["kind"] in ("samples", "span", "repeats", "covariate"):          # each error is read on `sets` sets of noise of their own, so the ratio of two moves by about one over the root of that number: 0.90 at 30 sets read 0.71 at 100 in a cold test (round 3). A noise lever scales the same draws, and a characterization reads the same sets, so theirs does not move
                row["known_within"] = row["error_over_that_under_the_records"] / math.sqrt(max(sets - 1, 1))
            if first not in row:
                row["why_not_measured"] = "the model shares no parameter and no term was proposed for it, so nothing measures the lever"
        experiment.append(row)
    ranked = sorted([r for r in experiment if r.get("error_over_that_under_the_records") is not None], key=lambda r: r["error_over_that_under_the_records"])
    first = chosen["term"] if chosen else None          # the term that stands in the twin of the levers, named from the same choice
    for r in model:
        r.pop("_term", None)
    return {"for_the_model": model, "for_the_experiment": ranked[:3], "every_lever": experiment, "width_of_the_curves": unit,
            "the_term_in_the_twin_of_the_levers": first, "the_noise_leads": leading[0]["answer"] if noisy else None,
            "under_the_records": None if here is None else {k: v for k, v in here.items() if k != "sets"},
            "fewer_than_three": [f"{len(model)} for the model: the report and the dictionary hold no more"] * (len(model) < 3)
            + [f"{len(ranked[:3])} for the experiment: the description holds {sum(not r['from_the_package'] for r in experiment)} levers and the package added {sum(r['from_the_package'] for r in experiment)}, of which {len(ranked)} were measured"] * (len(ranked) < 3),
            "kinds_proposed": sorted({r["his_kind"] for r in ranked[:3]})}
