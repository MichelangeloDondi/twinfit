"""The loop without an agent, the owner's message 140: "The non agentic version is the same, the only difference is that
the decision of which term of the model to add would be done by the senior researcher after the statistical machinery
has provided where to search new terms of the model".

    python -m twinfit loop DIR --n N [--sets S]
    run(DIR, n=N, sets=S)                    the same, from Python, with its summary returned
    run(DIR, n=N, sets=S, report=False)      only what keeps or refuses a term: no report of where to search and no
                                             proposals, for a count of many draws. Each model says that it was not read

DIR holds the description, `description.toml` (description.py), whose addenda name the file of the initial model, and
the terms that the researcher writes, one file each, `terms/01.py`, `terms/02.py` and so on. The file of the model
defines `model(x, shared, own, design)`, `SHARED`, `OWN`, `START_SHARED` and `start_own(record)`, and may define
`centre(shared, own, design)` with `CENTRE`, the name of the own parameter that places a curve, `bounds_own(record)`,
`BOUNDS_SHARED`, `NAME` and `HELD`, a dictionary of the values that the model reads and does not fit. An ansatz of the
description sets a value of `HELD`, or a covariate that every record gives at its value, and carries its uncertainty
into the bar of theory of every other value of the output. A term file defines `term`, a
`twinfit.mechanisms.Mechanism`, and `why`, the researcher's reason in a sentence of six words or more.

Each run reads every file again from the start and remakes the same analysis, so the folder alone holds it. Its record
is the notebook in `DIR/loop/notebook`, which is written again at each run. Each model is a rung of the notebook:

  1. the model is fitted, and its twin, whose prediction is stamped before it runs, admits the fit or does not
  2. the report of where to search and the six proposals, in `DIR/loop/model_N.json` and `DIR/loop/model_N.md`
  3. the next term file, if it is there: its term is fitted from the fit and read on synthetic sets that hold it, then
     kept when its fit ends inside its bounds, its Bayesian criterion is lower than that of the model before by more than
     best.TIED, and the twin admits the larger model. A term refused leaves the model as it was, with its reason, and the
     next file is read against the same model

The loop stops at the first term file that is not there, and says which file to write, or after `n` term files. Its
output is that of best.py on the models kept, in `DIR/loop/best.json` and `DIR/loop/best.md`.
"""
from __future__ import annotations

import fcntl
import importlib.util
import inspect
import json
import math
import os
import time
import re
import shutil
from pathlib import Path

import numpy as np

from .best import TIED, TOLD_APART, best, own_bars
from .candidates import ALONG
from .description import read as read_description
from .fit import fit
from .identify import identify
from .linear import Linear
from .noise import check as noise_check
from .notebook import FIT, MIN_SETS, Notebook
from .problem import Problem, Record
from .proposals import noise_leads, proposals
from .pool import asking, how_many
from .refine import carried, develop_one, held_on_the_twin, on_the_twin
from .twin import ENDED, admission, bands
from .where import TIE as TIE_PLACES, where

NEEDED = ("model", "SHARED", "OWN", "START_SHARED", "start_own")
GENERIC = Path(__file__).resolve().parent / "agent" / "model_template.py"          # the first model where the description names none: one peak or dip of one width, on its own height, place and baseline
ACTIONABLE = ("where", "shape", "what", "changes along", "matched")          # the questions whose answer points to where a term is missing: seen by says how a misfit is read, and the noise is a question of its own
NOISE_KINDS = {"a noise that grows with the signal": "Where the description gives no `[noise]`, the package reads its law in the height of the curve from neighbouring samples, the default, from = \"neighbours\". On a curve sharper than three samples read it from what the first model leaves, from = \"residuals\", and compare the chi square a degree of freedom of the two. "
                                                     "Where you know it, give it as a column of the curves, from = \"column\" and column = its name, or as a law of your own in model.py, from = \"model\". "
                                                     "For a counter it is the root of the fitted counts and not of the counts of the records, which move with their own noise: section 9 of the tutorial gives the five ways and their keys.",          # the page named a column and a counter alone, and the way that held on a trace of a photomultiplier was in the tutorial alone (the fresh-machine red team of 0.1.0rc14)
               "glitches": "The residuals hold samples far from the curve, more than the noise allows. A glitch is such a sample alone, in one record and among its neighbours: "
                           "where they stand at one abscissa in several records, or in one band of the height of the curve, they are no glitches, and the law of the noise there is too small or the model misses there, "
                           "which the check of the noise by bands of the height tells. Look at the samples of largest whitened residual, which the report of the model names with how many share their abscissa with another record, and where the instrument explains them, "
                           "give those samples a larger noise in a column of the curves, or leave them out of the files.",          # the report named as glitches 20 samples at the same abscissas in four records, where the fitted law of the noise dipped (the blind red team of 0.1.0rc13)
               "a noise correlated from sample to sample": "The package draws the noise of each sample apart from the others, so the bars of the fit are too small where it is "
                                                           "correlated, as a filter of the detector leaves it. Keep one sample in every few, as many as the correlation spans, "
                                                           "so that the noise of those kept is apart, or read the bars as too small."}          # the kinds of the report that are read in the noise and not in the mean
NO_BAR = ("the fit gives it no finite bar, as for a parameter at a bound of its search or one that the records cannot tell from another. "
          "Look at its bounds in the model, and at what else in the model moves the curve as it does")
UNTESTED = ("The twin draws its synthetic sets at the noise that the description gives, so it takes that noise as right: its coverage and its bias cannot see a noise "
            "stated wrong. The residuals against the height of the fitted curve say whether it is.")


def glitch_samples(problem, result, beyond: float = 5.0, most: int | None = 20) -> list:
    """The samples whose whitened residual, the residual over the noise of its sample, is beyond `beyond` in size, the
    largest first and at most `most`: each with its curve, its abscissa and that residual. The report named glitches
    and told the researcher to look at the samples of largest residual, and a blind agent found its 16 with numpy of
    its own (round 7 of the cold tests)."""
    out = []
    for k, rec in enumerate(problem.records):
        r = (rec.y - problem.predict(k, result.shared, result.own[k])) / rec.sd
        out += [(rec.name, float(x), float(v)) for x, v in zip(rec.x, r) if abs(v) > beyond]
    return sorted(out, key=lambda s: -abs(s[2]))[:most]          # `most` None: every one


def shared_abscissas(problem, far: list) -> int:
    """How many of the samples of `far`, as glitch_samples gives them, stand within one step of the samples of an
    abscissa where a sample of another record is in `far` too: a glitch is alone, and samples far from the curve at one
    place in several records say that the noise or the model misses there (the blind red team of 0.1.0rc13)."""
    step = min((float(np.median(np.diff(np.sort(rec.x)))) for rec in problem.records if len(rec.x) > 1), default=0.0)
    return sum(any(other != name and abs(x - y) <= step for other, y, _ in far) for name, x, _ in far)


def correlation_span(problem, result, under: float = 0.05, most: int = 200) -> tuple:
    """The correlation of neighbouring whitened residuals at one sample, and the number of samples at which it falls under
    `under` in size, each the median over the curves, the second None where it does not fall within `most`. The report
    said keep one sample in every few, as many as the correlation spans, and printed no span: a blind agent computed it
    from residuals_0.csv by numpy of its own (round 9 of the cold tests)."""
    first, spans = [], []
    for k, rec in enumerate(problem.records):
        r = (rec.y - problem.predict(k, result.shared, result.own[k])) / rec.sd
        r = r - r.mean()
        total = float(r @ r)
        if r.size < 3 or total <= 0.0:
            continue
        rho = [float(r[:-j] @ r[j:]) / total for j in range(1, min(most, r.size - 1) + 1)]
        first.append(rho[0])
        spans.append(next((j + 1 for j, c in enumerate(rho) if abs(c) < under), None))
    if not first:
        return None, None
    found = sorted(x for x in spans if x is not None)
    return float(np.median(first)), (found[len(found) // 2] if len(found) == len(spans) else None)




THIN_KEPT = 10          # the fewest samples of a curve that a copy of thin keeps: thin refuses fewer, and the report does not advise it
MISFIT = 2.0          # the chi square for each degree of freedom beyond which a correlation of neighbouring residuals is read as the misfit of the mean, and no thinning is advised


def thin_said(problem, span: int, how: str) -> str:
    """The clause of a page that advises thin over `span` samples, `how` its words, or why it does not: where one
    sample in every `span` keeps fewer than THIN_KEPT of the shortest curve, a copy cannot be fitted, and thin refuses
    it. A report advised `thin DIR 180` on curves of 1081 samples, which kept 6 (round 18 of the cold tests)."""
    shortest = min(len(r.x) for r in problem.records)
    if shortest // span < THIN_KEPT:
        return (f"falls under 0.05 at {span} samples, where one sample in every {span} keeps {shortest // span} of the {shortest} of the shortest curve, too few to fit, "
                "so thinning cannot tell a noise correlated that far from a term that the curves lack. ")
    return f"falls under 0.05 at {span} samples, {how.format(span=span)}. "


def misfit_first(result) -> str | None:
    """The sentence that a correlation of neighbouring residuals is mostly the misfit of the mean, where the chi square for
    each degree of freedom passes MISFIT, in the place of the advice to thin. At 40 a degree of freedom the report advised
    keeping one sample in every 180, six of a curve of 1081, for a correlation that a missing term of the mean made
    (round 18 of the cold tests, both the field and the blind agent)."""
    per = result.chi2 / max(result.dof, 1)
    if per <= MISFIT:
        return None
    return (f"At a chi square of {per:.3g} for each degree of freedom the correlation of neighbouring residuals is mostly the misfit of the mean, "
            "which a missing term makes: write the term of the mean that the report puts first, and read the correlation again before any thinning.")


def tied_noise(answers: list, problem, result) -> list:
    """The lines of a page where a noise correlated from sample to sample stands beyond the twin but does not lead: the
    span of its correlation and the way to thin. Only the lead was given them, and with the correlated kind tied the
    three proposals were terms that hold the noise (round 17 of the cold tests)."""
    if noise_leads(answers) or not any(a["question"] == "what" and a["answer"] == "a noise correlated from sample to sample" and a["beyond_the_twin"] for a in answers):
        return []
    rho, span = correlation_span(problem, result)
    if rho is None:
        return []
    if misfit_first(result):
        return ["", f"A noise correlated from sample to sample stands beyond the twin too, though not first. {misfit_first(result)}"]
    return ["", f"A noise correlated from sample to sample stands beyond the twin too, though not first: the correlation of neighbouring whitened residuals is {rho:.2f} at one sample, the median over the curves, and "
            + (thin_said(problem, span, "so `python -m twinfit thin DIR {span}` makes copies that keep one sample in every {span}") if span else "does not fall under 0.05 within 200 samples. ")
            + "A term proposed below that the thinned copies refuse holds the noise and not the curve."]


def noise_said(problem, result) -> dict:
    """The check of the noise that the twin cannot make (review 21): the mean squared whitened residual in six bands of
    the height of the fitted line, which is one under a right model with a right noise, and whether each band is one
    within three standard errors. Where it is not, the bars of the fit are wrong by about the root of that number. Where
    many samples share one height, as on the baseline of a curve whose wings are flat, some bands are empty and fewer are
    read: the sentence says how many. It said six where four were printed (the cold tests of row M169)."""
    rows = noise_check(problem, result)
    right = all(abs(r["mean_squared_residual"] - 1.0) <= 3.0 * r["standard_error"] for r in rows)
    fewer = "" if len(rows) == 6 else "Fewer than six bands are read where many samples lie at one height, as on a flat baseline. "
    said = (f"The mean squared whitened residual, the residual over the noise of its sample squared, in {len(rows)} bands of the height of the fitted curve above its lowest value, from low to high: "
            + ", ".join(f"{r['mean_squared_residual']:.2f}" for r in rows) + ", where 1 is right. " + fewer
            + ("Each is 1 within three standard errors." if right else "Not every band is 1 within three standard errors: there the noise that the description gives is not that of "
               "the records, or the model misses. Where the report finds an answer of the curve beyond the twin, a place, a kind or a setting, read it first as the model that misses, and read this check again once the term "
               "is in. Where the answer that leads is one of the noise, a noise that grows with the curve or glitches, the noise comes first, as the section of the proposals says. Where nothing stands beyond the twin and a band still is not 1, the noise is wrong there: a band above 1 says that the stated noise is too small at that height, "
               "and the bars of the fit too small by about the square root of the band, and a band under 1 the reverse. `[noise] from = \"residuals\"` in the description reads it from the second "
               "differences of what the first model leaves, which a term that the model lacks, smooth over three samples, does not reach, or give a law of your own."))          # the sentence asked a reading that no key gave (the blind red team of 0.1.0rc16, mislead 5)
    if not right:          # each band outside with its standard error: "best.md says how" stood over the general paragraph alone, and bands of 1.11 to 1.14 read as a failure with no measure of how far (the fresh red team of 0.1.0rc19, mislead 6)
        out = [r for r in rows if abs(r["mean_squared_residual"] - 1.0) > 3.0 * r["standard_error"]]
        said += " The bands outside: " + ", ".join(f"{r['mean_squared_residual']:.2f} +- {r['standard_error']:.2f}, {abs(r['mean_squared_residual'] - 1.0) / r['standard_error']:.1f} of its standard errors from 1" for r in out) + "."
    if not right:          # bands read as independent while the same report measured a correlation of 0.70 at one sample (round 15 of the cold tests)
        rho, span = correlation_span(problem, result)
        if rho is not None and abs(rho) > 0.3:
            said += (f" The whitened residuals of neighbouring samples correlate at {rho:.2f}, so the standard errors of these bands, which take the samples apart, are too small"
                     + (f" by about the root of {span}, the samples over which the correlation falls" if span else "") + ": read them again on the samples thinned, and not as a noise stated wrong.")
    return {"bands": rows, "right": right, "said": said}


def equation_of(source) -> str:
    """The equation of the first model: EQUATION of its file, in mathtext, or else the expression that its function
    `model` returns, as it is written. The visual log writes it at each turn, his message 309."""
    said = getattr(source, "EQUATION", "")
    if said:
        return str(said)
    try:
        lines = inspect.getsource(source.model).splitlines()
    except (OSError, TypeError, AttributeError):
        return ""
    back = [line.strip()[len("return "):] for line in lines if line.strip().startswith("return ")]
    return back[-1] if back else ""


def effects_of(iterations: list, index: int) -> tuple:
    """What each term of model `index` does to its curves at its fit: the model at the fit less the same model with the
    parameter of that term at 0, where every term of the package is absent, all else held. The first model has no term,
    and each of its shared parameters is moved by one bar of its fit instead. Returns the rows that name each effect and
    the arrays, x{r} the abscissa of record r and e{j}_{r} effect j on it, in the units of the values."""
    it = iterations[index]
    problem, result = it["problem"], it["result"]
    held = it.get("held") or {}          # a term held at its theory adds no shared parameter: its cell is set to 0 and back (the cold round 22: a held term had no panel)
    by = {row["model"]: row for row in iterations}
    chain, at = [], it
    while at.get("after"):
        if at.get("parameter") in problem.shared or (at.get("parameter") in held and callable(held[at["parameter"]][2])):
            chain.insert(0, at)
        at = by.get(at["after"], {})
    shared, own = np.asarray(result.shared, dtype=float), [np.asarray(row, dtype=float) for row in result.own]
    bars = result.sd
    moves = [(row["added"], row["parameter"], problem.shared.index(row["parameter"]) if row["parameter"] in problem.shared else None, "at zero") for row in chain] or \
            [(name, name, j, "one bar") for j, name in enumerate(problem.shared)]
    rows, arrays = [], {f"x{r}": np.asarray(rec.x, dtype=float) for r, rec in enumerate(problem.records)}

    def change(place, how, value):
        """The change of every curve when the parameter at `place` is moved: at zero, or by one bar."""
        moved = shared.copy()
        moved[place] = value
        out = {}
        for r in range(len(problem.records)):
            fitted = np.asarray(problem.predict(r, shared, own[r]), dtype=float)
            out[r] = fitted - np.asarray(problem.predict(r, moved, own[r]), dtype=float) if how == "at zero" else np.asarray(problem.predict(r, moved, own[r]), dtype=float) - fitted
        return out
    for j, (label, name, place, how) in enumerate(moves):
        try:
            if place is None:          # a held term: the cell of its theory at 0, then back at its centre
                centre, bar, setting = held[name]
                fitted = [np.asarray(problem.predict(r, shared, own[r]), dtype=float) for r in range(len(problem.records))]
                setting(0.0)
                try:
                    got = {r: fitted[r] - np.asarray(problem.predict(r, shared, own[r]), dtype=float) for r in range(len(problem.records))}
                finally:
                    setting(centre)
                value, spread, how = float(centre), float(bar), "held"
            else:
                value, spread = float(shared[place]), float(bars[place])
                got = change(place, how, 0.0 if how == "at zero" else value + spread)
                if how == "at zero" and not all(np.all(np.isfinite(v)) for v in got.values()):          # a term whose parameter at zero leaves no curve, as a tail of scale 0 that replaces the line: moved by one bar instead (the cold round 22: a panel of nan)
                    got, how = change(place, "one bar", value + spread), "one bar"
            if not all(np.all(np.isfinite(v)) for v in got.values()):
                continue
        except (ValueError, FloatingPointError, IndexError, ZeroDivisionError):
            continue
        arrays.update({f"e{j}_{r}": v for r, v in got.items()})
        rows.append({"name": label, "parameter": name, "value": value, "bar": spread, "how": how, "key": j})
    return rows, arrays


def write_equations(out_dir: Path, iterations: list) -> Path:
    """equations.json of a run: for each model of the loop, the model it extends, the term it adds, the formula of that
    term, mathtext with f the model before it and a its parameter, the names of its shared and own parameters, its
    records with the median noise of each, and what each of its terms does to its curves, whose arrays are
    effects_{k}.npz beside it. Written after each model, so that the visual log of a run that was stopped still has the
    models it reached."""
    out_dir = Path(out_dir)
    k = len(iterations) - 1
    rows_of_effects, arrays = effects_of(iterations, k)
    np.savez(out_dir / f"effects_{k}.npz", **arrays)
    path = out_dir / "equations.json"
    rows = json.loads(path.read_text()) if path.is_file() and k else []
    it = iterations[k]
    at = (it.get("held") or {}).get(it.get("parameter"))          # a term held at its theory: its centre and bar, for its line of the equation (the cold round 22)
    rows = rows[:k] + [{"model": it["model"], "after": it.get("after"), "added": it["added"], "formula": it.get("formula", ""), "parameter": it.get("parameter"),
                        "held_at": [float(at[0]), float(at[1])] if at is not None and it.get("parameter") not in it["problem"].shared else None,
                        "file": it.get("file"), "letters": it.get("letters") or {}, "shared": list(it["problem"].shared), "own": list(it["problem"].own),
                        "records": [rec.name for rec in it["problem"].records], "noise": [float(np.median(rec.sd)) for rec in it["problem"].records],
                        "effects": rows_of_effects}]
    path.write_text(json.dumps(rows, indent=1) + "\n")
    return path


def _module(path: Path, label: str, held: dict | None = None, names=None):
    """A file of the folder as a module. `held`: the dictionary HELD of the model file, which a term file then reads by
    the same name, so that an ansatz on a number of a term reaches the term as it reaches the model. A blind agent
    reached it through the globals of the model it wrapped (round 4 of the cold tests). `names`: the module of the model
    file, whose public names, its functions and its constants, a term file then reads by their names, with no import: a
    tester copied the response of his instrument into a term file, where `import model` was refused (round 14)."""
    spec = importlib.util.spec_from_file_location(label, path)
    module = importlib.util.module_from_spec(spec)
    if names is not None:
        module.__dict__.update({k: v for k, v in vars(names).items() if not k.startswith("_")})
    if held is not None:
        module.HELD = held
    try:
        spec.loader.exec_module(module)
    except (NameError, ImportError, SyntaxError, AttributeError, TypeError, IndexError) as err:          # a file that fails as it is read is refused in words, with no trace of the code
        hint = (". A term file reads HELD only where model.py defines it" if isinstance(err, NameError) and "HELD" in str(err)
                else ". A term file reads the names of model.py, its functions and its constants, by their names and with no import: write them so" if isinstance(err, ImportError) and getattr(err, "name", None) == "model" else "")
        raise ValueError(f"{Path(path).name} could not be read: {type(err).__name__}: {err}{hint}") from None
    return module


def list_of_terms(d, source=None) -> tuple:
    """The researcher's list of terms: the file that `list_of_terms` names under [addenda] of the description, whose
    `TERMS` is a list of terms written as a term file writes one, Mechanisms of the package or of one's own. `init`
    writes it empty, and the package holds none: with no list the report reads the residuals alone (the owner's message
    342: "For the list of terms, allow the user to define it, but keep it empty at the beginning"). A file that is
    missing, or that holds anything but terms, or two terms of one name, is refused."""
    named = d.addenda.get("list_of_terms")
    if not named:
        return ()
    path = d.home / str(named)
    if not path.is_file():
        raise ValueError(f"the description names {named} as the list of terms, and the folder holds no such file: write it with TERMS = [], or take the line out of [addenda]")
    from .mechanisms import Mechanism
    got = getattr(_module(path, "list_of_terms", getattr(source, "HELD", None), source), "TERMS", None)
    if not isinstance(got, (list, tuple)) or not all(isinstance(m, Mechanism) for m in got):
        raise ValueError(f"{named} must hold TERMS, a list of terms written as a term file writes one, Mechanisms of the package or of your own: it holds {type(got).__name__}")
    names = [m.name for m in got]
    if len(set(names)) < len(names):
        raise ValueError(f"{named} names {sorted({n for n in names if names.count(n) > 1})[0]} twice: the report tells the terms of the list by their names")
    return tuple(got)


def initial(folder: Path, with_module: bool = False):
    """The description of the folder, the problem of its initial model on the records that the description gives, and
    the name of the own parameter that places a curve, or None. With `with_module`, the module of the model file too:
    its `HELD` is the one that the model reads, whatever the model is, a function, a partial or an object (review 22,
    minor 19)."""
    d = read_description(folder / "description.toml")
    filled = [key for key in ("name", "what") if "TO FILL" in str(getattr(d, key, ""))]
    if filled:          # a project whose name and what still read TO FILL ran with no word, and every report quotes them (round 16 of the cold tests)
        raise ValueError(f"description.toml still reads TO FILL in {' and '.join(filled)}: the reports quote them, so write the problem in words there first")
    if not d.addenda.get("model") and (d.home / "model.py").is_file():          # a model of one's own left beside a description that does not name it would be passed over for the generic one
        raise ValueError("the folder holds model.py and the description names no first model: write model = \"model.py\" under [addenda], or remove model.py to start from the generic model of the package")
    path = d.home / d.addenda["model"] if d.addenda.get("model") else GENERIC          # the description alone, with no model: the generic model of the package, row M205 and his message 140
    m = _module(path, "initial_model")
    lacking = [name for name in NEEDED if not hasattr(m, name) and not (name == "START_SHARED" and callable(getattr(m, "start_shared", None)))]
    if lacking:
        raise ValueError(f"the file of the initial model defines no {', '.join(lacking)}")
    named, placed = getattr(m, "CENTRE", None), getattr(m, "centre", None)
    if named is not None and named not in m.OWN:
        raise ValueError(f"the model names {named!r} as its CENTRE, and its own parameters are {', '.join(m.OWN)}: CENTRE is the own parameter that places a curve")
    if (named is None) != (placed is None):
        raise ValueError("the model gives " + ("CENTRE and no function centre" if placed is None else "a function centre and no CENTRE") + ": the report reads the core and the wings of a curve from both, or from neither")
    bounds = getattr(m, "BOUNDS_SHARED", None)
    records = d.records()
    starts = m.start_shared(records) if callable(getattr(m, "start_shared", None)) else m.START_SHARED          # starts read from the curves, where the model reads them: a width started at 1.0 for a line 0.016 wide moved no curve (round 15)
    try:
        out = (d, Problem(records, m.model, m.SHARED, m.OWN, starts, m.start_own, getattr(m, "centre", None), getattr(m, "NAME", "the initial model"),
                          bounds, getattr(m, "bounds_own", None)), getattr(m, "CENTRE", None))
    except ValueError as err:          # "one start for each shared parameter" with no file and no count, after a pasted model kept the bounds of the one before (round 13 of the cold tests)
        if "start" not in str(err) and "bounds" not in str(err):
            raise
        held = "none" if bounds is None else " and ".join(str(len(part)) for part in bounds) if all(hasattr(part, "__len__") for part in bounds) else "not two lists"
        raise ValueError(f"{path.name}: {err}: SHARED names {len(m.SHARED)}, the starts of the shared parameters hold {len(starts)}, and the lower and upper bounds of BOUNDS_SHARED hold {held}. "
                         "Give one start and one lower and one upper bound for each shared parameter, or leave BOUNDS_SHARED out") from None
    return out + (m,) if with_module else out


def far_from_zero(problem) -> str:
    """Why a fit may not end, where the abscissa sits far from zero against the span of the curves: the fit takes the
    derivative by a parameter with a step that grows with the size of the parameter, so the place of a curve at 587.6
    nm moved by steps of 0.006 nm, a third of the width of a line 0.02 nm wide, and the first fit of a cold test did not
    end (round 5). An empty string where the abscissa is near zero."""
    lo, hi = min(float(np.min(r.x)) for r in problem.records), max(float(np.max(r.x)) for r in problem.records)
    middle, span = 0.5 * (lo + hi), (hi - lo) or 1.0
    if abs(middle) < 3.0 * span:          # it said fifty spans, and a first fit failed at 13.5 and held at 6.8 (round 11 of the cold tests)
        return ""
    return (f". The abscissa sits at {middle:.6g}, {abs(middle) / span:.0f} times its span from zero, the span from the least abscissa of the curves to the largest, and the fit then may move the place of a curve by steps too coarse "
            f"for its width: write abscissa_offset = {middle:.6g} under [problem] of the description, and the model then reads the abscissa less it")


AGAIN = 4          # the searches more that a first fit gets, each from where the last one stopped
SEARCH = 80          # the steps of one search, the default of the fit
STEPS = SEARCH * (AGAIN + 1)          # the steps that a synthetic set of a twin gets in one search, as many as the first fit of the records gets in its searches


def searched(problem):
    """The fit of a first model from its starts, and again from where each search stopped while it has not ended, up to
    AGAIN times. A search takes 80 steps at most: a first model whose one centre stands more than a width from the curves
    of the higher pressures came to its minimum in 161, and was refused as a fit that did not end, whatever its starts
    and its bounds (round 19 of the tests of the tutorial, the acceptance of release 0.0.5)."""
    result = fit(problem, max_iter=SEARCH)
    for _ in range(AGAIN):
        if result.converged:
            break
        result = fit(problem, shared0=result.shared, own0=result.own, max_iter=SEARCH)
    return result


def why_no_end(problem) -> str:
    """What may keep a first fit from ending, read from its starts and its records, beyond an abscissa far from zero: a
    start of the place of a curve near an end of its abscissa, as the template's start takes the end of a sloping
    baseline for a dip, and a curve so narrow in the unit of the abscissa that the search, which steps each own
    parameter by at least a hundred millionth of the largest own value of its curve, moves its place by a share of its
    width. A newcomer met "look at its starts and its bounds" twelve times, for a line 0.0015 nm wide of 3000 counts
    and for a start on the edge of a slope, and the refusal named neither (round 15 of the cold tests)."""
    said = []
    ends, narrow, outside = [], [], []
    start = np.asarray(problem.start_shared, dtype=float) if np.size(problem.start_shared) else np.zeros(0)
    for k, rec in enumerate(problem.records):
        x, y = np.asarray(rec.x, dtype=float), np.asarray(rec.y, dtype=float)
        try:
            own = problem.own_start(k)
            place = problem.centre_of(k, start, own) if getattr(problem, "centre", None) is not None else None
        except Exception:          # a start that cannot be read says nothing here: the refusal before names it
            continue
        span = float(np.ptp(x)) or 1.0
        if place is not None and not x.min() - span <= place <= x.max() + span:
            outside.append((rec.name, place, float(x.min()), float(x.max())))
        elif place is not None and min(abs(place - x.min()), abs(place - x.max())) < 0.05 * span:
            ends.append(rec.name)
        up = np.abs(y - np.median(y))
        wide = float(np.ptp(x[up >= 0.5 * up.max()])) if up.max() > 0 else span
        if own.size and wide > 0 and 1e-5 * 1e-3 * float(np.max(np.abs(own))) > 1e-3 * wide:
            narrow.append((rec.name, wide, float(np.max(np.abs(own)))))
    if outside:          # an offset declared and the starts left in the abscissa as recorded: "did not end in 5 searches" with no word of the frame (round 20 of the tests of the tutorial)
        name, place, lo, hi = outside[0]
        said.append(f"The starts put the place of {name} at {place:.6g}, outside its abscissa, which runs from {lo:.6g} to {hi:.6g}: where the description gives abscissa_offset, "
                    "the model reads the abscissa less it, so the starts and the bounds of a place are written less it too")
    if ends:
        said.append(f"The start that start_own gives puts the place of {', '.join(ends[:3])} near an end of its abscissa: on a sloping baseline the sample farthest from the median can be an end of the curve, "
                    "so start from the largest distance from a straight line through the two ends of each curve")
    if narrow:
        name, wide, top = narrow[0]
        said.append(f"The curve {name} is about {wide:.3g} wide in the unit of the abscissa, where its largest own value is {top:.3g}, and the search steps its place by at least {1e-8 * top:.3g}, a share of that width: "
                    "write the abscissa in a unit in which a curve is about 1 wide, or the values in one in which its height is about 1")
    return "".join(f". {s}" for s in said)


def held_of(ansatz: dict, problem, values: dict | None = None) -> dict:
    """The values that the description holds in its ansatz, as best.uncertainty takes them, {name: (centre, bar, at)}:
    at(value) is the problem of the model with that value in its place. It is a covariate that every record gives at the
    value of the ansatz, set to the new value in every record, or a value that the model file reads from its dictionary
    `HELD`, here `values`, which the model reads when it is computed. An ansatz on anything else is refused: it would
    carry no bar onto the other values, and the output would say that nothing is held."""
    out = {}
    for name, a in (ansatz or {}).items():
        centre, bar = float(a["value"]), float(a["uncertainty"])
        if all(name in rec.design for rec in problem.records):
            given = sorted({float(rec.design[name]) for rec in problem.records})
            if any(abs(v - centre) > 1e-9 * max(abs(centre), 1.0) for v in given):
                raise ValueError(f"the ansatz holds {name} at {centre:g}, and the records give {', '.join(f'{v:g}' for v in given)}: a covariate held by an ansatz has its value in every record")
            out[name] = (centre, bar, lambda v, name=name: problem.with_records([Record(r.x, r.y, r.sd, dict(r.design, **{name: v}), r.name) for r in problem.records]))
        elif values is not None and name in values:
            def at(v, name=name):
                values[name] = v
                return problem
            out[name] = (centre, bar, at)
        else:
            raise ValueError(f"the ansatz on {name!r} names neither a covariate of the records nor a value that the model file holds in its dictionary HELD: it would put no bar on the others")
    return out


def _criterion(problem, result) -> float:
    return float(result.chi2 + (problem.samples - result.dof) * math.log(problem.samples))


def held_term(problem, result, term):
    """The working model with a term held at the centre that the theory gives, `term.theory = (centre, bar)`, and its
    fit: no parameter is added. The value sits in a cell that `at` of the output moves one bar either way, through every
    model built on this one, so that the bar of the theory reaches each other value (rule 10, his messages 133 and 135).
    Returns the problem, its fit, and the entry of the held values, {parameter: (centre, bar, cell)}."""
    if term.theory is None:
        raise ValueError(f"{term.name}: a term held by the loop needs the centre and the bar that the theory gives, as `theory=(centre, bar)`")
    centre, bar = (float(v) for v in term.theory)
    cell = {"value": centre}
    n = len(problem.shared)
    wrapped = term.wrap(problem.model, problem, n)

    def model(x, shared, own, design):
        return wrapped(x, np.append(np.asarray(shared, dtype=float), cell["value"]), own, design)
    held = Problem(problem.records, model, problem.shared, problem.own, list(result.shared), [np.array(row, dtype=float) for row in result.own], problem.centre,
                   f"{problem.name} with {term.name} held at {centre:g}", problem._bounds_shared, problem._bounds_own)
    return held, fit(held, result.shared, result.own), {term.parameter: (centre, bar, cell)}


def held_terms(holding: dict, problem) -> dict:
    """The held terms as the output takes held values, {name: (centre, bar, at)}: `at(value)` sets the cell of the term
    and gives the problem of this model, whose model reads the cell."""
    return {name: (centre, bar, (lambda value, c=cell, p=problem: (c.update(value=float(value)), p)[1])) for name, (centre, bar, cell) in holding.items()}


def residuals(problem, result, path: Path) -> None:
    """The residuals of a fit, one line a sample: the curve, its abscissa, its value, the fitted value, the residual and
    its noise. A tester rebuilt the fitted curves by hand to estimate the noise from the wings (round 5 of the cold
    tests)."""
    import csv
    with open(path, "w", newline="") as f:
        out = csv.writer(f)
        out.writerow(["curve", "abscissa", "value", "fitted", "residual", "noise"])
        for k, rec in enumerate(problem.records):
            fitted = np.asarray(problem.predict(k, result.shared, result.own[k]), dtype=float)
            out.writerows([rec.name, f"{x:.10g}", f"{y:.10g}", f"{m:.10g}", f"{y - m:.6g}", f"{s:.6g}"] for x, y, m, s in zip(rec.x, rec.y, fitted, rec.sd))


def _rung(book, rung, by, problem, result, d, sets, folder, k, centre, report: bool = True, held: dict | None = None, exclude=(), refused=(), reuse=None, code: str = "", terms=(), dictionary: bool = False) -> dict:
    """One model as a rung of the notebook: its fit, its twin, its report of where to search and its six proposals. With
    `reuse`, a `twinfit.reuse.Reuse`, each costly step is read back where a run before computed it from the same `code`,
    the digest of the files that make the model, and the same records and fit."""
    def again(label, compute, *parts):
        return compute() if reuse is None else reuse(label, compute, code, problem, result, sets, *parts)
    residuals(problem, result, folder / f"residuals_{k}.csv")
    lin = Linear(problem, result)
    book.open_rung(rung, by=by)
    book.step("intake", rung, {"records": len(problem.records), "samples": problem.samples})
    book.step("identify", rung, identify(lin) if problem.shared else {"parameters": {}, "weakest": None, "why": "the model has no shared parameter: there is nothing to identify"})
    stamp = f"twin of {rung}"
    book.stamp(stamp, "The fit of this model returns, on its own synthetic sets at the noise of the records, the values that made them, "
                      "within the bands of the rule of the twin: its bias, its spread and its coverage.", rung)          # a name of the code in a page for the researcher (round 9)
    adm = again("the twin of a model", lambda: admission(problem, result, sets=sets, start="fit", workers=how_many(None, sets)))          # the sets of the twin side by side, as the proposals
    book.step("twin", rung, {"admission": adm}, stamp=stamp)
    book.grade(stamp, holds=adm["admitted"], what_came=adm["rule"] + (". Admitted" if adm["admitted"] else ". Not admitted"))
    book.admit(FIT, rung)
    book.step("fit", rung, {"shared": dict(zip(problem.shared, result.shared.tolist())), "standard_deviations": dict(zip(problem.shared, result.sd.tolist()))})
    rho, span = correlation_span(problem, result)          # the bar of the chi square where neighbouring samples correlate: the visual log gave that of independent samples alone, 0.015 where the guide gave 0.096 (the blind red teams of 0.1.0rc19, mislead 2 and mislead 1)
    book.step("read", rung, {"chi_square": result.chi2, "degrees_of_freedom": result.dof, "correlation_of_neighbours": rho, "samples_of_the_correlation": span})
    if not report:
        book.step("propose", rung, {"candidates": [], "leading_set": {}, "not_read": "the loop was run for a count, with no report of where to search and no proposals"})
        table = {"model": rung, "terms_by": by, "chi_square": result.chi2, "degrees_of_freedom": result.dof, "criterion": _criterion(problem, result), "admitted_on_the_twin": adm["admitted"],
                 "rule_of_the_twin": adm["rule"], "measured_on_the_twin": adm["measured"], "bands_are_wider_than_stated": adm["bands_are_wider_than_stated"], "noise": noise_said(problem, result),
                 "parameters": _values(problem, result)}
        book.step("report", rung, {"what_it_does_not_show": "where to search and the proposals: this run reads only what keeps or refuses a term", "default_choice": None})
        (folder / f"model_{k}.json").write_text(json.dumps(table, indent=1, default=float) + "\n")
        return table
    listed = "blind" if not terms and not dictionary else "the researcher's list: " + "; ".join(m.name for m in terms) + ("; the package's dictionary" if dictionary else "")          # the key of the reuse names the list, whose file the code digest holds too
    report = again("the report of a model", lambda: where(problem, result, sets=sets, fixed=d.fixed(), centre=centre, terms=terms, dictionary=dictionary), d.fixed(), centre, listed)
    six = again("the proposals of a model", lambda: proposals(problem, result, report, d.levers(), sets=sets, held=held, exclude=exclude, refused=refused, setup=d.setup, terms=terms),
                report, d.levers(), held, exclude, refused, d.setup, listed)          # no proposal of a term that the model holds, and one refused before says so
    book.step("propose", rung, {"candidates": six["for_the_model"], "leading_set": report["first_of_each_question"], "where": report, "proposals": six})
    book.step("report", rung, {"what_it_does_not_show": "the cause of what is missing: the report names a place, a class, a covariate and a family of statistics, "
                                                      "and the term of physical origin that fills them is the researcher's choice", "default_choice": None})
    table = {"model": rung, "terms_by": by, "chi_square": result.chi2, "degrees_of_freedom": result.dof, "criterion": _criterion(problem, result),
             "admitted_on_the_twin": adm["admitted"], "rule_of_the_twin": adm["rule"], "measured_on_the_twin": adm["measured"], "bands_are_wider_than_stated": adm["bands_are_wider_than_stated"], "noise": noise_said(problem, result),
             "where": report, "proposals": six, "parameters": _values(problem, result), "mode": report["mode"]}
    (folder / f"model_{k}.json").write_text(json.dumps(table, indent=1, default=float) + "\n")
    if "people" in d.output:
        (folder / f"model_{k}.md").write_text(_page(table, problem, result))
    if "agents" in d.output:          # the output for agents, his message 316: `[output] for` in the description
        (folder / f"brief_{k}.md").write_text(brief(table, problem, result, k, folder))
    return table


ASKS = {"where": "the core of the curves, their flanks, their wings or their asymmetry",
        "shape": "how the residuals are shaped: a slow trend over the span of each curve, or a periodic structure",
        "what": "glitches, a noise that grows with the signal, or a noise correlated from sample to sample",
        "changes along": "which part of the residuals changes from curve to curve along a setting, read on the residuals themselves",
        "seen by": "the family of statistics that sees the misfit most",
        "noise": "the shape of the noise: skewed to one side, or with heavy tails",
        "matched": "which term of your list accounts for most of what every family reads, asked only where the description names a list"}          # the questions of the report, as the page defines them (row M167). The report holds no list of terms (the owner's messages 337 and 340)


def noise_from_residuals(problem, result, passes: int = 2):
    """The first model fitted again with the noise of each sample from the law that `twinfit.noise.from_residuals` reads
    in what the fit leaves, twice, from the fit before: `[noise] from = "residuals"` in the description. The records of
    every later model carry that noise, and the twin draws its sets with it."""
    from .noise import from_residuals
    for _ in range(passes):
        curves = [np.asarray(problem.predict(k, result.shared, result.own[k]), dtype=float) for k in range(len(problem.records))]
        sds = from_residuals([rec.y - c for rec, c in zip(problem.records, curves)], curves)
        problem = problem.with_records([Record(rec.x, rec.y, sd, rec.design, rec.name) for rec, sd in zip(problem.records, sds)])
        result = fit(problem, result.shared, result.own)
        if not result.converged:
            raise ValueError("the fit of the first model with the noise read from its residuals did not end: give the noise as a column or a value")
    return problem, result


def after_the_run(folder, out: dict) -> Path | None:
    """The page for an agent of the last model, brief_K.md, made right at the end of the run: written when its model was
    fitted, it could not know that a term tried after it was refused, and it named the next file by the model's number,
    so that an agent that read it alone wrote again a file already refused (the blind red team of 0.1.0rc13). It now
    names each file refused after that model with its reason, and the next free number of terms/."""
    loop_dir = Path(folder) / "loop"
    briefs = sorted(loop_dir.glob("brief_*.md"), key=lambda p: int(p.stem.split("_")[1]))
    if not briefs:
        return None
    last = briefs[-1]
    taken = [int(p.stem) for p in (Path(folder) / "terms").glob("[0-9][0-9].py")] if (Path(folder) / "terms").is_dir() else []
    free = f"terms/{max(taken, default=0) + 1:02d}.py"
    text = re.sub(r"terms/\d\d\.py", free, last.read_text())
    gone = [r for r in out.get("refused") or [] if r.get("file")]
    if gone:
        said = "Tried after this model and refused, not to be written again: " + "; ".join(f"terms/{r['file']}, {r['why']}" for r in gone) + "."
        text = text.replace("\nNext:", "\n" + said + "\n\nNext:", 1)
    last.write_text(text)
    return last


def _values(problem, result) -> dict:
    """The shared parameters of a fit, each with its value and its bar, as model_K.json keeps them."""
    return {n: [float(v), float(s)] for n, v, s in zip(problem.shared, result.shared, result.sd)}


def moved(folder, k: int, names) -> list:
    """How each parameter of model `k` stood at every model before it, read from the model_J.json of `folder`: a value
    that falls where a term enters is a share that the term took. A fresh agent a turn read only the last model, and one
    Gaussian width of 2.7 took two widths of the model together where the first alone is 0.87 (the attempts
    of his message 316 of the repository it is made in). A line a parameter that some earlier model holds."""
    if folder is None:
        return []
    past = {}
    for j in range(k):
        path = Path(folder) / f"model_{j}.json"
        held = json.loads(path.read_text()).get("parameters") if path.exists() else None
        for name, (v, s) in (held or {}).items():
            past.setdefault(name, []).append((j, v, s))
    return [f"- {name}: " + ", ".join(f"{' ± '.join(_pm(v, s))} at model {j}" for j, v, s in past[name]) for name in names if name in past]


def chi_bar(chi: float, dof: int, rho, span) -> str:
    """The chi square for each degree of freedom with the bar of a right model, as the visual log gives it, `log.chi_said`:
    that of independent samples, and where neighbouring residuals correlate, that bar times the root of the samples over
    which they correlate (the blind red teams of 0.1.0rc19, mislead 2 and mislead 1)."""
    alone = math.sqrt(2.0 / max(dof, 1))
    wide = f", ± {alone * math.sqrt(span):.2g} if its noise correlates over {span} samples as these residuals do" if rho is not None and abs(rho) > 0.3 and span and span > 1 else ""
    return f"chi square per degree of freedom {chi / max(dof, 1):.3f}, a right model gives 1 ± {alone:.2g} if its samples are independent{wide}."          # one standard deviation, the ± of every bar, where best.md gives three: "within" read as a tolerance (the readers of the red team of his message 434, agents R10), and a longer line cut the slide's first sentence (review 141)


def noise_first(table: dict) -> str:
    """The sentence that says the noise comes first, or nothing, by the rule of `noise_leads`, at the head of where to
    search and in the brief. The page said it only under the list of every answer, after a table that put the core
    first, and the brief said write terms/01.py with proposal 1, while the proposals for the experiment read another
    rule, the first kind of the question what whatever the misfit of the mean (the blind red team of 0.1.0rc19, misleads
    7 and 3)."""
    answers = (table.get("where") or {}).get("answers", [])
    noise = noise_leads(answers)
    if not noise:
        return ""
    got = next((a for a in answers if a.get("question") == "what" and a.get("answer") == noise), None)
    excess = f", excess {got['excess']:.3g} above every answer of the mean," if got and got.get("excess") is not None else ","
    return f"The noise comes first: the first kind beyond the twin of the question what is {noise}{excess} which no term file holds."


def brief(table: dict, problem, result, k: int, folder=None) -> str:
    """The page of one model for an agent, brief_K.md, his messages 315 and 316 of the repository it is made in: an
    attempt of an agent read 9.1 million tokens, most of them its own conversation read again at each step, and model_K.md
    read whole was a sixth of them. Under 400 words: the fit, the twin, the parameters, the answers beyond the twin, the
    three proposals for the model with their term files, the three for the experiment, and the next step by the rule of
    the loop, with how each value stood at the models before, read from the model_J.json of `folder`. model_K.md stays
    the page for people, and model_K.json holds everything."""
    dof = max(int(table.get("degrees_of_freedom") or 1), 1)
    lines = [f"# Model {k}, for an agent", "",
             f"Fit: {chi_bar(table['chi_square'], dof, *correlation_span(problem, result))} Bayesian criterion {table['criterion']:.1f}, lower is better.",
             ("Twin: admitted on the synthetic sets of this report. best.md reads the best model again on other sets: quote a value only where best.md says that it may be quoted." if table.get("admitted_on_the_twin")
              else "Twin: NOT admitted on the synthetic sets of this report: quote no value of this model.")]          # "admitted" read as quotable, while best.md refused the same model on its own sets (the fresh red team of 0.1.0rc20, mislead 5)
    if problem.shared:
        lines.append("Parameters: " + "; ".join(f"{n} {' ± '.join(_pm(float(v), float(s)))}" for n, v, s in zip(problem.shared, result.shared, result.sd)) + ".")
        if history := moved(folder, k, problem.shared):
            lines += ["How each value stood at the models before, where a fall is a share that a later term took:"] + history
    where = table.get("where") or {}
    beyond = [a for a in where.get("answers", []) if a.get("beyond_the_twin")]
    lines.append("Answers beyond the twin: " + ("; ".join(f"{a['question']}: {a['answer']}, excess {a['excess']:.2g}" for a in beyond[:4]) if beyond
                                         else "none, so the records ask for no more term") + ".")
    tie = (where.get("tied_with_the_first") or {}).get("where") or []
    if tie:          # two places that the records do not tell apart, named as one answer (his messages 388 and 389)
        lines.append(f"The first place, {(where.get('first_of_each_question') or {}).get('where')}, is tied with {' and '.join(tie)}: the records do not tell which leads, so a term of either place is as asked.")
    lead = noise_first(table)
    if lead:
        kind = noise_leads(where.get("answers", []))
        lines.append(lead + " It is read in the noise of each sample, in a copy of the folder made without loop/ and agent/. " + NOISE_KINDS.get(kind, f"model_{k}.md says how under its list of every answer."))          # the page sent to model_K.md, which the card says is for people, for the key of the copy (the blind red team of 0.1.0rc20, mislead 4)
    props = table.get("proposals") or {}
    if props.get("for_the_model"):
        lines += ["", "Proposals for the model, each measured on the twin:"]
        for i, r in enumerate(props["for_the_model"][:3]):
            said = f"{i + 1}. {r['term']}: `{r.get('as_a_term_file') or r.get('no_term_of_the_package') or 'no term file'}`, {'beyond' if r.get('beyond_the_twin') else 'within'} the twin"
            if r.get("gain_in_chi_square") is not None:
                said += f", fall of chi square {r['gain_in_chi_square']:.1f}"
            if r.get("refused_before"):
                said += f", refused before: {r['refused_before']}"
            lines.append(said)
    if props.get("for_the_experiment"):
        lines += ["", "Proposals for the experiment, the error of the known values with the change over that under the records: "
                  + "; ".join(f"{r['lever']} {r['error_over_that_under_the_records']:.2f}" for r in props["for_the_experiment"][:3]) + "."]
    nxt = f"terms/{k + 1:02d}.py"
    lines += ["", ("Next: the noise first, in a copy of the folder, then " + nxt + " there from the proposals of its report." if lead and beyond else
                   "Next: write " + nxt + " with proposal 1, or the term of physics that it stands for, its reason in `why`, then run the loop again." if beyond else
                   "Next: stop. A term that the physics asks and the records cannot show is held at its theory in " + nxt + " with hold = True."),
              f"model_{k}.md is the page for people, model_{k}.json holds everything."]
    return "\n".join(lines) + "\n"


def _pm(value: float, bar: float) -> tuple:
    """A value and its bar as the page prints them: the bar at two significant digits, the value at the same decimals."""
    if not (bar > 0.0 and math.isfinite(bar)):
        return f"{value:.6g}", f"{bar:.2g}"
    d = max(0, 1 - int(math.floor(math.log10(bar))))
    return f"{value:.{d}f}", f"{bar:.{d}f}"


def _two(x: float) -> str:
    """A bar at two significant digits in plain decimals, 0.000085 and not 8.5e-05 among plain numbers (round 6 of the cold
    tests)."""
    if not (x > 0.0 and math.isfinite(x)):
        return f"{x:g}"
    d = max(0, 1 - int(math.floor(math.log10(x))))
    return f"{x:.{d}f}"


TERM_WORDS = (("second component", "- **A second component**: a copy of the curve displaced along the abscissa, of a size of its own, a term that a term file may write. A second shape of its own width or of its own law, as a second decay of another lifetime, is written in `model.py` or as a Mechanism, as sections 10 and 11 of the tutorial show."),
              ("size", "- **The size** of a term: the one parameter that it adds, which `python -m twinfit terms` defines for each. For a term that adds a part of the curve it is the area of that part over that of the curve without it, so that its share of the whole is a / (1 + a), and never a distance."),
              ("of the level", "- **The level** of a curve is its baseline."))          # the words of the terms that a researcher writes, defined on a page that uses them: the report names none (the owner's message 340), and a kept term may


def used_words(lines: list) -> list:
    """The page with the words of TERM_WORDS that it uses defined in its list of words, after the others: a kept second
    component read as a second lifetime, and its size read as the distance it was written at (rounds 12 and 13 of the
    cold tests), while a page that holds no term defines none of them."""
    at = next((k for k, line in enumerate(lines) if line == "## The words of this page"), None)
    if at is None:
        return lines
    end = next((k for k in range(at + 2, len(lines)) if not lines[k].startswith("- **") and not lines[k].startswith("  ")), len(lines))
    body = "\n".join(lines[end:]).lower()
    return lines[:end] + [entry for word, entry in TERM_WORDS if word in body] + lines[end:]


def _words(sets: int, beyond: float) -> list:
    """The words of the page, defined before it uses them (row M167): a researcher who opens the report cold reads what
    a twin, a bar, an excess, a level, beyond the twin and a scatter are before the first number that uses them."""
    return ["## The words of this page", "",
            "- **The twin**: synthetic records drawn from this fitted model with the noise that the description gives, each fitted again as the records were. "
            "What the records show beyond what the twin shows is what the model lacks.",
            "- **The bar of the fit**: the standard deviation of a parameter that the fit gives from the curvature of its chi square, at first order. The first order bound is that bar.",
            "- **A scatter**: the standard deviation of a parameter over the synthetic sets of the twin, which the rule of the twin calls the spread of the estimate. A bias of 3 scatters is a value three of its own spreads from the truth of the twin.",
            "- **The excess** of an answer: its statistic on the records less the median of the same statistic on the synthetic sets, over the span from that median to their 95th percentile. "
            "An excess of 1 is the 95th percentile of the twin, and an excess under 0 an answer that the records show less than the median of the twin does.",
            f"- **The chance** of an answer: the share of the synthetic sets on which its statistic is as large as on the records. With {sets} sets the smallest chance is {1.0 / (sets + 1):.3f}, "
            "and the answers that share it are ordered by their excess.",
            "- **The level** of a curve, in a polynomial of the level, is its baseline.",
            "- **A search**: the fit of one synthetic set, which ends when its steps no longer lower its chi square. A **Mechanism** is a term written in Python, as section 11 of the tutorial shows.",
            "- **Seen by**: the family of statistics that sees the misfit most. It says how the misfit is read, and not which term to add.",
            "- **The known parameters**: the shared parameters of the model, those of every curve. **The rule of the twin**: the bands within which the fit of this model must give back, "
            "on its synthetic sets, the values that made them, its searches ending on nearly all sets, its bias under a share of its spread and its spread near its bar, for its bars to be read.",
            f"- **Beyond the twin**: an excess larger than {beyond:.2f}, the largest excess that a synthetic set of a right model reaches on 95 of 100 of them. "
            "Some answer is beyond it on 5 of 100 sets where the model lacks nothing. The line is read on the twin of this model, so it differs a little from one model to the next.", ""]


def _span_of(problem, result):
    """The samples over which neighbouring whitened residuals correlate, where they correlate at more than 0.3, else None,
    and None where the records give no residuals to read, as a page made from a table alone."""
    try:
        rho, span = correlation_span(problem, result)
    except (AttributeError, TypeError, IndexError):
        return None
    return span if rho is not None and abs(rho) > 0.3 else None


def per_degree(chi2: float, dof: int, noise_right: bool = True, span=None) -> str:
    """The chi square of a fit per degree of freedom, with what it says: near 1 within three of its standard deviations,
    the root of 2 over the degrees of freedom, where the model and the stated noise are right. The page gave the chi
    square and its degrees of freedom alone, and 18881 on 3592 read as approval beside "the twin admits the fit" (round 3
    of the cold tests)."""
    ratio = chi2 / max(dof, 1)
    sd = math.sqrt(2.0 / max(dof, 1))
    each = f"{sd:.2g} each" + (f", or {sd * math.sqrt(span):.2g} where the noise correlates over {span} samples as the residuals do" if span and span > 1 else "")          # the bar of correlated neighbours, as the visual log and the page for an agent give it (the blind red teams of 0.1.0rc19, mislead 2 and mislead 1)
    if abs(ratio - 1.0) <= 3.0 * sd:          # a sum near 1 where the bands of the noise below are not: a noise stated too large at some heights and too small at others balances in the sum (round 6 of the cold tests)
        return (f"{ratio:.3g} a degree of freedom, 1 within three of its standard deviations, {each}: " + ("the model and the stated noise agree with the records in the sum, which sums away what the table of where to search reads, so a chi square near 1 does not say that the model lacks nothing" if noise_right else
                "the sum agrees with the records, but the check of the noise below does not: the noise is stated wrong at some heights, and the errors balance in the sum"))
    if ratio > 100.0:          # 9130 a degree of freedom from a start seven widths from the line read as a term missing (round 9 of the cold tests)
        return f"{ratio:.3g} a degree of freedom, above 1 by more than three of its standard deviations, {each}: the model lacks something large, or its fit began far from the curves, so look at the fitted curves beside the records before a term is written"          # starts first misled where a shape was missing whole, 475 a degree of freedom (round 10 of the cold tests, before its launch)
    if ratio > 1.0:
        return f"{ratio:.3g} a degree of freedom, above 1 by more than three of its standard deviations, {each}: the model lacks something, or the noise is stated too small"
    return f"{ratio:.3g} a degree of freedom, below 1 by more than three of its standard deviations, {each}: the noise is stated too large"


def _scatters(part: dict) -> str:
    """A bias in scatters, or the bias itself where the twin's scatter is zero: the page ended on "division by zero" for a
    term whose fits on its twin all gave one value (round 11 of the cold tests)."""
    if not part.get("scatter"):
        return f"{part['bias']:.3g}, its scatter on the twin being zero,"
    return f"{_one(part['bias'] / part['scatter'])} scatters, {part['bias']:.3g} in its own unit,"          # 2.3 scatters of a wrong model read as 2.3 bars, where the scatter was 70 bars and the fit moved by 0.108 (round 20 of the tests of the tutorial)


def _one(v: float) -> str:
    """A number at one decimal, with 0.0 where it rounds to minus zero: the reports printed -0.0 (round 7 of the cold
    tests)."""
    said = f"{v:.1f}"
    return "0.0" if said == "-0.0" else said


def _page(table, problem, result) -> str:
    """The report of one model for the researcher, in words and tables. It defines its words before it uses them, gives
    the first answer of each question before the list of every answer, says why a table holds less than it asks, and
    names the kind of each proposal for the experiment (row M167)."""
    rep = table["where"]
    lines = [f"# {table['model']}", ""] + _words(rep["sets"], rep["largest_excess_on_the_twin"])
    lines += ["## The fit", "",
              f"- **Chi square**: {table['chi_square']:.1f} on {table['degrees_of_freedom']} degrees of freedom, {per_degree(table['chi_square'], table['degrees_of_freedom'], table['noise'].get('right', True), _span_of(problem, result))}.",
              f"- **The twin {'admits' if table['admitted_on_the_twin'] else 'does not admit'} the fit**: the fit of this model "
              + ("gives back, on synthetic sets drawn from it, the values that made them, within the bands of its rule. That says the fit can be read, and not that the model is right."
                 if table['admitted_on_the_twin'] else "does not give back, on synthetic sets drawn from it, the values that made them within the bands of its rule: its bars are not to be read."),
              f"- **What the twin measured**, each bias in spreads of the estimate on its own sets: {table['measured_on_the_twin']}.",
              f"- **The rule of the twin**: {table['rule_of_the_twin']}." + (f" The verdict above used the wider bands, since the twin read {rep['sets']} sets." if table.get("bands_are_wider_than_stated") else ""),
              f"- **The noise**: {UNTESTED} {table['noise']['said']}", ""]          # labelled lines: one paragraph of 270 words read as one block (round 3 of the cold tests)
    if problem.shared:
        lines += ["| parameter | value | bar of the fit |", "|---|---|---|"] + [f"| {n} | {' | '.join(_pm(v, s))} |" for n, v, s in zip(problem.shared, result.shared, result.sd)]
    else:
        lines += ["The model shares no parameter among the records: each record has its own alone, so no value is known across them, "
                  "and no answer along a setting can be read. The values of each curve are in best.md."]
    unit = (table.get("proposals") or {}).get("width_of_the_curves")
    if getattr(problem, "centre", None) is not None and unit and table["chi_square"] > 10.0 * max(table["degrees_of_freedom"], 1):
        places = [problem.centre_of(k, result.shared, result.own[k]) for k in range(len(problem.records))]
        if max(places) - min(places) > unit:          # each spectrum of a doublet started on a line that its noise chose, and no page said so (round 12 of the cold tests)
            lines += [f"The places of the curves differ by {max(places) - min(places):.3g}, more than a width of the curves, while the fit leaves {table['chi_square'] / max(table['degrees_of_freedom'], 1):.3g} of chi square "
                      "a degree of freedom: the starts may sit on different lines of a multiplet, or on different features. Start every curve on the same feature, in `start_own`, before a term is read."]
    lines += ["", "## Where to search", ""] + ([noise_first(table) + " The table orders the answers by their excess, and a term of the curve is read once the noise holds it, as the end of this section says.", ""] if noise_first(table) else []) + [f"Read on {rep['sets']} synthetic sets. The first answer beyond the twin of each question, beyond meaning an excess above {rep['largest_excess_on_the_twin']:.2f} here:", "",
              "| question | what it asks | first answer beyond the twin | its excess |", "|---|---|---|---|"]
    kinds = sorted([a for a in rep["answers"] if a["question"] == "shape" and a["beyond_the_twin"]], key=lambda a: -a["excess"])          # the shapes beyond the twin: a slow trend and two cycles over the span can tie
    tied = [a for a in kinds if a["excess"] >= kinds[0]["excess"] / 1.5] if kinds else []
    six = table["proposals"]
    for q, asks in ASKS.items():
        if q == "matched" and q not in rep["first_of_each_question"]:
            continue          # no list of terms: the question is not asked
        mine = [a for a in rep["answers"] if a["question"] == q]
        first = rep["first_of_each_question"].get(q) or next((a["answer"] for a in sorted(mine, key=lambda a: -a["excess"]) if a["beyond_the_twin"]), None)          # the row and the list of every answer say one verdict (round 4 of the cold tests)
        if first is not None:
            got = next((a for a in mine if a["answer"] == first), None)
            said = first
            if q == "where" and rep.get("where_read_along"):          # the place of a part that changes along a setting, where no place stands beyond the twin alone (row M485)
                said, got = f"{first}, read along a setting, {rep['where_read_along']}: no place stands beyond the twin alone", next(a for a in rep["answers"] if a["question"] == "changes along" and a["answer"] == rep["where_read_along"])
            places = (rep.get("tied_with_the_first") or {}).get("where") or []
            if q == "where" and places and not rep.get("where_read_along"):          # a hidden term's misfit read 21.5 in the flanks and 21.1 in the core, and the row named the core alone (his messages 388 and 389)
                said = f"{first}, tied with {' and '.join(places)}: their excesses differ by less than {TIE_PLACES:g} standard deviations of their difference as the records scatter them, so the records do not tell which of these places leads"
            if q == "shape" and len(tied) > 1:          # the row named one shape and the proposals put another first (round 3 of the cold tests)
                said = (f"{first}, tied with {', '.join(a['answer'] for a in tied if a['answer'] != first)} within a factor of 1.5 of its excess: the proposals below put the kind of the largest excess first, then the others by the gain of their term")
            lines.append(f"| {q} | {asks} | {said} | {_one(got['excess'])} |")
        elif mine:
            top = max(mine, key=lambda a: a["excess"])
            lines.append(f"| {q} | {asks} | none beyond the twin. The largest: {top['answer']} | {_one(top['excess'])} |")
        else:
            lines.append(f"| {q} | {asks} | none: no answer was read for this question | none |")
    if len(tied) > 1:          # an odd misfit, as a slope, lights several kinds at once: their order then says little, and the first kept was an impostor (the cold tests of row M169)
        lines += ["", f"{len(tied)} kinds stand beyond the twin within a factor of 1.5 of the largest excess: {', '.join(a['answer'] for a in tied)}. One misfit can light several kinds at once, "
                  "and their order then says little: try each, with `python -m twinfit scan` where a number of the term is chosen, and read what each does to the known parameters."]
    if not any(a["beyond_the_twin"] for a in rep["answers"] if a["question"] in ACTIONABLE):
        unread = [q for q in ACTIONABLE if not any(a["question"] == q for a in rep["answers"])]
        lines += ["", "No answer of where, shape, what or changes along, nor of a term of your list, stands beyond the twin: on these records the model leaves nothing that its own synthetic sets do not leave"
                  + (f", though no answer of {' or '.join(unread)} was read" if unread else "")
                  + (", and the families of statistics that stand past their line name no term" if any(a["beyond_the_twin"] for a in rep["answers"]) else "")
                  + ". A term much smaller than the noise can still be missing, and more sets, or records of more samples, read smaller ones."]
    lines += [f"- {said}" for said in rep["cannot_be_read"]]
    lines += ["", "Every answer, the first twelve by their chance and then their excess:", "",
              "| question | answer | chance | excess | beyond the twin |", "|---|---|---|---|---|"]
    lines += [f"| {a['question']} | {a['answer']} | {a['level']:.3f} | {_one(a['excess'])} | {'yes' if a['beyond_the_twin'] else 'no'} |" for a in rep["answers"][:12]]
    noisy = [a for a in [noise_leads(rep["answers"])] if a]
    lines += tied_noise(rep["answers"], problem, result)
    if noisy:          # an answer that no term of the mean fills was met with proposals of the mean alone (round 3 of the cold tests)
        lines += ["", f"The first kind beyond the twin is {noisy[0]}, which no term of the mean holds: it is read in the noise of each sample. {NOISE_KINDS[noisy[0]]} "
                  "Do it in a copy of this folder, so that the curves as they were taken stay as they are, and run the loop on the copy: a term file cannot hold it."]
        if noisy[0] == "a noise correlated from sample to sample":          # the span that "one sample in every few" asks for, and that a term kept on one choice of samples may be refused on another (round 9)
            rho, span = correlation_span(problem, result)
            if rho is not None and misfit_first(result):
                lines += ["", misfit_first(result)]
            elif rho is not None:
                lines += ["", f"The correlation of neighbouring whitened residuals is {rho:.2f} at one sample, the median over the curves, and "
                          + (thin_said(problem, span, "so keep one sample in every {span}, which `python -m twinfit thin DIR {span}` does in three copies") if span else "does not fall under 0.05 within 200 samples. ")
                          + "Run the loop on several such copies, each from another first sample: a term that one keeps and the others refuse is not to be believed."]
        every = glitch_samples(problem, result, most=None) if noisy[0] == "glitches" else []          # every sample beyond 5, of which 20 are shown: the page gave the number shown for the number beyond, 20 for 46 (the blind red team of 0.1.0rc18, mislead 10)
        far = every[:20]
        if far:          # the samples named, which a blind agent found by numpy of its own (round 7)
            lines += ["", f"The samples whose whitened residual, the residual over the noise of its sample, is beyond 5 in size: {len(every)}, the largest first, of which {len(far)} are shown.", "",
                      "| curve | abscissa | whitened residual |", "|---|---|---|"]
            lines += [f"| {name} | {x + 0.0:.6g} | {v:.1f} |" for name, x, v in far]          # -0 for an abscissa of zero (round 9)
            shared = shared_abscissas(problem, every)
            lines += ["", f"Of these {len(every)}, {shared} stand at an abscissa where a sample of another record is far from the curve too, within one step of the samples: "
                      + ("no glitches, but a place where the law of the noise is too small or the model misses." if shared else "each alone, as a glitch is.")]
    lines += ["", "## Three proposals for the model", ""]
    if noisy:          # a cosine of period 4.68 that fitted 16 spikes was ranked first and kept (round 7 of the cold tests)
        lines += [f"{noisy[0][0].upper() + noisy[0][1:]} leads, and the proposals below are read on residuals that it leads: a term kept on them may fit it and not a part of the curve, "
                  "and is not to be believed before the noise holds it. So here the noise comes first, and a term only where an answer of the curve leads. The proposals for the experiment are measured with such a term too.", ""]          # one order where the page gave two (the blind red team of 0.1.0rc13)
    lines += ["Each is a direction read on the residuals, for an answer of where, of shape or of what changes along a setting: "
              "those beyond the twin first, then by the gain in chi square that the shape that the residuals hold there brings when it is fitted beside the model. "
              "The package names no term: which term of physical origin has that place, shape and setting is the researcher's to say, in a term file."
              + (" On its twin, synthetic curves drawn with the term at its fitted size are fitted without it and with it: the bias without it is what the term would cost the known "
                 "parameters if it were real, and says nothing of whether it is." if problem.shared else ""), ""]
    if six.get("width_of_the_curves") is not None:
        lines += [f"In the names of the proposals a width is the full width of the fitted curves at half their height above their median value, here {six['width_of_the_curves']:.3g} in the unit of the abscissa: "
                  "where the wings of a curve fill its window, it is narrower than the width at half its peak, and it is not a parameter of the model. "
                  "Above the centre is toward larger abscissa.", ""]
    for m in six["for_the_model"]:
        lines.append(f"- {m['term']}, for {m['stands_for']}, {'beyond' if m['beyond_the_twin'] else 'within'} the twin at an excess of {m['excess_in_the_report']:.1f}"
                     + ((f": a gain in chi square of {m['gain_in_chi_square']:.1f}" if m["gain_in_chi_square"] >= 0 else f": at its rough number the shape raises the chi square by {-m['gain_in_chi_square']:.1f}, so its number is far off: scan it")
                        if "gain_in_chi_square" in m else f": {m.get('why_not_measured', '')}"))          # a gain of -1089128.4 for the right term, unread (round 11 of the cold tests)
        if m.get("as_a_term_file") and " changes along " in m["term"]:          # the square of the quencher proposed for a lifetime of 1/(1 + K Q), 36 bars off (round 14 of the cold tests)
            lines.append("  - a polynomial of the setting about its mean, which no law of the physics need be: where the physics gives one, as 1/(1 + K Q) or the root of 1 + P/P0, write that law in `model.py`, as section 9 of the tutorial shows")
        ripple = m.get("ripple")
        if ripple:          # a period of 2.82, the width of the curves, for a ripple of 0.35, and no word of its frame (round 14 of the cold tests)
            lines.append(f"  - its period, {ripple['period']:.4g} in the unit of the abscissa, is read from the residuals, a periodogram summed over the curves: scan about it to fix it")
            about = ripple.get("coherence_about_the_centres")
            if about is not None and ripple["coherence_about_the_abscissa"] > about + 0.2:
                lines.append(f"  - the curves agree at that period at equal abscissa, a coherence of {ripple['coherence_about_the_abscissa']:.2f} of 1, more than about their centres, {about:.2f}: "
                             "a ripple fixed to the abscissa, as of the instrument, which a term placed about the centre of each curve does not hold. Write it as a Mechanism with x in the place of x - centre, as section 11 of the tutorial shows")
        if m.get("as_a_term_file"):
            lines.append(f"  - as a term file: `{m['as_a_term_file']}`")
            pk = m.get("the_term_of_the_package")
            if pk and pk["stands"]:
                lines.append(f"  - that term, fitted beside the model: its size, {pk['parameter']}, is {' ± '.join(_pm(pk['size'], pk['bar']))}, a gain in chi square of {pk['gain_in_chi_square']:.1f}"
                             + (": a negative size takes the shape away from the curve, as wings that fall faster than the model's, and may stand for another term" if pk["size"] < 0 and pk["parameter"].startswith(("pedestal", "second component", "tail")) else ""))          # a negative pedestal ranked with no word on its sign (round 6)
            elif pk:
                lines.append(f"  - that term, fitted beside the model, does not stand: {pk['why_not']}")
            if m.get("no_term_of_the_package"):          # a line with a word on what it does not hold, as a background of one size where each record may carry its own (round 6)
                lines.append(f"  - {m['no_term_of_the_package'][0].upper() + m['no_term_of_the_package'][1:]}")
        elif m.get("no_term_of_the_package") and not m["no_term_of_the_package"].startswith("the package names no term"):          # the line above the list says it once for every direction
            lines.append(f"  - as a term file: none, {m['no_term_of_the_package']}")
        if m.get("refused_before"):          # level(1) proposed again after its refusal, with no word (round 8 of the cold tests)
            lines.append(f"  - refused before, in {m['refused_before']}")
        if m.get("follows"):          # a step in proportion to the height was proposed as a level along the angle that the height followed (round 8)
            f = m["follows"]
            lines.append(f"  - its size in each record is {' ± '.join(_pm(f['share'], f['bar']))} times {f['own']}, an own value of the model, one share in every record within its bars, a chi square of "
                         f"{f['chi_square']:.1f} on {f['degrees_of_freedom']}: a term whose size is a share of {f['own']} holds it, a Mechanism of your own that reads own[{f['index']}], "
                         f"as section 11 of the tutorial shows. A size that follows {f['own']} follows every setting that {f['own']} follows: read this before a term along a setting")
        for name, v in (m.get("known_parameters") or {}).items():
            lines.append(f"  - if the term were real at that size, {name} would be biased on its twin by {_scatters(v['without_the_term'])} without it and by "
                         f"{_scatters(v['with_the_term'])} with it: what the term would cost, and not a sign that it is real. A scatter of a model that lacks a term can be many bars of its fit"
                         + (f". But the working model reads {name} on the records {abs(v['the_records_in_scatters']):.0f} scatters from what it reads on this twin, so the term in this shape and at this size "
                            "does not make the records, and what it would cost there may be far larger: fit it, and read how far the values move" if abs(v.get("the_records_in_scatters") or 0.0) > 3.0 else "")
                         + (". This shape is read on the residuals, which hold none of what the fit took into its own parameters, so this cost is a lower bound: "
                            "a term of the physics in this place may cost far more" if m.get("cost_is_a_lower_bound") else ""))          # read as an argument to add a term that moved a shared width by 15 bars (round 7)
    lines += [f"- Fewer than three: {f}" for f in six["fewer_than_three"] if f.split(" ")[2:4] == ["the", "model:"]]
    if rep["first_of_each_question"].get("where") == "the asymmetry" and any(a["question"] == "where" and a["answer"] == "the asymmetry" and a["beyond_the_twin"] for a in rep["answers"]):
        lines += ["", "The asymmetry leads. Where no proposal above holds it once fitted, a shape that no term of the package writes may be missing, asymmetric about the place of the line: "
                  "an exponential fall of the line on one side, its height times x over one plus x squared, x the distance from its place in its widths, "
                  "or the derivative of the line. A Mechanism of your own writes it, as section 11 of the tutorial shows, and the physics of the records says which."]          # no cause named, which the package never names (the readers of the red team of his message 434, tutorial R3)          # a Fano line led by the asymmetry, and no proposal named a shape for it (round 10 of the cold tests)
    term = six.get("the_term_in_the_twin_of_the_levers")
    lines += ["", "## Three proposals for the experiment", "",
              "Each is of one of three kinds: a better setup, a better characterization of the setup, or a better data collection. Each is measured on the twin"
              + (f", with the first proposal for the model that stands beyond the twin in it, {term}" if term
                 else f", with the model as it is, since the largest answer beyond the twin of the question what is {six['the_noise_leads']}, a noise of the samples that no term holds: read them again once the noise holds it"          # it said the report leads with it, under a table that put the core first (the fresh red team of 0.1.0rc18, mislead 3)
                 + (", and more samples over the same span buy less than the twin says, since it draws them apart where the records hold them correlated" if six["the_noise_leads"] == "a noise correlated from sample to sample" else "")
                 if six.get("the_noise_leads") else ", with the model as it is, since no proposal for the model stands beyond the twin")
              + ("" if problem.shared else ". The model shares no parameter, so what a change buys is read on the size of that term") + ". The error of a value is the root of its squared bias and its squared scatter "
              "there, over the same under the records as they are: under 1, the change helps.", ""]
    for e in six["for_the_experiment"]:
        lines.append(f"- {e['his_kind'][0].upper() + e['his_kind'][1:]}, {e['lever']}: {e['words']}. In the twin, {e['how_it_enters_the_twin']}. "
                     f"The error of {'the known parameters' if problem.shared else 'the size of that term'} is {e['error_over_that_under_the_records']:.2f} of that under the records"
                     + (f", known within about {e['known_within']:.2f} on {table['where']['sets']} sets, since the change draws noise of its own" if e.get("known_within") else ""))
        if len(e.get("of_each") or {}) > 1:          # each parameter with its own spread: the headline, their geometric mean, hid a width that grew worse (round 6 of the cold tests)
            lines.append("  - of each parameter, the number above being their geometric mean, the n-th root of the product of the n of them: " + ", ".join(f"{n} {v:.2f}" + (f" ± {v * e['known_within'] / e['error_over_that_under_the_records']:.2f}" if e.get("known_within") else "")
                                                                                         for n, v in e["of_each"].items()) + (". Over 1, the change makes that parameter worse" if any(v > 1.0 for v in e["of_each"].values()) else ""))
    shown = {e["lever"] for e in six["for_the_experiment"]}
    lines += [f"- Not measured, {e['his_kind']}, {e['lever']}: {e['why_not_measured']}" for e in six.get("every_lever", []) if e["lever"] not in shown and e.get("why_not_measured")]
    lines += [f"- Also measured, not among the three, {e['his_kind']}, {e['lever']}: {e['words']}. The error is {e['error_over_that_under_the_records']:.2f} of that under the records"
              + (f", known within about {e['known_within']:.2f}" if e.get("known_within") else "") for e in six.get("every_lever", [])
              if e["lever"] not in shown and e.get("error_over_that_under_the_records") is not None]          # a lever of the researcher that ranks out is said, not dropped (round 2)
    lines += [f"- Fewer than three: {f}" for f in six["fewer_than_three"] if f.split(" ")[2:4] == ["the", "experiment:"]]
    return "\n".join(used_words(lines)) + "\n"


def noise_settled(d, problem, result, source):
    """The first model fitted again with the noise that the description asks: read from its residuals, `[noise] from =
    "residuals"`, or from the researcher's law of model.py, `from = "model"`, and as it stands otherwise. The loop and
    `scan` and `matrix` read the first model through it alike: `model_of` left the noise as the records gave it, so that
    on a folder whose noise is read from the residuals `scan` gave a fall of 21.7 to a term that the loop refused, and
    636 to a pedestal that the loop kept at a fall of 222 (the blind red team of 0.1.0rc19, mislead 5)."""
    if d.noise["from"] == "residuals":          # the law of the noise read from what the first model leaves, his message 316: the agent of a blind run built it by hand
        return noise_from_residuals(problem, result)
    if d.noise["from"] == "model":          # the researcher's law, `variance(mean, record, own)` of model.py, settled at the fit of the first model, his message 316
        if not callable(getattr(source, "variance", None)):
            raise ValueError("the description gives the noise from the model, and model.py defines no variance(mean, record, own): the variance of each sample from the fitted mean of its curve, the curve and its own parameters")
        from .noise import with_law
        for _ in range(2):          # from the fit that the search found, as twinfit.noise.settle does from a fit of its own
            problem = with_law(problem, result, source.variance)
            result = fit(problem, result.shared, result.own)
        if not result.converged:
            raise ValueError("the fit of the first model with the noise of its law did not end: look at the law, or give the noise as a column")
    return problem, result


def model_of(folder, skip=None) -> tuple:
    """The model that the term files of a folder give, fitted as `scan` and `matrix` read it: the first model with the noise
    that the description asks, `noise_settled`, then each term file in turn, kept by the Bayesian criterion alone, a term that gives its theory and says `hold` held at it,
    and the file `skip` left out. No twin is run. Returns the problem, the fit, the name of the centre, the files used,
    the dictionary HELD of model.py with the ansatz of the description in it, and model.py's module."""
    if not (Path(folder) / "description.toml").is_file():          # a folder that is not a project ended in a traceback
        raise ValueError(f"{folder} holds no description.toml: give the folder of a project, as `python -m twinfit init DIR` writes it")
    d, problem, centre, source = initial(Path(folder), with_module=True)
    held = getattr(source, "HELD", None)
    for name, a in (d.addenda.get("ansatz") or {}).items():          # the values of the ansatz in HELD, as the loop sets them
        if held is not None and name in held:
            held[name] = float(a["value"])
    result = searched(problem)
    problem, result = noise_settled(d, problem, result, source)          # the noise of the loop: without it `scan` read another criterion than the loop's (the blind red team of 0.1.0rc19, mislead 5)
    used = []
    for path in sorted((Path(folder) / "terms").glob("[0-9][0-9].py")):          # the model that the term files give, by the criterion alone
        if skip is not None and path.resolve() == Path(skip).resolve():
            continue
        term = getattr(_module(path, f"term_{path.stem}", held, source), "term", None)
        if term is None:
            continue
        if term.theory is not None and getattr(_module(path, f"term_{path.stem}", held, source), "hold", False):
            problem, result = held_term(problem, result, term)[:2]
            used.append(f"{path.name}, held")
            continue
        try:
            got = develop_one(problem, result, term)
        except ValueError:          # a term that builds on one that the criterion left out, which stopped the scan whole (round 3 of the cold tests)
            continue
        if got["ended"] and not got["at_its_bound"] and _criterion(problem, result) - _criterion(got["_problem"], got["_fit"]) > TIED:
            problem, result = got["_problem"], got["_fit"]
            used.append(path.name)
    return problem, result, centre, used, held, source


def matrix(folder, top: int = 5) -> dict:
    """What every family of statistics says on the records of a folder, for a researcher at the command line (the owner's
    messages 262 to 267: the clues of the generalised moments, the windows, the transformations of the value and the
    chi functions as the pillars of the decision). On the model that the term files give, fitted as `scan` fits it: the
    reading of each family on the residuals, a chi square with its degrees of freedom and its level under a right model
    at first order, the shape of the noise, `where.noise_shape`, and the noise matrix, `families.noise_matrix`, which
    statistic sees which law of the noise by simulation. No list of terms is read (the owner's messages 337 and 340). No
    twin is run: the loop's report gives the twin's levels, and a reading here is the evidence beside them."""
    from .families import catalogue, family, noise_matrix
    from .where import noise_shape
    problem, result, centre, used, _, _ = model_of(folder)
    lin = Linear(problem, result)
    readings = {}
    for name in catalogue(centre):
        try:
            readings[name] = family(lin, name, centre).reading()
        except (ValueError, np.linalg.LinAlgError) as why:
            readings[name] = {"refused": str(why)}
    r = np.asarray(lin.r, dtype=float)
    ordered = [r[lin.rows_of(k)][np.argsort(np.asarray(rec.x, dtype=float), kind="stable")] for k, rec in enumerate(problem.records)]
    return {"folder": str(folder), "beside": used, "readings": readings, "noise": noise_shape(ordered), "noises": noise_matrix(lin, centre)}


def matrix_page(out: dict) -> str:
    """The page that `python -m twinfit matrix DIR` prints from `matrix`: what each family reads, the shape of the noise,
    and which statistic sees which law of the noise."""
    beside = f"with {', '.join(out['beside'])}" if out["beside"] else "its first model"
    lines = [f"The families of statistics on the model of {out['folder']}, {beside}, each term file kept by the Bayesian criterion alone. No twin is run: a level here is that of a chi square "
             "under a right model at first order, which the families of many statistics do not reach, and the report of the loop gives the levels of the twin.", "",
             "## What each family reads", "",
             "Each family turns the residuals into a few statistics. Its reading is their chi square in the covariance that they have under a right model, "
             "with as many degrees of freedom as the family has directions that the model does not hold. The smallest levels first.", "",
             "| family | chi square | degrees of freedom | level |", "|---|---|---|---|"]
    read = {n: r for n, r in out["readings"].items() if "refused" not in r}
    for name, r in sorted(read.items(), key=lambda kv: kv[1]["level"]):
        lines.append(f"| {name} | {r['chi_square']:.1f} | {r['degrees_of_freedom']} | {r['level']:.2g} |")
    lines += [f"- Not read: {name}, {r['refused']}" for name, r in out["readings"].items() if "refused" in r]
    lines += ["", "## The shape of the noise", "",
              "The whitened residuals of every curve as they stand, their skewness over its standard error under a normal law and their excess kurtosis over its own, "
              "each a standard normal deviate at first order where the noise is normal. A misfit of the mean leaves tails of its own, and the report reads them once the shape "
              "that leads is taken out: " + ", ".join(f"{name} {v:.2f}" for name, v in out["noise"].items()) + "."]
    nm = out.get("noises")
    if nm:
        from .families import NOISE_WORDS
        kinds = list(nm["chance"])
        lines += ["", "## Which statistic sees which law of the noise", "",
                  f"Each law of the noise is drawn {nm['draws']} times on these records, whitened and of unit variance on every sample, so that the noise that the description states "
                  "stays right on average and only its law is wrong, and taken through the fit as the residuals take it. The line of each statistic is where it stands on 95 of 100 "
                  f"draws of a normal noise of the same design, so each cell is a chance by simulation at the level {nm['level']:g}, and owes nothing to first order. "
                  "After the families come the five readings of the noise that the report gives, the records together, and two opponents that practitioners read, the runs test and Durbin and Watson, which the report does not read. The last row reads every other statistic at one line, "
                  "as the report reads its answers: the chance that the largest excess of them all passes where a normal noise stands on 95 of 100 draws. The laws:", ""]
        lines += [f"- {kind}: {NOISE_WORDS[kind]}" for kind in kinds]
        lines += ["", "| statistic | " + " | ".join(kinds) + " |", "|---|" + "---|" * len(kinds)]
        for name in nm["statistics"]:
            lines.append(f"| {name} | " + " | ".join(f"{nm['chance'][k][name]:.2f}" for k in kinds) + " |")
    return "\n".join(lines) + "\n"


def scan(folder, expression: str, values) -> dict:
    """A term whose argument is a choice and not a fitted parameter, the distance of `second`, the stretch of `pedestal`
    or of `tail`, the period of `periodic`, tried at each of `values`: `expression` is the term with x in the place of
    that argument, as "second(x)" or "tail(x, 'upper')", read among the terms of mechanisms.py. Each is fitted beside
    the model that the term files of the folder give, each kept there by the Bayesian criterion alone, and the row says
    whether its fit ended inside its bounds, the fall of the criterion, and the size with its bar. No twin is run: the
    scan says where to look, and the loop, with the term written at the chosen value, decides whether it is kept (row
    M168: the report names a second component of the chromatograms and its first order guess of the distance, 0.35
    widths, is not the 1.4 seconds that drew them)."""
    from . import mechanisms
    folder = Path(folder)
    own = (folder / expression if (folder / expression).is_file() or not Path(expression).is_file() else Path(expression)) if expression.strip().endswith(".py") else None          # a path from DIR or from where the command runs: my_nmr/terms/04.py was refused as no term_at (round 6)          # a term file of the researcher's own, which defines term_at(x): a blind agent could not scan the scale of its own tail (round 4 of the cold tests)
    if own is None and "x" not in expression:
        raise ValueError(f"the term {expression!r} holds no x: write x where the scanned value goes, as second(x)" + (". along has no number to choose: write it in a term file, and the loop fits its size" if expression.strip().startswith("along(") else ""))
    if own is not None and not own.is_file():
        raise ValueError(f"{expression}: no such file, from {folder} or from where the command runs")
    problem, result, centre, used, held, source = model_of(folder, skip=own)
    made = None          # the term of a file of one's own at the number x, where one is scanned
    if own is not None:
        made = getattr(_module(own, "term_scanned", held, source), "term_at", None)
        if made is None:
            raise ValueError(f"{expression}: a term file of your own is scanned when it defines term_at(x), the term at the number x, beside term = term_at(...) for the loop")
    rows = []

    def line_of(v: float) -> str:
        return f"term = term_at({v:g})" if own is not None else "term = mechanisms." + re.sub(r"\bx\b", f"{v:g}", expression)
    for v in values:
        try:
            term = _term_at(made, own, expression, mechanisms, v)          # a whole number as one: level(2.0) was named "level power 2.0" and fitted beside the kept level(2) (round 6 of the cold tests)
        except NameError as why:
            raise ValueError(f"{expression}: {why}. A term is one of `python -m twinfit terms`, and a word in it is written in quotes, as tail(x, 'upper')") from None
        except SyntaxError:
            raise ValueError(f"{expression} is not a term with x in it, as second(x) or tail(x, 'upper')") from None
        if term.parameter in problem.shared:          # a term that the model holds already was scanned beside itself, and gave a size of -0 with a bar of 1992 (round 3)
            if all(_held_at(made, own, expression, mechanisms, w, problem) for w in values):
                raise ValueError(f"{term.parameter} is a parameter of the model that the term files give already: scan a term that it does not hold, or scan it in a folder without its term file")
            rows.append({"x": float(v), "term": term.name, "line": line_of(float(v)), "how": "held already by a term file", "ended": False, "fall_of_the_criterion": None, "size": None, "bar": None, "shared": {}})
            continue          # the other values are scanned: the whole call was refused for the one that a file held (round 8 of the cold tests)
        got = develop_one(problem, result, term)
        line = line_of(float(v))          # the line of a term file that holds this term
        how = "ended" if got["ended"] and not got["at_its_bound"] else ("at a bound of its size" if got["ended"] else "did not end")
        rows.append({"x": float(v), "term": term.name, "line": line, "how": how, "ended": bool(got["ended"] and not got["at_its_bound"]),
                     "fall_of_the_criterion": float(_criterion(problem, result) - _criterion(got["_problem"], got["_fit"])) if "_fit" in got else None,
                     "size": got.get("parameter", {}).get("value"), "bar": got.get("parameter", {}).get("bar"),
                     "shared": dict(zip(problem.shared, (float(s) for s in got["_fit"].shared))) if "_fit" in got else {}})          # how the chosen x moves the other values: a width quoted at one distance moved 1.2 bars across the range (round 5 of the cold tests)
    top, within, edge, parabola, beside = profile_top(rows, values, line_of)
    return {"expression": expression, "beside": used, "rows": rows, "best": top, "within_one": within, "at_the_edge": edge, "from_a_parabola": parabola, "next_to_no_end": beside}


def profile_top(rows: list, values, line_of) -> tuple:
    """The top of the profile of a scan: the row of the largest fall of the criterion among the fits that ended, the
    values within 1 of it, whether they reach an end of the values scanned, whether a parabola through the three falls
    about the largest gives the top, and the values beside the largest whose fit did not end. The parabola is taken only
    through three neighbouring values of the scan that all ended, and its top only between the outer two: through the
    ended values on either side of a gap where the fits did not end it gave a top at 42 where the largest fall stood at
    200 (round 17 of the cold tests)."""
    ended = [r for r in rows if r["ended"] and r["fall_of_the_criterion"] is not None]
    top = max(ended, key=lambda r: r["fall_of_the_criterion"]) if ended else None
    near = sorted(r["x"] for r in ended if top is not None and r["fall_of_the_criterion"] >= top["fall_of_the_criterion"] - 1.0)          # the profile at one standard deviation: the criterion of models of one size differs as their chi square
    edge = bool(near) and (near[0] == min(float(v) for v in values) or near[-1] == max(float(v) for v in values))
    within, parabola = ((near[0], near[-1]) if near else None), False
    every = sorted(rows, key=lambda r: r["x"])
    at = next((i for i, r in enumerate(every) if r is top), None)
    beside = [every[i]["x"] for i in ((at - 1, at + 1) if at is not None else ()) if 0 <= i < len(every) and not every[i]["ended"]]          # with no top, nothing beside it: `at - 1` on None stopped a scan of a pedestal bounded at zero whose every fit stood at its bound (review 108, moderate 1)
    if near and near[0] == near[-1] and at is not None and 0 < at < len(every) - 1 and not beside:          # the range under the step of the values: a parabola through the three falls about the largest gives it, where a finer scan asked for a finer one again (round 5 of the cold tests)
        xs, fs = [every[i]["x"] for i in (at - 1, at, at + 1)], [every[i]["fall_of_the_criterion"] for i in (at - 1, at, at + 1)]
        a, b, _ = np.polyfit(xs, fs, 2)
        if a < 0.0 and xs[0] <= -b / (2.0 * a) <= xs[2]:
            centre, half = -b / (2.0 * a), 1.0 / math.sqrt(-a)
            within, parabola = (centre - half, centre + half), True
            digits = max(0, 1 - int(math.floor(math.log10(half)))) if half > 0 else 3
            top = dict(top, x=round(centre, digits), line=line_of(round(centre, digits)), scanned=top["x"])          # the value to write is the top of the parabola: the scan recommended a value outside the range that it printed (round 6 of the cold tests)
    return top, within, edge, parabola, beside



def _term_at(made, own, expression: str, mechanisms, v):
    """The term of a scan at the value v: the term file of one's own at v, or the expression with v in the place of x,
    a whole number passed as one, so that a copy of a kept `level(2)` is refused (round 6 of the cold tests)."""
    return made(float(v)) if own is not None else eval(expression, {"__builtins__": {}}, {**{k: getattr(mechanisms, k) for k in dir(mechanisms) if not k.startswith("_")}, "x": int(v) if float(v).is_integer() else float(v)})


def _held_at(made, own, expression: str, mechanisms, v, problem) -> bool:
    """Whether the term of a scan at the value v adds a parameter that the model holds already."""
    return _term_at(made, own, expression, mechanisms, v).parameter in problem.shared


def _named(path: Path, err: Exception) -> ValueError:
    """A refusal met while a term file is fitted or read on its twin, with the file named: the loop stopped on "record '':
    a value or an abscissa is not a number" and named no file, where a term of one's own overflowed (round 10 of the cold
    tests)."""
    if isinstance(err, IndexError):
        return ValueError(f"{path.name}: {err}. Its wrap reads own[k] or shared[k] past those that the model holds: own[k] in the order of OWN, shared[:n] those before the term and shared[n] its own")
    if "not a number" not in str(err) and "overflow" not in str(err):          # a misspelt setting was followed by a sentence on values that are not numbers (round 13 of the cold tests)
        return ValueError(f"{path.name}: {err}")
    return ValueError(f"{path.name}: {err}. The term gives a value that is not a number somewhere, as an exponential that overflows: "
                      "compute its wrap on the curves alone, with the starts of the model")


TALK = None          # a function that the command sets to say each step of the loop as it starts: 8.5 minutes went by with no word, and its cost grew with each term (round 13 of the cold tests)


_STATUS = None          # the file of the progress of the run in course, and its start


def _say(text: str) -> None:
    """Says a step of the loop where the command asked for it, and nothing elsewhere."""
    if TALK is not None:
        TALK(text)
    if _STATUS is not None:
        try:
            _STATUS[0].write_text(json.dumps({"at": time.strftime("%Y-%m-%d %H:%M:%S"), "seconds": round(time.time() - _STATUS[1]), "said": text}) + "\n")
        except OSError:
            pass


def by_a_hair(adm: dict) -> bool:
    """Whether the twin refused a fit only by statistics that each miss their band by less than their own standard error
    on the sets read: a bias over the spread within one over the root of the sets, a spread within one over the root of
    twice the sets, a coverage within the root of 0.95 times 0.05 over the sets. Every search must have ended as the rule
    asks, and nothing else may be named. A right term was refused at 100 sets alone, its coverage 0.88 against a band
    from 0.885, and kept at 30, 50, 150, 200 and 500 (round 13 of the cold tests)."""
    b, n = adm.get("bands"), int(round(adm.get("ended", 0.0) * adm.get("sets", 0)))
    if adm.get("admitted") or not b or n < 30 or adm.get("ended", 0.0) < ENDED or adm.get("why_not"):
        return False
    misses = []
    for c in adm.get("parameters", {}).values():
        if c is None:
            return False
        if c.get("bias_over_spread") is not None and not abs(c["bias_over_spread"]) < b["bias_over_spread_under"]:
            misses.append((abs(c["bias_over_spread"]) - b["bias_over_spread_under"]) * math.sqrt(n))
        lo, hi = b["spread_over_bound"]
        if c.get("spread_over_bound") is not None and not lo <= c["spread_over_bound"] <= hi:
            misses.append(max(lo - c["spread_over_bound"], c["spread_over_bound"] - hi) * math.sqrt(2.0 * n))
        lo, hi = b["coverage"]
        if c.get("coverage") is not None and not lo <= c["coverage"] <= hi:
            misses.append(max(lo - c["coverage"], c["coverage"] - hi) / math.sqrt(0.95 * 0.05 / n))
    return bool(misses) and max(misses) <= 1.0


def too_few(twin: dict) -> str:
    """Why a term whose twin ended on too few sets is refused: the fits with it ended on fewer synthetic sets than a
    reading needs, which is how a term that the records cannot hold apart from those before it shows itself."""
    if not twin["read"]:          # "0 of 30", then "run with --sets 100", then "0 of 100", for a right term whose sets the twin read in one search of 80 steps (round 20 of the tests of the tutorial)
        return (f"on its twin the fits with it ended on none of its {twin['sets']} synthetic sets, each searched in {STEPS} steps, so more sets will not read it: "
                "the records do not hold it apart from the terms before it, or its search does not end, as where its curve jumps as a parameter moves. "
                "A term of one's own is written smooth in its parameters. For a term of the package, write the law that it stands for in `model.py` in its place, "
                "as the gallery's clock writes the shift of its line along the power")
    return (f"on its twin the fits with it ended on {twin['read']} of {twin['sets']} synthetic sets, and a term is read on {MIN_SETS} synthetic sets at least: "
            "the records do not hold it apart from the terms before it, or its search does not end, as where its curve jumps as a parameter moves: write it smooth in its parameters")          # "read on 30 at least" read beside "ended on 20 of 30", and a step at the place kept searches from ending (round 11)


def term_file(path: Path, m) -> None:
    """Refuses a term file that does not hold a term and its reason, naming which of the three it lacks: one message
    for the three left a researcher to guess (the cold tests of row M169)."""
    from .mechanisms import Mechanism
    if not hasattr(m, "term"):
        raise ValueError(f"{path.name} defines no `term`: write term = mechanisms.second(1.0), or another of `python -m twinfit terms`, or a Mechanism of your own")
    if not isinstance(m.term, Mechanism):
        raise ValueError(f"{path.name}: `term` is a {type(m.term).__name__}, and a term is a twinfit Mechanism, as `python -m twinfit terms` lists them")
    import inspect
    try:
        takes = len(inspect.signature(m.term.wrap).parameters)
    except (TypeError, ValueError):
        takes = 3
    if takes != 3:          # a wrap of two arguments ended the loop on a trace of the code (round 7 of the cold tests)
        raise ValueError(f"{path.name}: the wrap of its Mechanism takes {takes} arguments, and it takes three, wrap(model, problem, n): section 11 of the tutorial gives one")
    if not hasattr(m, "why"):          # "holds 0 words" where the file wrote no why (round 7 of the cold tests)
        raise ValueError(f"{path.name} defines no `why`: write why = \"...\", the reason of the researcher for the term, in six words or more")
    words = len(str(m.why).split())
    if words < 6:
        raise ValueError(f"{path.name}: `why`, the reason of the researcher for the term, holds {words} words, and it takes six or more")


def run(folder, n: int = 10, sets: int = 30, report: bool = True, max_seconds: float | None = None, workers: int | None = None, dictionary: bool = False) -> dict:
    """The loop on a folder, with the output of the last run kept until this one ends: a run that stops on a fault puts
    the last `loop/` back, so that its best.md is not lost (the cold tests of row M169). One run at a time in a folder,
    by a lock, and the last whole run kept through a kill too: two runs in one folder destroyed each other's output, and
    a run killed before it ended left the last good one where the next run removed it (the red team of the package, his
    message 181). `max_seconds` ends the run between two of its steps once that many seconds have passed. `workers`: the
    twins of each model's proposals side by side on that many cores, every core where it is None, one with
    TWINFIT_WORKERS=1, and the same numbers whatever it is (`twinfit.pool.side_by_side`). See `_run`."""
    folder = Path(folder)
    if not folder.is_dir():
        near = sorted(p.name for p in folder.parent.glob(folder.name.split("_thin_")[0] + "_thin_*") if p.is_dir()) if "_thin_" in folder.name and folder.parent.is_dir() else []
        if near:          # a copy of thin named from memory, `_2` for `_1`, was sent to init (round 18 of the cold tests)
            raise ValueError(f"{folder} is no folder: the copies that thin wrote beside it are {', '.join(near)}")
        raise ValueError(f"{folder} is no folder: `python -m twinfit init {folder}` writes one to fill, and `python -m twinfit example nmr {folder}` a worked one")
    lock = (folder / ".twinfit.lock").open("a+")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.seek(0)
        held = lock.read().strip()
        lock.close()
        raise ValueError(f"another run of the loop works in {folder}" + (f", process {held}" if held else "") + ": wait for it to end, or run in a copy of the folder, since two runs in one folder destroy each other's output") from None
    lock.seek(0)
    lock.truncate()
    lock.write(f"{os.getpid()}\n")
    lock.flush()
    began, cpu, ended = time.time(), _cpu_seconds(), "a fault"
    try:
        with asking(workers):
            out = _kept_through(folder, n, sets, report, max_seconds, dictionary)
            ended = "ended"
            after_the_run(folder, out)
            return out
    finally:
        try:
            cost(folder, began, _cpu_seconds() - cpu, how_many(workers, 10 ** 6), ended)
        finally:
            lock.close()


COSTS = "loop.costs.tsv"          # beside loop/: a line a run of the loop, what it cost


def _cpu_seconds() -> float:
    """The CPU of this process and of the workers it has waited for, in seconds."""
    import resource
    return sum(getattr(resource.getrusage(who), f) for who in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN) for f in ("ru_utime", "ru_stime"))


def cost(folder: Path, began: float, cpu: float, workers: int, ended: str) -> Path:
    """A line of `loop.costs.tsv` in the folder for a run of the loop: when it began, its seconds, the CPU of the run and
    of its workers in seconds, the workers it could use, the version of the package and how it ended. The cost of an
    analysis is then read from the folder, whoever started the run: an agent that starts the loop in the background
    leaves its CPU out of its own count (the red team of his message 319, S3)."""
    from .version import __version__
    path = Path(folder) / COSTS
    new = not path.exists()
    with path.open("a") as f:
        if new:
            f.write("began\tseconds\tcpu_seconds\tworkers\tversion\tended\n")
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(began))}\t{time.time() - began:.1f}\t{cpu:.1f}\t{workers}\t{__version__}\t{ended}\n")
    return path


def _kept_through(folder: Path, n: int, sets: int, report: bool, max_seconds: float | None, dictionary: bool = False) -> dict:
    """run, in a folder that this process holds alone."""
    out_dir, kept = folder / "loop", folder / "loop.before"
    if kept.exists():
        if (kept / "best.md").exists() and not (out_dir / "best.md").exists():          # a run killed before it ended: the last whole run comes back, and what the killed one left is set aside
            if out_dir.exists():
                if (folder / "loop.killed").exists():
                    shutil.rmtree(folder / "loop.killed")
                out_dir.rename(folder / "loop.killed")
            kept.rename(out_dir)
        else:
            shutil.rmtree(kept)
    if out_dir.exists():
        out_dir.rename(kept)
    try:
        try:
            out = _run(folder, n, sets, report, max_seconds, dictionary)
        except ValueError as err:          # any refusal of a run on an abscissa far from zero names that cause too: the twin's fits fail on it as the first fit does (round 5)
            try:
                hint = "" if "abscissa_offset" in str(err) else far_from_zero(initial(folder)[1])
            except Exception:          # a description or a model that cannot be read says so in the refusal itself
                hint = ""
            if hint:
                raise ValueError(str(err) + hint) from None
            raise
    except BaseException as err:
        if out_dir.exists():
            shutil.rmtree(out_dir)
        if kept.exists():
            kept.rename(out_dir)
            stale(out_dir, err)
        raise
    if kept.exists():
        shutil.rmtree(kept)
    return out


def seen_by_beyond(answers: list) -> list:
    """The families of statistics that stand beyond the twin in a report, its answers of seen by, which name no term: a
    noise of the root of the counts of the records stood beyond it in three families while the command said that the
    records ask for no more (round 15 of the cold tests)."""
    return [a["answer"] for a in answers if a["question"] == "seen by" and a["beyond_the_twin"]]


def asks_no_term(answers: list) -> bool:
    """Whether a report asks for no term: no answer of where, shape, what or changes along stands beyond the twin. A
    family of statistics past its line, an answer of seen by, names no term, and the loop said write the next file where
    only those stood past it (round 4 of the cold tests)."""
    return not any(a["beyond_the_twin"] for a in answers if a["question"] in ACTIONABLE)


def stale(out_dir: Path, err: BaseException) -> None:
    """Marks the output that a refused run puts back as that of the run before it: best.md was read as the answer to a
    run that had been refused (round 3 of the cold tests). REFUSED.md says why, and the first line of best.md says that
    the page is older than the last run. The next run that ends writes both again."""
    said = f"The run of {time.strftime('%Y-%m-%d %H:%M:%S')} was " + ("interrupted" if isinstance(err, KeyboardInterrupt) else f"refused: {err or type(err).__name__}")          # an interrupt read "was refused: ." (the red team of the package, his message 181)
    (out_dir / "REFUSED.md").write_text(f"# The last run was refused\n\n{said}\n\nThe other files of this folder are those of the run before it, which ended.\n")
    page = out_dir / "best.md"
    if page.exists() and not page.read_text().startswith("**Not the last run**"):
        page.write_text(f"**Not the last run**: {said}. This page is of the run before it: REFUSED.md says why.\n\n" + page.read_text())
    numbers = out_dir / "best.json"
    if numbers.exists():          # a program that read best.json got the numbers of the run before with no word (the red team of the package, his message 181)
        try:
            held = json.loads(numbers.read_text())
            held["stale"] = said
            numbers.write_text(json.dumps(held, indent=1, default=float) + "\n")
        except (ValueError, OSError):
            pass


def builds_on(term, problem, refused: list) -> str | None:
    """Why a term is not read, where it changes a parameter that the model does not hold: the term file that adds it was
    refused, or it is named wrong. The second file of a recipe of two stopped the whole loop when the first had been
    refused (round 3 of the cold tests)."""
    named = term.name.split(ALONG.rstrip())[0] if term.family == "grows with a covariate" and ALONG.rstrip() in term.name else None          # along(p, key) names itself "p changes along key"
    if named is None or named in problem.shared:
        return None
    by = [r for r in refused if r.get("parameter") == named]
    if by:
        return f"it changes {named!r}, the parameter of {by[-1]['term']} in {by[-1]['file']}, which was refused: a term that builds on a refused term is not fitted"
    return None


def in_place(iterations: list, refused: list, made: dict, kept: list) -> list:
    """The models in which a term refused by the Bayesian criterion takes the place of the term kept just before it,
    each fitted on the model before that kept term and set aside: best() holds it in the spread over the models where
    its criterion is within TIED of the model of the kept term, and names it where it is lower by more. Two causes of one
    misfit fit alike, and the one written first was kept while the other, fitted on top of it, added nothing: a dead time
    and a width that grows with the intensity tied within 0.1 of the criterion, eight bars of the width apart, and the
    page said that the records ask for no more (round 16 of the cold tests). Each refused entry that was fitted so gets
    `in_place`, with the file it stands for and its criterion less that of the model of the kept term."""
    by = {it["model"]: it for it in iterations}
    fitted = {k["term"]: k["file"] for k in kept if "held_at" not in k}
    out = []
    for r in refused:
        on = by.get(r.get("on"))
        if not r["why"].startswith(("the Bayesian criterion", "the twin does not admit")) or on is None or not on.get("after") or on.get("added") not in fitted or r["file"] not in made:          # a right term refused by the twin behind a wrong one kept first, told to move a bound that a term of the package does not take (round 17 of the cold tests)
            continue
        base = by[on["after"]]
        try:
            developed = develop_one(base["problem"], base["result"], made[r["file"]])
        except (ValueError, IndexError):
            continue
        if not developed["ended"] or developed["at_its_bound"]:
            continue
        less = _criterion(on["problem"], on["result"]) - _criterion(developed["_problem"], developed["_fit"])
        r["in_place"] = {"of": fitted[on["added"]], "term": on["added"], "criterion_less": float(less)}
        out.append({"model": f"{on['model']} with {r['file']} in the place of {fitted[on['added']]}", "problem": developed["_problem"], "result": developed["_fit"],
                    "added": r["term"], "after": base["model"], "aside": True, "tried": [], "refused": [], "held": base.get("held")})
    return out


def dropped(problem, result, larger, bigger) -> list:
    """The own parameters that moved a curve under the model before a term and move none under the model with it: a
    term that replaces the line and returns a curve of its own drops what it does not add again, as a baseline, and the
    loop kept such a tail 16 bars from its truth, its baseline printed with no bar (round 16 of the cold tests)."""
    out = []
    for j, name in enumerate(problem.own):
        if j >= len(larger.own) or larger.own[j] != name:
            continue
        moves = []
        for p, r in ((problem, result), (larger, bigger)):
            own = np.array(r.own[0], dtype=float)
            step = 1e-3 * (abs(own[j]) + 1.0)
            up = own.copy()
            up[j] += step
            base = p.predict(0, r.shared, own)
            moves.append(bool(np.any(np.abs(p.predict(0, r.shared, up) - base) > 1e-9 * (np.max(np.abs(base)) + 1e-300))))
        if moves[0] and not moves[1]:
            out.append(name)
    return out


def _run(folder, n: int = 10, sets: int = 30, report: bool = True, max_seconds: float | None = None, dictionary: bool = False) -> dict:
    """The loop on a folder, up to `n` term files. Returns what was kept, what was refused, the file to write next, and
    the output of best.py. Fewer than 30 synthetic sets are refused: the notebook asks that many of a twin."""
    if sets < 30:
        raise ValueError(f"the twin of each model needs 30 synthetic sets at least, so that its coverage and its bias can be read: {sets} were asked")
    folder = Path(folder)
    if not (folder / "description.toml").is_file():          # refused before anything is written into the folder: a run on --help made a folder of that name (the cold tests of row M169)
        raise ValueError(f"{folder} holds no description.toml: `python -m twinfit init {folder}` writes one to fill, and `python -m twinfit example nmr {folder}` a worked one")
    out_dir = folder / "loop"          # run() moved the last one aside: each run remakes the analysis from the files of the folder
    out_dir.mkdir(parents=True)
    began = time.time()
    global _STATUS
    _STATUS = (out_dir / "STATUS.json", began)          # what the run is doing, after each of its steps, for whoever waits on it (the red team of the package, his message 181)

    from .reuse import Reuse, code_of
    reuse = Reuse(folder)          # what the runs before computed from the same files, read back (the red team of his message 327)
    code = code_of(folder)

    def budget(step: str) -> None:          # a run stops at its budget only once it has computed a step anew, so the next run reads it back and goes further
        if max_seconds is not None and time.time() - began > max_seconds and reuse.computed:
            raise ValueError(f"the run stopped at its budget of {max_seconds:g} seconds, {step}, after {time.time() - began:.0f}: a larger `--max-seconds`, fewer term files or fewer sets end it")
    book = Notebook(out_dir / "notebook")
    d, problem, centre, source = initial(folder, with_module=True)
    terms = list_of_terms(d, source)          # the researcher's, empty by default (his message 342)
    ansatz = d.addenda.get("ansatz") or {}
    if not d.addenda.get("model"):          # row M205: the description alone, as his message 140 gives the input
        _say("the description names no first model: the loop starts from the generic one of the package, one peak or dip of one width on its own height, place and baseline, "
             "twinfit/agent/model_template.py. Where the curves are no such peak, write model.py from the physics of the description, as the paragraph of the guide before its first iteration says, and name it under [addenda]")
    for k in range(1, n + 1):          # every term file read before the first fit: a reason of two words refused after three models were fitted (round 16 of the cold tests)
        early = folder / "terms" / f"{k:02d}.py"
        if not early.is_file():
            break
        m_early = _module(early, f"early_{k:02d}", getattr(source, "HELD", None), source)
        term_file(early, m_early)
        try:
            m_early.term.wrap(problem.model, problem, len(problem.shared))
        except ValueError as err:          # a setting that the records do not carry, said after "is fitted" (round 16 of the cold tests)
            if "carry no setting" in str(err):
                raise ValueError(f"{early.name}: {err}") from None
        except Exception:          # a term that reads what a kept term adds is wrapped on that model only, later
            pass
    values = getattr(source, "HELD", None)          # the values that the model file holds, which an ansatz sets
    for name, a in ansatz.items():
        if values is not None and name in values:
            values[name] = float(a["value"])
    held_of(ansatz, problem, values)          # an ansatz that would put no bar is refused before the first fit
    try:
        result = searched(problem)
    except KeyError as err:          # the model reads a setting by another name: the refusal gave the key alone (round 5 of the cold tests)
        raise ValueError(f"the model reads {err} from the settings of a curve, `design`, and the settings that the description gives are {', '.join(problem.records[0].design) or 'none'}") from None
    except ValueError as err:
        if "moves no curve" in str(err):          # the message pointed at the start of the shared parameter alone, and the fault was in start_own (round 3 of the cold tests)
            one = [key for key in problem.records[0].design if len({str(r.design.get(key)) for r in problem.records}) == 1 and len(problem.records) > 1] if problem.records else []
            raise ValueError(f"{err}, and at start_own: a curve started away from its peak or its dip, or at no height, is moved by no shared parameter. "
                             "A start far from the value that the records ask moves no curve either, as a width started at 1.0 for a line 0.016 wide: give START_SHARED the scale of the curves"
                             + (f". A parameter that multiplies a setting minus a value that every curve takes moves no curve: {', '.join(one)} takes one value over all the curves" if one else "")
                             + (". A parameter that multiplies a setting minus its one value moves no curve where there is one curve alone" if len(problem.records) == 1 else "")) from None          # a shift per bar on one curve, whose pressure less 1 was 0 (round 13 of the cold tests)
        raise
    if not result.converged:
        raise ValueError(f"the fit of the initial model did not end in {AGAIN + 1} searches of {SEARCH} steps at most"
                         + (f", the last of which stopped after {result.iterations} step{'s' * (result.iterations != 1)}, where no step lowered its chi square any more" if result.iterations < SEARCH else "")
                         + f", its chi square {result.chi2 / max(result.dof, 1):.4g} for each degree of freedom: look at its starts and its bounds before the loop" + far_from_zero(problem) + why_no_end(problem))
    problem, result = noise_settled(d, problem, result, source)
    rung = "model 0"
    _say(f"model 0: the first model is fitted, and its twin and its report are read on {sets} synthetic sets")
    _rung(book, rung, "the initial model of the description, from its addenda", problem, result, d, sets, out_dir, 0, centre, report, held_of(ansatz, problem, values), reuse=reuse, code=code, terms=terms, dictionary=dictionary)
    budget(f"after {rung}")
    iterations = [{"model": rung, "problem": problem, "result": result, "added": problem.name, "tried": [], "refused": [], "held": held_of(ansatz, problem, values), "formula": equation_of(source),
                   "letters": {str(k): str(v) for k, v in (getattr(source, "LETTERS", None) or {}).items()}}]          # each letter of EQUATION with its meaning, which the visual log writes under it (the blind red team of 0.1.0rc12)
    write_equations(out_dir, iterations)
    kept, refused, nxt, holding, made = [], [], None, {}, {}
    for i in range(1, n + 1):
        path = folder / "terms" / f"{i:02d}.py"
        if not path.is_file():
            nxt = str(path.relative_to(folder))
            break
        budget(f"before {path.name}")
        code = code_of(folder, *(folder / "terms" / f"{j:02d}.py" for j in range(1, i + 1)))          # the term files read so far make every model from here
        m = _module(path, f"term_{i:02d}", values, source)
        term_file(path, m)
        term = m.term
        made[path.name] = term
        iterations[-1]["tried"].append(term.name)
        built = builds_on(term, problem, refused)
        if built:
            refused.append({"term": term.name, "file": path.name, "why": built, "parameter": term.parameter})
            iterations[-1]["refused"].append({"term": term.name, "why": built})
            continue
        if getattr(m, "hold", False):          # a term that stands on theory is held at it, not fitted: rule 10, his messages 133 and 135
            told = carried(problem, result, term)          # the fits at the centre of the theory and one bar either side
            larger, bigger, entry = held_term(problem, result, term)
            if not (told["ended"] and bigger.converged):
                refused.append({"term": term.name, "file": path.name, "why": "the fit with it held at the centre that the theory gives, or one bar from it, did not end"})
                iterations[-1]["refused"].append({"term": term.name, "why": refused[-1]["why"]})
                continue
            stamp = f"{term.name} held on the twin of {rung}"
            book.stamp(stamp, f"Synthetic sets drawn with {term.name} at sizes from its theory, its centre and its bar, and fitted with it held at the centre, give back the other values.", rung)
            twin = reuse("the twin of a held term", lambda: held_on_the_twin(problem, result, told, sets=sets, max_iter=STEPS, table_of=None, workers=how_many(None, sets)),
                         code, problem, result, told, sets, report)
            budget(f"after the twin of {path.name}")
            twin.setdefault("share_of_sets_in_which_something_is_seen", None)          # not read in a count, as for a fitted term
            if twin["read"] < MIN_SETS:
                book.note(stamp, too_few(twin))
                refused.append({"term": term.name, "file": path.name, "why": too_few(twin)})
                iterations[-1]["refused"].append({"term": term.name, "why": too_few(twin)})
                continue
            book.step("candidate on the twin", rung, dict(twin, term_file=path.name), stamp=stamp)
            book.grade(stamp, holds=twin["read"] >= 2, what_came=f"{twin['read']} of {twin['sets']} sets ended")
            holding.update(entry)
            book.choose(rung, term.name, "held", by=str(m.why))
            problem, result = larger, bigger
            before, rung = rung, f"model {len(iterations)}"
            now = held_of(ansatz, problem, values) | held_terms(holding, problem)
            _rung(book, rung, str(m.why), problem, result, d, sets, out_dir, len(iterations), centre, report, now, [k["term"] for k in kept] + [term.name], refused, reuse=reuse, code=code, terms=terms, dictionary=dictionary)
            budget(f"after {rung}")
            iterations.append({"model": rung, "problem": problem, "result": result, "added": term.name, "after": before, "tried": [], "refused": [], "held": now, "formula": term.formula, "parameter": term.parameter, "file": path.name})
            write_equations(out_dir, iterations)
            kept.append({"term": term.name, "file": path.name, "model": rung, "held_at": entry[term.parameter][0], "bar_of_the_theory": entry[term.parameter][1], "why": str(m.why)})
            continue
        _say(f"{path.name}: {term.name} is fitted, and read on {sets} synthetic sets: each term costs more than the one before it, since the model grows")
        try:
            developed = develop_one(problem, result, term)
        except (ValueError, IndexError) as err:          # REFUSED: record '': a value or an abscissa is not a number, with no file named (round 10), and a trace of the code for own[k] out of range (round 11)
            raise _named(path, err) from None
        why_not, hair, fall = None, False, None
        if not developed["ended"] or developed["at_its_bound"]:
            got = developed.get("parameter") or {}          # the value, its start and the bounds: a width refused so at its bound gave none of them, and a researcher read the source to find that its setting was taken raw (the blind red team of 0.1.0rc18, block 1)
            said = f", {got['name']} at {got['value']:.4g} from its start at {got['started_at']:.4g}, its bounds {term.lower:.4g} and {term.upper:.4g}" if "value" in got else ""
            why_not = developed.get("why_not") or f"its fit ended at the bound of its own parameter{said}, where the term is absent or as large as it may be"
            if "moves no curve" in why_not and kept:          # a term that replaces the line, after a kept term, was sent to its starts and its bounds (round 10)
                why_not += ". A term that replaces the line and calls no `model` drops the terms kept before it, whose parameters then move no curve: such a term goes first, as 01.py, in a copy without the files before it"
            if "does not give a number at the start" in why_not:          # a tail of 1/t started at 0, refused with no word of its start (round 14 of the cold tests)
                why_not += ". The start of the term, or one of the model, gives a value that is not a number: start a parameter that the curve divides by away from zero, and compute the wrap on the curves alone, with the starts of the model, before the loop"
            if "no step lowered" in why_not:          # a shift of the line along the power, on a model whose place is own to each curve, was refused for its steps and not for its cause (round 9 of the cold tests)
                why_not += ". A term that moves the curves as an own parameter of each curve already moves them gives its search no step: read whether one holds the same freedom, and write it in `model.py` in its place"
        else:
            stamp = f"{term.name} on the twin of {rung}"
            book.stamp(stamp, f"Synthetic sets drawn from the model with {term.name} at its fitted size give back that size within its bar, and the working model is moved by it.", rung)
            try:
                twin = reuse("the twin of a term", lambda: on_the_twin(problem, result, developed, sets=sets, max_iter=STEPS, table_of=None, workers=how_many(None, sets)), code, problem, result, developed, sets, report)          # the twin of a term on every core: one of 30 sets ran 58 minutes on one core in a blind run of 0.1.0rc12
            except ValueError as err:          # the value that is not a number met where the twin drew its sets (round 10)
                raise _named(path, err) from None
            twin.setdefault("share_of_sets_in_which_something_is_seen", None)          # not read in a count, where no table of candidates is made on its sets
            if twin["read"] < MIN_SETS:          # the notebook refuses a candidate read on fewer, and the loop stopped there with a trace of the code (the cold tests of row M169)
                book.note(stamp, too_few(twin))
                refused.append({"term": term.name, "file": path.name, "why": too_few(twin)})
                iterations[-1]["refused"].append({"term": term.name, "why": too_few(twin)})
                continue
            book.step("candidate on the twin", rung, dict(twin, term_file=path.name), stamp=stamp)
            back = twin.get("what_the_whole_procedure_returns", {}).get(term.parameter, {})
            book.grade(stamp, holds=bool(back) and abs(back.get("bias_in_scatters", 99.0)) < 1.0, what_came=f"its size comes back {back.get('bias_in_scatters', float('nan')):.2f} scatters from what was put in")
            larger, bigger = developed["_problem"], developed["_fit"]
            gain = _criterion(problem, result) - _criterion(larger, bigger)
            adm = reuse("the twin of a larger model", lambda: admission(larger, bigger, sets=sets, start="fit"), code, larger, bigger, sets)
            budget(f"after the twin of {path.name}")
            lost = dropped(problem, result, larger, bigger)
            if lost:
                why_not = (f"with it, the curve reads no {', '.join(lost)}, which the model before it read: a term that returns a curve of its own and does not call `model` "
                           "drops what it does not add again, so add it to the curve that its wrap returns, as `+ own[k]` for a baseline, section 11 of the tutorial")
            elif gain <= TIED:
                moved = f"rises by {-gain:.2f}" if gain < 0 else f"falls by only {gain:.2f}"          # the opening, "the Bayesian criterion", is what toy/chosen.py reads
                why_not = f"the Bayesian criterion {moved} with it, where it must fall by more than {TIED:g}: these records do not tell the model with it from the model without it"
            elif not adm["admitted"]:
                why_not = "the twin does not admit the fit of the model with it: " + "; ".join(adm.get("why_not", []) + adm["failed"] or [adm["rule"]])
                hair, fall = by_a_hair(adm), float(gain)          # what the records say of the model that the twin alone refused, which best.md names (the blind red team of 0.1.0rc18, mislead 13)
                got = developed.get("parameter") or {}
                why_not += f". With it the Bayesian criterion falls by {gain:.1f}" + (f", and {got['name']} is fitted at {' ± '.join(_pm(got['value'], got['bar']))}" if math.isfinite(got.get("bar", math.nan)) and "name" in got else "")          # a refusal with no size and no fall, so that nothing said whether 60 minutes at 100 sets were worth it (the blind red team of 0.1.0rc20)
                if hair:          # a right term refused at 100 sets alone by 0.005 of coverage, and the command said write the next term (round 13 of the cold tests)
                    why_not += (f". Each statistic that fails misses its band by less than its own standard error on these {sets} sets, so the refusal is a draw of the twin as much as a verdict on the term: "
                                "run the loop again on another number of sets, 200 or 500, before the term is set aside")          # the statistic that failed, with its value and its band
                table_before = out_dir / f"model_{len(iterations) - 1}.json"
                if table_before.is_file() and not json.loads(table_before.read_text())["admitted_on_the_twin"]:          # a slope refused at 100 sets alone, while the model without it missed the same bands (round 12 of the cold tests)
                    why_not += (f". The model without it is not admitted by its twin either, so this refusal does not choose between the two, while the Bayesian criterion falls by {gain:.1f} with it. "
                                f"A bias read on {sets} sets has a standard error of {1.0 / math.sqrt(sets):.2g} of a spread, and a band missed by less than that can be met at another number of sets: "
                                "read both models on 200 sets, and quote neither before its twin admits it")
                if sets < 100:          # a true step refused at 30 sets and kept at 50, 70 and 100, with no word (round 8 of the cold tests)
                    why_not += f". Read on {sets} sets, a refusal by the twin can turn on more: `--sets 100` reads it again before the term is set aside"
        if why_not is not None:
            refused.append({"term": term.name, "file": path.name, "why": why_not, "parameter": term.parameter, "by_a_hair": hair, "on": rung} | ({"criterion_falls": fall} if fall is not None else {})
                           | ({"size": developed["parameter"]["value"], "bar": developed["parameter"]["bar"]}
                              if developed.get("ended") and not developed.get("at_its_bound") and math.isfinite(developed.get("parameter", {}).get("bar", math.nan)) else {}))          # what the records say of its size (round 3)
            iterations[-1]["refused"].append({"term": term.name, "why": why_not, "by_a_hair": hair})          # into best.json: the gallery said "refused" of a term refused by a hair on 30 sets and kept on 100 (round 19)
            continue
        book.choose(rung, term.name, "fitted", by=str(m.why))
        report_before = out_dir / f"model_{len(iterations) - 1}.json"
        asked = (not asks_no_term(json.loads(report_before.read_text())["where"]["answers"])) if report and report_before.is_file() else None          # a cosine that the report put within its twin, kept with no word (round 12 of the cold tests)
        problem, result = developed["_problem"], developed["_fit"]
        before, rung = rung, f"model {len(iterations)}"
        now = held_of(ansatz, problem, values) | held_terms(holding, problem)
        _rung(book, rung, str(m.why), problem, result, d, sets, out_dir, len(iterations), centre, report, now, [k["term"] for k in kept] + [term.name], refused, reuse=reuse, code=code, terms=terms, dictionary=dictionary)
        budget(f"after {rung}")
        iterations.append({"model": rung, "problem": problem, "result": result, "added": term.name, "after": before, "tried": [], "refused": [], "held": now, "formula": term.formula, "parameter": problem.shared[-1], "file": path.name})
        write_equations(out_dir, iterations)
        kept.append({"term": term.name, "file": path.name, "model": rung, "parameter": term.parameter, "why": str(m.why), "asked": asked})
    try:
        held_models = {k["model"] for k in kept if "held_at" in k}          # a term held at its theory stays in the model: the best was once chosen among models that lacked it (round 5 of the cold tests)
        first = max((i for i, it in enumerate(iterations) if it["model"] in held_models), default=0)
        budget("before best.md")
        _say(f"best.md: the five parts of each value of the best model, on {max(sets, 30)} synthetic sets")
        swapped = [s for s in in_place(iterations, refused, made, kept) if s["after"] in {it["model"] for it in iterations[first:]}]
        output = reuse("the best model", lambda: best(iterations[first:] + swapped, sets=max(sets, 30), start="fit", workers=how_many(None, max(sets, 30))),
                       code, [(it["model"], it["problem"], it["result"], it.get("held"), it.get("added"), it.get("after"), it.get("aside")) for it in iterations[first:] + swapped], max(sets, 30))          # each held value carries its bar of theory: the fit made again at the value less and more one bar
        by, at, chain = {it["model"]: it for it in iterations}, iterations[first].get("after"), []
        while at:          # the terms of the models before the held one: "Its terms" named the held term alone (round 11 of the cold tests)
            chain.insert(0, by[at]["added"])
            at = by[at].get("after")
        output["terms"] = chain + output["terms"]
    finally:
        for name, a in ansatz.items():
            if values is not None and name in values:
                values[name] = float(a["value"])
    chosen = next(it for it in iterations if it["model"] == output["best"])
    chosen_problem, chosen_result = chosen["problem"], chosen["result"]
    plain = {k: v for k, v in output.items() if k != "iterations"} | {"iterations": output["iterations"]}
    noise = noise_said(chosen_problem, chosen_result)
    plain["noise"] = noise
    plain["own_parameters"] = [{"record": rec.name, "settings": dict(rec.design), **{name: {"value": v, "bar_of_the_fit": b} for name, (v, b) in row.items()}}
                               for rec, row in zip(chosen_problem.records, own_bars(chosen_problem, chosen_result))]          # each curve's own values, as the time of each eclipse (the cold tests of row M169)
    plain["own_trends"] = own_trends(plain["own_parameters"])
    lo, hi = min(float(np.min(r.x)) for r in chosen_problem.records), max(float(np.max(r.x)) for r in chosen_problem.records)
    plain["level_unit"] = (0.5 * (lo + hi), 0.5 * (hi - lo) or 1.0)          # as mechanisms.level takes it
    plain["per_unit"] = per_unit(chosen_problem, chosen_result)
    plain["abscissa_offset"] = d.offset
    plain["places"] = places(chosen_problem, chosen_result)
    plain["abscissa_unit"] = d.abscissa[1]
    plain["setup"] = d.setup          # the instruments and settings that the description gives, row M216
    plain["generic"] = not d.addenda.get("model")
    plain["order"] = d.order
    plain["chi_square_per_degree"] = float(chosen_result.chi2) / max(int(chosen_result.dof), 1)
    held = sorted({name for c in output["parameters"].values() for name in c["bar_of_theory"]["of_each"]})
    plain["held_that_moves_nothing"] = [name for name in held if all(not c["bar_of_theory"]["of_each"].get(name) for c in output["parameters"].values())]          # review 22, finding 5
    last = json.loads((out_dir / f"model_{len(iterations) - 1}.json").read_text()) if report else None
    quiet = bool(last) and asks_no_term(last["where"]["answers"])
    leads = noise_leads(last["where"]["answers"]) if last else None
    past = [a for a in last["where"]["answers"] if a["beyond_the_twin"] and a["question"] in ACTIONABLE] if last else []
    weak = len(past) == 1 and past[0]["excess"] < 1.5 * last["where"]["largest_excess_on_the_twin"]
    plain["quoted_with_that_said"] = with_caution(weak, sets, leads)          # "these values are quoted with that said" over lines that each said not yet to be quoted (round 16 of the cold tests)          # one answer just past its line: the guide stops there, and the command said write the next file (round 2)
    first_past = None if quiet or not last else next((f"{a['question']}: {a['answer']}" for a in sorted(last["where"]["answers"], key=lambda a: -a["excess"]) if a["beyond_the_twin"] and a["question"] in ACTIONABLE), None)
    lacks = (f"the last report leads with {leads}, which the noise of the samples must hold before this value is read" if leads
             else (f"the last report has one answer just past its line, {first_past}, as a model that lacks nothing shows on some of its synthetic sets: `--sets 100` reads it again before these values are quoted" if sets < 100
                   else f"the last report has one answer just past its line on {sets} sets, {first_past}, as a model that lacks nothing shows on some of its synthetic sets: `--sets 500` is the check that may follow, and where it stays past there, the model may lack a term, which its proposals name" if sets < 500
                   else f"the last report has one answer just past its line on {sets} sets, {first_past}, which a model that lacks nothing shows on few of them: the model may lack a term, which its proposals name")
                   if weak and first_past
             else f"the last report still has an answer beyond its twin, {first_past}, so the model may lack a term, and a term kept in its place can move this value" if first_past else None)          # the command said the loop stops and best.md said write the next file and quote nothing (round 9)
    plain["one_answer_just_past"] = weak
    plain["no_term_worth"] = no_term_worth(last, kept, refused)
    plain["seen_by_beyond"] = seen_by_beyond(last["where"]["answers"] if last else [])          # a noise of the root of the counts of the records, read beyond the twin by three families, and the command said the records ask for no more (round 15 of the cold tests)
    plain["noise_leads"], plain["the_last_report_asks"] = leads, lacks
    stop = stop_sentence(leads, weak, plain["no_term_worth"], quiet, f"terms/{len(iterations):02d}.py", sets)
    if report and stop:          # the brief of the last model says of the end of the path what best.md and the console say (the blind red teams of 0.1.0rc22)
        agree_brief(out_dir / f"brief_{len(iterations) - 1}.md", stop)
    plain["noise_from"] = d.noise.get("from")          # where the description takes the noise: the command named the differences of neighbouring samples beside a noise given as a value (round 19)
    plain["noise_from_a_column"] = d.noise.get("from") == "column"          # a column of the root of the observed counts raised families of statistics to hundreds beyond the twin while the check of the noise passed (round 17 of the cold tests)
    hidden = sorted(p.name for p in (folder / "terms").glob("[0-9][0-9].py") if nxt and p.name > Path(nxt).name)          # files after the first number missing, which the loop does not read
    beyond = sorted(p.name for p in (folder / "terms").glob("[0-9][0-9].py") if not nxt and int(p.stem) > n)          # files past --n, which the loop left in silence (round 3)
    plain["not_read"] = ([f"{', '.join(hidden)}, since {Path(nxt).name} is missing and the loop reads the term files in their order up to the first that is"] if hidden else []) \
        + ([f"{', '.join(beyond)}, since --n {n} reads {n} term files"] if beyond else [])          # said by the command alone, and best.md held no word of it (round 8 of the cold tests)
    best_table = out_dir / f"{str(output['best']).replace(' ', '_')}.json"          # the table of the best model, whose twin its report reads
    by_its_report = best_table.is_file() and json.loads(best_table.read_text()).get("admitted_on_the_twin") is False
    if by_its_report or output.get("admitted_on_the_twin") is False:          # a bar quoted where the report of the same model said that its twin does not admit the fit (round 18 of the cold tests)
        output["parameters"] = {name: dict(c, admitted_on_the_twin=False, refused_by="its report" if by_its_report else "this page") for name, c in output["parameters"].items()}
    plain["quotable"] = {name: not not_quoted(dict(c, sets=max(sets, 30)), noise["right"], None if plain["quoted_with_that_said"] else lacks) for name, c in output["parameters"].items()}          # a verdict that a program reads, where it parsed the prose of best.md (the red team of the package, his message 181)
    (out_dir / "best.json").write_text(json.dumps(plain, indent=1, default=float) + "\n")
    page = _best_page(output, noise, plain, refused, nxt, max(sets, 30), kept, quiet)
    if (Path(folder) / "literature.md").is_file():          # the search of the literature that [literature] asks, physics beside the data (his message 316)
        page += "\nThe literature: `literature.md` beside the description holds the terms that the papers give, each with its size and its paper, from the search that `[literature]` asks. A term of it enters as a term file, its paper in `why`.\n"
    (out_dir / "best.md").write_text(page)
    misnamed = sorted(p.name for p in (folder / "terms").glob("*.py") if not (len(p.stem) == 2 and p.stem.isdigit()))
    reuse.prune()          # a whole run: what the files of the folder no longer make goes
    return {"steps_read_back": reuse.read, "steps_computed": reuse.computed, "kept": kept, "refused": refused, "write_next": nxt, "best": output["best"], "terms": output["terms"], "output": output, "noise": noise, "notebook": str(book.folder), "noise_from": plain["noise_from"],
            "not_read": hidden, "not_read_past_n": beyond, "n": n, "nothing_beyond_the_twin": quiet, "one_answer_just_past": weak, "misnamed": misnamed, "noise_leads": leads,
            "the_last_report_asks": lacks, "report": f"model_{len(iterations) - 1}.md", "sets": sets, "own_trends": plain["own_trends"], "seen_by_beyond": plain["seen_by_beyond"],
            "noise_from_a_column": plain["noise_from_a_column"], "no_term_worth": plain["no_term_worth"],
            "ties": [f"{r['file']} in the place of {r['in_place']['of']}" for r in refused if r.get("in_place") and abs(r["in_place"]["criterion_less"]) < TIED],
            "better_set_aside": [f"{r['file']} in the place of {r['in_place']['of']}" for r in refused if r.get("in_place") and r["in_place"]["criterion_less"] >= TIED]}


def stop_sentence(leads: str | None, weak: bool, no_term_worth: bool, quiet: bool, nxt: str, sets: int = 30) -> str | None:
    """The one sentence of the end of a path, for the brief of the last model, in the order of best.md: the noise first,
    one answer just past its line, no answer left worth a term, nothing beyond the twin. None where the path goes on.
    One answer just past its line is read again on 100 sets below 100, and on 500 at 100 or more, where the sentence
    advised `--sets 100` after a run of 100 (the readers of the red team of his message 434, agents R3)."""
    if leads:
        return f"Next: stop, the noise first. The first kind beyond the twin is {leads}, which a term file cannot hold: a copy of the folder whose noise holds it, then the loop on the copy."
    if weak and sets < 100:
        return "Next: stop. One answer stands just past its line in the last report, as a model that lacks nothing shows on some of its synthetic sets: `--sets 100` reads it again before a term is written."
    if weak and sets < 500:
        return f"Next: stop. One answer stands just past its line in the last report, read on {sets} sets, as a model that lacks nothing shows on some of its synthetic sets: `--sets 500` is the check that may follow, and where it stays past there, the model may lack a term, which its proposals name."          # from 100 sets the next step is named (the readers of the red team of his message 436)
    if weak:
        return f"Next: stop. One answer stands just past its line in the last report, read on {sets} sets, which a model that lacks nothing shows on few of them: the model may lack a term, which its proposals name, and a term is written from its physics."          # at 500 sets the sentence sent to 500 sets (review 141)
    if no_term_worth:
        return f"Next: stop. The last term tried was refused, and no proposal beyond the twin gains in chi square what a term of one parameter pays to be kept: no answer left is worth a term. A term that the physics asks for is held at its theory, in {nxt} with hold = True."
    if quiet:
        return f"Next: stop. No answer of where, what or along stands beyond the twin: the records ask for no more term. A term that the physics asks for is held at its theory, in {nxt} with hold = True."
    return None


def agree_brief(path, sentence: str) -> bool:
    """Writes `sentence` in the place of the line of the brief that begins with Next, and says whether it did."""
    path = Path(path)
    if not path.is_file():
        return False
    lines = path.read_text().split("\n")
    at = [k for k, line in enumerate(lines) if line.startswith("Next: ")]
    if not at:
        return False
    lines[at[-1]] = sentence
    path.write_text("\n".join(lines))
    return True


def no_term_worth(last: dict | None, kept: list, refused: list) -> bool:
    """Whether no answer left is worth a term: the last term file tried was refused, and no proposal for the model of the
    last report that stands beyond its twin gains in chi square more than a term of one parameter must pay to be kept,
    TIED and the logarithm of the samples, which the degrees of freedom bound from below. The command asked for a next
    term there (the blind red team of 0.1.0rc16, mislead 6), and asked it four times where the proposals beyond the twin
    gained 0.2 and 0.1 and one within it 9.7 (the blind red teams of 0.1.0rc18, mislead 4, and of 0.1.0rc19, mislead 7).
    A proposal of an older report, with no word of its twin, is read as beyond it."""
    if not last or not refused or (kept and kept[-1]["file"] > refused[-1]["file"]):
        return False
    props = (last.get("proposals") or {}).get("for_the_model") or []
    dof = last.get("degrees_of_freedom")
    cost = TIED + (math.log(int(dof) + 1) if dof else 0.0)
    gains = [p.get("gain_in_chi_square") for p in props if p.get("beyond_the_twin", True)]
    return bool(props) and all(g is not None and g < cost for g in gains)


def own_trends(own_rows: list) -> list:
    """Each own value of the model that follows a setting along a line: a constant through its values across the curves
    lies beyond three standard deviations of its chi square, and a line through them within them. Each with the slope
    per unit of the setting and its bar, the chi square about the line and the number of curves. A line through the
    place of each curve, 400 of its bars across five powers, was left to the reader as "read them down each column", and
    the report, which reads the shared parameters alone, said that the records ask for no more (round 9 of the cold
    tests). `own_rows` is the table of best.json: one row a curve, its settings and each own value with its bar."""
    out = []
    if len(own_rows) < 3:
        return out
    names = [k for k in own_rows[0] if k not in ("record", "settings")]
    keys = [k for k, v in (own_rows[0].get("settings") or {}).items() if isinstance(v, (int, float)) and not isinstance(v, bool)]
    for name in names:
        v = np.array([float(r[name]["value"]) for r in own_rows])
        s = np.array([float(r[name]["bar_of_the_fit"] or 0.0) for r in own_rows])
        if not np.all(s > 0):
            continue
        w = 1.0 / s ** 2
        flat = float(np.sum(w * (v - np.sum(w * v) / np.sum(w)) ** 2))
        n = v.size
        if flat <= (n - 1) + 3.0 * math.sqrt(2.0 * (n - 1)):
            continue          # one value for every curve holds
        found = len(out)
        for key in keys:
            z = np.array([float(r["settings"][key]) for r in own_rows])
            if np.ptp(z) == 0.0:
                continue
            A = np.column_stack([np.ones(n), z]) * np.sqrt(w)[:, None]
            coef, *_ = np.linalg.lstsq(A, v * np.sqrt(w), rcond=None)
            chi = float(np.sum(w * (v - coef[0] - coef[1] * z) ** 2))
            if n > 2 and chi <= (n - 2) + 3.0 * math.sqrt(2.0 * (n - 2)):
                cov = np.linalg.inv(A.T @ A)
                out.append({"own": name, "setting": key, "slope": float(coef[1]), "bar": float(math.sqrt(cov[1, 1])), "chi_square": chi, "degrees_of_freedom": n - 2, "curves": n})
        if len(out) == found:          # it moves beyond its bars along no line of a setting: a place that stepped 40 bars between two temperatures was not named (round 11 of the cold tests)
            out.append({"own": name, "setting": None, "chi_square": flat, "degrees_of_freedom": n - 1, "curves": n})
    return out


def per_unit(problem, result) -> list:
    """For each shared parameter that changes along a setting by a term of `along` of power 1, p (1 + a z), z the setting
    less its mean over the curves, over its standard deviation with their number: p is the value at the mean of the
    setting, and its growth per unit of the setting is a p over that deviation, with the bar that the covariance of the
    fit gives the two. A tester converted it by hand and quoted the width without the setting at which it holds (round 5
    of the cold tests)."""
    out = []
    names = list(problem.shared)
    for j, name in enumerate(names):
        if " along " not in name or name.endswith("^2"):
            continue
        p, key = name.rsplit(" along ", 1)
        if p not in names or not all(key in r.design for r in problem.records):
            continue
        i = names.index(p)
        z = np.array([float(r.design[key]) for r in problem.records])
        mean, sd = float(z.mean()), float(z.std())
        if sd == 0.0:
            continue
        a, v = float(result.shared[j]), float(result.shared[i])
        grad = np.array([a / sd, v / sd])          # by p, then by a
        C = np.asarray(result.cov, dtype=float)[np.ix_([i, j], [i, j])]
        var = float(grad @ C @ grad)
        chain = [name]          # the terms stacked on its size, each the size of the one before along the same setting: p (1 + a z (1 + b z (1 + ...)))
        while f"{chain[-1]} along {key}" in names:
            chain.append(f"{chain[-1]} along {key}")
        where = [i] + [names.index(c) for c in chain]
        values = np.array([float(result.shared[k]) for k in where])
        Cc = np.asarray(result.cov, dtype=float)[np.ix_(where, where)]

        def along(q, zc):          # the parameter at a standardized setting, through the whole chain: the first term alone gave a width at 0 K of 0.711 where the model gave 1.154 (round 17 of the cold tests)
            inner = 1.0
            for size in q[:0:-1]:
                inner = 1.0 + size * zc * inner
            return q[0] * inner

        def value_and_bar(zc):
            steps = 1e-6 * (np.abs(values) + 1e-6)
            g = np.array([(along(values + np.eye(len(values))[k] * steps[k], zc) - along(values - np.eye(len(values))[k] * steps[k], zc)) / (2.0 * steps[k]) for k in range(len(values))])
            w = float(g @ Cc @ g)
            return along(values, zc), math.sqrt(w) if w > 0 else None
        at = []          # the value at each value of the setting, with the bar that the covariance gives: a tester propagated the size at each temperature by hand, taking the two as independent (round 6)
        for value in sorted(set(z.tolist())):
            at.append((value,) + value_and_bar((value - mean) / sd))
        z0 = (0.0 - mean) / sd
        zero = value_and_bar(z0)          # the value where the setting is zero, by the line, with the bar of the covariance: a tester extrapolated it by hand (round 9 of the cold tests)
        near = z.min() - np.ptp(z) <= 0.0 <= z.max() + np.ptp(z)          # zero within one span of the values taken, as 0 mW beside 2 to 10, and not 0 K beside 290 to 320
        out.append({"parameter": p, "setting": key, "mean": mean, "deviation": sd, "per_unit": a * v / sd, "bar": math.sqrt(var) if var > 0 else None, "at": at,
                    "at_zero": zero if near else None, "chain": chain})
    return out


def numbered(name: str, others: list) -> bool:
    """Whether the name of a kept term holds a number of its own, chosen by hand, once the names of the other parameters
    are taken out of it: "f0 along pressure" was said to hold a chosen number, the 0 of f0 (round 20 of the tests of the
    tutorial)."""
    rest = name
    for other in sorted(others, key=len, reverse=True):
        if other != name:
            rest = re.sub(rf"(?<!\w){re.escape(other)}(?!\w)", " ", rest)
    return re.search(r"(?<![A-Za-z_])\d", rest) is not None


def places(problem, result) -> list:
    """The shared parameters that are a place: each moves the place of every curve, as the model's `centre` reads it, by
    as much as it moves. The model reads the abscissa less its offset, and a researcher copied a resonance at 75.00009 kHz
    as "f0: 0.00009 ± 0.00017", with the offset in a sentence after the list (round 20 of the tests of the tutorial)."""
    if problem.centre is None:
        return []
    out = []
    for j, name in enumerate(problem.shared):
        base = np.asarray(result.shared, dtype=float)
        step = 1e-6 * max(abs(float(base[j])), 1.0)
        moved = base.copy()
        moved[j] += step
        try:
            ratios = [(problem.centre_of(k, moved, result.own[k]) - problem.centre_of(k, base, result.own[k])) / step for k in range(len(problem.records))]
        except Exception:          # a centre that cannot be read at the fit says nothing of a place
            continue
        if ratios and all(abs(r - 1.0) < 1e-4 for r in ratios):
            out.append(name)
    return out


def back_from_offset(c: dict, plain: dict) -> str:
    """A place read less the offset of the abscissa, given again with the offset added back and the unit of the abscissa,
    or nothing where no offset is given or the place has no bar."""
    offset = plain.get("abscissa_offset") or 0.0
    if not offset or c.get("bar_of_the_fit") is None or not math.isfinite(c["bar_of_the_fit"]):
        return ""
    value, bar = _pm(c["value"] + offset, c["bar_of_the_fit"])
    unit = f" {plain['abscissa_unit']}" if plain.get("abscissa_unit") else ""
    return f", that is {value} ± {bar}{unit} with the offset of the abscissa, {offset:g}, added back"


def to_quote(name: str, c: dict, noise_right: bool, at: str = "", lacks: str | None = None) -> str:
    """The line of one parameter to quote, from its five parts: the value with the bar of the fit, then the bar of theory
    and the half range over the models where each is not zero, as separate terms, or why it is not yet to be quoted. A
    newcomer met five parts and no rule to make one number of them (round 4 of the cold tests)."""
    if c["bar_of_the_fit"] is None:
        return f"- {name}: not to be quoted alone, since it is not told apart from {', '.join(c['not_told_apart_from'])}: quote the combination that the records fix"
    if not math.isfinite(c["bar_of_the_fit"]):          # -5 ± nan of the fit (round 8 of the cold tests)
        return f"- {name}: {c['value']:.6g}, not to be quoted: {NO_BAR}"
    value, bar = _pm(c["value"], c["bar_of_the_fit"])
    parts = [f"{value} ± {bar} of the fit"]
    why = not_quoted(c, noise_right, lacks)
    if c["bar_of_theory"]["of_each"] and c["bar_of_theory"]["together"] > 1e-3 * c["bar_of_the_fit"]:          # a held value that moves a parameter by the arithmetic alone moves it by none: 2.0e-10 was printed (round 5)
        each = {k: v for k, v in c["bar_of_theory"]["of_each"].items() if v and v > 1e-3 * c["bar_of_the_fit"]}
        parts.append(f"± {_two(c['bar_of_theory']['together'])} of theory" + (" (" + ", ".join(f"{_two(v)} from {k}" for k, v in each.items()) + ")" if len(each) > 1 else ""))          # which held value moves it, where several do (round 6)
    if c["spread_over_the_models"]:
        parts.append(f"± {_two(c['spread_over_the_models'])} over the models")
    return f"- {name}: " + " ".join(parts) + at + (f". Not yet to be quoted: {', and '.join(why)}" if why else "")


def setup_said(setup: dict) -> list:
    """The section of best.md that gives the setup as the description does: each instrument or setting with its value,
    its uncertainty and its unit where they are known, and how it is known. An item assumed or not recorded is where a
    better characterization of the setup starts, and a term of a physical origin starts from the setup (his message
    184, row M216)."""
    if not setup:
        return []
    def value(item):
        if "value" not in item:
            return "no value given"
        number = " ± ".join(_pm(item["value"], item["uncertainty"])) if item.get("uncertainty") else f"{item['value']:g}"
        return number + (f" {item['unit']}" if item.get("unit") else "")
    return ["", "## The setup, as the description gives it", "", "Each item with how it is known. One that is assumed or not recorded is where a better characterization of the setup starts, "
            "and a lever of kind setup on it is measured on the twin only where the model holds a value of its name.", ""] + [f"- {name}: {item['words']}, {value(item)}, {item.get('known', 'measured')}" for name, item in setup.items()]


def not_read(output: dict) -> str:
    """Why a coverage and a bias are not read in best.md: how many searches of the synthetic sets of its table ended, of
    how many, under the 30 that they need. Those sets are not the sets of the report of the model, which said that its
    searches ended on 1.00 of them while best.md said fewer than 30 ended, and nothing said why (round 18 of the cold
    tests)."""
    ended, sets = output.get("searches_that_ended"), output.get("sets_of_the_twin")
    counted = f"{ended} of the {sets} searches" if ended is not None and sets else "fewer than 30 searches"
    return (f"Not read: {counted} of the synthetic sets of this table ended, under the 30 that a coverage and a bias need, so the twin gives neither. "
            "These sets are drawn for this table and are not those of the report of the model, whose count is in its rule of the twin. More sets, `--sets 100`, read them.")


def own_to_quote(own: list, shared: bool) -> list:
    """The lines of "To quote" for the parameters of each curve. A model whose width is own to each curve, at each
    temperature of an annealing, gave a "To quote" with nothing in it, and a researcher no value to quote (round 18 of the
    cold tests): the values are in the section of the parameters of each curve, each with its bar of the fit, and the
    twin, which reads the shared parameters alone, gives them no coverage and no bias."""
    if not own:
        return []
    names = [k for k in own[0] if k not in ("record", "settings")]
    return [f"- {', '.join(names)}, of each curve: " + ("the model shares no parameter, so these are the values to quote. " if not shared else "")
            + "Each is in the section of the parameters of each curve with its bar of the fit, which counts what the shared parameters leave uncertain. "
            "The twin reads no coverage and no bias of them: quote each with its bar of the fit alone, and with the settings of its curve"]


def not_quoted(c: dict, noise_right: bool, lacks: str | None = None) -> list:
    """Why a parameter is not yet to be quoted, from its five parts, or nothing: one place for every line of "To quote",
    since the lines of a growth along a setting gave values with no caution where the noise failed its check (the red
    team of the package, his message 181)."""
    if c.get("bar_of_the_fit") is None:
        return [f"it is not told apart from {', '.join(c.get('not_told_apart_from') or [])}"]
    if not math.isfinite(c["bar_of_the_fit"]):
        return [NO_BAR]
    why = []
    if c.get("admitted_on_the_twin") is False:
        why.append("the twin does not admit the fit of this model, as its report says, so its bars are not to be read: where a statistic misses its band by little, more sets, `--sets 500`, read it again" if c.get("refused_by", "its report") == "its report"
                   else "the twin of this page does not admit the fit of this model on its own synthetic sets, where the twin of its report did on its sets: read the coverage of each parameter below, and run more sets to tell the two apart")          # until round 19 it said "as its report says" beside a report that said the twin admits the fit
    if not noise_right:
        why.append("the noise fails its check, which a missing term causes as well as a noise stated wrong")
    if lacks:          # a shared width 15 bars from its truth, under a wrong kept term, was given with no word while the last report still asked for a term (round 7)
        why.append("the last report asks for more, as the paragraph above says")
    lo, hi = bands(int(c.get("sets", 30)))["coverage"]
    if c["coverage_on_the_twin"] is not None and not lo <= c["coverage_on_the_twin"] <= hi:          # one rule with the twin: two standard deviations refused at 500 sets a coverage of 0.93 that the band of the twin, 0.92 to 0.98, admitted (round 14 of the cold tests)
        why.append(f"its coverage, {c['coverage_on_the_twin']:.2f}, is outside {lo:.3g} to {hi:.3g}, the band of the rule of the twin on {int(c.get('sets', 30))} sets")
    if c["bias_on_the_twin_in_bars"] is not None and abs(c["bias_on_the_twin_in_bars"]) > 1.0:
        why.append(f"its bias on the twin is {c['bias_on_the_twin_in_bars']:.2f} bars, and best.json gives the value less it")
    if c.get("absent_from"):
        why.append(f"{', '.join(c['absent_from'])}, which the records do not tell from this model, has no such parameter, so the records do not say that its term is there")
    return why


def with_caution(weak: bool, sets: int, leads) -> bool:
    """Whether the values of best.md are quoted with the answer just past its line said beside them: one answer just
    past its line, read on 100 sets or more, where the noise does not lead. The loop stops there, and best.md quoted
    nothing before 500 sets, hours more, so that a researcher who followed the rule of stopping ended with no number to
    write (the fresh red team of 0.1.0rc20, mislead 6)."""
    return bool(weak and sets >= 100 and not leads)


def twin_alone(refused) -> list:
    """The refused terms that the twin alone refused while the Bayesian criterion falls by more than TIED with them, fitted
    on the model before them or in the place of the term kept before them: the records prefer the model with such a term."""
    return [r for r in refused if str(r.get("why", "")).startswith("the twin does not admit")
            and max(r.get("criterion_falls") or 0.0, (r.get("in_place") or {}).get("criterion_less", 0.0)) > TIED]


def _alone_said(r) -> str:
    """A term that the twin alone refused, in words: how far the criterion falls with it, and the statistic of its twin."""
    falls, place = r.get("criterion_falls") or 0.0, r.get("in_place") or {}
    how = ([f"lowers the Bayesian criterion of {r['on']} by {falls:.1f}"] if falls > TIED else []) + ([f"lowers it by {place['criterion_less']:.1f} in the place of {place['of']}"] if place.get("criterion_less", 0.0) > TIED else [])
    reason = r["why"].split(": ", 1)[-1].split(". ", 1)[0]
    return f"{r['term']}, of {r['file']}, " + " and ".join(how) + (", refused by a hair" if r.get("by_a_hair") else "") + f", where its twin reads {reason}"


def _best_page(output, noise, plain, refused, nxt, sets: int = 30, kept=(), quiet: bool = False) -> str:
    """The output of the loop for the researcher: the best model, the words of its table defined before the table, and
    the five parts of the uncertainty of each parameter, never folded into one (row M168: a newcomer met five headings
    that no line defined)."""
    lines = [f"# The best model: {output['best']}", "", "Its terms: " + ", ".join(output["terms"]) + "."]
    if plain.get("generic"):
        lines += ["", "The description names no first model, so the loop started from the generic one of the package, one peak or dip of one width on its own height, place and baseline, `twinfit/agent/model_template.py`: "
                  "where the physics gives another curve, write it as `model.py`, name it under [addenda] and run the loop again."]
    lines += [""] + [f"- {k['term']}, in {k['file']}" + (f", {k['parameter']}" if k.get("parameter") else "") + f": {k['why']}" for k in kept if k.get("why")] if any(k.get("why") for k in kept) else []          # a term of one's own was quoted with no word of what it is: its reason travels with it (round 10 of the cold tests)
    alone = twin_alone(refused)
    if alone:          # the best model named 4347 of the criterion above a rival of its own table, which the twin refused by a hair, and the page said to run more sets, which did not end in 1 h 40 (the blind red team of 0.1.0rc18, mislead 13)
        lines += ["", "The records prefer a model that its twin refused, so the model named here is the best of those that their twin admits, and not the one that fits best: "
                  + ". ".join(_alone_said(r) for r in alone) + ". "
                  "Where each statistic that fails misses its band by less than its own standard error, the twin does not choose between the two models and the criterion does. "
                  "Where it misses by more, the bars of the model with the term are too narrow or its values biased on its own twin. "
                  "Quote no value below as that of the line without saying beside it that a model which fits better was refused by the twin alone."]
    left = plain.get("chi_square_per_degree")
    if left is not None and left > 10.0:          # 299 a degree of freedom after a kept tail of the package, and no word that a term of one's own was needed (round 10 of the cold tests)
        lines += ["", f"The best model leaves {left:.3g} of chi square a degree of freedom. Where the kept terms of the package do not bring it near 1, the shape may be one that no term of the package writes, "
                  "and a Mechanism of your own may be needed, as section 11 of the tutorial shows: what the terms kept here hold is then not to be believed."]
    order = plain.get("order")
    unknown = bool(order) and str(order).strip().lower().startswith(("not known", "unknown"))          # "not known" printed as an order that the description gives (round 13 of the cold tests)
    lines += ["", f"The description says that the order in which the curves were taken is not known: {order}. A setting taken in the order of the time cannot be told from a drift." if unknown
              else f"The curves were taken in this order, as the description says: {order}. Where a setting was taken in the order of the time, it cannot be told from a drift." if order
              else "The description does not say in which order the curves were taken, `order` under [problem]: a setting taken in the order of the time cannot be told from a drift."]
    lines += ["",
             "The best model is the one of the smallest Bayesian criterion, the chi square plus the number of fitted parameters times the logarithm of the number of samples. "
             f"The models whose criterion is within {TIED:g} of it are those that the records cannot tell apart from it.", "",
             "| model | its term | chi square | fitted parameters | Bayesian criterion | less that of the best |", "|---|---|---|---|---|---|"]
    lines += [f"| {row['model']} | {row['added'] if row.get('after') else 'the first model'} | {row['chi_square']:.2f} | {row['fitted_parameters']} | {row['criterion']:.2f} | {row['criterion_less_that_of_the_best']:.2f} |"
              for row in output.get("iterations", []) if row["ended"]]          # the criterion of a model was in model_K.json alone (round 16 of the cold tests)
    ties = [r for r in refused if r.get("in_place") and abs(r["in_place"]["criterion_less"]) < TIED]
    if ties:          # two causes of one misfit tied within 0.1, and the page said the records ask for no more (round 16 of the cold tests)
        lines += ["", "Two causes fit alike: " + ". ".join(f"{r['term']}, of {r['file']}, fitted in the place of {r['in_place']['term']}, of {r['in_place']['of']}, gives a criterion {abs(r['in_place']['criterion_less']):.2f} "
                                                        + ("lower" if r["in_place"]["criterion_less"] > 0 else "higher") + f" than the model with {r['in_place']['term']}" for r in ties)
                  + ". The records do not tell them apart: each model is in the spread over the models below, and the physics, or a design that moves one cause and not the other, chooses between them."]
    better = [r for r in refused if r.get("in_place") and r["in_place"]["criterion_less"] >= TIED and r not in alone]          # a term that the twin alone refused is said above: in a copy its twin reads it again
    if better:
        lines += ["", "A refused term is better in the place of a kept one: " + ". ".join(f"{r['term']}, of {r['file']}, fitted in the place of {r['in_place']['term']}, of {r['in_place']['of']}, lowers the criterion by {r['in_place']['criterion_less']:.2f}" for r in better)
                  + ". Take the kept term out in a copy of the folder and write this one as its file, then run the loop on the copy."]
    lines += ["", "## The words of the table", "",
             "- **The value**: the fitted value of the parameter in the best model. A parameter that a term adds is defined by its term: `python -m twinfit terms` lists them.",
             "- **The bar of the fit**: the standard deviation of the parameter that the fit gives at first order, from the curvature of its chi square.",
             "- **The coverage on the twin**: the share of synthetic sets of the best model on which the interval of 1.96 bars of the fit holds the value that made them. "
             f"0.95 is right, and with the {sets} sets read here a coverage is known within {math.sqrt(0.95 * 0.05 / sets):.2f}, one standard deviation. Well under 0.95, the bar of the fit is too small, and well over it, too large. "
             "The sets are drawn from fixed seeds, so a run gives the same numbers again, and two projects of one design and close values read close coverages: that coverage is one draw of its uncertainty, which more sets narrow.",          # 0.87 on three projects of one design read as a property of the data (round 12 of the cold tests)
             "- **The bias on the twin, in bars**: the mean of the values fitted on those sets less the value that made them, in bars of the fit. It is given beside the value and not taken out of it: `best.json` gives the value less it too.",
             "- **The bar of theory**: where the model holds values of theory, half the move of the parameter when each held value is set one bar of theory below and above its centre. "
             "None held: the description holds no value of theory, which an ansatz in its addenda would, `[addenda.ansatz.NAME]` with its value and its uncertainty.",
             "- **The spread over the models**: half the range of the parameter over the best model and the models that the loop fitted and that the records cannot tell apart from it. "
             "A term that was never written into a term file is in no model, so two terms that fit alike and were not both tried leave it at 0.",
             f"- **Not told apart**: a parameter correlated with another beyond {TOLD_APART:g} gets no bar of the fit, since the records fix a combination of the two and not each.", "",
             "## The parameters", "", f"{UNTESTED} {noise['said']}", "",
             "| parameter | value | bar of the fit | coverage on the twin | bias on the twin, in bars | bar of theory | spread over the models |", "|---|---|---|---|---|---|---|"]
    for name, c in output["parameters"].items():
        value, bar = ((f"{c['value']:.6g}", "none: not told apart from " + ", ".join(c["not_told_apart_from"])) if c["bar_of_the_fit"] is None
                      else (f"{c['value']:.6g}", "none: " + NO_BAR) if not math.isfinite(c["bar_of_the_fit"]) else _pm(c["value"], c["bar_of_the_fit"]))          # -5 ± nan for a parameter at its bound (round 8)
        twin = ("not read", "not read") if c["coverage_on_the_twin"] is None else (f"{c['coverage_on_the_twin']:.2f}", "not read" if c["bias_on_the_twin_in_bars"] is None else f"{c['bias_on_the_twin_in_bars']:.2f}".replace("-0.00", "0.00"))          # -0.00 in the column of the bias, and a bias of no bar that ended the page in a trace of the code (round 8 of the cold tests)
        if c["coverage_on_the_twin"] is not None and not noise["right"]:          # the twin takes the noise of the description as right: where it fails its check, these are not to be read (review 21)
            twin = tuple(f"{v}, not to be read" for v in twin)
        theory = (("not moved by the held values" if not c["bar_of_theory"]["together"] else "0") if c["bar_of_the_fit"] and c["bar_of_theory"]["together"] <= 1e-3 * c["bar_of_the_fit"] else _two(c["bar_of_theory"]["together"])) if c["bar_of_theory"]["of_each"] else "none held"          # a 0 read as exactly zero where the model did not read the held value (round 6)          # two significant digits, 0.00080 and not 0.0008 (round 3), and 0 for a move of the arithmetic alone (round 5)
        spread = "0" if c["spread_over_the_models"] == 0 else _two(c["spread_over_the_models"])
        lines.append(f"| {name} | {value} | {bar} | {twin[0]} | {twin[1]} | {theory} | {spread} |")
    lines += [f"- {name}: {', '.join(c['absent_from'])} has no such parameter, and the records do not tell that model from this one, so its spread over the models reads this model alone, and the records do not say that its term is there."
              for name, c in output["parameters"].items() if c.get("absent_from")]
    if not noise["right"]:
        lines += ["", "The noise fails its check: the coverage and the bias on the twin are not to be read, since the twin draws its sets at that noise."]
    if any(c["coverage_on_the_twin"] is None for c in output["parameters"].values()):          # "not read", defined nowhere, beside a report that read 0.97 on other sets (round 15 of the cold tests)
        lines += ["", not_read(output)]
    if sets < 100:          # the table of 30 sets read as settled (round 4 of the cold tests)
        lines += ["", f"Read on {sets} synthetic sets, the coverage and the bias are provisional: `--sets 100` reads each coverage within {math.sqrt(0.95 * 0.05 / 100):.2f}."]
    lines += ["", "## To quote", "", "Each parameter as a value with the bar of the fit, then, as separate terms, the bar of theory and the half range over the models where they are not zero. "
              "Where one bar is asked, add them in squares, since each comes from another source. "
              "The bias on the twin is not taken out while it is under one bar. Where the noise, the coverage or the bias fails, the line says so."]          # separate or in squares, which a newcomer asked (round 9)
    grows = {g["parameter"]: g for g in plain.get("per_unit", [])}
    caution = bool(plain.get("quoted_with_that_said"))
    if plain.get("the_last_report_asks"):
        lines += ["", f"The values below are quoted with this said beside them: {plain['the_last_report_asks']}." if caution else f"No value below is to be quoted yet: {plain['the_last_report_asks']}."]
    share = {f"{g['parameter']} along {g['setting']}": g for g in grows.values()}          # a size along a setting read as a growth per unit, or with no unit (round 9 of the cold tests)
    lines += [""] + [to_quote(name, dict(c, sets=sets), noise["right"], f", at the mean of {grows[name]['setting']} over the curves, {grows[name]['mean']:g}" + (back_from_offset(c, plain) if name in plain.get("places", []) else "") if name in grows
                              else f", a share of {share[name]['parameter']} for each deviation of {share[name]['setting']} over the curves, {share[name]['deviation']:.4g}: the line below gives it per unit" if name in share
                              else back_from_offset(c, plain) if name in plain.get("places", []) else "",
                              None if caution else plain.get("the_last_report_asks"))
                     for name, c in output["parameters"].items()]
    lines += own_to_quote(plain.get("own_parameters") or [], bool(output["parameters"]))
    lines += [f"- {g['parameter']} {'falls' if g['per_unit'] < 0 else 'grows'} by {' ± '.join(_pm(abs(g['per_unit']), g['bar'])) if g['bar'] else format(abs(g['per_unit']), '.6g')} per unit of {g['setting']}: "
              f"the size of its term times {g['parameter']} over {g['deviation']:.6g}, the deviation of {g['setting']} over the curves, with the bar that the covariance of the two gives"
              + (f", at the mean of {g['setting']}. Its values below are those of the {len(g['chain'])} terms stacked along {g['setting']} together, {', '.join(g['chain'])}" if len(g.get("chain") or []) > 1 else "")
              + (". At each value of " + g["setting"] + ": " + ", ".join(f"{s:g}, " + (" ± ".join(_pm(v, b)) if b else f"{v:.6g}") for s, v, b in g.get("at", [])) if g.get("at") else "")
              + (f". At {g['setting']} 0, by the line: " + (" ± ".join(_pm(*g["at_zero"])) if g["at_zero"][1] else f"{g['at_zero'][0]:.6g}")
                 + (back_from_offset({"value": g["at_zero"][0], "bar_of_the_fit": g["at_zero"][1]}, plain) if g["parameter"] in plain.get("places", []) and g["at_zero"][1] else "")
                 + ", past the settings of the records: it holds as far as the line holds there, another form of the term gives another value, and its bar carries none of that, so compare the forms in copies and quote the range of their values"          # a width at zero power of 1.0676 ± 0.0025 by the line, where the law of the physics gave 1.0011 (round 13 of the cold tests)
                 if g.get("at_zero") and 0.0 not in [s for s, _, _ in g.get("at", [])] else "")
              + (f". Not yet to be quoted: {', and '.join(cautions)}" if g["parameter"] in output["parameters"]
                 and (cautions := not_quoted(dict(output["parameters"][g["parameter"]], sets=sets), noise["right"], None if caution else plain.get("the_last_report_asks"))) else "") for g in grows.values()]          # a growth and its values at each setting with no caution where the noise failed (the red team of the package, his message 181)
    if plain.get("abscissa_offset"):
        lines += [""]
        lines += [f"The abscissa is read less its offset, {plain['abscissa_offset']:g}: every place is from it."]
    held = [k for k in kept if "held_at" in k]
    lines += ["", "## Terms held at their value of theory", ""] + [f"- {k['term']}, in {k['file']}: held at {k['held_at']:g} ± {k['bar_of_the_theory']:g}, which the bar of theory of each other value carries" for k in held] if held else []          # a held term was in the list of terms and nowhere else (round 4)
    small = [k for k in kept if "held_at" not in k and k.get("parameter") in output["parameters"] and output["parameters"][k["parameter"]]["bar_of_the_fit"]
             and abs(output["parameters"][k["parameter"]]["value"]) < 2.0 * output["parameters"][k["parameter"]]["bar_of_the_fit"]]
    if small:          # a second component kept first, at 0.0113 ± 0.0086 once the right term came, with no word (round 9 of the cold tests)
        lines += ["", "## Kept terms now within two bars of none", ""] + [f"- {k['parameter']}, kept in {k['file']}: a term kept before another is not looked at again when the next arrives. "
                                                                          "Try the model without it, in a copy, as the section of the refused terms says, and keep the one of the smaller Bayesian criterion." for k in small]
    unasked = [k for k in kept if k.get("asked") is False and "held_at" not in k and k["term"] in output["terms"]]
    if unasked:          # a cosine that the report put at an excess of 0.7 was kept, and its size quoted with no caution (round 12 of the cold tests)
        lines += ["", "## Kept terms that no report asked for", ""] + [f"- {k['term']}, kept in {k['file']}: the report of the model before it put no answer of where, what, along, changes along or matched beyond its twin, "
                                                                      "so the records did not ask for it, and the fall of the Bayesian criterion alone kept it. A term chosen after looking at such a report fits what the twin also shows on some of its sets: "
                                                                      "quote the values of the model without it beside these, and its size as what the records allow and not as a finding." for k in unasked]
    lines += setup_said(plain.get("setup") or {})
    lines += ["", "## Notes", "", "- The bias on the twin here is in bars of the fit, on sets of its own: the report of each model reads it in spreads of the estimate, on other sets, so the two differ.",
              f"- The coverage on the twin here is read on sets of its own too, and a coverage counted on {sets} sets scatters by about {math.sqrt(0.95 * 0.05 / sets):.2f} from one draw of the sets to another: "
              "the coverage of a parameter in the report of its model can differ from this one by two or three times that, with no fault in either."]          # 0.97 in model_1.md and 0.87 in best.md, with the bias explained and the coverage not (round 18 of the cold tests)
    if plain.get("seen_by_beyond") and plain.get("noise_from_a_column"):
        lines += [f"- {', '.join(plain['seen_by_beyond'][:3])} read beyond the twin, and the noise of each sample is a column of the curves. Where that column comes from the values of the records, "
                  "as the root of the observed counts, it moves with their own noise and these families read it, while the check of the noise can pass: give it from the fitted values of the model, "
                  "section 9 of the tutorial, before a term is written for what they read."]
    if kept:          # three forms of one law gave widths at zero power 0.036 apart, against a bar of 0.0011 and a spread of 0 (round 13 of the cold tests)
        lines += ["- The spread over the models reads the models of this folder alone. A form of a term tried in a copy, as a law of the physics in `model.py` in the place of a term of the package, is not in it: "
                  "compare such forms by their Bayesian criterion and quote the range of their values beside the five parts."]          # the notes ran into one paragraph, a level's unit glued to the line before (round 8 of the cold tests)
    area = [name for name, c in output["parameters"].items() if name.startswith(("tail ", "pedestal ", "second component at ")) and c["value"] < 0]
    lines += [f"- {name} is negative: a term of negative area takes something away from the curve, and may stand in for another term. Read it against the report before it is believed." for name in area]
    lines += [f"- {name}: the size a of a ((x - {plain['level_unit'][0]:.6g}) / {plain['level_unit'][1]:.6g})^{name.split()[-1]} under every curve, in the unit of the values, "
              "the middle and the half span of the abscissas of all the curves together" for name in output["parameters"] if re.fullmatch(r"level power \d+", name) and plain.get("level_unit")]          # the unit of a level, which no page said (round 3)
    chosen_number = [name for name in output["parameters"] if name.startswith(("tail ", "pedestal ", "second component at ", "cosine of period ", "sine of period ")) and re.search(r"\d", name)]          # "tail share" and "tail length" of model.py, both fitted, were said chosen by their names (round 15 of the cold tests)
    chosen_number += [k["parameter"] for k in kept if k.get("parameter") in output["parameters"] and k["parameter"] not in chosen_number and "held_at" not in k
                      and numbered(str(k["parameter"]), list(output["parameters"]))]          # an outflow at a scanned distance and a width chosen by hand, its bar 2.2 times too small, with no word (round 13 of the cold tests)
    if chosen_number:          # the distance of a second component had no line to quote (round 11 of the cold tests)
        lines += [f"- {', '.join(chosen_number)}: the number in its name is chosen and not fitted, so it has no bar of the fit, where the row above gives the size of the term and its bar: quote that number with the range that `python -m twinfit scan` gives about it, "
                  "and hold it with an ansatz of half that range where its uncertainty is to reach the other values, or write it in `model.py` as a shared parameter, which the fit then gives a bar and the twin a coverage."
                  + (" The size of a second component is its area over that of the curve, and the place stays that of the first line, so the middle of a doublet is the place plus half its distance."
                     if any(name.startswith("second component at ") for name in chosen_number) else "")]          # a coupling read from the label, and the place of a doublet read as its middle (round 12 of the cold tests)
        lines += [f"- {', '.join(chosen_number)}: its size holds at the number written in its term file, or in `HELD` of the model where the term reads it there, a distance, a stretch or a period. "
                  "`python -m twinfit scan` says how far the records fix that number, and the size moves with it."]
    lines += [f"- The held value {name} moves no parameter of this model: the model does not read it where a parameter is fitted, so its bar of theory of zero says nothing of {name}." for name in plain["held_that_moves_nothing"]]
    own = plain.get("own_parameters") or []
    if own:
        names = [k for k in own[0] if k not in ("record", "settings")]
        keys = list(own[0].get("settings") or {})          # beside its settings, so that a trend of an own value along one reads at a glance: a blind agent found one by hand (round 4 of the cold tests)
        lines += ["", "## The parameters of each curve", "", "Each with its bar of the fit, which counts what the shared parameters leave uncertain: the twin reads the shared ones alone. "
                  "An own value that moves with a setting is a term of its own: read them down each column.", ""]
        for t in [t for t in plain.get("own_trends") or [] if t["setting"] is None]:
            lines += [f"- {t['own']} moves from curve to curve beyond its bars, a chi square of {t['chi_square']:.1f} on {t['degrees_of_freedom']} about one value for all, along no line of a setting: "
                      "read it down its column beside the order in which the curves were taken, since a step between two curves can be a drift of the instrument as well as the physics. "
                      "Where it follows a setting along a curve and not a line, as a height that saturates with the power, write the curve that the physics gives as shared parameters, as for a line."]          # a height that saturated, with a remedy for a line alone (round 13 of the cold tests)
        for t in [t for t in plain.get("own_trends") or [] if t["setting"] is not None]:          # the place of each curve, 400 bars across five powers, with no number beside it (round 9 of the cold tests)
            lines += [f"- {t['own']} moves with {t['setting']}: a line through its values gives {' ± '.join(_pm(t['slope'], t['bar']))} per unit of {t['setting']}, "
                      f"{abs(t['slope']) / t['bar']:.0f} bars from none, a chi square of {t['chi_square']:.1f} on {t['degrees_of_freedom']} about the line. While {t['own']} is own to each curve, "
                      f"no question of the report reads it: write it in `model.py` as shared parameters that read `design[\"{t['setting']}\"]`, as section 9 of the tutorial shows, and the report reads what is left. "
                      f"Across {t['curves']} curves a line can also hold by chance, so ask whether the physics gives it."]
        lines += [""] if plain.get("own_trends") else []
        lines += ["| curve | " + " | ".join(keys + names) + " |", "|---|" + "---|" * len(keys + names)]
        lines += [f"| {row['record']} | " + " | ".join([f"{row['settings'][k]:g}" if isinstance(row["settings"].get(k), (int, float)) else str(row["settings"].get(k)) for k in keys]
                                                     + [" ± ".join(_pm(row[n]["value"], row[n]["bar_of_the_fit"])) if row[n]["bar_of_the_fit"] else f"{row[n]['value']:.6g}" for n in names]) + " |" for row in own]
        if any(not row[n]["bar_of_the_fit"] for row in own for n in names):          # own values printed with no bar under a heading that promised one (round 16 of the cold tests)
            lines += ["", "A value with no ± has no bar of the fit: the records do not fix it apart from the others, or the curve of the model does not read it."]
    by_the_sets = lambda r: sets < 100 and ("the twin does not admit" in r["why"] or "synthetic sets, and a term is read on" in r["why"] or "its search stopped after" in r["why"]) and "ended on none of its" not in r["why"]          # noqa: E731
    lines += ["", "Tried and refused:" + ("" if refused else " none.")] + [f"- {r['term']}, in {r['file']}: {r['why']}"
                                                                         + (f". Its size, fitted: {' ± '.join(_pm(r['size'], r['bar']))}" if r.get("size") is not None else "")
                                                                         + (f". Fitted alone in the place of {r['in_place']['of']}, its criterion is {abs(r['in_place']['criterion_less']):.2f} " + ("lower" if r["in_place"]["criterion_less"] > 0 else "higher")
                                                                            + f" than that of the model with {r['in_place']['of']}" if r.get("in_place") else "")
                                                                         + (f". Read on {sets} synthetic sets, this refusal can turn on more: `--sets 100` reads it again before the term is set aside" if by_the_sets(r) else "")
                                                                         for r in refused]          # right terms refused at 30 sets for searches of the twin, and the advice of best.md was to write them smooth (round 15 of the cold tests)          # an empty heading read as a list that was lost (review 21)
    after = [r for r in refused if kept and r["file"] > min(k["file"] for k in kept)]
    if after:          # the true gain was refused behind two wrong kept terms, and no page said how to take one out (round 7 of the cold tests)
        lines += ["", "A term refused after others were kept may be refused because a kept term stands in its place. To try it alone, copy this folder, keep in the copy's terms/ "
                  "only that file, renamed 01.py, and run the loop on the copy. To take a kept term out, copy the folder, leave out its file and rename those after it so that the numbers run on."]
    lines += ([""] + [f"Not read: {said}." for said in plain["not_read"]]) if plain.get("not_read") else []
    leads = plain.get("noise_leads")
    if nxt and leads:          # the page said write the next term file where the first kind beyond the twin was glitches (round 7 of the cold tests)
        lines += ["", f"The loop stops here: the first kind beyond the twin in the last report is {leads}, which a term file cannot hold. {NOISE_KINDS[leads]} "
                  "Do it in a copy of this folder and run the loop on the copy. Until then, a term kept here may fit it and not a part of the curve."]
    elif nxt and plain.get("one_answer_just_past"):          # the command said the loop stops, the page said write the next file (round 9 of the cold tests)
        lines += ["", "The loop stops here: one answer of where, what, along, changes along or matched, the questions that name a term, stands just past its line in the last report, as a model that lacks nothing shows on some of its synthetic sets. "
                  + ("`--sets 100` reads it again before a term is written." if sets < 100 else f"It was read on {sets} sets: `--sets 500` is the check that may follow, and where it stays past there, the model may lack a term, which its proposals name." if sets < 500
                     else f"It was read on {sets} sets: the model may lack a term, which its proposals name, and a term is written from its physics.")]
    elif nxt and plain.get("no_term_worth"):          # the page asked for a next term where no answer left could be kept (the blind red team of 0.1.0rc16, mislead 6)
        lines += ["", f"The loop stops here: the last term tried was refused, and no proposal of the last report beyond its twin gains in chi square what a term of one parameter pays to be kept, {TIED:g}, the margin of the Bayesian criterion, and the logarithm of the samples, so no answer left is worth a term. "
                  f"A term that the physics asks for is held at its value of theory, in {nxt}, with hold = True."]
    elif nxt and quiet:          # the page said write the next file where the command said the records ask for no more (round 4 of the cold tests)
        lines += ["", f"The loop stops here: no answer of where, what or along stands beyond the twin, so the records ask for no more term. A term that the physics asks for is held at its value of theory, in {nxt}, with hold = True."]
    elif nxt:
        lines += ["", f"The loop stops here: write {nxt} with the next term, then run it again."]
    return "\n".join(lines) + "\n"
