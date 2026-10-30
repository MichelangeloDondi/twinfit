"""The readings of the package counted on real draws (the red team of 2026-10-13, rows M43 to M45 of the plan).

Until this program, a page said "the chi square would see the missing term on 11 of 100 draws, the best family on 99".
No draw stood behind it: it was a first order chance computed on one fit, and the best family was picked with the truth
in hand. Each of the 16 readings was tested at 5 in 100 on its own, with no level for all of them together. And the
opponent was the chi square of the fit alone. Here every share is a count of draws, with its standard error:

    null       draws of the 8 synthetic traces of a right model, under that model. How often each reading reports, how
               often any of them does, and the one level at which any of them reports on 5 draws in 100
    lacking    draws of the traces of a detail under the model that lacks its term, at the size of config.toml and at
               the size that the recorded traces ask. How often the misfit is reported by the chi square of the fit, by
               the usual checks of an analyst, by each family at the level of the null, by the test of the named term,
               and by the fit made again with the term at the weights of the first fit, with AIC and BIC

    python3 tools/autoqueue.py add toy/power.py
                                        the counts of the table, through the queue (his message 147): the count in
                                        which nobody names the term, `HIDDEN` draws of each case whose term is hidden
                                        and twice as many of each right model, after the one level for all the readings
                                        of the case and its licence of first order, counted on synthetic sets of its
                                        first draw. The counts of the named terms, `null` and `lacking`, were made
                                        before theory was held: they are not in the table, and `KINDS` names them when
                                        they are made again.
                                        And the cut window: `CUT_DRAWS` draws of the first model, each with its window
                                        whole, at a half and at a quarter, the width by the fit, by the moments of the
                                        package and by two usual rules of moments. The queue asks `units`, runs `slice`
                                        and reads `progress`, and a count outside it is refused
    python toy/power.py table           counts the kept rows again and writes toy/power.json
    python toy/power.py --check         ends without error when toy/power.json holds counts of rows made from the files
                                        as they are, which give the same shares when they are counted again

A row is kept in runs/ as it is made, so an interruption loses one row. The first draws of a right model are the draws
of `toy/draws.json`. No number of the physics is written here.
"""
from __future__ import annotations

import copy
import fcntl
import hashlib
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2_law

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import twinfit  # noqa: E402
from twinfit import criteria, families, propose, standard  # noqa: E402
from twinfit import unsuspected as U  # noqa: E402

from toy import analysis, draws, fits, model  # noqa: E402

LEVEL = 0.05          # the level of a test on its own, and the share of right draws on which all the readings together may report
NULLS = 150          # draws of a right model: a false alarm of 5 of 100 is told from 10 of 100 on 150, row M204 (400 until version 0.0.5)
LACKING = 30          # 30 a hidden case, the first stage of row M474 at 0.1.0rc16, its right model kept at `RIGHT`. Until then 100: draws of a model that lacks a term: a share near 70 of 100 is read within 9 points on 100, row M204 (200 until version 0.0.5)
WORKERS = 9          # processes side by side: the machine has ten cores
FIT = analysis.SHORT["every sample: the fit"]          # the reading that is the chi square of the fit
AGAIN = 40          # the steps of a fit made again: it starts at the first fit, and a term that the traces do not hold walks to its bound and stays
NOTHING = "nothing is named"          # the kind of the first guess on a draw on which no term leads
SMALLER = (0.5, 0.25)          # the parts of the size of the twin at which a missing term is drawn too: at the whole size every opponent that is told the term sees it on every draw, and a comparison at the ceiling ranks nothing (review 10, moderate 3)
SETS = 200          # the synthetic sets on which the level for all the readings of a case is counted by the package
LICENCE = 100          # the synthetic sets on which first order is set beside the fit made again, for the licence of a case
HALF = 0.5          # the part of the size of the twin at which a hidden term is drawn too
RIGHT = 50          # the draws of each of the four right models, 200 pooled, whose false alarm C5b of toy/thresholds.toml reads on at least 180
HIDDEN = LACKING          # the draws of a case whose term is hidden, and twice as many of a model that is right, which `units` makes: the sizes of row M204. Until 2026-10-17 they were 200 and 400, constants of their own, and the sizes that row M204 set reached no count
APART = 0.10          # a size that the recorded traces ask is drawn as a second size when it differs from the one of config.toml by more than this share


def said(value: float) -> float:
    """A number as a row keeps it: six significant digits, which a count of draws does not need more of."""
    return float(f"{float(value):.6g}")


def with_term(d: int, term: str, const: dict, config: dict, free: list):
    """The model of detail `d` with one further term of the toy, as the package asks a model: the shared parameters are
    those of `free`, then the size of the term. Every other term that the model of detail `d` does not hold stays absent,
    at zero, as in `analysis.of_the_terms`."""
    last = max(int(k) for k in config["detail"])
    inside = {name for name, e in config["parameter"].items() if e["enters"] <= d}
    base = dict(model.truth(config), **{name: 0.0 for name in config["parameter"] if name not in inside})
    names = list(free) + [term]

    def of(x, shared, own, design):
        values = dict(base, **{name: float(v) for name, v in zip(names, shared)})
        u = np.asarray(x, dtype=float) - own[1]
        top = float(model.shape(last, 0.0, design["T"], values, const))
        return own[0] * model.shape(last, u, design["T"], values, const) / top + own[2] + own[3] * u
    return of


def again(problem, result, d: int, term: str, const: dict, config: dict) -> dict:
    """The fit made again with one term more, at the weights of the first fit: the fall of the chi square is then, at
    first order, the likelihood ratio of the two models. With weights of its own the larger model weighs the samples by
    its own line, and the difference of the two chi squares is no such ratio. A term whose lower bound is not zero
    leaves the two models not quite nested, and the fall may be a little under zero. Returns the fall, the fitted size
    of the term with its uncertainty, whether the search ended away from a bound, and whether the larger fit says that
    two of its parameters are one: a size is then a size of both, and no size of the term (review 9, moderate 4)."""
    free, entry = list(problem.shared), config["parameter"][term]
    entries = [config["parameter"][name] for name in free] + [entry]
    larger = problem.with_model(with_term(d, term, const, config, free), tuple(free) + (term,), problem.own, list(result.shared) + [entry["start"]], list(result.own),
                                name=f"detail {d} with {term}", bounds_shared=([e["lower"] for e in entries], [e["upper"] for e in entries]))
    got = twinfit.fit(larger, list(result.shared) + [entry["start"]], result.own, max_iter=AGAIN)
    finite = bool(np.all(np.isfinite(got.cov)))
    bars = np.sqrt(np.clip(np.diag(got.cov), 0.0, None)) if finite else np.full(len(free) + 1, np.nan)
    ended = bool(got.converged and not got.at_bound and np.isfinite(bars[-1]) and bars[-1] > 0.0)
    return {"gain": said(result.chi2 - got.chi2), "value": said(got.shared[-1]), "bar": said(bars[-1]) if ended else None, "ended": ended, "tied": bool(finite and fits.tied(got.cov))}


def guessed(lin, terms: list) -> dict:
    """What the table of the package names with no hint but the names of the terms of the toy: the first candidate, how
    many lead, how many are left when the ones that do the work of an earlier one are taken out, and for each term of
    the toy its rank, whether it leads, and whether it stands with the first: it is the first, or it is like it, or it
    does its work. The package marks the terms that are alike, and a guess can only name a class of them."""
    table = propose.propose(lin, [c for _, _, c in terms] + [c for c in U.held_by_the_residuals(lin) if c.family == "trend"])          # beside the terms, the one basis fixed before the residuals are read, as `unsuspected.read` ranks by default: a shape built from the residuals scores their own size and leads whatever is missing
    ranked = [c for c in table["candidates"] if c["status"] == "ranked"]
    leading = [c for c in ranked if c["name"] in table["leading_set"]]
    kept = []
    for c in leading:          # in the order of their rank: one that is like a kept one, or does the work of the first, names nothing new
        if not (c.get("does_the_work_of_the_first") and kept) and not any(c["name"] in k.get("like", []) for k in kept):
            kept.append(c)
    first = ranked[0] if ranked else None
    names = {cand.name: name for name, _, cand in terms}
    of_terms = {}
    for rank, c in enumerate(ranked):
        if c["name"] in names:
            with_first = bool(first is not None and (c is first or c["name"] in first.get("like", []) or (c["name"] in table["leading_set"] and c.get("does_the_work_of_the_first"))))
            of_terms[names[c["name"]]] = {"rank": rank + 1, "leads": bool(c["name"] in table["leading_set"]), "with_first": with_first}
    return {"looked_at": int(table["candidates_looked_at"]), "first": first["name"] if first else None, "first_kind": first["family"] if first else None,
            "first_level": said(first["family_level"]) if first else None, "leading": len(leading), "classes": len(kept), "terms": of_terms}


def read_size(amplitude: float, its_bar: float, added: float, entry: dict) -> dict:
    """The size that a test of a term reads, at first order, from the amplitude of the term against the size that was
    added. Where the line changes as the square of the size, the amplitude is the square of the size over the one that
    was added, and a negative one is the sign that the term cannot have: no size is read then."""
    if entry.get("effect") != "square":
        return {"size": said(amplitude * added), "bar": said(its_bar * added)}
    if not amplitude > 0.0:
        return {"size": None, "bar": None}
    root = math.sqrt(amplitude)
    return {"size": said(root * added), "bar": said(its_bar * added / (2.0 * root))}


def row(d: int, free: list | None, rows: list, const: dict, config: dict) -> dict:
    """Everything that one draw says under one model: the level of each of the 16 readings, the usual checks, the test of
    each term of the toy that the model lacks, the fit made again with each, and the guess of the package."""
    problem, result, lin = analysis.at(d, rows, const, config, free)
    told = bool(result.converged and not result.at_bound and not fits.tied(result.cov))
    out = {"told": told, "values": {name: said(v) for name, v in zip(problem.shared, result.shared)}, "chi2": said(result.chi2), "dof": int(result.dof), "samples": int(problem.samples),
           "readings": {analysis.SHORT[name]: said(families.family(lin, name, analysis.CENTRE).reading()["level"]) for name in analysis.SHORT},
           "usual": {"runs": said(standard.runs(lin)["level"]), "bins": said(standard.binned(lin)["level"]), "neighbours": said(standard.serial(lin)["level"])}}
    terms = analysis.of_the_terms(problem, result, const, config, small=True, detail=d)
    out["named"], out["again"] = {}, {}
    for name, size, cand in terms:
        test = standard.against(lin, cand.columns[0], name)
        out["named"][name] = dict({"score": said(test["score"]), "level": said(test["level"]), "visible": said(test["visible_share"])},          # the share of the term that the fit leaves to be seen: under `analysis.HIDDEN` no size is read (review 9, moderate 4)
                                  **read_size(test["size"], test["its_standard_deviation"], size, config["parameter"][name])) if test["seen"] else None
        out["again"][name] = again(problem, result, d, name, const, config)
    out["guess"] = guessed(lin, terms)
    return out


def right(config: dict) -> list:
    """Each right model with its own traces, as (its name in the tables, the detail of the traces, the detail of the
    model, its free parameters, nothing planted)."""
    return [(draws.key(d, free, config), d, d, free, None) for d, free in draws.every(config)]


def asked(detail: int, root: Path = ROOT) -> dict:
    """The size that the recorded traces ask of the term of a detail: what the joint fit of all 8 gives it under the model
    with that term, as the built folder holds it. Empty when the folder holds no such fit."""
    _, config = model.load(Path(root) / "toy", root)
    record = Path(root) / config["folder"]["name"] / "BUILD.json"
    if not record.exists():
        return {}
    entry = json.loads(record.read_text())["details"].get(str(detail))
    held = Path(root) / config["folder"]["name"] / entry["folder"] / "numbers.json" if entry else None
    if held is None or not held.exists():
        return {}
    fit = (json.loads(held.read_text()).get("recorded") or {}).get("with", {}).get("all")
    return {name: float(fit["values"][name]) for name, e in config["parameter"].items() if fit and e["enters"] == detail and name in fit["values"]}


def lacking(config: dict, root: Path = ROOT) -> list:
    """Each case in which a model lacks the term of the next detail, as (its name, the detail of the traces, the detail of
    the model, its free parameters, the sizes at which the term is drawn, or None for those of config.toml)."""
    out = []
    for made in config["folder"]["details"][1:]:
        d, free = made - 1, fits.lacking(made, config)
        term = [name for name, e in config["parameter"].items() if e["enters"] == made]
        out.append((f"{made} under {draws.key(d, free, config)}", made, d, free, None))
        other = {name: size for name, size in asked(made, root).items() if abs(size / float(config["parameter"][name]["true"]) - 1.0) > APART}
        if other and set(other) == set(term):
            at = " and ".join(f"{said(size):g}" for size in other.values())          # the size is in the name, and so in the file of the rows: rows drawn at one size are never counted as rows at another
            out.append((f"{made} under {draws.key(d, free, config)}, at {at} as the recorded traces ask", made, d, free, other))
        for part in SMALLER:          # and at parts of the size of the twin, where an opponent that is told the term no longer sees it on every draw
            out.append((f"{made} under {draws.key(d, free, config)}, at {part:g} of the size of the twin", made, d, free, {name: said(part * float(config["parameter"][name]["true"])) for name in term}))
    return out


def hidden(config: dict) -> list:
    """The cases of the count in which nobody names the term (CLAUDE.md, rule 2d, and his message 133: "the cases when nobody
    can name a terms"). For each detail after the first, the model that lacks its term: on its own traces, where it is
    right, and on the traces of that detail, whose term is hidden from it, at the size of the twin and at a half of it.
    As (its name, the detail of the traces, the detail of the model, its free parameters, the sizes or None)."""
    out = []
    for made in config["folder"]["details"][1:]:
        d, free = made - 1, fits.lacking(made, config)
        under = draws.key(d, free, config)
        half = {name: said(HALF * float(e["true"])) for name, e in config["parameter"].items() if e["enters"] == made}
        out += [(f"{under}, right", d, d, free, None), (f"{made} hidden under {under}", made, d, free, None), (f"{made} hidden under {under}, at {HALF:g} of the size of the twin", made, d, free, half)]
    return out


def drawn_with(config: dict, sizes: dict | None) -> dict:
    """The values that the traces of a case are drawn with: those of config.toml, with the sizes of the case in place of
    the true ones."""
    if not sizes:
        return config
    drawn = copy.deepcopy(config)
    for term, size in sizes.items():
        drawn["parameter"][term]["true"] = float(size)
    return drawn


def level_of(case: tuple, const: dict, config: dict, root: Path = ROOT, workers: int = WORKERS) -> dict:
    """What the package counts once for a case, on its first draw, as a researcher would on his one set of records: the
    one level for all the readings, on `SETS` synthetic sets of the fit. And where a generic term is fitted, the licence
    of first order, on `LICENCE` synthetic sets of the model with that term. Kept in runs/, under what a row depends on."""
    name, made, d, free, sizes = case
    path = Path(root) / "runs" / f"power_level_{name.replace(' ', '_').replace(',', '')}_{digest(root)[:12]}.json"
    if path.exists():
        return json.loads(path.read_text())
    problem, result, lin = analysis.at(d, fits.synthetic(made, const, drawn_with(config, sizes), 1), const, config, free)
    level = U.level_for_all(problem, result, list(analysis.SHORT), analysis.CENTRE, sets=SETS, workers=workers)
    out = {"level": level, "licence": None, "no_licence": "the case is a right model: no term is fitted"}
    if made != d:
        known = U.known_parameters(problem, result, lin, names=list(analysis.COMBINED), centre=analysis.CENTRE)
        try:
            out.update(licence=U.licence(problem, result, known, sets=LICENCE, workers=workers), no_licence=None)
        except ValueError as why:
            out["no_licence"] = str(why)
    path.parent.mkdir(exist_ok=True)
    aside = path.with_suffix(f".{os.getpid()}.part")
    aside.write_text(json.dumps(out))
    os.replace(aside, path)
    return out


def unsuspected(d: int, free: list | None, rows: list, const: dict, config: dict, level: dict) -> dict:
    """Everything that one draw says under one model when nobody names the term: the four answers of the package
    (twinfit/unsuspected.py), with no term of the toy in its list. Is something missing, at the one level of the case.
    Of which kind, by the classes of the generic dictionary. And each known parameter by five estimates: the fit, the
    generic term of the largest gain, fitted and at first order, and the combination of the families that is free of
    that term, or of the leading class. Beside them the usual checks of an analyst, which need no name either."""
    problem, result, lin = analysis.at(d, rows, const, config, free)
    told = bool(result.converged and not result.at_bound and not fits.tied(result.cov))
    got = U.read(problem, result, level, names=list(analysis.COMBINED))
    kind, known = got["leading_class"], got["known"]
    pair = lambda e: None if e is None or e.get("value") is None else [said(e["value"]), said(e["bar"])]          # noqa: E731
    return {"told": told, "chi2": said(result.chi2), "dof": int(result.dof),
            "readings": {analysis.SHORT[name]: said(p) for name, p in got["missing"]["levels"].items()},
            "usual": {"runs": said(standard.runs(lin)["level"]), "bins": said(standard.binned(lin)["level"]), "neighbours": said(standard.serial(lin)["level"])},
            "missing": bool(got["missing"]["something_is_missing"]),
            "kind": None if kind is None else {"families": kind["families"], "members": len(kind["members"]),
                                               "bias": {name: [said(span["lowest"]), said(span["highest"])] for name, span in kind["bias_if_a_member_is_what_is_missing_in_standard_deviations"].items()}},
            "term": None if known["term"] is None else {"name": known["term"]["name"], "family": known["term"]["family"], "gain": said(known["term"]["gain"])},
            "tried": {t["term"]: (None if t["gain"] is None else said(t["gain"])) for t in known["terms_tried"]},
            "known": {name: {"fit": pair(e["the_fit"]), "free": pair(e["free_of_the_class"]), "generic": pair(e["with_a_generic_term"]), "first_order": pair(e["with_that_term_at_first_order"]), "free_term": pair(e["free_of_that_term"])}
                      for name, e in known["parameters"].items()}}


def one(job: tuple) -> tuple:
    """One row, for a pool of processes: its place and what it holds. A draw on which the row cannot be made is kept as
    a row that says why: one failure in a pool threw away every row that was still queued (review 9, moderate 3). With
    the level of its case, the row is that of the case in which nobody names the term."""
    name, made, d, free, sizes, draw, const, config = job[:8]
    level = job[8] if len(job) > 8 else None
    try:
        traces_of = fits.synthetic(made, const, drawn_with(config, sizes), draw)
        return name, draw, row(d, free, traces_of, const, config) if level is None else unsuspected(d, free, traces_of, const, config, level)
    except (Exception, SystemExit) as why:
        return name, draw, {"failed": f"{type(why).__name__}: {why}"[:300]}


MAKES = ("said", "with_term", "again", "guessed", "read_size", "row", "one", "drawn_with", "level_of", "unsuspected")          # the functions of this program that make a row
MADE_WITH = ("AGAIN", "SETS", "LICENCE")          # and the constants that they read
RUNS = ("toy/analysis.py", "toy/power.py")          # the programs whose imports say which files of the package a row runs
READS = ("at", "of_the_terms", "SHORT", "CENTRE", "COMBINED", "TERM", "CONVERGED", "HALVINGS")          # what a row takes from the analysis: its other functions draw and combine, and move no row


def sources(text: str, names: tuple) -> bytes:
    """The text of the functions and of the constants that are named, as a program holds them at its top: what a digest
    reads of a program of which only a part makes a row."""
    import ast
    out, found = b"", set()
    for node in ast.parse(text).body:
        named = [node.name] if isinstance(node, ast.FunctionDef) else [target.id for target in node.targets if isinstance(target, ast.Name)] if isinstance(node, ast.Assign) else []
        if set(named) & set(names):
            out += ast.get_source_segment(text, node).encode() + b"\0"
            found |= set(named) & set(names)
    if found != set(names):          # a name that was changed would leave its part out of the digest, and the rows would go stale unseen (review 10, minor 10)
        raise SystemExit(f"REFUSED: {sorted(set(names) - found)} are not at the top of the program that a row is made with")
    return out


def digest(root: Path = ROOT) -> str:
    """What a row depends on: every program and every value that a kept fit of the last detail depends on, of the
    analysis what a row takes from it, and of this program the functions that make a row with the constants that they
    read. The last detail, whatever the detail of the traces: a row tries every term that its model lacks, the later
    ones too, so a change of a later term moves the rows of the models before it, as it moves the readings of a
    subfolder (review 9, moderate 1). What counts the rows is not in it: a repair of a count or of the table leaves the
    rows as they are, which are hours of fits. Of the package, the files that a row runs: a file of it that no row
    reaches moves none."""
    root = Path(root)
    h = hashlib.sha256(draws.digest(root).encode())
    h.update(sources((root / "toy" / "analysis.py").read_text(), READS))
    h.update(sources((root / "toy" / "power.py").read_text(), MAKES + MADE_WITH))
    draws.hashed(h, root, draws.package_of(root, RUNS))          # the files of the package that a row runs: a kept fit runs fewer
    return h.hexdigest()


def kept(root: Path, kind: str, name: str) -> Path:
    """The file of the rows of one case, under what its rows depend on."""
    return Path(root) / "runs" / f"power_{kind}_{name.replace(' ', '_').replace(',', '')}_{digest(root)[:12]}.jsonl"


def read(path: Path) -> dict:
    """The rows that a file keeps, by their draw. A line that is no row, as a power cut leaves the last one, is set
    aside into a file beside it and the rest is read: the draw is then made again. A draw that the file holds twice is
    refused: two runs wrote into it, and a count would take both (review 9, moderate 2 and minor 20)."""
    out, torn = {}, []
    for line in (path.read_text().splitlines() if path.exists() else []):
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            torn.append(line)
            continue
        if not isinstance(r, dict) or "draw" not in r:          # a line that is no row is set aside as a torn one is (review 10, minor 10)
            torn.append(line)
            continue
        if r["draw"] in out:
            raise SystemExit(f"REFUSED: {path.name} holds draw {r['draw']} twice: two runs wrote into it. Keep one line a draw before it is counted")
        out[r["draw"]] = r
    if torn:
        with path.with_suffix(".torn").open("a") as aside:          # a second tear does not write over the first
            aside.write("\n".join(torn) + "\n")
        path.write_text("".join(json.dumps(out[k]) + "\n" for k in sorted(out)))
    return out


def run(kind: str, n: int, root: Path = ROOT, workers: int = WORKERS, last: int | None = None) -> dict:
    """The rows of every case of a kind, `n` draws each, side by side. A row that a run before this one kept is not made
    again. `last`: only the cases whose traces are of a detail up to this one, since the first tiers of the model come
    first (the owner's message 124). With no workers the rows are made one after the other, in this process. One run
    of a kind at a time: the kept rows are read once the lock is held, so that a second run waits and then finds what
    the first made. Each row is kept as it comes back, a failed one too. Returns for each case the rows in the order of
    their draws."""
    from concurrent.futures import ProcessPoolExecutor, as_completed
    root = Path(root)
    const, config = model.load(root / "toy", root)
    cases = [c for c in {"null": lambda: right(config), "lacking": lambda: lacking(config, root), "unsuspected": lambda: hidden(config)}[kind]() if last is None or c[1] <= last]
    (root / "runs").mkdir(exist_ok=True)
    paths = {name: kept(root, kind, name) for name, *_ in cases}
    draws_of = lambda c: 2 * n if kind == "unsuspected" and c[1] == c[2] else n          # noqa: E731  a right model of that count has twice the draws of a case: its shares are near 5 in 100
    with (root / "runs" / f"power_{kind}.lock").open("w") as mine:
        fcntl.flock(mine, fcntl.LOCK_EX)
        got = {name: read(paths[name]) for name in paths}
        levels = {c[0]: (level_of(c, const, config, root, max(workers, 1))["level"],) for c in cases} if kind == "unsuspected" else {c[0]: () for c in cases}          # counted once a case, before its rows
        jobs = [(c[0], c[1], c[2], c[3], c[4], 1 + k, const, config) + levels[c[0]] for c in cases for k in range(draws_of(c)) if 1 + k not in got[c[0]]]

        def keep(made: tuple) -> None:
            name, draw, made_row = made
            got[name][draw] = dict(made_row, draw=draw)
            with paths[name].open("a") as fh:
                fh.write(json.dumps(got[name][draw]) + "\n")
        if workers < 1:
            for job in jobs:
                keep(one(job))
        else:
            with ProcessPoolExecutor(min(workers, max(len(jobs), 1))) as pool:
                for done in as_completed([pool.submit(one, job) for job in jobs]):
                    keep(done.result())
    most = {c[0]: draws_of(c) for c in cases}
    return {name: [got[name][k] for k in sorted(got[name]) if k <= most[name]] for name in got}


KINDS = {"unsuspected": HIDDEN}          # the counts that the queue makes, with the draws of a case: those of the named terms, `null` and `lacking`, were made before theory was held and are not in the table


def level_path(root: Path, name: str) -> Path:
    """Where `level_of` keeps the level of a case: the queue asks whether it is made before it starts the rows of the case."""
    return Path(root) / "runs" / f"power_level_{name.replace(' ', '_').replace(',', '')}_{digest(root)[:12]}.json"


def cases_of(kind: str, config: dict, root: Path = ROOT) -> list:
    return {"null": lambda: right(config), "lacking": lambda: lacking(config, root), "unsuspected": lambda: hidden(config)}[kind]()


def units(root: Path = ROOT) -> list:
    """What the queue asks first (his message 147): each case of each kind with its draws not yet kept, the level of each
    case of the count in which nobody names the term as a unit that comes before its draws, and the draws of the cut
    window. A right model of that count has twice the draws of a case."""
    root = Path(root)
    const, config = model.load(root / "toy", root)
    out = []
    for kind, n in KINDS.items():
        for name, made, d, free, sizes in cases_of(kind, config, root):
            after = []
            if kind == "unsuspected":
                ready = level_path(root, name).exists()
                out.append({"condition": f"level|{name}", "units": [] if ready else [1], "after": [], "kind": "level"})
                after = [] if ready else [f"level|{name}"]
            most = RIGHT if kind == "unsuspected" and made == d else n
            kept_rows = read(kept(root, kind, name))
            out.append({"condition": f"{kind}|{name}", "units": [k for k in range(1, most + 1) if k not in kept_rows], "after": after, "kind": kind})
    window = read(cut_kept(root))
    out.append({"condition": "cut_window", "units": [k for k in range(1, CUT_DRAWS + 1) if k not in window], "after": [], "kind": "cut_window"})
    return out


def slice_of(condition: str, first: int, last: int, limit: float = float("inf"), root: Path = ROOT) -> dict:
    """The units of a condition from `first` to `last` not yet made, each kept as it ends, until `limit` seconds have
    passed: the level of a case, made on one core, the draws of a case, or the draws of the cut window. The rows are
    made by the functions that `run` and `run_cut` use, so that a row of the queue is the row of a run."""
    import resource
    import time
    root = Path(root)
    const, config = model.load(root / "toy", root)
    kind, _, name = condition.partition("|")
    case = next((c for c in cases_of("unsuspected" if kind == "level" else kind, config, root) if c[0] == name), None) if name else None
    began, ids, seconds = time.time(), [], []
    for k in range(first, last + 1):
        t = time.time()
        if kind == "level":
            if level_path(root, name).exists():
                continue
            level_of(case, const, config, root, 1)
        else:
            path = cut_kept(root) if kind == "cut_window" else kept(root, kind, name)
            if k in read(path):
                continue
            if kind == "cut_window":
                draw, made_row = cut_one((k, const, config))
            else:
                extra = (level_of(case, const, config, root, 1)["level"],) if kind == "unsuspected" else ()
                _, draw, made_row = one(tuple(case) + (k, const, config) + extra)
            with path.open("a") as fh:          # one line a draw, written whole: waves of one case append side by side
                fh.write(json.dumps(dict(made_row, draw=draw)) + "\n")
        ids.append(k)
        seconds.append(time.time() - t)
        if time.time() - began > limit:
            break
    return {"done": ids, "seconds": seconds, "rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


def progress(root: Path = ROOT) -> str:
    """The preliminary result after a wave: for each case with rows, the rows kept, those that failed, and in the count in
    which nobody names the term the share of draws on which the package says that something is missing."""
    root = Path(root)
    _, config = model.load(root / "toy", root)
    said_ = []
    for kind in KINDS:
        for name, *_ in cases_of(kind, config, root):
            rows = list(read(kept(root, kind, name)).values())
            if not rows:
                continue
            good = [r for r in rows if "failed" not in r]
            line = f"{kind} {name}: {len(rows)} rows, {len(rows) - len(good)} failed"
            if kind == "unsuspected" and good:
                line += f", something missing on {sum(bool(r['missing']) for r in good)}"
            said_.append(line)
    return "; ".join(said_) or "no row yet"


def share(flags) -> dict:
    """A share of draws with its standard error, and the draws that it counts."""
    flags = [bool(f) for f in flags]
    n = len(flags)
    p = sum(flags) / n if n else float("nan")
    return {"share": p, "standard_error": math.sqrt(p * (1.0 - p) / n) if n else float("nan"), "of": n}


def lowest(r: dict, without_the_fit: bool = False) -> float:
    """The lowest level among the readings of a row: that of the reading that reports first."""
    return min(v for k, v in r["readings"].items() if not (without_the_fit and k == FIT))


def of_the_null(rows: list) -> dict:
    """What the draws of a right model say. For each reading the share of draws on which it reports at `LEVEL` on its own.
    The share on which any of them does. And the one level for all of them together: the level under which the lowest
    of the 16 falls on `LEVEL` of the draws. The same for the usual checks, for the test of each term of the toy and for
    the fit made again with it. Each of those two has its own line here: the level of the test under which it falls on
    `LEVEL` of these draws, and the 95th centile of the fall of the chi square. A term whose size is held away from
    zero by a bound follows neither the chi square law (review 9, minor 14). Every draw is counted, also one whose fit
    does not tell its parameters apart, as an analyst meets it: `told` says how many do. Three more, for the opponents
    that are told nothing (review 10): `any_usual`, how often any of the three usual checks reports at `LEVEL`,
    `level_for_the_usual_together`, the one level for the three, and `level_of_the_fit_at_the_5th_centile`, the level
    under which the chi square of the fit alone falls on `LEVEL` of these draws."""
    out = {"draws": len(rows), "ended": int(sum(r["told"] for r in rows)), "told": share(r["told"] for r in rows),
           "readings": {name: share(r["readings"][name] < LEVEL for r in rows) for name in rows[0]["readings"]},
           "any_reading": share(lowest(r) < LEVEL for r in rows),
           "level_for_all_together": float(np.quantile([lowest(r) for r in rows], LEVEL)),
           "usual": {name: share(r["usual"][name] < LEVEL for r in rows) for name in rows[0]["usual"]},
           "guess_names_a_term": share(r["guess"]["leading"] > 0 for r in rows), "terms": {}}
    out["as_many_tests_on_their_own"] = math.log(1.0 - LEVEL) / math.log(1.0 - out["level_for_all_together"]) if 0.0 < out["level_for_all_together"] < 1.0 else float("nan")
    out["any_usual"] = share(min(r["usual"].values()) < LEVEL for r in rows)          # the three checks of an analyst are an opponent too, and they get what the readings get: one level for the three together (review 10, moderate 1)
    out["level_for_the_usual_together"] = float(np.quantile([min(r["usual"].values()) for r in rows], LEVEL))
    out["level_of_the_fit_at_the_5th_centile"] = float(np.quantile([r["readings"][FIT] for r in rows], LEVEL))          # and the chi square of the fit a level of its own: at 5 in 100 it does not report on 5 draws in 100 of every right model (review 10, minor 8)
    half = [lowest(r) for r in rows[::2]], [lowest(r) for r in rows[1::2]]          # the level found on one half of the draws, tried on the other half
    out["share_at_that_level_on_other_draws"] = share([v < float(np.quantile(half[0], LEVEL)) for v in half[1]] + [v < float(np.quantile(half[1], LEVEL)) for v in half[0]])
    line = float(chi2_law.isf(LEVEL, 1))
    for name in rows[0]["named"]:
        gains = np.array([r["again"][name]["gain"] for r in rows])
        levels = [r["named"][name]["level"] if r["named"][name] is not None else 1.0 for r in rows]          # a term that the fit takes whole has no test: it never reports
        samples = rows[0]["samples"]
        out["terms"][name] = {"test_of_the_term": share(level < LEVEL for level in levels), "level_of_the_test_at_the_5th_centile": float(np.quantile(levels, LEVEL)),
                              "fit_again": share(gains > line), "gain_at_the_95th_centile": float(np.quantile(gains, 1.0 - LEVEL)),
                              "AIC": share(criteria.criteria(g, 1, samples, 8, samples - rows[0]["dof"])["AIC"] > 0.0 for g in gains),
                              "BIC": share(criteria.criteria(g, 1, samples, 8, samples - rows[0]["dof"])["BIC_n_the_samples"] > 0.0 for g in gains)}
    return out


def of_the_lacking(rows: list, null: dict, term: list, true: dict) -> dict:
    """What the draws of a model that lacks a term say, beside the null of the same model on its own traces. How often
    the misfit is reported: by the chi square of the fit and by the usual checks at `LEVEL`, by each family and by any
    of them at the level for all together, by the test of the named term at `LEVEL` and at its own level of the null,
    by the fit made again with it, by AIC and by BIC. Then which term is named: how often the missing one has the
    largest score among the terms of the toy, and the largest gain when each is fitted in turn, which is what a choice
    by AIC or BIC names first. The chi square of the fit and the usual checks are counted at `LEVEL` and at a level of
    their own on the right model too, as the readings are. And the guess, on the draws on which it names anything: how often the missing term is
    the first that the package names, how often it stands with the first, and the kind of the first. A draw on which
    no term leads is no hit, whatever ranks first there (review 9, moderate 5). Then the size, with its bias and the
    share of its intervals that hold the truth: the fitted one of the fit made again, where its search ended and its
    parameters are told apart, and the first order one of the test, where the fit leaves `analysis.HIDDEN` of the term
    to be seen (review 9, moderate 4). And `others`, the control of the guess: for each term of the toy that the model
    lacks and the traces do not hold, how often it is the first guess and how often it stands with the first."""
    together, line, samples = null["level_for_all_together"], float(chi2_law.isf(LEVEL, 1)), rows[0]["samples"]
    names = lambda r: r["guess"]["leading"] > 0          # noqa: E731  the guess names something on this draw
    out = {"draws": len(rows), "level_for_all_together": together, "told": share(r["told"] for r in rows),
           "chi_square_of_the_fit": share(r["readings"][FIT] < LEVEL for r in rows),
           "chi_square_of_the_fit_at_its_level": share(r["readings"][FIT] < null["level_of_the_fit_at_the_5th_centile"] for r in rows),
           "usual": {name: share(r["usual"][name] < LEVEL for r in rows) for name in rows[0]["usual"]},
           "any_usual_at_their_level": share(min(r["usual"].values()) < null["level_for_the_usual_together"] for r in rows),
           "readings_on_their_own": {name: share(r["readings"][name] < LEVEL for r in rows) for name in rows[0]["readings"]},
           "readings_at_the_level_for_all": {name: share(r["readings"][name] < together for r in rows) for name in rows[0]["readings"]},
           "any_reading_at_the_level_for_all": share(lowest(r) < together for r in rows),
           "guess_names_a_term": share(names(r) for r in rows), "terms": {}}
    kinds = {}
    for r in rows:
        kind = str(r["guess"]["first_kind"]) if names(r) else NOTHING
        kinds[kind] = kinds.get(kind, 0) + 1
    out["kind_of_the_first"] = {k: v / len(rows) for k, v in sorted(kinds.items(), key=lambda kv: -kv[1])}
    for name in term:
        gains = np.array([r["again"][name]["gain"] for r in rows])
        levels = [r["named"][name]["level"] if r["named"][name] is not None else 1.0 for r in rows]
        mine = lambda r: r["guess"]["terms"].get(name, {})          # noqa: E731
        fitted = [(r["again"][name]["value"], r["again"][name]["bar"]) for r in rows if r["again"][name]["ended"] and not r["again"][name]["tied"]]
        first = [(r["named"][name]["size"], r["named"][name]["bar"]) for r in rows if r["named"][name] is not None and r["named"][name]["size"] is not None and r["named"][name]["visible"] >= analysis.HIDDEN]
        sizes = {}
        for label, pairs in (("fit_again", fitted), ("test_of_the_term", first)):
            v, b = np.array([p[0] for p in pairs], dtype=float), np.array([p[1] for p in pairs], dtype=float)
            sizes[label] = {"of": len(pairs), "mean": float(v.mean()) if len(v) else None, "spread": float(v.std(ddof=1)) if len(v) > 2 else None,
                            "holds_the_truth": float(np.mean(np.abs(v - true[name]) <= b)) if len(v) else None}
        out["terms"][name] = {"true": true[name], "test_of_the_term": share(level < LEVEL for level in levels),
                              "test_of_the_term_at_the_level_of_the_null": share(level < null["terms"][name]["level_of_the_test_at_the_5th_centile"] for level in levels),
                              "fit_again": share(gains > line), "fit_again_at_the_centile_of_the_null": share(gains > null["terms"][name]["gain_at_the_95th_centile"]),
                              "AIC": share(criteria.criteria(g, 1, samples, 8, samples - rows[0]["dof"])["AIC"] > 0.0 for g in gains),
                              "BIC": share(criteria.criteria(g, 1, samples, 8, samples - rows[0]["dof"])["BIC_n_the_samples"] > 0.0 for g in gains),
                              "largest_score_among_the_terms": share(max((n for n in r["named"] if r["named"][n] is not None), key=lambda n: r["named"][n]["score"], default=None) == name for r in rows),
                              "largest_gain_among_the_terms": share(max(r["again"], key=lambda n: r["again"][n]["gain"]) == name for r in rows),
                              "ranks_first": share(mine(r).get("rank") == 1 for r in rows),
                              "is_the_first": share(names(r) and mine(r).get("rank") == 1 and mine(r).get("leads", False) for r in rows),
                              "stands_with_the_first": share(names(r) and mine(r).get("with_first", False) for r in rows), "size": sizes}
    of = lambda r, name: r["guess"]["terms"].get(name, {})          # noqa: E731
    out["others"] = {name: {"is_the_first": share(names(r) and of(r, name).get("rank") == 1 and of(r, name).get("leads", False) for r in rows),
                            "stands_with_the_first": share(names(r) and of(r, name).get("with_first", False) for r in rows)}
                     for name in rows[0]["again"] if name not in term}          # the control: the terms of the toy that the model lacks and the traces do not hold. Where one of them stands with the first guess as often as the missing term, the guess names a class and not a term (review 10, moderate 2)
    return out


ESTIMATES = ("fit", "generic", "first_order", "free_term", "free")          # the five estimates of a known parameter that a row of the unsuspected case keeps, by their names in the row: the fit, with a generic term fitted and at first order, and the combination of the families that is free of that term, or of the leading class


def of_the_known(rows: list, true: dict) -> dict:
    """What each estimate of each known parameter does over the draws of a case, against its true value: its bias in
    spreads of its own estimates and in its own bars, the bias in its own bar on the median draw, which a few wild draws
    do not drive as they drive the mean (review 15), and the share of the draws on which its interval of one bar holds
    the truth, where 0.68 is right. Where the package gives no such estimate on a draw, since nothing was reported or
    no combination was of use, a researcher keeps the fit: so the fit stands in on that draw, and `given` says on how
    many draws the estimate itself was made. `where_given` counts those draws alone."""
    out = {}
    for name in rows[0]["known"]:
        out[name] = {"true": float(true[name]), "estimates": {}}
        for which in ESTIMATES:
            given = [r["known"][name][which] for r in rows]
            used = np.array([g if g is not None else r["known"][name]["fit"] for g, r in zip(given, rows)], dtype=float)
            off, spread = used[:, 0] - float(true[name]), float(used[:, 0].std(ddof=1)) if len(used) > 1 else float("nan")
            own = np.array([g for g in given if g is not None], dtype=float).reshape(-1, 2)
            out[name]["estimates"][which] = {
                "given": share(g is not None for g in given), "mean": float(used[:, 0].mean()), "spread": spread, "mean_bar": float(used[:, 1].mean()),
                "bias_in_spreads": float(off.mean() / spread) if spread > 0 else float("nan"), "bias_in_bars": float(np.mean(off / used[:, 1])),
                "median_bias_in_bars": float(np.median(off / used[:, 1])), "holds_the_truth": share(np.abs(off) <= used[:, 1]),
                "where_given": {"bias_in_bars": float(np.mean((own[:, 0] - float(true[name])) / own[:, 1])) if len(own) else None, "holds_the_truth": share(np.abs(own[:, 0] - float(true[name])) <= own[:, 1])}}
    return out


def of_the_right(rows: list, true: dict) -> dict:
    """What the draws of a right model say in the count where nobody names a term. How often the package says that
    something is missing, at the one level that it counted itself for the case: near `LEVEL` when that level is
    right. How often any reading reports at `LEVEL` on its own, and the level for all that these draws give, to set
    beside the one of the package. The levels at which the opponents that need no name report on `LEVEL` of these
    draws: the chi square of the fit, and the three usual checks together. And what each estimate of the known
    parameters does where nothing is missing: what the package costs when it is not needed."""
    return {"draws": len(rows), "told": share(r["told"] for r in rows), "missing": share(r["missing"] for r in rows),
            "any_reading": share(lowest(r) < LEVEL for r in rows), "level_for_all_on_these_draws": float(np.quantile([lowest(r) for r in rows], LEVEL)),
            "chi_square_of_the_fit": share(r["readings"][FIT] < LEVEL for r in rows), "level_of_the_fit_at_the_5th_centile": float(np.quantile([r["readings"][FIT] for r in rows], LEVEL)),
            "any_usual": share(min(r["usual"].values()) < LEVEL for r in rows), "level_for_the_usual_together": float(np.quantile([min(r["usual"].values()) for r in rows], LEVEL)),
            "a_class_leads": share(r["kind"] is not None for r in rows), "a_term_is_fitted": share(r["term"] is not None for r in rows), "known": of_the_known(rows, true)}


def of_the_hidden(rows: list, right_: dict, true: dict) -> dict:
    """What the draws of a case say when its term is hidden and nobody names it. How often the package says that
    something is missing, beside the two opponents that need no name either, each at the level at which it reports
    the right model on `LEVEL` of its draws. Of which kind: how often a class of generic terms leads, the families of
    that class, and the generic term of the largest gain. What each estimate of the known parameters does. And on how
    many draws first order gives the same value as the fit made again, within `U.SAME` of the bar of that fit: what
    the licence of the package says beforehand, counted afterwards."""
    kinds, terms = {}, {}
    for r in rows:
        kind = "no class leads" if r["kind"] is None else ", ".join(r["kind"]["families"])
        term = "no term is fitted" if r["term"] is None else r["term"]["name"]
        kinds[kind], terms[term] = kinds.get(kind, 0) + 1, terms.get(term, 0) + 1
    same = {}
    for name in rows[0]["known"]:
        both = [(r["known"][name]["first_order"], r["known"][name]["generic"]) for r in rows if r["known"][name]["first_order"] is not None and r["known"][name]["generic"] is not None]
        same[name] = share(abs(first[0] - again_[0]) <= U.SAME * again_[1] for first, again_ in both)
    return {"draws": len(rows), "told": share(r["told"] for r in rows), "missing": share(r["missing"] for r in rows),
            "chi_square_of_the_fit_at_its_level": share(r["readings"][FIT] < right_["level_of_the_fit_at_the_5th_centile"] for r in rows),
            "any_usual_at_their_level": share(min(r["usual"].values()) < right_["level_for_the_usual_together"] for r in rows),
            "a_class_leads": share(r["kind"] is not None for r in rows), "families_of_the_leading_class": {k: v / len(rows) for k, v in sorted(kinds.items(), key=lambda item: -item[1])},
            "term_of_the_largest_gain": {k: v / len(rows) for k, v in sorted(terms.items(), key=lambda item: -item[1])},
            "known": of_the_known(rows, true), "first_order_gives_the_same_value": same}


WINDOWS = (1.0, 0.5, 0.25)          # the parts of the window of a trace that are kept about the centre of the line, for the count of a window that cuts the curve
MOMENTS = "moments: powers 0 to 4"          # the family of the package whose own estimate is set beside the usual moments
EDGE = 0.05          # the share of the samples at each end of a cut trace that an analyst takes for the baseline
CUT_DRAWS = 30          # 30, the first stage of row M474, until then 200: the draws of the right model on which the cut window is counted
CUTS = ("cut", "usual_widths", "cut_window_row", "WINDOWS", "MOMENTS", "EDGE")          # what a row of the cut window is made with, in this program


def cut(rows: list, part: float, const: dict) -> list:
    """The traces with their window cut to a part of the one that was recorded, about the centre of the line: every
    sample within that part of the half window, and each array of a trace cut alike."""
    half = float(part) * float(const["axis"]["half_window_MHz"])
    out = []
    for one_row in rows:
        keep = np.abs(np.asarray(one_row["nu"], dtype=float)) <= half
        out.append({k: (np.asarray(v)[keep] if isinstance(v, np.ndarray) and v.shape == keep.shape else v) for k, v in one_row.items()})
    return out


def usual_widths(one_row: dict) -> dict:
    """The full width at half maximum of the line of one trace by two usual rules of moments, each of which takes the
    window for the whole curve (his message 133: "when the current methods are affected by window truncation
    approximations"). The baseline is the mean of the samples at the two ends, as an analyst reads it. By the area:
    the area of a Lorentzian is pi / 2 times its height times its width. By the second moment: the width of a
    Gaussian is 2 sqrt(2 ln 2) times the root of its variance. A window that cuts the wings takes area from the
    first, sets the baseline too high for both, and bounds a variance that has no bound on a Lorentzian."""
    nu, y = np.asarray(one_row["nu"], dtype=float), np.asarray(one_row["y"], dtype=float)
    ends = max(int(EDGE * nu.size), 2)
    z = y - float(np.mean(np.concatenate([y[:ends], y[-ends:]])))
    area = float(np.trapezoid(z, nu))
    top = float(np.convolve(z, np.ones(9) / 9.0, mode="same").max())          # the top as it is read on a trace with noise: a mean over nine samples
    centre = float(np.trapezoid(nu * z, nu)) / area
    variance = float(np.trapezoid((nu - centre) ** 2 * z, nu)) / area
    return {"by_area": 2.0 * area / (math.pi * top), "by_variance": 2.0 * math.sqrt(2.0 * math.log(2.0)) * math.sqrt(max(variance, 0.0))}


def cut_window_row(draw: int, const: dict, config: dict) -> dict:
    """One draw of the count of a window that cuts the curve: the 8 synthetic traces of the first detail, a Lorentzian
    of the natural width, under their own model, with the window whole and cut to each part of `WINDOWS`. At each,
    the natural width by four estimates: the fit, the family of moments of the package, which sets the moments of the
    cut traces against the moments of the model over the same window, and the two usual rules, as their mean over the
    8 traces with the standard error of that mean."""
    out = {}
    for part in WINDOWS:
        some = cut(fits.synthetic(1, const, config, draw), part, const)
        problem, result, lin = analysis.at(1, some, const, config)
        got = families.family(lin, MOMENTS, analysis.CENTRE).estimate()
        usual = [usual_widths(one_row) for one_row in some]
        mean_of = lambda key: [said(np.mean([u[key] for u in usual])), said(np.std([u[key] for u in usual], ddof=1) / math.sqrt(len(usual)))]          # noqa: E731
        out[f"{part:g}"] = {"fit": [said(result.shared[0]), said(result.sd[0])], "package": [said(got.shared[0]), said(got.sd[0])] if got.converged else None,
                            "by_area": mean_of("by_area"), "by_variance": mean_of("by_variance")}
    return out


def cut_one(job: tuple) -> tuple:
    """One row of the cut window, for a pool of processes. A draw on which it cannot be made is kept as a row that says why."""
    draw, const, config = job
    try:
        return draw, cut_window_row(draw, const, config)
    except (Exception, SystemExit) as why:
        return draw, {"failed": f"{type(why).__name__}: {why}"[:300]}


def cut_kept(root: Path) -> Path:
    """The file of the rows of the cut window, under what they depend on: what a kept fit of the first detail depends
    on, the files of the package that the analysis runs, and the functions of this program that make such a row."""
    root = Path(root)
    h = hashlib.sha256(draws.digest(root, 1).encode())
    h.update(sources((root / "toy" / "analysis.py").read_text(), ("at", "CENTRE")))
    h.update(sources((root / "toy" / "power.py").read_text(), CUTS))
    draws.hashed(h, root, draws.package_of(root, RUNS))
    return root / "runs" / f"power_cut_window_{h.hexdigest()[:12]}.jsonl"


def run_cut(n: int = CUT_DRAWS, root: Path = ROOT, workers: int = WORKERS) -> list:
    """The rows of the cut window, `n` draws, side by side, each kept as it is made. A row that was kept is not made again."""
    from concurrent.futures import ProcessPoolExecutor, as_completed
    root = Path(root)
    const, config = model.load(root / "toy", root)
    path = cut_kept(root)
    path.parent.mkdir(exist_ok=True)
    got = read(path)

    def keep(made: tuple) -> None:
        got[made[0]] = dict(made[1], draw=made[0])
        with path.open("a") as fh:
            fh.write(json.dumps(got[made[0]]) + "\n")
    jobs = [(1 + k, const, config) for k in range(n) if 1 + k not in got]
    if workers < 1:
        for job in jobs:
            keep(cut_one(job))
    else:
        with ProcessPoolExecutor(min(workers, max(len(jobs), 1))) as pool:
            for done in as_completed([pool.submit(cut_one, job) for job in jobs]):
                keep(done.result())
    return [got[k] for k in sorted(got) if k <= n]


def of_the_cut_window(rows: list, true: float) -> dict:
    """What the draws say of each estimate of the width at each window: its mean, the share of the true width that it
    is off by, its bias in its own bars, and the share of the draws on which its interval of one bar holds the true
    width, where 0.68 is right. The bar of a usual rule is the standard error of its mean over the 8 traces."""
    out = {}
    for part in rows[0]:
        if part == "draw":
            continue
        out[part] = {}
        for which in ("fit", "package", "by_area", "by_variance"):
            got = np.array([r[part][which] for r in rows if r[part][which] is not None], dtype=float).reshape(-1, 2)
            off = got[:, 0] - float(true)
            out[part][which] = {"given": len(got), "mean": float(got[:, 0].mean()), "off_by_a_share_of_the_true_width": float(off.mean() / float(true)), "spread": float(got[:, 0].std(ddof=1)),
                                "mean_bar": float(got[:, 1].mean()), "bias_in_bars": float(np.mean(off / got[:, 1])), "holds_the_truth": share(np.abs(off) <= got[:, 1])}
    return out


def could_hide(d: int, free: list | None, made: int, const: dict, config: dict, level: float) -> dict:
    """What the 8 traces could hide of the term of detail `made` under the model that lacks it (row M68 of the plan),
    in the linear picture of that model at its true values, with no noise drawn. The best that any reading of these
    traces can do is the statistic that is aimed at the term itself, which only one who names the term can make: the
    size at which that statistic sees the term on half of the sets, at `level` and at `LEVEL`. Under that size the
    records hide the term from every reading. And what each known parameter is off by, in its bar, at that size and
    at the size of the twin: what the model absorbed of the term, which no residual shows. The size at which the
    statistic sees the term is first order, no draw: the counts of the hidden cases stand beside it. What the fit
    leaves to be seen and what each known parameter is off by at the size of the twin are given by the fit made again,
    with no noise drawn, and at first order beside it: the counts refuse first order its licence for the laser noise
    and the transit time (CLAUDE.md, rule 5b, his message 135)."""
    from toy import ceiling, traces
    from twinfit.ledger import size_at
    problem, own = ceiling.at_truth(d, traces.conditions(const), const, config, free=free)
    result = twinfit.fit(problem, shared0=problem.start_shared, own0=own)
    lin = twinfit.Linear(problem, result)
    name, added, cand = next(t for t in analysis.of_the_terms(problem, result, const, config, small=True, detail=d) if config["parameter"][t[0]]["enters"] == made)
    entry, column = config["parameter"][name], cand.columns[0]
    seen = lin.visible(column)
    absorbed = np.asarray(lin.absorbed(column), dtype=float)
    twin = (float(entry["true"]) / added) ** 2 if entry.get("effect") == "square" else float(entry["true"]) / added          # the amplitude of the column at the size of the twin

    def at(alpha: float) -> dict:
        amplitude = size_at(float(seen @ seen), 1, alpha, 0.5)
        size = read_size(amplitude, 0.0, added, entry)["size"] if math.isfinite(amplitude) else None
        return {"level": float(alpha), "size": size, "over_the_size_of_the_twin": None if size is None else said(size / float(entry["true"])),
                "parameters_off_in_bars": {p: said(absorbed[j] * amplitude / result.sd[j]) for j, p in enumerate(problem.shared)} if math.isfinite(amplitude) else None}
    with_it, own_with = ceiling.at_truth(d, traces.conditions(const), const, config, traces_of=made, free=free)          # traces with no noise drawn that hold the term at the size of the twin
    again = twinfit.fit(with_it, shared0=with_it.start_shared, own0=own_with)          # the model that lacks the term, fitted again to them
    norm = sum(float(np.sum(((a.y - b.y) / a.sd) ** 2)) for a, b in zip(with_it.records, problem.records))          # the whole term, in the noise of the traces that hold it
    return {"term": name, "unit": entry["unit"], "size_of_the_twin": float(entry["true"]),
            "share_that_the_fit_leaves_visible": said(math.sqrt(float(again.chi2) / norm)) if norm > 0 else None,          # of the norm, as the package's visibility is
            "share_that_the_fit_leaves_visible_at_first_order": said(lin.visibility(column)),
            "seen_on_half_of_the_sets_by_the_statistic_aimed_at_it": {"at_the_level_for_all": at(level), "at_the_level_of_one_test": at(LEVEL)},
            "parameters_off_in_bars_at_the_size_of_the_twin": {p: said((float(again.shared[j]) - float(problem.start_shared[j])) / result.sd[j]) for j, p in enumerate(problem.shared)},
            "parameters_off_in_bars_at_the_size_of_the_twin_at_first_order": {p: said(absorbed[j] * twin / result.sd[j]) for j, p in enumerate(problem.shared)}}


def counted(root: Path = ROOT, details: list | None = None) -> dict:
    """The table of the counts, from the rows that were kept: for each right model its null, and for each case in which
    a model lacks a term what the draws say beside the null of the same model. `details`: the details whose cases are
    counted, all of them when None: the null of every right model with rows is counted whatever they are, since a case
    stands beside the null of the detail before its own (review 9, minor 12). A case with no kept row is left out. A
    draw on which no row could be made is in the rows and in no share: `failed` says how many there are."""
    root = Path(root)
    const, config = model.load(root / "toy", root)
    out = {"level": LEVEL, "null": {}, "lacking": {}}
    made_from = digest(root)

    def rows_of(kind: str, name: str) -> tuple:
        rows = [r for _, r in sorted(read(kept(root, kind, name)).items())]
        return rows, [r for r in rows if "failed" not in r]
    for name, made, d, free, _ in right(config):
        rows, good = rows_of("null", name)
        if good:
            out["null"][name] = dict(of_the_null(good), failed=len(rows) - len(good), detail=d, free=free, made_from=made_from, rows=rows)
    for name, made, d, free, sizes in lacking(config, root):
        rows, good = rows_of("lacking", name)
        null = out["null"].get(draws.key(d, free, config))
        if good and null and (details is None or made in details):
            term = [n for n, e in config["parameter"].items() if e["enters"] == made]
            true = dict(model.truth(config), **(sizes or {}))
            out["lacking"][name] = dict(of_the_lacking(good, null, term, true), failed=len(rows) - len(good), traces_of_detail=made, detail=d, free=free, sizes=sizes, made_from=made_from, rows=rows)
    out["map"] = {name: map_of(d, free, const, config) for name, made, d, free, _ in right(config) if name in out["null"]}          # for each right model that was counted, what the 8 traces can tell apart under it
    out["unsuspected"] = {}
    cases = hidden(config)
    for case in sorted(cases, key=lambda c: c[1] != c[2]):          # the right models first: a case stands beside the right model of its own
        name, made, d, free, sizes = case
        rows, good = rows_of("unsuspected", name)
        beside = out["unsuspected"].get(f"{draws.key(d, free, config)}, right")
        if not good or (made != d and beside is None) or not (details is None or made in details):
            continue
        true = dict(model.truth(config), **(sizes or {}))
        once = level_of(case, const, config, root)
        entry = of_the_right(good, true) if made == d else of_the_hidden(good, beside, true)
        out["unsuspected"][name] = dict(entry, failed=len(rows) - len(good), traces_of_detail=made, detail=d, free=free, sizes=sizes, level_of_the_package={k: v for k, v in once["level"].items() if k != "readings"},
                                        licence=once["licence"], no_licence=once["no_licence"], made_from=made_from, rows=rows)
    out["unsuspected"] = {c[0]: out["unsuspected"][c[0]] for c in cases if c[0] in out["unsuspected"]}          # in the order of the cases
    for name, made, d, free, sizes in cases:          # beside each case at the size of the twin: what the traces could hide of its term, at first order
        if made != d and sizes is None and name in out["unsuspected"]:
            out["unsuspected"][name]["could_hide"] = could_hide(d, free, made, const, config, out["unsuspected"][name]["level_of_the_package"]["level_for_all"])
    window = read(cut_kept(root))
    good = [r for _, r in sorted(window.items()) if "failed" not in r]
    if good:
        out["cut_window"] = dict(windows=of_the_cut_window(good, model.truth(config)["natural_width"]), draws=len(good), failed=len(window) - len(good), true=model.truth(config)["natural_width"],
                                 half_window_MHz=const["axis"]["half_window_MHz"], rows=[r for _, r in sorted(window.items())])
    return out


def text_of(table: dict) -> str:
    """The table of the counts as the text of its file: the counts on lines of their own, which a reader can open, and
    each row of a draw on one line, so that the file of some thousands of rows stays small and a later count changes
    only the lines of its own rows."""
    light, rows = json.loads(json.dumps(table)), {}
    for kind in ("null", "lacking", "unsuspected"):
        for name, entry in light.get(kind, {}).items():
            mark = f"@ the rows of {kind} {name} @"
            rows[json.dumps(mark)] = "[\n" + ",\n".join(json.dumps(r, separators=(",", ":")) for r in entry["rows"]) + "\n]"
            entry["rows"] = mark
    if "cut_window" in light:
        mark = "@ the rows of the cut window @"
        rows[json.dumps(mark)] = "[\n" + ",\n".join(json.dumps(r, separators=(",", ":")) for r in light["cut_window"]["rows"]) + "\n]"
        light["cut_window"]["rows"] = mark
    text = json.dumps(light, indent=1)
    for mark, lines in rows.items():
        text = text.replace(mark, lines)
    return text + "\n"


def queue_of(root: Path, program: str, made_from: str) -> dict:
    """The waves that the queue spent on the rows of a count, which its file holds beside them: `spent` of
    tools/autoqueue.py, loaded by its path. No check of a count reads it, since the ledger of a copy is not where the
    count is checked."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("toy_autoqueue", ROOT / "tools" / "autoqueue.py")
    queue = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(queue)
    return queue.spent(Path(root), program, made_from)


def write(root: Path = ROOT, details: list | None = None) -> Path:
    """Writes toy/power.json: whole, or the one before, with the waves that the queue spent on its rows."""
    path = Path(root) / "toy" / "power.json"
    part = path.with_suffix(".json.part")
    part.write_text(text_of(dict(counted(root, details), package=twinfit.__version__, queue=queue_of(root, "toy/power.py", digest(root)))))          # the version under the table, which tools/met.py reads against the release
    os.replace(part, path)
    return path


def alike(a, b) -> bool:
    """Whether two counts are the same count: every key and every item, a number within a part in a thousand million,
    and two numbers that are no number alike."""
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(alike(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(alike(x, y) for x, y in zip(a, b))
    if isinstance(a, bool) or isinstance(b, bool) or not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
        return a == b
    return (math.isnan(a) and math.isnan(b)) or abs(a - b) <= 1e-9 * max(1.0, abs(a), abs(b))


def check(root: Path = ROOT) -> list:
    """What is wrong with toy/power.json: an empty list when every count that it holds was made from the files as they
    are and comes out the same, whole, when its rows are counted again (review 9, moderate 6: three shares of a table
    were counted again, and a hand could write the rest). A table whose shares a hand wrote differs here, and so does
    one that holds a draw twice or says another number of draws than it holds. A model or a case that the files do not
    give is refused. A model or a case that is not there is not: a table may hold fewer tiers than the files give, and
    the tests of the counts hold which ones it has to hold (review 11, minor 8)."""
    root = Path(root)
    path = root / "toy" / "power.json"
    if not path.exists():
        return ["toy/power.json is not there"]
    held, wrong = json.loads(path.read_text()), []
    _, config = model.load(root / "toy", root)
    truth, made_from = model.truth(config), digest(root)

    def rows_of(name: str, entry: dict) -> list | None:
        good = [r for r in entry["rows"] if "failed" not in r]
        if len({r["draw"] for r in entry["rows"]}) != len(entry["rows"]) or entry["draws"] != len(good) or entry.get("failed") != len(entry["rows"]) - len(good):
            wrong.append(f"{name} holds a draw twice, or says another number of draws than it holds")
            return None
        return good
    const, _ = model.load(root / "toy", root)
    models, cases = {c[0]: c for c in right(config)}, {c[0] for c in lacking(config, root)}
    if held.get("level") != LEVEL or not set(held["null"]) <= set(models) or not set(held["lacking"]) <= cases:
        wrong.append("the table holds another level than the program, or a model or a case that the files do not give")          # review 10, minor 9
    for name in sorted(set(held.get("map", {})) & set(models)):
        if not alike(held["map"][name], json.loads(json.dumps(map_of(models[name][2], models[name][3], const, config)))):
            wrong.append(f"the map of the model {name} is not the one that the files give")
    for name, entry in held["null"].items():
        if entry["made_from"] != made_from:
            wrong.append(f"the null of the model {name} was made from other files than those on disk")
        good = rows_of(f"the null of the model {name}", entry)
        again_ = of_the_null(good) if good else None
        if again_ is not None and not alike(again_, {k: entry.get(k) for k in again_}):
            wrong.append(f"the null of the model {name} does not give its shares when its rows are counted again")
    if sorted(held.get("map", {})) != sorted(held["null"]):
        wrong.append("the map of what the traces can tell apart is not there for every right model that was counted")
    for name, entry in held["lacking"].items():
        if entry["made_from"] != made_from:
            wrong.append(f"the case {name} was made from other files than those on disk")
        null = held["null"].get(draws.key(entry["detail"], entry["free"], config))
        if null is None:
            wrong.append(f"the case {name} has no null of its model")
            continue
        good = rows_of(f"the case {name}", entry)
        term = [n for n, e in config["parameter"].items() if e["enters"] == entry["traces_of_detail"]]
        again_ = of_the_lacking(good, null, term, dict(truth, **(entry["sizes"] or {}))) if good else None
        if again_ is not None and not alike(again_, {k: entry.get(k) for k in again_}):
            wrong.append(f"the case {name} does not give its shares when its rows are counted again")
    mine = held.get("unsuspected", {})
    if not set(mine) <= {c[0] for c in hidden(config)}:
        wrong.append("the table holds a case of a hidden term that the files do not give")
    for name, entry in mine.items():
        if entry["made_from"] != made_from:
            wrong.append(f"the case {name}, where nobody names the term, was made from other files than those on disk")
        good = rows_of(f"the case {name}", entry)
        true = dict(truth, **(entry["sizes"] or {}))
        beside = mine.get(f"{draws.key(entry['detail'], entry['free'], config)}, right")
        if entry["traces_of_detail"] != entry["detail"] and beside is None:
            wrong.append(f"the case {name} has no right model beside it")
            continue
        again_ = (of_the_right(good, true) if entry["traces_of_detail"] == entry["detail"] else of_the_hidden(good, beside, true)) if good else None
        if again_ is not None and not alike(json.loads(json.dumps(again_)), {k: entry.get(k) for k in again_}):
            wrong.append(f"the case {name} does not give its shares when its rows are counted again")
        if "could_hide" in entry and not alike(entry["could_hide"], json.loads(json.dumps(could_hide(entry["detail"], entry["free"], entry["traces_of_detail"], const, config, entry["level_of_the_package"]["level_for_all"])))):
            wrong.append(f"what the traces could hide in the case {name} is not what the files give")
    window = held.get("cut_window")
    if window is not None:
        good = [r for r in window["rows"] if "failed" not in r]
        if len({r["draw"] for r in window["rows"]}) != len(window["rows"]) or window["draws"] != len(good) or window["true"] != truth["natural_width"]:
            wrong.append("the cut window holds a draw twice, another number of draws than it says, or another true width than the files")
        elif not alike(json.loads(json.dumps(of_the_cut_window(good, truth["natural_width"]))), window["windows"]):
            wrong.append("the cut window does not give its numbers when its rows are counted again")
    return wrong


ALIKE = 0.9          # two terms whose visible parts have a cosine beyond this are alike, as in the table of the package


def map_of(d: int, free: list | None, const: dict, config: dict) -> dict:
    """The map of what 8 traces can tell apart under one model: at its true values with no noise drawn, the cosine
    between what the fit leaves visible of two terms, for each pair of terms of the toy that the model lacks, and for
    each such term the generic candidate of the package that is nearest to it, with those that are alike it. At 1 two
    terms are one to these traces, and the residuals cannot say which is missing: the guess then names a class. A term
    of which the fit leaves less than `analysis.HIDDEN` to be seen is taken by the model: it has no shape of its own
    left, and it is listed apart. A candidate with several columns is left out, since it has no single shape."""
    from itertools import combinations
    from toy import ceiling, traces
    problem, own = ceiling.at_truth(d, traces.conditions(const), const, config, free=free)
    result = twinfit.fit(problem, shared0=problem.start_shared, own0=own)
    lin = twinfit.Linear(problem, result)
    terms = analysis.of_the_terms(problem, result, const, config, small=True, detail=d)
    units, seen = {}, {}
    from toy import bench          # the toy's own list of named shapes: the map asks which of them each term of the toy resembles, a question of the toy and not of the package, which holds no list (his message 340)
    from twinfit.candidates import covariates_of
    for name, cand, mine in [(name, cand, True) for name, _, cand in terms] + [(c.name, c, False) for c in bench.candidates(lin, covariates=covariates_of(problem))]:
        if len(cand.columns) != 1:
            continue
        v = lin.visible(cand.columns[0])
        size = float(np.sqrt(v @ v))
        if mine:
            seen[name] = said(lin.visibility(cand.columns[0]))
        if size > 0.0 and (not mine or seen[name] >= analysis.HIDDEN):
            units[name] = v / size
    named = [name for name in seen if name in units]
    generic = [name for name in units if name not in seen]
    return {"detail": d, "free": list(problem.shared), "visible_share": seen, "taken_by_the_model": [name for name in seen if name not in units], "terms": named,
            "pairs": {f"{a} and {b}": said(abs(float(units[a] @ units[b]))) for a, b in combinations(named, 2)},
            "alike": {a: sorted(name for name in generic if abs(float(units[a] @ units[name])) >= ALIKE) for a in named},
            "nearest": {a: max(((name, said(abs(float(units[a] @ units[name])))) for name in generic), key=lambda pair: pair[1], default=None) for a in named}}


def main(argv: list) -> int:
    if argv == ["units"]:
        print(json.dumps(units(ROOT)))
        return 0
    if argv[:1] == ["slice"] and len(argv) >= 4:
        if not os.environ.get("TOY_QUEUE"):
            print("REFUSED: a slice is a wave of the queue: `python3 tools/autoqueue.py add toy/power.py` (his message 147)", file=sys.stderr)
            return 1
        print(json.dumps(slice_of(argv[1], int(argv[2]), int(argv[3]), float(argv[5]) if argv[4:5] == ["--for"] else float("inf"), ROOT)))
        return 0
    if argv == ["progress"]:
        print(progress(ROOT))
        return 0
    if argv == ["digest"]:          # what the queue judges a wave by: a change that moves no unit makes none stale
        print(digest(ROOT))
        return 0
    if argv[:1] in (["null"], ["lacking"], ["unsuspected"], ["cut_window"]):
        print("REFUSED: the counts run through the queue: `python3 tools/autoqueue.py add toy/power.py` (his message 147)", file=sys.stderr)
        return 1
    if argv == ["table"]:
        print(write(ROOT))
        return 0
    if argv == ["--check"]:
        wrong = check(ROOT)
        for line in wrong:
            print("REFUSED: " + line)
        print("the counts on real draws: " + ("every one was made from the files as they are and counts again to the same shares" if not wrong else f"{len(wrong)} thing(s) wrong"))
        return 1 if wrong else 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
