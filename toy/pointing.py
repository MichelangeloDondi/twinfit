"""Where the report of the package points, counted on real draws of the toy (row M86 of the plan, his messages 138 and 139).

That a term is missing says nothing in a real case (message 138): the point is to identify it. Here each term of details
2 to 5 is hidden in turn from the model that lacks it, at the size of the twin and at a half of it, and the report of
where to search (twinfit/where.py) is read on each draw, with no name of a term of the toy in it. What the physics of
each term says, written in `EXPECTED` before any count, is set beside where the report points: the place, and the
covariate along which the term changes. On a right model the report should point nowhere, and the share of its draws on
which some answer stands beyond the twin is its rate of false alarms.

    python3 tools/autoqueue.py add toy/pointing.py
                                        the count, through the queue (his message 147): the twin of each case in `PARTS`
                                        units, then `DRAWS` draws of each case whose term is hidden and twice as many of
                                        each right model, each draw a unit read against the twin of its case. The queue
                                        asks `units`, runs `slice` and reads `progress`, and a slice outside it is refused
    python toy/pointing.py table        counts the kept rows and writes toy/pointing.json
    python toy/pointing.py --check      ends without error when toy/pointing.json holds counts of rows made from the
                                        files as they are, which give the same shares when they are counted again

The twin of a case is read once, on `SETS` synthetic sets of the fit of its first draw, and every draw of the case is
read against it. A researcher reads the twin of the fit of his own records, which differs from draw to draw: row M104
counts the right models against the twin of each draw (review 21, finding 9). A row is kept in runs/ as it is made. This
program is apart from toy/power.py so that its rows and those of the counts there never make each other stale.
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

from twinfit import twin, where  # noqa: E402

from toy import analysis, draws, fits, model, power  # noqa: E402

DRAWS = 30          # 30 a hidden case and 60 a right model, the first stage of rows M474 and M481, 100 where a share lies within two standard errors of its threshold. Until then 100: 100 a hidden case and 200 a right model at 0.1.0rc6, sized by the precision of a share within about 0.10 in 95 cases of 100, as row M204 asks, where its count at 0.0.8 was stale against the package he will receive (the red team of his message 303, severe 1). 400 since the slap of 2026-10-17 14:09 until then
SECOND, LINE = 100, 0.60          # the second stage of the rule above: a hidden case whose share of the place or of its setting lies within two standard errors of the line of C1 and C4 is counted to 100 draws. At 0.1.0rc22 the place of the collisions, 0.567 +- 0.090, and the setting of the collisions and of the transit, 0.667 +- 0.086 and 0.633 +- 0.088 (his message 402)
SETS = 100          # the synthetic sets of the twin of a case, read once
OWN = 50          # 50 a right model, each draw read against its own twin at the loop's OWN_SETS sets, as a researcher reads one record (the red team of his message 308, severe 3): at 10 a rate of 5 to 14 in 100 held 0 or 1 hit. 10 at 0.1.0rc6 until then: a draw read against its own twin builds that twin, 5 to 17 minutes of a core a draw at 0.0.8, so 50 a right model was 28 hours of the 33 that the count would take, where the false alarms of the red team of message 303 come from the 200 draws of each right model against the twin of its case. 50 since his message 196 until then
OWN_SETS = 30          # the sets of the twin of a draw read against its own fit: those of `twinfit.loop.run` and of the command line, so the false alarm counted is the one that a researcher meets (the red team of his message 308, severe 3)
PARTS = 10          # the twin of a case is made in this many units of the queue, each of SETS / PARTS sets, so that no unit outlasts a pre-wave
EXPECTED = {          # what the physics of each term of the toy says of where it shows, written before any count. The place is the one that the package's own report gives on noiseless traces, toy/places.json, rows M203 and M223: the transit time was written as the wings, a guess, which his message 184 questioned
    "collisions": {"where": "the core", "setting": "T", "physics": "a Lorentzian width proportional to the density of the vapour, which grows steeply with the temperature and is one at 130 C whatever the power"},
    "gas_width": {"where": "the core", "setting": None, "physics": "one Lorentzian width, the same in every trace"},
    "laser_noise": {"where": "the core", "setting": None, "physics": "one Gaussian width, the same in every trace"},
    "transit": {"where": "the flanks", "setting": "T", "physics": "the crossing of the beam by the atoms, a width that grows as the root of the temperature, with wings of its own"},
}


WRITTEN = {          # the answers of shape, of the noise and of what changes along a setting that the physics of each term admits, written from the physics alone before any count of the blind report (the red team of his message 343, S1), outside MAKES so that a grade moves no row. Each term of the toy is a symmetric width of the line: no slow trend and no periodic structure over the span, and no law of the noise, so the shape and the noise are answered by none of them, and a width along the temperature changes the core, the flanks or the wings along it
    "collisions": {"shape": None, "what": None, "changes along": ("the core along T", "the flanks along T", "the wings along T")},
    "gas_width": {"shape": None, "what": None, "changes along": None},
    "laser_noise": {"shape": None, "what": None, "changes along": None},
    "transit": {"shape": None, "what": None, "changes along": ("the core along T", "the flanks along T", "the wings along T")},
}
QUESTIONS = (("place", "where"), ("shape", "shape"), ("what", "what"), ("changes_along", "changes along"))          # the four questions of the report, his messages 202 and 308, by the name of their column and of their question
BEFORE_AT = "859fb07"          # the commit of the work copy whose toy/pointing.json holds the shares at 0.0.8, set beside those of 0.1.0rc6 (the red team of his message 303, severe 1)


def cases(config: dict) -> list:
    """The cases of the count, as toy/power.py makes those in which nobody names the term: for each detail after the first,
    the model that lacks its term on its own traces, where it is right, and on the traces of that detail at the size of
    the twin and at a half of it. Each as (its name, the detail of the traces, the detail of the model, its free
    parameters, the sizes or None, the term that is hidden or None)."""
    out = []
    for name, made, d, free, sizes in power.hidden(config):
        term = None if made == d else next(t for t, e in config["parameter"].items() if e["enters"] == made)
        out.append((name, made, d, free, sizes, term))
    last = max(int(k) for k in config["detail"])
    if not any(c[1] == c[2] == last for c in out):          # the last model right on its own traces, the one the twin of his line is: its false alarms were never counted (the red team of his message 303, severe 1)
        out.append((f"{last}, right", last, last, list(config["detail"][str(last)]["free"]), None, None))
    return out


def twin_path(root: Path, name: str, part: int | None = None) -> Path:
    tail = "" if part is None else f"_part{part}"
    return Path(root) / "runs" / f"pointing_null_{name.replace(' ', '_').replace(',', '')}_{digest(root)[:12]}{tail}.json"


def twin_part(case: tuple, part: int, const: dict, config: dict, root: Path = ROOT) -> list:
    """One unit of the twin of a case: the statistics of the report on its sets from SETS / PARTS * (part - 1), each fitted
    again as the first draw was, as `where.twin_sets` makes them. Kept in runs/ as it ends."""
    name, made, d, free, sizes, _ = case
    problem, result, _ = analysis.at(d, fits.synthetic(made, const, power.drawn_with(config, sizes), 1), const, config, free)
    sets = twin.replicas(problem, result, SETS, seed=23)
    size = SETS // PARTS
    got = [where._one_set((sets[i], result, 80, analysis.CENTRE)) for i in range(size * (part - 1), size * part)]
    path = twin_path(root, name, part)
    path.parent.mkdir(exist_ok=True)
    aside = path.with_suffix(".part")
    aside.write_text(json.dumps([None if s is None else {"|".join(k): v for k, v in s.items()} for s in got]))
    aside.replace(path)
    return got


def null_of(case: tuple, const: dict, config: dict, root: Path = ROOT) -> list:
    """The twin of a case, read once: the statistics of the report on `SETS` synthetic sets of the fit of its first draw,
    from its parts where the queue made them. Kept in runs/, under what a row depends on."""
    name = case[0]
    path = twin_path(root, name)
    if not path.exists():
        parts = [twin_path(root, name, k) for k in range(1, PARTS + 1)]
        for k, part in enumerate(parts, start=1):
            if not part.exists():          # a part that the queue did not make is made here
                twin_part(case, k, const, config, root)
        made = [json.loads(part.read_text()) for part in parts]
        aside = path.with_suffix(".whole")
        aside.write_text(json.dumps([s for part in made for s in part if s is not None]))
        aside.replace(path)
    return [{tuple(k.split("|", 1)): v for k, v in s.items()} for s in json.loads(path.read_text())]


def one(job: tuple) -> tuple:
    """One row: where the report points on one draw of a case, against the twin of the case."""
    name, made, d, free, sizes, draw, const, config, null = job
    try:
        problem, result, _ = analysis.at(d, fits.synthetic(made, const, power.drawn_with(config, sizes), draw), const, config, free)
        report = where.where(problem, result, sets=OWN_SETS, null=null, centre=analysis.CENTRE)          # with no twin given, the twin of this draw's own fit at the loop's sets
        beyond = [[a["question"], a["answer"], power.said(a["excess"])] for a in report["answers"] if a["beyond_the_twin"]]
        return name, draw, {"first": report["first_of_each_question"], "tied": report.get("tied_with_the_first") or {}, "beyond": beyond, "threshold": power.said(report["largest_excess_on_the_twin"])}
    except (Exception, SystemExit) as why:
        return name, draw, {"failed": f"{type(why).__name__}: {why}"[:300]}


MAKES = ("cases", "twin_part", "null_of", "one", "EXPECTED", "SETS", "OWN_SETS", "PARTS")          # what a row is made with in this program
RUNS = ("toy/analysis.py", "toy/pointing.py")          # the programs whose imports say which files of the package a row runs


def digest(root: Path = ROOT) -> str:
    """What a row depends on: the kept fits of the last detail, what a row takes from the analysis, the parts of this
    program that make a row, and the files of the package that it runs."""
    root = Path(root)
    h = hashlib.sha256(draws.digest(root).encode())
    h.update(power.sources((root / "toy" / "analysis.py").read_text(), ("at", "CENTRE")))
    h.update(power.sources((root / "toy" / "pointing.py").read_text(), MAKES))
    draws.hashed(h, root, draws.package_of(root, RUNS))
    return h.hexdigest()


def kept(root: Path, name: str) -> Path:
    return Path(root) / "runs" / f"pointing_{name.replace(' ', '_').replace(',', '')}_{digest(root)[:12]}.jsonl"


def written(term: str) -> dict:
    """The answers that the physics of a term admits to each of the four questions of the blind report: the place of
    `EXPECTED`, the shape, the noise and what changes along a setting of `WRITTEN`. None where nothing of the term
    answers it."""
    e, w = EXPECTED[term], WRITTEN[term]
    return {"where": None if e["where"] is None else (e["where"],), "shape": w["shape"], "what": w["what"], "changes along": w["changes along"]}


def graded(rows: list, question: str, admitted) -> dict:
    """One question on the rows of a hidden term, as the red team of his message 308 asked it: the share of draws on
    which the first answer is one that the physics admits, the share on which such an answer stands among those beyond
    the twin, and the share on which the report is silent on the question. Where the physics admits none, as for a
    covariate of a term that changes along nothing, a draw is right only where the report spoke and gave this question no
    answer: silence counts apart and never as right."""
    spoke = [bool(r["beyond"]) for r in rows]
    beyond = [{a for q, a, *_ in r["beyond"] if q == question} for r in rows]          # a row of the table keeps the question and the answer, a unit of toy/planted.py its excess too
    silent = [r["first"].get(question) is None for r in rows]
    if admitted is None:
        first = [s and not b for s, b in zip(spoke, beyond)]
        among = first
    else:
        first = [r["first"].get(question) in admitted for r in rows]
        among = [bool(b & set(admitted)) for b in beyond]
    return {"admitted": None if admitted is None else list(admitted), "first": power.share(first), "among_beyond": power.share(among), "silent": power.share(silent)}


def shares(rows: list, term: str | None) -> dict:
    """The shares of the rows of one case: on a hidden term, how often the report puts first the place that its physics
    says, how often what changes along the setting that it says, and both, with how often each answer comes first, and each of the four questions
    graded by `graded`. On every case, how often some answer stands beyond the twin, which on a right model is the rate
    of its false alarms."""
    cell = {"draws": len(rows), "something_beyond_the_twin": power.share([bool(r["beyond"]) for r in rows])}
    if term is not None:
        e = EXPECTED[term]
        place = [r["first"].get("where") == e["where"] or e["where"] in (r.get("tied") or {}).get("where", []) for r in rows]          # a tie that holds the place of the physics is right, his decisions of messages 388 and 389
        setting = [(r["first"].get("changes along") or "").endswith(f" along {e['setting']}") if e["setting"] else r["first"].get("changes along") is None for r in rows]          # the setting read from what changes along it, the blind report holding no question along
        cell.update(expected={"where": e["where"], "setting": e["setting"]}, place=power.share(place), setting=power.share(setting), both=power.share([a and b for a, b in zip(place, setting)]))
        firsts = {}
        for r in rows:
            for q, a in r["first"].items():
                firsts.setdefault(q, {}).setdefault(str(a), 0)
                firsts[q][str(a)] += 1
        cell["first_of_each_question"] = {q: dict(sorted(v.items(), key=lambda kv: -kv[1])) for q, v in firsts.items()}
        admitted = written(term)
        cell["questions"] = {name: graded(rows, question, admitted[question]) for name, question in QUESTIONS}
    return cell


def small(r: dict) -> dict:
    """A row kept small in the table: its first answers, the places tied with the first, and the question and answer of
    each that stands beyond the twin, so that a grade can be read again (the red team of his message 308, severe 2).
    Until 2026-10-26 it left the tie out, and C1 of 0.1.0rc22 read 0.475, the first place alone, on rows that held it."""
    return {"draw": r["draw"], "first": r["first"], "tied": r.get("tied") or {}, "beyond": [[q, a] for q, a, *_ in r["beyond"]]}


def own_kept(root: Path, name: str) -> Path:
    """The rows of a right model whose draws are each read against the twin of their own fit, as a researcher reads one
    record (row M104), apart from the rows read against the twin of the case."""
    return Path(root) / "runs" / f"pointing_own_{name.replace(' ', '_').replace(',', '')}_{digest(root)[:12]}.jsonl"


def near(case: dict | None) -> bool:
    """Whether a hidden case of the table, counted on its first stage, has its share of the place, or of the setting
    where its term changes along one, within two standard errors of `LINE`: it is then counted to `SECOND` draws."""
    if not case or case.get("term") is None or case.get("draws", 0) < DRAWS:
        return False
    shares = [case["place"]] + ([case["setting"]] if (case.get("expected") or {}).get("setting") else [])
    return any(abs(s["share"] - LINE) < 2.0 * s["standard_error"] for s in shares)


def units(root: Path = ROOT) -> list:
    """What the queue asks first: for each case the parts of its twin not yet made, and its draws not yet kept, which
    come after the twin of their case."""
    root = Path(root)
    const, config = model.load(root / "toy", root)
    out = []
    table = root / "toy" / "pointing.json"
    counted_cases = json.loads(table.read_text()).get("cases", {}) if table.exists() else {}
    for c in cases(config):
        whole = twin_path(root, c[0]).exists()
        parts = [] if whole else [k for k in range(1, PARTS + 1) if not twin_path(root, c[0], k).exists()]
        out.append({"condition": "twin|" + c[0], "units": parts, "after": [], "kind": "twin"})
        kept_rows = power.read(kept(root, c[0]))
        most = 2 * DRAWS if c[1] == c[2] else (SECOND if near(counted_cases.get(c[0])) else DRAWS)
        out.append({"condition": c[0], "units": [k for k in range(1, most + 1) if k not in kept_rows], "after": ["twin|" + c[0]] if parts else [], "kind": "right" if c[5] is None else "hidden"})
        if c[5] is None:          # a right model: its draws read each against its own twin too
            own = power.read(own_kept(root, c[0]))
            out.append({"condition": "own|" + c[0], "units": [k for k in range(1, OWN + 1) if k not in own], "after": [], "kind": "own"})
    return out


def slice_of(condition: str, first: int, last: int, limit: float = float("inf"), root: Path = ROOT) -> dict:
    """The units of a condition from `first` to `last` not yet made, each kept as it ends, until `limit` seconds have
    passed: parts of the twin of a case, or draws of a case read against its twin."""
    import resource
    root = Path(root)
    const, config = model.load(root / "toy", root)
    case = next(c for c in cases(config) if c[0] == condition.split("|", 1)[-1])
    began, ids, seconds = time.time(), [], []
    for k in range(first, last + 1):
        t = time.time()
        if condition.startswith("twin|"):
            if twin_path(root, case[0], k).exists():
                continue
            twin_part(case, k, const, config, root)
        elif condition.startswith("own|"):
            path = own_kept(root, case[0])
            if k in power.read(path):
                continue
            name, draw, row = one((case[0], case[1], case[2], case[3], case[4], k, const, config, None))          # with no twin given, the report makes the twin of this draw's own fit
            with path.open("a") as fh:
                fh.write(json.dumps(dict(row, draw=draw)) + "\n")
        else:
            path = kept(root, case[0])
            if k in power.read(path):
                continue
            name, draw, row = one((case[0], case[1], case[2], case[3], case[4], k, const, config, null_of(case, const, config, root)))
            with path.open("a") as fh:          # one line a draw, written whole: waves of one case append side by side
                fh.write(json.dumps(dict(row, draw=draw)) + "\n")
        ids.append(k)
        seconds.append(time.time() - t)
        if time.time() - began > limit:
            break
    return {"done": ids, "seconds": seconds, "rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


def progress(root: Path = ROOT) -> str:
    """The preliminary result after a wave: for each case with rows, how often the report puts first the place that its
    physics says, or on a right model how often something stands beyond the twin, so far."""
    root = Path(root)
    _, config = model.load(root / "toy", root)
    said = []
    for name, made, d, free, sizes, term in cases(config):
        rows = [r for r in power.read(kept(root, name)).values() if "failed" not in r]
        if not rows:
            continue
        if term is None:
            said.append(f"{name}: beyond the twin on {sum(bool(r['beyond']) for r in rows)} of {len(rows)}")
            own = [r for r in power.read(own_kept(root, name)).values() if "failed" not in r]
            if own:
                said.append(f"{name}, each draw against its own twin: beyond it on {sum(bool(r['beyond']) for r in own)} of {len(own)}")
        else:
            said.append(f"{name}: {EXPECTED[term]['where']} first on {sum(r['first']['where'] == EXPECTED[term]['where'] for r in rows)} of {len(rows)}")
    return "; ".join(said) or "no row yet"


def counted(root: Path = ROOT) -> dict:
    """The table of every case from the rows kept in runs/: its shares, and its rows, each kept small, its first answers
    and how many answers stand beyond the twin, so that a copy of the tree with no runs/ counts them again."""
    root = Path(root)
    const, config = model.load(root / "toy", root)
    from twinfit import version
    out = {"digest": digest(root), "package": version.__version__, "expected": EXPECTED, "written": {t: {k: (list(v) if v else None) for k, v in w.items()} for t, w in WRITTEN.items()}, "cases": {}}
    before = before_of(root)
    for name, made, d, free, sizes, term in cases(config):
        kept_rows = power.read(kept(root, name)).values()
        rows = [small(r) for r in kept_rows if "failed" not in r]
        if rows:
            out["cases"][name] = {"term": term, "failed": sum(1 for r in kept_rows if "failed" in r), **shares(rows, term), "before": before.get(name), "rows": rows}
        own_rows = power.read(own_kept(root, name)).values() if term is None else []
        rows = [small(r) for r in own_rows if "failed" not in r]
        if len(own_rows) >= OWN:          # row M104: the same right model, each draw against its own twin, once every draw of it is made
            out.setdefault("own_twin", {})[name] = {"term": None, "failed": sum(1 for r in own_rows if "failed" in r), "sets": OWN_SETS, **shares(rows, None), "rows": rows}
    return out


def before_of(root: Path = ROOT, held: dict | None = None) -> dict:
    """The shares of each case at 0.0.8, from toy/pointing.json at `BEFORE_AT`, or from `held`: the place, the covariate,
    both, and something beyond the twin, set beside those of this count. Empty where git does not hold that commit."""
    if held is None:
        import subprocess
        got = subprocess.run(["git", "-C", str(root), "show", f"{BEFORE_AT}:toy/pointing.json"], capture_output=True, text=True)
        if got.returncode != 0:
            return {}
        held = json.loads(got.stdout)
    return {name: {"package": "0.0.8", "draws": c["draws"], **{k: c[k] for k in ("place", "along", "both", "something_beyond_the_twin") if k in c}} for name, c in held.get("cases", {}).items()}


def write(root: Path = ROOT) -> Path:
    path = Path(root) / "toy" / "pointing.json"
    path.write_text(json.dumps(dict(counted(root), queue=power.queue_of(root, "toy/pointing.py", digest(root))), indent=1) + "\n")          # with the waves that the queue spent on its rows (review 22, minor 18)
    return path


CONFUSED = 0.5          # a pair of terms that the report names the same on more than this share of paired draws, all four questions at once and one of them named, is a pair that the tiers cannot tell apart


def confusion(root: Path = ROOT, rows_of=None) -> dict:
    """Row M406, the red team of his messages 303, moderate 2, and 308, moderate 4: the confusion of the four terms from
    the kept rows of the count, at no compute. For each pair of hidden terms at the size of the twin, and at a half of it,
    draw k of the one set beside draw k of the other: the share of the pairs on which the report gives both the same
    answer to each question of the blind report, the place, the shape, the noise and what changes along a setting, to
    all four, to all four with one of them named, both silent on all four, and some answer apart, each with its
    standard error. The pairs that the report names the same on more
    than CONFUSED of the draws are named: the five tiers cannot tell them apart. Two terms on which the report is silent
    are not told apart either, but the report did not confuse them: it said nothing (review 84, moderate 3). `rows_of(name)`: the rows of a case, read from runs/ where it is not given."""
    from itertools import combinations
    from twinfit import version
    _, config = model.load(Path(root) / "toy", root)
    rows_of = rows_of or (lambda name: [r for r in power.read(kept(root, name)).values() if "failed" not in r])
    hidden = [(name, term, sizes) for name, made, d, free, sizes, term in cases(config) if term is not None]
    out = {"package": version.__version__, "digest": digest(root), "pairs": [], "cannot_tell_apart": []}
    questions = (("place", "where"), ("shape", "shape"), ("noise", "what"), ("change", "changes along"))          # the questions of the blind report (the red team of his message 343, S1)
    for half in (False, True):
        of = {term: {r["draw"]: r["first"] for r in rows_of(name)} for name, term, sizes in hidden if bool(sizes) == half}
        for a, b in combinations(sorted(of), 2):
            draws = sorted(set(of[a]) & set(of[b]))
            same = {f"same_{label}": power.share([of[a][k].get(q) == of[b][k].get(q) for k in draws]) for label, q in questions}
            every = [all(of[a][k].get(q) == of[b][k].get(q) for _, q in questions) for k in draws]
            silent = [all(of[a][k].get(q) is None and of[b][k].get(q) is None for _, q in questions) for k in draws]
            pair = {"terms": [a, b], "size": "half of the twin" if half else "the twin", "draws": len(draws), **same, "same_all": power.share(every),
                    "same_named": power.share([e and not s for e, s in zip(every, silent)]), "both_silent": power.share(silent), "different": power.share([not e for e in every])}
            pair["standard_error"] = pair["same_all"]["standard_error"]
            (out["pairs"] if not half else out.setdefault("pairs_at_half", [])).append(pair)
            if draws and pair["same_named"]["share"] > CONFUSED:
                out["cannot_tell_apart"].append(f"{a} and {b}, at {pair['size']}")
    return out


def check(root: Path = ROOT) -> list:
    """What is wrong with toy/pointing.json: rows made from other files, or shares that its own rows no longer give."""
    path = Path(root) / "toy" / "pointing.json"
    if not path.exists():
        return ["toy/pointing.json is not written: add `toy/pointing.py` to the queue, then `python toy/pointing.py table`"]
    held = json.loads(path.read_text())
    if held.get("digest") != digest(root):
        return ["toy/pointing.json was counted from rows made from other files than those on disk now: run the count again"]
    wrong = [name for part in ("cases", "own_twin") for name, cell in held.get(part, {}).items() if {k: v for k, v in cell.items() if k not in ("term", "failed", "rows", "before", "sets")} != shares(cell["rows"], cell["term"])]
    return [f"toy/pointing.json does not hold the shares that its own rows give, for {', '.join(wrong)}: write the table again"] if wrong else []


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["units"]:
        print(json.dumps(units()))
    elif args[:1] == ["slice"] and len(args) >= 4:
        if not os.environ.get("TOY_QUEUE"):
            sys.exit("REFUSED: a slice is a wave of the queue: `python3 tools/autoqueue.py add toy/pointing.py` (his message 147)")
        print(json.dumps(slice_of(args[1], int(args[2]), int(args[3]), float(args[5]) if args[4:5] == ["--for"] else float("inf"))))
    elif args == ["progress"]:
        print(progress())
    elif args == ["digest"]:          # what the queue judges a wave by: a change that moves no unit makes none stale
        print(digest())
    elif args[:1] == ["run"]:
        sys.exit("REFUSED: the count runs through the queue: `python3 tools/autoqueue.py add toy/pointing.py` (his message 147)")
    elif args == ["confusion"]:          # row M406: from the kept rows, at no compute
        got = confusion()
        (ROOT / "toy" / "confusion.json").write_text(json.dumps(got, indent=1) + "\n")
        print("written: toy/confusion.json. The tiers cannot tell apart: " + ("; ".join(got["cannot_tell_apart"]) or "no pair"))
    elif args[:1] == ["table"]:
        print(f"written: {write().relative_to(ROOT)}")
    elif args[:1] == ["--check"]:
        wrong = check()
        print("\n".join(wrong) or "where the report points: counted from the files as they are, and the kept rows give the same shares")
        sys.exit(1 if wrong else 0)
    else:
        print(__doc__)
        sys.exit(1)
