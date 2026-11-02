"""The small term in the core: which statistic of the whole line sees the collisions at the size of the twin (row M504,
toward 0.1.1, his message 402 and `runs/reports/release-0.1.0rc22/C1_SIZES.md`).

The collisions are a Lorentzian width that grows with the density, steep along the temperature, and the model that lacks
them holds one width for every trace. The fit takes the mean of that broadening into its width, and what it leaves is
read by the report in windows along the abscissa, which no window isolates: at the size of the twin the report is silent
on 31 draws of 100. This program reads two statistics of the whole line that no window gives, blind to the term, for
each shared parameter p of the user's model and each numeric setting of the design: `p along SETTING`, the excess of p
that each record asks on its own, read through the derivative of the mean by p and set against the setting over the
records, and `p across the records`, the spread of those excesses beyond their noise, whatever the setting. Each is one
number that grows where the fit is worse, a standard normal or a chi square at first order where nothing is missing,
and its line is the 95th percentile of itself on the twin of the case, the report's own construction. Counted on the
same draws as the pointing count, the collisions hidden at a half, one, two and three times the size of the twin and
the model right on its own traces, with the report's share from toy/pointing.json beside each.

    python3 tools/autoqueue.py add toy/broad.py
    python toy/broad.py table
    python toy/broad.py --check

Since rows M507 and M508 the count holds two more readings: each draw of a right model against the twin of its own fit at
the loop's default of OWN_SETS sets, the false alarm a researcher sees (`own|` rows, the table's `own_twin`), and each
hidden case against the twin that carries the misfit, the sets of the first draw's fit with what the fit leaves added to
their mean, as the tie's twin is made (`misfit|` and `vs_misfit|` rows, the table's `against_the_misfit`).
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from twinfit import twin, where  # noqa: E402
from twinfit.candidates import covariates_of  # noqa: E402
from twinfit.fit import fit  # noqa: E402
from twinfit.linear import Linear  # noqa: E402

from toy import analysis, draws, fits, model, pointing, power  # noqa: E402

TERM = "collisions"          # the term the count was made for, row M504: the one the report misses at the size of the twin. Since row M506 every term is counted
LINE = 0.95          # the quantile of a statistic on the sets of the twin that is its line, as the report's (where.where)
SETS = pointing.SETS          # the twin of a case: as many sets as the pointing count reads
PARTS = pointing.PARTS          # made in as many units of the queue as the pointing count makes its twin
MAX_ITER = 80          # the iterations of the fit of a set of the twin, as where._one_set gives them
OWN = pointing.OWN          # the draws of a right model read each against the twin of its own fit, row M507
OWN_SETS = pointing.OWN_SETS          # at the loop's default of sets, the false alarm a researcher sees
MISFIT_SEED = 24          # the seed of the twin that carries the misfit, row M508, apart from the twin of the case
OWN_MORE_SETS = pointing.SETS          # row M513: the hidden cases of the release against the twin of each draw's own fit at the sets of a case too, as `OWN_MORE_SETS` of the work copy's toy/pointing.py named it, which the release copy's does not hold: this program came from that copy and was never run in the release copy (2026-11-01)
STORED_SIZES = (2.0, 3.0)          # the sizes of the twin at which each hidden term is drawn for 0.1.1, his message 402, as the work copy's toy/pointing.py named them
STORED_DRAWS = pointing.SECOND          # a stored case counted to the second stage, his message 402, as the work copy's toy/pointing.py counted it



def stored_cases(config: dict) -> list:
    """The cases kept for 0.1.1 and read by no table of the release: each hidden term at each of STORED_SIZES times the
    size of the twin, in the form of `pointing.cases`, made from its case at HALF, as the work copy's toy/pointing.py
    made them."""
    out = []
    for name, made, d, free, sizes, term in pointing.cases(config):
        if sizes:
            for m in STORED_SIZES:
                out.append((name.replace(f"at {power.HALF:g} of the size", f"at {m:g} of the size"), made, d, free, {k: power.said(m / power.HALF * v) for k, v in sizes.items()}, term))
    return out


def stored_case(name: str) -> bool:
    """Whether a case is one of `stored_cases`, by its name."""
    return any(f", at {m:g} of the size of the twin" in name for m in STORED_SIZES)

def statistics(lin) -> dict:
    """For each shared parameter and each numeric setting, `p along SETTING`: the sum over the records of the centred
    setting times the record's own projection of its residuals on the derivative of the mean by p, over the root of its
    variance, a standard normal number at first order where the parameter is one value across the records. And `p
    across the records`: the spread of the records' excesses of p beyond their noise, the chi square of the weighted
    deviations from their weighted mean, with one degree of freedom less than the records. Both whitened and freed of
    the own parameters of each record, as the linear picture gives them."""
    problem = lin.problem
    keys = covariates_of(problem)
    out = {}
    for j, name in enumerate(problem.shared):
        weight, excess = [], []
        for k in range(len(problem.records)):
            a = lin.J[lin.rows_of(k), j]
            w = float(a @ a)
            weight.append(w)
            excess.append(float(a @ lin.r[lin.rows_of(k)]) / w if w > 0.0 else 0.0)
        weight, excess = np.array(weight), np.array(excess)
        held = weight > 0.0
        if held.sum() < 3:
            continue          # a parameter that fewer than three records read says nothing of a change across them
        for key in keys:
            z = np.array([float(rec.design[key]) for rec in problem.records])[held]
            z = z - z.mean()
            if not np.any(z):
                continue
            out[f"{name} along {key}"] = float(abs(np.sum(z * weight[held] * excess[held]) / np.sqrt(np.sum(z * z * weight[held]))))
        mean = float(np.sum(weight[held] * excess[held]) / np.sum(weight[held]))
        out[f"{name} across the records"] = float(np.sum(weight[held] * (excess[held] - mean) ** 2))
    return out


def cases(config: dict) -> list:
    """The cases of this count: every case of the pointing count, kept and stored, each hidden term at every size and
    each model right on its own traces (row M506, the small term on every term). Until 2026-10-27 the cases of TERM alone
    and the right model that lacks it, row M504, whose rows stand under the same names."""
    return pointing.cases(config) + stored_cases(config)


def null_path(root: Path, name: str, part: int | None = None) -> Path:
    tail = "" if part is None else f"_part{part}"
    return Path(root) / "runs" / f"broad_null_{name.replace(' ', '_').replace(',', '')}_{digest(root)[:12]}{tail}.json"


def null_part(case: tuple, part: int, const: dict, config: dict, root: Path = ROOT) -> list:
    """One unit of the twin of a case: the statistics on its sets from SETS / PARTS * (part - 1), each drawn from the fit
    of the first draw with its noise and fitted again, as `where.twin_sets` makes them, with the twin's own seed. Kept
    in runs/ as it ends. A set whose fit does not end is None."""
    name, made, d, free, sizes, _ = case
    problem, result, _ = analysis.at(d, fits.synthetic(made, const, power.drawn_with(config, sizes), 1), const, config, free)
    sets = twin.replicas(problem, result, SETS, seed=23)
    size = SETS // PARTS
    got = []
    for i in range(size * (part - 1), size * part):
        try:
            again = fit(sets[i], result.shared, result.own, max_iter=MAX_ITER)
            got.append(statistics(Linear(sets[i], again)) if again.converged else None)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            got.append(None)
    path = null_path(root, name, part)
    path.parent.mkdir(exist_ok=True)
    aside = path.with_suffix(".part")
    aside.write_text(json.dumps(got))
    aside.replace(path)
    return got


def null_of(case: tuple, const: dict, config: dict, root: Path = ROOT) -> list:
    """The twin of a case, read once: the statistics on its SETS sets, from its parts where the queue made them."""
    path = null_path(root, case[0])
    if not path.exists():
        parts = [null_path(root, case[0], k) for k in range(1, PARTS + 1)]
        for k, part in enumerate(parts, start=1):
            if not part.exists():
                null_part(case, k, const, config, root)
        made = [json.loads(part.read_text()) for part in parts]
        aside = path.with_suffix(".whole")
        aside.write_text(json.dumps([s for part in made for s in part if s is not None]))
        aside.replace(path)
    return json.loads(path.read_text())


def own_kept(root: Path, name: str) -> Path:
    return Path(root) / "runs" / f"broad_own_{name.replace(' ', '_').replace(',', '')}_{digest(root)[:12]}.jsonl"


def own_more_kept(root: Path, name: str) -> Path:
    """The rows of a hidden case of the release whose draws are each read against the twin of their own fit at
    OWN_MORE_SETS sets, row M513, kept apart from those at the default, with the sets in the name."""
    return Path(root) / "runs" / f"broad_own{OWN_MORE_SETS}_{name.replace(' ', '_').replace(',', '')}_{digest(root)[:12]}.jsonl"


def one_own_more(job: tuple) -> tuple:
    """One row of row M513: a draw read against the twin of its own fit at OWN_MORE_SETS sets. The body of `one_own` with
    the sets changed, and not `one_own` with a parameter, because `one_own` is in MAKES and the digest of every row kept,
    so a change of it would move the name of every row of every case, as `pointing.one_more` says of `pointing.one`."""
    name, made, d, free, sizes, draw, const, config = job
    try:
        problem, result, lin = analysis.at(d, fits.synthetic(made, const, power.drawn_with(config, sizes), draw), const, config, free)
        got = statistics(lin)
        sets = twin.replicas(problem, result, OWN_MORE_SETS, seed=23)          # the first OWN_SETS of them those of `one_own`
        null = []
        for s in sets:
            try:
                again = fit(s, result.shared, result.own, max_iter=MAX_ITER)
                null.append(statistics(Linear(s, again)) if again.converged else None)
            except (np.linalg.LinAlgError, FloatingPointError, ValueError):
                null.append(None)
        line = thresholds(null)
        return name, draw, {"statistics": got, "line": line, "sets": OWN_MORE_SETS, "beyond": {k: bool(got[k] > line[k]) for k in line if k in got}}
    except (Exception, SystemExit) as why:
        return name, draw, {"failed": f"{type(why).__name__}: {why}"[:300]}


def one_own(job: tuple) -> tuple:
    """One row of row M507: a draw of a right model read against the twin of its own fit at the loop's default of OWN_SETS
    sets, as a researcher runs the loop: the statistics of the draw against their line on those sets."""
    name, made, d, free, sizes, draw, const, config = job
    try:
        problem, result, lin = analysis.at(d, fits.synthetic(made, const, power.drawn_with(config, sizes), draw), const, config, free)
        got = statistics(lin)
        sets = twin.replicas(problem, result, OWN_SETS, seed=23)
        null = []
        for s in sets:
            try:
                again = fit(s, result.shared, result.own, max_iter=MAX_ITER)
                null.append(statistics(Linear(s, again)) if again.converged else None)
            except (np.linalg.LinAlgError, FloatingPointError, ValueError):
                null.append(None)
        line = thresholds(null)
        return name, draw, {"statistics": got, "line": line, "sets": OWN_SETS, "beyond": {k: bool(got[k] > line[k]) for k in line if k in got}}
    except (Exception, SystemExit) as why:
        return name, draw, {"failed": f"{type(why).__name__}: {why}"[:300]}


def misfit_path(root: Path, name: str, part: int | None = None) -> Path:
    tail = "" if part is None else f"_part{part}"
    return Path(root) / "runs" / f"broad_misfit_null_{name.replace(' ', '_').replace(',', '')}_{digest(root)[:12]}{tail}.json"


def misfit_part(case: tuple, part: int, const: dict, config: dict, root: Path = ROOT) -> list:
    """One unit of the twin that carries the misfit, row M508: the sets of the first draw's fit with what the fit leaves,
    averaged over neighbouring samples, added to their mean, as `where.misfit` gives it and the tie's twin is made, each
    fitted again, and the statistics on each. Kept in runs/ as it ends. A set whose fit does not end is None."""
    name, made, d, free, sizes, _ = case
    problem, result, _ = analysis.at(d, fits.synthetic(made, const, power.drawn_with(config, sizes), 1), const, config, free)
    sets = twin.replicas(problem, result, SETS, seed=MISFIT_SEED, planted=where.misfit(problem, result))
    size = SETS // PARTS
    got = []
    for i in range(size * (part - 1), size * part):
        try:
            again = fit(sets[i], result.shared, result.own, max_iter=MAX_ITER)
            got.append(statistics(Linear(sets[i], again)) if again.converged else None)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            got.append(None)
    path = misfit_path(root, name, part)
    path.parent.mkdir(exist_ok=True)
    aside = path.with_suffix(".part")
    aside.write_text(json.dumps(got))
    aside.replace(path)
    return got


def misfit_of(case: tuple, const: dict, config: dict, root: Path = ROOT) -> list:
    """The twin that carries the misfit of a case, read once, from its parts where the queue made them."""
    path = misfit_path(root, case[0])
    if not path.exists():
        parts = [misfit_path(root, case[0], k) for k in range(1, PARTS + 1)]
        for k, part in enumerate(parts, start=1):
            if not part.exists():
                misfit_part(case, k, const, config, root)
        made = [json.loads(part.read_text()) for part in parts]
        aside = path.with_suffix(".whole")
        aside.write_text(json.dumps([s for part in made for s in part if s is not None]))          # the sets that carry the misfit and whose fit ended
        aside.replace(path)
    return json.loads(path.read_text())


def vs_misfit_kept(root: Path, name: str) -> Path:
    return Path(root) / "runs" / f"broad_vs_misfit_{name.replace(' ', '_').replace(',', '')}_{digest(root)[:12]}.jsonl"


def thresholds(null: list) -> dict:
    """The line of each statistic that the sets of the twin give: its LINE quantile over the sets on which it was read,
    20 at least, as the report asks of its sets."""
    keys = sorted({k for s in null if s is not None for k in s})
    out = {}
    for key in keys:
        values = [s[key] for s in null if s is not None and key in s]
        if len(values) >= 20:
            out[key] = float(np.quantile(values, LINE))
    if not out:
        raise ValueError(f"no statistic is read on 20 sets of the twin of {len(null)}: a line needs 20 at least")
    return out


def one(job: tuple) -> tuple:
    """One row: the statistics on the fit of one draw of a case, each against its line on the twin of the case."""
    name, made, d, free, sizes, draw, const, config, null = job
    try:
        _, _, lin = analysis.at(d, fits.synthetic(made, const, power.drawn_with(config, sizes), draw), const, config, free)
        got = statistics(lin)
        line = thresholds(null)
        return name, draw, {"statistics": got, "line": line, "beyond": {k: bool(got[k] > line[k]) for k in line if k in got}}
    except (Exception, SystemExit) as why:
        return name, draw, {"failed": f"{type(why).__name__}: {why}"[:300]}


MAKES = ("statistics", "cases", "null_part", "null_of", "thresholds", "one", "one_own", "misfit_part", "misfit_of", "TERM", "LINE", "SETS", "PARTS", "MAX_ITER", "OWN", "OWN_SETS", "MISFIT_SEED")          # what a row is made with in this program
RUNS = ("toy/analysis.py", "toy/broad.py")          # the programs whose imports say which files of the package a row runs


def digest(root: Path = ROOT) -> str:
    """What a row depends on: the kept fits of the last detail, what a row takes from the analysis, the parts of this
    program that make a row, and the files of the package that it runs."""
    root = Path(root)
    h = hashlib.sha256(draws.digest(root).encode())
    h.update(power.sources((root / "toy" / "analysis.py").read_text(), ("at", "CENTRE")))
    h.update(power.sources((root / "toy" / "broad.py").read_text(), MAKES))
    draws.hashed(h, root, draws.package_of(root, RUNS))
    return h.hexdigest()


def kept(root: Path, name: str) -> Path:
    return Path(root) / "runs" / f"broad_{name.replace(' ', '_').replace(',', '')}_{digest(root)[:12]}.jsonl"


def most_of(case: tuple) -> int:
    """The draws of a case that the pointing count reads, so that the two counts stand on the same draws."""
    name, made, d = case[0], case[1], case[2]
    return 2 * pointing.DRAWS if made == d else STORED_DRAWS if stored_case(name) else pointing.SECOND


def units(root: Path = ROOT) -> list:
    """What the queue asks first: for each case the parts of its twin not yet made, and its draws not yet kept, which
    come after the twin of their case."""
    root = Path(root)
    const, config = model.load(root / "toy", root)
    out = []
    for c in cases(config):
        whole = null_path(root, c[0]).exists()
        parts = [] if whole else [k for k in range(1, PARTS + 1) if not null_path(root, c[0], k).exists()]
        out.append({"condition": "twin|" + c[0], "units": parts, "after": [], "kind": "twin"})
        kept_rows = power.read(kept(root, c[0]))
        out.append({"condition": c[0], "units": [k for k in range(1, most_of(c) + 1) if k not in kept_rows], "after": ["twin|" + c[0]] if parts else [], "kind": "right" if c[5] is None else "hidden"})
        if c[5] is None:          # a right model: its draws read each against the twin of its own fit at the default too, row M507
            own = power.read(own_kept(root, c[0]))
            out.append({"condition": "own|" + c[0], "units": [k for k in range(1, OWN + 1) if k not in own], "after": [], "kind": "own"})
            more = power.read(own_more_kept(root, c[0]))          # row M514: and at the sets of a case, the false alarm of the statistic that the one level of 0.1.1 reads at decision 13's sets
            out.append({"condition": f"own{OWN_MORE_SETS}|" + c[0], "units": [k for k in range(1, OWN + 1) if k not in more], "after": [], "kind": "own_more"})
        else:          # a hidden case: the twin that carries the misfit, and its draws read against that line, row M508
            whole_m = misfit_path(root, c[0]).exists()
            parts_m = [] if whole_m else [k for k in range(1, PARTS + 1) if not misfit_path(root, c[0], k).exists()]
            out.append({"condition": "misfit|" + c[0], "units": parts_m, "after": [], "kind": "twin"})
            vs = power.read(vs_misfit_kept(root, c[0]))
            out.append({"condition": "vs_misfit|" + c[0], "units": [k for k in range(1, most_of(c) + 1) if k not in vs], "after": ["misfit|" + c[0]] if parts_m else [], "kind": "hidden"})
            size = "half" if "at 0.5" in c[0] else "2" if "at 2 of" in c[0] else "3" if "at 3 of" in c[0] else "1"          # row M513: a hidden case against the twin of each draw's own fit, at the default and at the sets of a case, the reading a researcher makes with the statistic, which row M507 made of the right models alone, and row M517 the stored cases at 2 and 3 times the size. A kind per term and size (lesson 71)
            draws = 2 * OWN if c[5] == TERM and not c[4] else OWN          # the collisions at the size of the twin on the 100 draws C1's line of 0.60 is read on, which decides whether 0.1.1 holds the small term (review 128, moderate 4)
            own = power.read(own_kept(root, c[0]))
            out.append({"condition": "own|" + c[0], "units": [k for k in range(1, draws + 1) if k not in own], "after": [], "kind": f"own_hidden_{c[5]}_{size}"})
            more = power.read(own_more_kept(root, c[0]))
            out.append({"condition": f"own{OWN_MORE_SETS}|" + c[0], "units": [k for k in range(1, draws + 1) if k not in more], "after": [], "kind": f"own_more_{c[5]}_{size}"})
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
            if null_path(root, case[0], k).exists():
                continue
            null_part(case, k, const, config, root)
        elif condition.startswith("misfit|"):
            if misfit_path(root, case[0], k).exists():
                continue
            misfit_part(case, k, const, config, root)
        elif condition.startswith("own|"):
            path = own_kept(root, case[0])
            if k in power.read(path):
                continue
            name, draw, row = one_own((case[0], case[1], case[2], case[3], case[4], k, const, config))
            with path.open("a") as fh:
                fh.write(json.dumps(dict(row, draw=draw)) + "\n")
        elif condition.startswith(f"own{OWN_MORE_SETS}|"):          # row M513: the same at the sets of a case
            path = own_more_kept(root, case[0])
            if k in power.read(path):
                continue
            name, draw, row = one_own_more((case[0], case[1], case[2], case[3], case[4], k, const, config))
            with path.open("a") as fh:
                fh.write(json.dumps(dict(row, draw=draw)) + "\n")
        elif condition.startswith("vs_misfit|"):
            path = vs_misfit_kept(root, case[0])
            if k in power.read(path):
                continue
            name, draw, row = one((case[0], case[1], case[2], case[3], case[4], k, const, config, misfit_of(case, const, config, root)))
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


def shares(rows: list) -> dict:
    """For each statistic, the share of the draws on which it stands beyond its line, with its standard error, over the
    draws on which it was read."""
    keys = sorted({k for r in rows for k in r["beyond"]})
    return {key: power.share([r["beyond"][key] for r in rows if key in r["beyond"]]) for key in keys}


def beside(root: Path, name: str, term) -> dict | None:
    """The report's share on the same case from toy/pointing.json, where it is written: the place first on a hidden
    term, something beyond the twin on a right model."""
    path = Path(root) / "toy" / "pointing.json"
    if not path.exists():
        return None
    held = json.loads(path.read_text())
    cell = held.get("cases", {}).get(name) or held.get("stored", {}).get(name)
    if cell is None:
        return None
    key = "something_beyond_the_twin" if term is None else "place"
    return {"package": held.get("package"), "reads": key, **{k: cell[key][k] for k in ("share", "standard_error", "of") if key in cell and k in cell[key]}}


def progress(root: Path = ROOT) -> str:
    """The preliminary result after a wave: for each case with rows, the share beyond the line of each statistic so far."""
    root = Path(root)
    _, config = model.load(root / "toy", root)
    said = []
    for name, made, d, free, sizes, term in cases(config):
        rows = [r for r in power.read(kept(root, name)).values() if "failed" not in r]
        if rows:
            got = shares(rows)
            said.append(f"{name}: " + ", ".join(f"{k} {got[k]['share']:.2f} of {got[k]['of']}" for k in got))
        for label, path in (("against its own twin at the default", own_kept(root, name)), ("against the twin that carries the misfit", vs_misfit_kept(root, name))):
            rows = [r for r in power.read(path).values() if "failed" not in r]
            if rows:
                got = shares(rows)
                said.append(f"{name} {label}: " + ", ".join(f"{k} {got[k]['share']:.2f} of {got[k]['of']}" for k in got))
    return "; ".join(said) or "no row yet"


def counted(root: Path = ROOT) -> dict:
    """The table of every case from the rows kept in runs/: the share of each statistic beyond its line, the report's
    share beside it, and the rows, each kept small, so that a copy of the tree with no runs/ counts them again."""
    root = Path(root)
    _, config = model.load(root / "toy", root)
    from twinfit import version
    out = {"digest": digest(root), "package": version.__version__, "term": TERM, "line": LINE, "cases": {}}
    for name, made, d, free, sizes, term in cases(config):
        kept_rows = [r for _, r in sorted(power.read(kept(root, name)).items())]
        rows = [{"draw": r["draw"], "beyond": r["beyond"]} for r in kept_rows if "failed" not in r]
        if rows:
            out["cases"][name] = {"term": term, "failed": sum(1 for r in kept_rows if "failed" in r), "statistics": shares(rows), "package": beside(root, name, term), "rows": rows}
        if term is None:          # row M507: the right model's draws against the twin of its own fit at the default, and row M514 at the sets of a case
            for part, sets, path in (("own_twin", OWN_SETS, own_kept(root, name)), ("own_twin_more", OWN_MORE_SETS, own_more_kept(root, name))):
                own_rows = [r for _, r in sorted(power.read(path).items())]
                rows = [{"draw": r["draw"], "beyond": r["beyond"]} for r in own_rows if "failed" not in r]
                if rows:
                    out.setdefault(part, {})[name] = {"term": None, "sets": sets, "failed": sum(1 for r in own_rows if "failed" in r), "statistics": shares(rows), "rows": rows}
        else:          # row M508: the hidden case's draws against the twin that carries the misfit
            vs_rows = [r for _, r in sorted(power.read(vs_misfit_kept(root, name)).items())]
            rows = [{"draw": r["draw"], "beyond": r["beyond"]} for r in vs_rows if "failed" not in r]
            if rows:
                out.setdefault("against_the_misfit", {})[name] = {"term": term, "failed": sum(1 for r in vs_rows if "failed" in r), "statistics": shares(rows), "rows": rows}
            for part, path, sets in (("own_twin", own_kept(root, name), OWN_SETS), ("own_twin_more", own_more_kept(root, name), OWN_MORE_SETS)):          # rows M513 and M517: every hidden case against the twin of each draw's own fit, at the default and at the sets of a case
                got_rows = [r for _, r in sorted(power.read(path).items())]
                rows = [{"draw": r["draw"], "beyond": r["beyond"]} for r in got_rows if "failed" not in r]
                if rows:
                    out.setdefault(part, {})[name] = {"term": term, "sets": sets, "failed": sum(1 for r in got_rows if "failed" in r), "statistics": shares(rows), "rows": rows}
    return out


def write(root: Path = ROOT) -> Path:
    path = Path(root) / "toy" / "broad.json"
    path.write_text(json.dumps(dict(counted(root), queue=power.queue_of(root, "toy/broad.py", digest(root))), indent=1) + "\n")
    return path


def check(root: Path = ROOT) -> list:
    """What is wrong with toy/broad.json: rows made from other files, or shares that its own rows no longer give."""
    path = Path(root) / "toy" / "broad.json"
    if not path.exists():
        return ["toy/broad.json is not written: add `toy/broad.py` to the queue, then `python toy/broad.py table`"]
    held = json.loads(path.read_text())
    if held.get("digest") != digest(root):
        return ["toy/broad.json was counted from rows made from other files than those on disk now: run the count again"]
    wrong = [name for part in ("cases", "own_twin", "own_twin_more", "against_the_misfit") for name, cell in held.get(part, {}).items() if cell["statistics"] != shares(cell["rows"])]
    return [f"toy/broad.json does not hold the shares that its own rows give, for {', '.join(wrong)}: write the table again"] if wrong else []


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["units"]:
        print(json.dumps(units()))
    elif args[:1] == ["slice"] and len(args) >= 4:
        if not os.environ.get("TOY_QUEUE"):
            sys.exit("REFUSED: a slice is a wave of the queue: `python3 tools/autoqueue.py add toy/broad.py` (his message 147)")
        print(json.dumps(slice_of(args[1], int(args[2]), int(args[3]), float(args[5]) if args[4:5] == ["--for"] else float("inf"))))
    elif args == ["progress"]:
        print(progress())
    elif args == ["digest"]:
        print(digest())
    elif args[:1] == ["run"]:
        sys.exit("REFUSED: the count runs through the queue: `python3 tools/autoqueue.py add toy/broad.py` (his message 147)")
    elif args[:1] == ["table"]:
        print(f"written: {write().relative_to(ROOT)}")
    elif args[:1] == ["--check"]:
        wrong = check()
        print("\n".join(wrong) or "the small term in the core: counted from the files as they are, and the kept rows give the same shares")
        sys.exit(1 if wrong else 0)
    else:
        print(__doc__)
        sys.exit(1)
