"""Where to search a new term of the model: one report, read from the whole machinery (the owner's messages 139 to 142).

A fit of a working model leaves residuals. The report asks of them the questions of a senior researcher, and gives each
answer its level, counted on synthetic sets of the fitted model that are fitted again as the records were:

    where    the core of the curves, within half a width of their centre, their flanks, from half a width to two, their
             wings, beyond two widths, or their asymmetry
    what     a polynomial of the level, a periodic function, a second component of the curve, a wider pedestal under
             it, a tail on one side, glitches, a noise that grows with the signal, or a noise correlated from sample
             to sample. A class names the classes whose leading shapes these records cannot tell from its own
    along    each covariate of the design whose value changes from record to record. One that takes one value in every
             record is named among those that cannot be read: nothing along it shows in these records
    seen by  each family of statistics of the package: the generalized moments, of the curve, of its derivative and of
             its running integral, their ratios, the windows along the abscissa and along the value, the crossings, the
             transformations of the value, and every sample, which is the chi square of the fit. The family that sees most
             says by which statistics the missing term is read: a periodic function by the cosines, glitches by every
             sample and by no family of fewer statistics

    report = where(problem, result, sets=100, centre="position")      the answers, each with its statistic, its level,
                                                                       its excess, and the median and 95th percentile of
                                                                       the twin from which the excess is counted
    statistics(problem, result, centre="position")                    the statistics alone, on the records as they are
    null = twin_sets(problem, result, sets=100)                        the statistics on the sets of the twin, once, for
    where(other, its_fit, null=null)                                   records of the same design read against them, as
                                                                       a count of many draws reads one twin of its case

`centre` names the own parameter that places a curve, as `families.catalogue` takes it. A family that the package
refuses on the records, or on more than a few sets of the twin, has no answer.

Each statistic is larger when the residuals hold more of what it reads, and each grows as the square of the residuals:
records much cleaner than the twin's put every answer beyond every set, where only the excesses order them, and a
statistic that grew as their fourth power led by arithmetic alone. On a clean ripple of the toy, at a thousandth of its
noise, the two statistics of the law of the noise stood at 1e12 and 2.5e10 against 2.3e6 for the best shape.

Its level is the share of the synthetic sets on which it is as large, and its excess says how far beyond them it stands: the statistic less the median of the sets, over
the span from that median to their 95th percentile, so that an excess of one is the 95th percentile of the twin. Where
several answers lie beyond every set their levels are all the smallest, and their excess orders them.

A dozen answers are read at once, so some answer stands beyond its own 95th percentile on most sets of a right model.
The twin says how large the largest excess is when nothing is missing: `largest_excess_on_the_twin` is the 95th
percentile, over the sets, of the largest excess of each set against the others. An answer beyond it is beyond the
twin, which some answer is on 5 of 100 of the twin's own sets. On other records of a right model it can be more often:
the twin is made from the one fit of these records, so a count on draws of a right model says how often, as the
problem's own twin makes them. The first answer of each question is the first of it that stands beyond the twin, or
none. A covariate along which no candidate can be made, since the model has no shared parameter to change along it, has
no answer: it is named among what cannot be read. The glitches and the law of the noise are read once the
shape that leads among the classes is taken out of the residuals: a shape that the model lacks is not a law of the noise.

The package names a place, a class and a covariate, never a cause. Which term of physical origin fills them is the
choice of the researcher, or of the agent that stands in the researcher's place. Nothing of any one field is read here.
"""
from __future__ import annotations

import numpy as np

from .candidates import ALONG, ALONG_SQUARE, Candidate, covariates_of, dictionary
from .families import catalogue, family
from .fit import fit
from .linear import Linear
from .propose import ABSORBED, RANK, _one
from .report import by_position, by_position_and_record
from .twin import replicas

CORE, WINGS = 0.5, 2.0          # the core is within half a width of the centre of a curve, the flanks from there to two widths, the wings beyond: a term that acts between the core and the wings, had no place of its own (the owner's messages 184 and 187)
PLACES = (("the core", lambda r: r["to"] <= CORE), ("the flanks", lambda r: r["from"] >= CORE and r["to"] <= WINGS), ("the wings", lambda r: r["from"] >= WINGS))          # each band of report.BANDS in one place
CLASSES = (("baseline", "a polynomial of the level"), ("periodic", "a periodic function"), ("second component", "a second component of the curve"),
           ("pedestal", "a wider pedestal under the curve"), ("tail", "a tail on one side"))
SIMILAR = 0.9        # two classes whose leading shapes, as the records see them, have a cosine beyond this are not told apart


def statistics(problem, result, centre: str | None = None) -> dict:
    """The statistics of the report on the records as they are: for each answer a number that grows with what the
    residuals hold of it. A whitened residual is a residual over its standard deviation. The glitches and the law of
    the noise are read on the residuals once the shape that leads among the classes is taken out of them."""
    return _read(problem, result, centre)[0]


def changes(problem, result) -> dict:
    """How the structure of what the fit leaves changes along each setting of the design, the owner's message 170: "the
    whole point of the package is to look at the structure (and how it's varying with the parameters!) of the residuals".
    For each covariate that takes more than two values across the records, and each place, the core, the flanks, the wings and the
    asymmetry, the readings of `report.by_position_and_record` of each band of that place are set against the
    covariate, centred and of unit length over the records: the sum over the records of covariate times reading, a
    standard normal number at first order where nothing changes along it. The statistic of the place is the sum of the
    squares of those numbers over its bands, so that it grows as the square of the residuals, as every other statistic
    of the report. A covariate with one value, or two, says nothing of a trend and is not read."""
    rows = by_position_and_record(problem, result)
    out = {}
    for key in covariates_of(problem):
        z = np.array([float(rec.design[key]) for rec in problem.records])
        z = z - z.mean()
        if len(set(np.round(z, 12))) <= 2 or not np.any(z):
            continue
        z = z / np.sqrt(z @ z)
        for place, pick, side in [(name, pick, "both_sides") for name, pick in PLACES] + [("the asymmetry", lambda r: True, "one_side_less_the_other")]:
            trend = [float(np.nansum(z * np.asarray(r[side], dtype=float))) for r in rows if pick(r)]
            if trend:
                out[("changes along", f"{place} along {key}")] = float(sum(t * t for t in trend))
    return out


def _read(problem, result, centre: str | None = None) -> tuple:
    """The statistics, and the leading shape of each class as the records see it, of length one."""
    lin = Linear(problem, result)
    out = {}
    rows = by_position(problem, result)
    for place, pick in PLACES:
        out[("where", place)] = float(sum(r["both_sides"] ** 2 for r in rows if pick(r)))
    out[("where", "the asymmetry")] = float(sum(r["one_side_less_the_other"] ** 2 for r in rows))
    out.update(changes(problem, result))
    scores, leading, shapes = {}, None, {}
    keys = covariates_of(problem)
    for cand in dictionary(lin, covariates=keys):
        entry = _one(lin, cand)
        if entry["status"] == "ranked":
            scores.setdefault(cand.family, []).append((cand.name, entry["score"]))
            if leading is None or entry["score"] > leading["score"]:
                leading = entry
            if cand.family not in shapes or entry["score"] > shapes[cand.family][0]:
                shapes[cand.family] = (entry["score"], entry["unit"])
    dense = periodic_best(lin)
    if dense is not None:          # the period that the dense scan finds, scored as every candidate is
        entry = _one(lin, dense[0])
        if entry["status"] == "ranked":
            scores.setdefault("periodic", []).append((dense[0].name, entry["score"]))
            if leading is None or entry["score"] > leading["score"]:
                leading = entry
            if "periodic" not in shapes or entry["score"] > shapes["periodic"][0]:
                shapes["periodic"] = (entry["score"], entry["unit"])
    for kind, said in CLASSES:
        out[("what", said)] = float(max((score for _, score in scores.get(kind, [])), default=0.0))
    r = np.asarray(lin.r, dtype=float)
    if leading is not None:
        r = r - leading["W"] @ (leading["W"].T @ r)          # the noise is read once the shape that leads is taken out: a shape the model lacks is not a law of the noise
    rows = [lin.rows_of(k) for k in range(len(problem.records))]
    ordered = [r[i][np.argsort(np.asarray(rec.x, dtype=float), kind="stable")] for i, rec in zip(rows, problem.records)]
    out[("what", "glitches")] = float(max(np.max(o * o) for o in ordered))          # the largest squared whitened residual
    heights = np.concatenate([np.asarray(problem.predict(k, result.shared, result.own[k]), dtype=float) for k in range(len(problem.records))])
    r = np.concatenate([r[i] for i in rows])
    h = (heights - heights.mean()) / (heights.std() or 1.0)
    out[("what", "a noise that grows with the signal")] = float(abs(h @ (r * r - 1.0)) / np.sqrt(max(h @ h * 2.0, 1e-300)))          # the trend of the squared whitened residuals with the height, a standard normal deviate on the twin
    pairs = sum(o.size - 1 for o in ordered)
    out[("what", "a noise correlated from sample to sample")] = float(abs(sum(float(o[:-1] @ o[1:]) for o in ordered)) / np.sqrt(max(pairs, 1)))          # the products of neighbours, a standard normal deviate on the twin
    for key in keys:
        along = [score for name, score in scores.get("grows with a covariate", []) if name.endswith(ALONG + key) or name.endswith(ALONG_SQUARE + key)]
        if along:          # with no candidate along it there is no answer: a zero would read as nothing that changes along it (review 21)
            out[("along", key)] = float(max(along))
    for name in catalogue(centre):
        try:
            out[("seen by", name)] = float(family(lin, name, centre).reading()["chi_square"])
        except (ValueError, np.linalg.LinAlgError):
            pass          # a family that the package refuses on these records reads nothing here
    return out, {said: shapes[kind][1] for kind, said in CLASSES if kind in shapes}


def periodic_best(lin, chunk: int = 64):
    """The periodic candidate whose score is largest over a dense ladder: a cosine and a sine for every period from the
    span of the records to four of their steps, half a Fourier step apart, so that the nearest keeps a correlation of
    0.9 or more with any sine over the span. The ladder of the dictionary, five octaves of the width, named a ripple of
    1.3 widths a second component on clean traces of the toy, and the dictionary keeps it, since the loop scores every
    candidate with all that `propose` reads. Here each period is scored as `propose._one` scores it, the residuals on
    the part of its columns that the fit cannot absorb, all at once and with nothing else. The candidate and its score by
    the scan, or None when no period fits."""
    pr = lin.problem
    centre = [lin.centre_of(k) for k in range(len(pr.records))]
    span = float(np.median([rec.x[-1] - rec.x[0] for rec in pr.records]))
    step = float(np.median([np.median(np.diff(rec.x)) for rec in pr.records]))
    cycles = np.arange(1.0, span / (4.0 * step) + 1e-9, 0.5)
    best, at = 0.0, None
    for lo in range(0, cycles.size, chunk):
        f = cycles[lo:lo + chunk] / span
        m = f.size
        W, total = np.zeros((lin.n, 2 * m)), np.zeros(2 * m)
        for k, rec in enumerate(pr.records):
            sd, Q = lin.parts[k]
            phase = 2.0 * np.pi * np.outer(np.asarray(rec.x, dtype=float) - centre[k], f)
            block = np.hstack([np.cos(phase), np.sin(phase)]) / sd[:, None]
            total += (block * block).sum(axis=0)
            W[lin.rows_of(k)] = block - Q @ (Q.T @ block)          # freed of the own parameters, as `Linear.reduce` does
        W -= lin.ledger.Q @ (lin.ledger.Q.T @ W)          # and of the shared ones, as `Ledger.visible` does: what the residuals keep
        b = W.T @ lin.r
        cc, ss, cs = (W[:, :m] ** 2).sum(axis=0), (W[:, m:] ** 2).sum(axis=0), (W[:, :m] * W[:, m:]).sum(axis=0)
        seen = np.sqrt(np.maximum(np.stack([cc, ss]), 0.0) / np.maximum(np.stack([total[:m], total[m:]]), 1e-300)).max(axis=0) >= ABSORBED
        for j in np.flatnonzero(seen):
            w, U = np.linalg.eigh(np.array([[cc[j], cs[j]], [cs[j], ss[j]]]))
            keep = w > RANK * w.max()
            c = U[:, keep].T @ np.array([b[j], b[m + j]])
            score = float(np.sum(c * c / w[keep]))
            if score > best:
                best, at = score, float(cycles[lo + j])
    if at is None:
        return None
    P = span / at
    return Candidate(f"a periodic term of {at:g} cycles over the span", "periodic",
                     [[np.cos(2.0 * np.pi * (rec.x - centre[k]) / P) for k, rec in enumerate(pr.records)],
                      [np.sin(2.0 * np.pi * (rec.x - centre[k]) / P) for k, rec in enumerate(pr.records)]], "any"), best


def _one_set(args):
    twin, result, max_iter, centre = args
    try:
        got = fit(twin, result.shared, result.own, max_iter=max_iter)
        return statistics(twin, got, centre) if got.converged else None
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return None


def twin_sets(problem, result, sets: int = 100, seed: int = 23, centre: str | None = None, max_iter: int = 80) -> list:
    """The statistics of the report on `sets` synthetic sets of the fitted model, each fitted again as the records were.
    The sets whose fit did not end are left out."""
    return [s for s in (_one_set((twin, result, max_iter, centre)) for twin in replicas(problem, result, sets, seed=seed)) if s is not None]


def where(problem, result, sets: int = 100, seed: int = 23, fixed=(), centre: str | None = None, max_iter: int = 80, null: list | None = None) -> dict:
    """The report of where to search: each answer with its statistic on the records, its level and its excess on
    `sets` synthetic sets of the fitted model, and whether it stands beyond the twin, first the smallest level and then
    the largest excess. `first_of_each_question` names the answer that leads each question. `fixed`: the covariates that
    the description says take one value in every record, named as what cannot be read. `null`: the statistics of the
    sets of a twin already read, by `twin_sets`, in place of new ones. Fewer than 20 sets that end are refused: a level
    and an excess need them."""
    seen, shapes = _read(problem, result, centre)
    null = twin_sets(problem, result, sets, seed, centre, max_iter) if null is None else null
    if len(null) < 20:
        raise ValueError(f"only {len(null)} synthetic sets of the fitted model ended: a level needs 20 at least")
    answers = []
    read = [key for key in seen if sum(key in t for t in null) >= 0.9 * len(null)]          # a family refused on more than a tenth of the sets of the twin is read on none
    for key in read:
        value, on_twin = seen[key], [t[key] for t in null if key in t]
        typical, high = np.quantile(on_twin, [0.5, 0.95])
        answers.append({"question": key[0], "answer": key[1], "statistic": value, "level": float((1 + sum(v >= value for v in on_twin)) / (len(on_twin) + 1)),
                        "excess": _excess(value, on_twin), "on_the_twin": {"median": float(typical), "95th percentile": float(high)}})
        if key[1] in shapes:
            answers[-1]["not_told_apart_from"] = [other for other, unit in shapes.items() if other != key[1] and abs(float(unit @ shapes[key[1]])) >= SIMILAR]
    largest = [max(_excess(t[key], [u[key] for u in null if u is not t and key in u]) for key in read if key in t) for t in null]
    beyond = float(np.quantile(largest, 0.95))
    for a in answers:
        a["beyond_the_twin"] = a["excess"] > beyond
    answers.sort(key=lambda a: (a["level"], -a["excess"]))
    first = {q: next((a for a in answers if a["question"] == q and a["beyond_the_twin"]), None) for q in ("where", "what", "along", "changes along", "seen by")}          # the first that is beyond: one before it with a smaller level may stand within (review 21)
    unread = [key for key in covariates_of(problem) if ("along", key) not in seen]
    return {"answers": answers, "first_of_each_question": {q: None if a is None else a["answer"] for q, a in first.items()},
            "sets": len(null), "largest_excess_on_the_twin": beyond,
            "cannot_be_read": [f"along {name}: it takes one value in every record" for name in fixed]
            + [f"along {name}: the model has no shared parameter that could change along it" for name in unread]}


def _excess(value: float, on_twin) -> float:
    """How far beyond the sets of the twin a statistic stands: one at their 95th percentile, zero at their median."""
    typical, high = np.quantile(on_twin, [0.5, 0.95])
    return float((value - typical) / max(high - typical, 1e-12))
