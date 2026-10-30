"""Where to search a new term of the model: one report, read from the whole machinery (the owner's messages 139 to 142).

A fit of a working model leaves residuals. The report asks of them the questions of a senior researcher, and gives each
answer its level, counted on synthetic sets of the fitted model that are fitted again as the records were:

    where    the core of the curves, within half a width of their centre, their flanks, from half a width to two, their
             wings, beyond two widths, or their asymmetry
    shape    how the residuals are shaped, read from them alone: a slow trend over the span of each curve, by the
             Legendre polynomials of orders 1 to 3, or a periodic structure along the abscissa, by a dense scan of
             periods. Each names the shapes that these records cannot tell from it
    what     glitches, a noise that grows with the signal, or a noise correlated from sample to sample
    changes along  which place of the residuals changes along each covariate of the design whose value changes from
             record to record. One that takes one value in every record is named among those that cannot be read
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
none. The glitches and the law of the noise are read once the shape of the residuals that leads is taken out of them: a
shape that the model lacks is not a law of the noise.

The report is blind (the owner's messages 337 and 340): it holds no list of terms and names none, and every answer is a
reading of the residuals. The package names a place, a shape and a covariate, never a cause. Which term of physical origin fills them is the
choice of the researcher, or of the agent that stands in the researcher's place. Nothing of any one field is read here.
"""
from __future__ import annotations

import numpy as np

from .candidates import Candidate, covariates_of, from_mechanism
from .families import catalogue, family
from .fit import fit
from .linear import Linear
from .propose import ABSORBED, RANK, _one
from .report import by_position, by_position_and_record
from .twin import replicas

CORE, WINGS = 0.5, 2.0          # the core is within half a width of the centre of a curve, the flanks from there to two widths, the wings beyond: a term that acts between the core and the wings, had no place of its own (the owner's messages 184 and 187)
PLACES = (("the core", lambda r: r["to"] <= CORE), ("the flanks", lambda r: r["from"] >= CORE and r["to"] <= WINGS), ("the wings", lambda r: r["from"] >= WINGS))          # each band of report.BANDS in one place
SHAPES_READ = (("trend", "a slow trend over the span of each curve"), ("periodic", "a periodic structure along the abscissa"))          # the answers of the question shape, each a reading of the residuals
TREND_CYCLES = 2.0          # a structure of fewer cycles over the span of a curve than this is a slow trend, of more a periodic one: the Legendre polynomials of orders 1 to 3 of `slow_trend` cross zero up to three times
SIMILAR = 0.9        # two classes whose leading shapes, as the records see them, have a cosine beyond this are not told apart
TIE = 1.96          # two places whose excesses differ by less than this many standard deviations of their difference on the records are named together, his decision of message 388
MISFIT_SETS = 30          # the sets of the twin that carries the misfit, which reads the scatter of an excess on the records, his message 389
SMOOTH = 9          # the samples over which a residual is averaged into the misfit that the twin carries, as the noise of neighbouring samples is read


def statistics(problem, result, centre: str | None = None, terms=(), dictionary: bool = False) -> dict:
    """The statistics of the report on the records as they are: for each answer a number that grows with what the
    residuals hold of it. A whitened residual is a residual over its standard deviation. The glitches and the law of
    the noise are read on the residuals once the shape that leads is taken out of them. `terms`: the researcher's list
    of terms, Mechanisms, empty by default (the owner's message 342). `dictionary`: the package's own list beside it,
    only where the researcher asks it (337 to 343)."""
    return _read(problem, result, centre, terms, dictionary)[0]


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


def slow_trend(lin) -> Candidate:
    """The Legendre polynomials of orders 1 to 3 over the span of each curve, as the columns of one reading: how much
    of the residuals is a slow trend over the span. A reading of the residuals and no term."""
    cols = []
    for order in (1, 2, 3):
        col = []
        for rec in lin.problem.records:
            x = np.asarray(rec.x, dtype=float)
            u = 2.0 * (x - x.min()) / max(float(np.ptp(x)), 1e-300) - 1.0
            col.append(np.polynomial.legendre.legval(u, [0.0] * order + [1.0]))
        cols.append(col)
    return Candidate("a slow trend over the span of each curve", "trend", cols, "any")


def matched_terms(lin, terms, centre: str | None = None, dictionary: bool = False) -> dict:
    """The researcher's list of terms set against the residuals: for each term, the part of what every family reads
    that its shape accounts for, `matched_families`, as the statistic of the question matched. A term asked at a sign
    that it cannot have scores zero. The package holds no such list, and an empty one reads nothing (the owner's message
    342: "For the list of terms, allow the user to define it, but keep it empty at the beginning"). With `dictionary`,
    the package's generic list stands beside it, only where the researcher asks it with a flag (337 to 343)."""
    if not terms and not dictionary:
        return {}
    cands = [from_mechanism(lin, m, float(m.start) if m.start else 0.1) for m in terms]
    if dictionary:
        from .dictionary import candidates as generic
        cands += generic(lin, covariates=covariates_of(lin.problem))
    try:
        got = matched_families(lin, cands, centre)
    except ValueError:
        return {}          # no family reads these records, so nothing is matched
    out = {("matched", c.get("key") or c["name"]): float(c["score"]) for c in got["candidates"]}          # by its key: a period in the name moves with each fit, and a twin's sets would read another name than the records'
    out.update({("matched", r.get("key") or r["name"]): 0.0 for r in got["of_the_sign_they_cannot_have"]})
    return out


def _read(problem, result, centre: str | None = None, terms=(), dictionary: bool = False) -> tuple:
    """The statistics, each a reading of the residuals alone, and the leading direction of each shape as the records
    see it, of length one. No list of terms of the package is read (the owner's messages 337 and 340): `terms` is the
    researcher's own, empty by default (342)."""
    lin = Linear(problem, result)
    out = {}
    rows = by_position(problem, result)
    for place, pick in PLACES:
        out[("where", place)] = float(sum(r["both_sides"] ** 2 for r in rows if pick(r)))
    out[("where", "the asymmetry")] = float(sum(r["one_side_less_the_other"] ** 2 for r in rows))
    out.update(changes(problem, result))
    read, leading = {}, None
    dense = periodic_best(lin)
    for kind, said in SHAPES_READ:
        cand = slow_trend(lin) if kind == "trend" else (dense[0] if dense is not None else None)
        entry = _one(lin, cand) if cand is not None else {"status": "none"}
        if entry["status"] == "ranked":
            out[("shape", said)] = float(entry["score"])
            read[said] = entry["unit"]
            if leading is None or entry["score"] > leading["score"]:
                leading = entry
    r = np.asarray(lin.r, dtype=float)
    if leading is not None:
        r = r - leading["W"] @ (leading["W"].T @ r)          # the noise is read once the shape that leads is taken out: a shape the model lacks is not a law of the noise
    for said, value in noise_reader(lin)(r).items():          # the skewness and the tails of the residuals, and the glitches, the noise that grows and the correlated noise (the owner's message 266)
        out[("noise" if said in NOISE_SHAPES else "what", said)] = value
    for name in catalogue(centre):
        try:
            out[("seen by", name)] = float(family(lin, name, centre).reading()["chi_square"])
        except (ValueError, np.linalg.LinAlgError):
            pass          # a family that the package refuses on these records reads nothing here
    for said, value in beyond_reader(lin)(lin.r).items():          # the records together, which no family gives (the owner's message 269)
        out[("seen by", said)] = value
    out.update(matched_terms(lin, terms, centre, dictionary))
    return out, read


from .families import NOISE_SHAPES, beyond_reader, catalogue, family, noise_reader, noise_shape  # noqa: E402,F401  the readings of the noise live with the families, so that no kept reading depends on this file through them (review 64, minor 2)


def matched_families(lin, candidates, centre: str | None = None, names=None) -> dict:
    """The owner's message 262: the clues of every family of statistics as the pillars of the decision of a turn. For each
    family of the catalogue but the fit of every sample, the fingerprint of the residuals, `Family.fingerprint`, and for
    each candidate the fingerprint of what the residuals keep of it, `Linear.visible`, both standardised in the
    covariance of the family's statistics under a right model. The fingerprints are joined over the families, and a
    candidate is scored by the part of the residuals' fingerprint that its own accounts for, the square of their product
    over the squared length of its own, with the amplitude that gives, and what each family reads once the candidate is
    taken out of the residuals at that amplitude. A family that the package refuses on the records is left out and
    named. A candidate that no family sees is left out. The score is a chi square of one degree of freedom under a right
    model only where the families' statistics are independent, which they are not: the twin gives its level, as it
    gives every level of the report. Returns the families read, the reading of each before any term, and the candidates
    ranked by their score, each with its amplitude and the reading of each family after it."""
    fams, left = {}, []
    for name in (list(names) if names is not None else [n for n in catalogue(centre) if n != "every sample: the fit"]):
        try:
            fams[name] = family(lin, name, centre)
        except (ValueError, np.linalg.LinAlgError) as err:
            left.append({"family": name, "why": str(err)})
    if not fams:
        raise ValueError("no family of the catalogue can be read on these records: " + "; ".join(f"{x['family']}: {x['why']}" for x in left))
    prints = {}
    for name, fam in list(fams.items()):          # a family read on its construction can still refuse its fingerprint, as one of more than LARGEST statistics on a hundred records: it is left out and named, and the match goes on (the code review of 0.1.0rc2, moderate 2)
        try:
            prints[name] = fam.fingerprint(lin.r)
        except (ValueError, np.linalg.LinAlgError) as err:
            left.append({"family": name, "why": str(err)})
            del fams[name]
    if not fams:
        raise ValueError("no family of the catalogue can be read on these records: " + "; ".join(f"{x['family']}: {x['why']}" for x in left))
    Z = np.concatenate(list(prints.values()))
    refused = []
    before = {name: float(p @ p) for name, p in prints.items()}
    ranked = []
    for cand in candidates:
        D = np.column_stack([lin.reduce(col) for col in cand.columns])
        V = np.column_stack([lin.ledger.visible(D[:, j]) for j in range(D.shape[1])])
        total = np.array([lin.norm2(col) for col in cand.columns])
        if total.max() == 0.0 or float(np.max(np.sqrt(np.maximum((V * V).sum(axis=0), 0.0) / np.maximum(total, 1e-300)))) < ABSORBED:          # absorbed whole: no statistic of these records sees it, as `propose._one` reads it
            continue
        T = np.column_stack([np.concatenate([fams[name].fingerprint(V[:, j]) for name in fams]) for j in range(V.shape[1])])
        G = T.T @ T
        w, U = np.linalg.eigh(G)
        keep = w > 1e-9 * max(float(w.max()), 1e-300)
        if not keep.any():
            continue
        amp = (U[:, keep] / w[keep]) @ (U[:, keep].T @ (T.T @ Z))          # the amplitude of each column that matches the fingerprint best
        if getattr(cand, "sign", 0) and amp.size == 1 and float(amp[0]) * cand.sign < 0.0:          # a width asked to be negative is no width, and its score ties with the kernel of the other kind (the owner's message 269)
            refused.append({"name": cand.name, "key": cand.key or cand.name, "amplitude": float(amp[0])})
            continue
        fitted = T @ amp
        score = float(fitted @ fitted)
        at, after = 0, {}
        for name, p in prints.items():
            m = p.size
            left_over = p - fitted[at:at + m]
            after[name] = float(left_over @ left_over)
            at += m
        ranked.append({"name": cand.name, "key": cand.key or cand.name, "family": cand.family, "parity": cand.parity, "score": score, "amplitude": amp.tolist(), "after": after})
    ranked.sort(key=lambda e: -e["score"])
    return {"families": list(fams), "left_out": left, "readings": before, "candidates": ranked, "of_the_sign_they_cannot_have": refused}


def periodic_best(lin, chunk: int = 64):
    """The periodic candidate whose score is largest over a dense ladder: a cosine and a sine for every period from half
    the span of the records, `TREND_CYCLES`, to four of their steps, half a Fourier step apart, so that the nearest keeps a correlation of
    0.9 or more with any sine over the span. Each period is scored as `propose._one` scores it, the residuals on
    the part of its columns that the fit cannot absorb, all at once and with nothing else. The candidate and its score by
    the scan, or None when no period fits."""
    pr = lin.problem
    centre = [lin.centre_of(k) for k in range(len(pr.records))]
    span = float(np.median([rec.x[-1] - rec.x[0] for rec in pr.records]))
    step = float(np.median([np.median(np.diff(rec.x)) for rec in pr.records]))
    cycles = np.arange(TREND_CYCLES, span / (4.0 * step) + 1e-9, 0.5)          # fewer cycles are a slow trend, which `slow_trend` reads: from one, a polynomial of the level was named a periodic structure of 1.5 cycles
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
    return Candidate(f"a periodic structure of {at:g} cycles over the span", "periodic",
                     [[np.cos(2.0 * np.pi * (rec.x - centre[k]) / P) for k, rec in enumerate(pr.records)],
                      [np.sin(2.0 * np.pi * (rec.x - centre[k]) / P) for k, rec in enumerate(pr.records)]], "any"), best


def _one_set(args):
    twin, result, max_iter, centre = args[:4]
    terms = args[4] if len(args) > 4 else ()
    dictionary = args[5] if len(args) > 5 else False
    try:
        got = fit(twin, result.shared, result.own, max_iter=max_iter)
        return (statistics(twin, got, centre, terms, dictionary) if terms or dictionary else statistics(twin, got, centre)) if got.converged else None          # the list passed only where there is one
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return None


def twin_sets(problem, result, sets: int = 100, seed: int = 23, centre: str | None = None, max_iter: int = 80, terms=(), dictionary: bool = False, planted=None) -> list:
    """The statistics of the report on `sets` synthetic sets of the fitted model, each fitted again as the records were.
    The sets whose fit did not end are left out. `planted`: one array a record added to the mean of every set."""
    return [s for s in (_one_set((twin, result, max_iter, centre, terms, dictionary)) for twin in replicas(problem, result, sets, seed=seed, planted=planted)) if s is not None]


def misfit(problem, result, span: int = SMOOTH) -> list:
    """What the fit leaves in each record, averaged over `span` neighbouring samples: the misfit that the records hold,
    with the noise averaged down, which the twin that carries the misfit adds to the mean of each of its sets."""
    out = []
    for k, rec in enumerate(problem.records):
        r = np.asarray(rec.y, dtype=float) - np.asarray(problem.predict(k, result.shared, result.own[k]), dtype=float)
        n = max(1, min(span, r.size)); n -= 1 - n % 2
        out.append(np.convolve(np.pad(r, n // 2, mode="edge"), np.ones(n) / n, mode="valid"))
    return out


def on_the_records(problem, result, answers: list, seed: int, centre, max_iter: int) -> None:
    """The scatter of the excess of each place beyond the twin as the records hold it, `scatter_on_the_records`, read on
    MISFIT_SETS sets of the twin that carries the misfit, each fitted again, in the units of the excess. A twin with no
    misfit scatters an excess by 0.3 to 0.6, while on draws of the toy with a term hidden the same excess scattered 2.5
    to 11 times more, and a tie read on the first named the flanks of the fifth hidden term with its core on 8 draws of 30 where the
    records do not tell them apart (his message 389). Read only where two places stand beyond the twin."""
    places = [a for a in answers if a["question"] == "where" and a.get("beyond_the_twin")]
    if len(places) < 2:
        return
    held = twin_sets(problem, result, MISFIT_SETS, seed + 1, centre, max_iter, planted=misfit(problem, result))
    for a in places:
        on = [t[("where", a["answer"])] for t in held if ("where", a["answer"]) in t]
        span = a["on_the_twin"]["95th percentile"] - a["on_the_twin"]["median"]
        a["scatter_on_the_records"] = float(np.std(on) / max(span, 1e-12)) if len(on) >= 10 else None


def where(problem, result, sets: int = 100, seed: int = 23, fixed=(), centre: str | None = None, max_iter: int = 80, null: list | None = None, terms=(), dictionary: bool = False) -> dict:
    """The report of where to search: each answer with its statistic on the records, its level and its excess on
    `sets` synthetic sets of the fitted model, and whether it stands beyond the twin, first the smallest level and then
    the largest excess. `first_of_each_question` names the answer that leads each question. `fixed`: the covariates that
    the description says take one value in every record, named as what cannot be read. `null`: the statistics of the
    sets of a twin already read, by `twin_sets`, in place of new ones. `terms`: the researcher's list of terms, which
    adds the question matched, empty by default. Fewer than 20 sets that end are refused: a level and an excess need
    them."""
    seen, leading_shapes = _read(problem, result, centre, terms, dictionary)
    null = twin_sets(problem, result, sets, seed, centre, max_iter, terms, dictionary) if null is None else null
    if len(null) < 20:
        raise ValueError(f"only {len(null)} synthetic sets of the fitted model ended: a level needs 20 at least")
    answers = []
    read = [key for key in seen if sum(key in t for t in null) >= 0.9 * len(null)]          # a family refused on more than a tenth of the sets of the twin is read on none
    for key in read:
        value, on_twin = seen[key], [t[key] for t in null if key in t]
        typical, high = np.quantile(on_twin, [0.5, 0.95])
        answers.append({"question": key[0], "answer": key[1], "statistic": value, "level": float((1 + sum(v >= value for v in on_twin)) / (len(on_twin) + 1)),
                        "excess": _excess(value, on_twin), "scatter": _scatter(on_twin), "on_the_twin": {"median": float(typical), "95th percentile": float(high)}})
        if key[1] in leading_shapes:
            answers[-1]["not_told_apart_from"] = [other for other, unit in leading_shapes.items() if other != key[1] and abs(float(unit @ leading_shapes[key[1]])) >= SIMILAR]
    largest = [max(_excess(t[key], [u[key] for u in null if u is not t and key in u]) for key in read if key in t) for t in null]
    beyond = float(np.quantile(largest, 0.95))
    for a in answers:
        a["beyond_the_twin"] = a["excess"] > beyond
    on_the_records(problem, result, answers, seed, centre, max_iter)
    answers.sort(key=lambda a: (a["question"] == "noise", a["level"], -a["excess"]))          # the shape of the noise after the answers on the mean
    first = {q: next((a for a in answers if a["question"] == q and a["beyond_the_twin"]), None) for q in ("where", "shape", "what", "changes along", "seen by", "noise") + (("matched",) if terms or dictionary else ())}
    along = None
    if first["where"] is None and first["changes along"] is not None:          # no place beyond the twin alone, and one that changes along a setting: its place, read along it. The report named no place on 26 draws of 30 with the self-broadening hidden, which it read along the temperature (row M485)
        along = first["changes along"]["answer"]
        place = along.split(" along ", 1)[0]
        first["where"] = next((a for a in answers if a["question"] == "where" and a["answer"] == place), {"question": "where", "answer": place})
    varies = {key for key in covariates_of(problem) if any(k[0] == "changes along" and k[1].endswith(f" along {key}") for k in seen)}
    return {"answers": answers, "first_of_each_question": {q: None if a is None else a["answer"] for q, a in first.items()}, "where_read_along": along,
            "tied_with_the_first": {"where": tied_with(first["where"], answers)},
            "sets": len(null), "largest_excess_on_the_twin": beyond, "mode": "default" if not terms and not dictionary else ", ".join(([f"the researcher's list of {len(terms)} term{'s' if len(terms) != 1 else ''}"] if terms else []) + (["the package's dictionary, asked by its flag"] if dictionary else [])),
            "cannot_be_read": [f"along {name}: it takes one value in every record" for name in fixed]
            + [f"along {name}: no place of the residuals is read along it" for name in covariates_of(problem) if name not in varies and name not in fixed]}


def _scatter(on_twin) -> float:
    """The standard deviation of an excess on the sets of the twin, in the units of the excess: one in about 1.6 for a
    statistic near a normal law."""
    typical, high = np.quantile(on_twin, [0.5, 0.95])
    return float(np.std(on_twin) / max(high - typical, 1e-12))


def tied_with(lead, answers: list) -> list:
    """The places that the records do not tell from the first, `lead`: beyond the twin at its level, their excess within
    TIE standard deviations of the difference of the two as the records scatter them, `on_the_records`. On traces with no
    noise the misfit of a hidden term stands at 21.5 in the flanks and 21.1 in the core, and the report put the core first on 24
    draws of 30 where the flanks are the place of its physics: the two are named together, his decisions of messages 388
    and 389. A place whose scatter on the records was not read ties with none."""
    if not lead or lead.get("scatter_on_the_records") is None:
        return []
    return [a["answer"] for a in answers if a is not lead and a["question"] == lead["question"] and a.get("beyond_the_twin") and a["level"] == lead["level"]
            and a.get("scatter_on_the_records") is not None and abs(a["excess"] - lead["excess"]) <= TIE * float(np.hypot(a["scatter_on_the_records"], lead["scatter_on_the_records"]))]


def _excess(value: float, on_twin) -> float:
    """How far beyond the sets of the twin a statistic stands: one at their 95th percentile, zero at their median."""
    typical, high = np.quantile(on_twin, [0.5, 0.95])
    return float((value - typical) / max(high - typical, 1e-12))
