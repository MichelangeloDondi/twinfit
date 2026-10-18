"""Synthetic records from a model: the twin. It gives every reading its law before the reading is used on the records.

A replica is a set of synthetic records with the abscissas, the design and the noise of the recorded ones, and a mean
that is known: the fitted working model, with or without a term added to it, with or without nuisances drawn. The same
fit and the same readings are run on each replica. What comes back, replica by replica, answers three questions that
first order alone cannot:

    admission        does the fit return the truth at this level of noise: its searches end, its bias is under a quarter
                     of its spread, its spread is that of the bound, its interval holds the truth as often as it says
    false alarms     under a right model, how often does the table of candidates pass its own family threshold
    a planted term   with a term in the records that the model lacks: how often is it seen, how often is it in the
                     leading set, and what do the estimates return. This is where a guessed term is tried before it
                     is ever added to the model

    check = noiseless(problem, result)                                   # first: with no noise the fit returns what was put in
    draws = replicas(problem, result, sets=200, level=1.0, seed=1, planted=None, nuisances=None)
    table = admission(problem, result, sets=200, level=1.0)
    table = planted(problem, result, candidate_columns, amplitude, sets=200)

The level multiplies the standard deviation of every sample: the ladder of noise is climbed by calling the same
function at 0.1, 0.3, 1, 3 and 10. A nuisance that is known with an uncertainty is drawn in every replica and held at
its central value in the fit: `nuisances(rng, record)` returns the design of the record with the drawn values in it.
"""
from __future__ import annotations

import numpy as np

from .candidates import dictionary
from .fit import fit
from .linear import Linear
from .pool import Lazy, elsewhere, in_order
from .problem import Record
from .propose import propose

ENDED, BIAS, COVERAGE, SPREAD = 0.98, 0.25, (0.92, 0.98), (0.85, 1.15)
HELD = 0.05              # a parameter held at a bound on more than one synthetic set in twenty is not admitted: until review 43 one held set refused a fit
Z = 1.959964
NOMINAL = 0.95                  # the coverage that an interval of Z standard deviations says it has


def bands(sets: int) -> dict:
    """The three bands of the rule of admission when `sets` sets were read. Each band is its stated width, or three
    standard errors of its statistic at that number of sets when that is wider: a bias over its spread has a standard
    error of one over the root of the sets, a ratio of spreads of one over the root of twice the sets, a coverage the
    binomial one. Stated: a bias under 0.25, a spread of 0.85 to 1.15 of the bound, a coverage of 0.92 to 0.98. They are
    held as stated from 144, 201 and 475 sets on. With 200 sets the stated band of the coverage is 1.9 standard errors
    each side: it refuses a right fit on 0.033 of the parameters at each level of noise, and on 0.28 of the fits with two
    parameters over five levels. The first trial of the package by an agent met it (trials/01). What the wider bands
    cost is in the table: with few sets a fit a little outside the stated bands is admitted too."""
    n = max(int(sets), 2)
    bias = max(BIAS, 3.0 / n ** 0.5)
    spread = max((SPREAD[1] - SPREAD[0]) / 2.0, 3.0 / (2.0 * (n - 1)) ** 0.5)
    cover = max((COVERAGE[1] - COVERAGE[0]) / 2.0, 3.0 * (NOMINAL * (1.0 - NOMINAL) / n) ** 0.5)
    return {"bias": bias, "spread": (max(0.0, 1.0 - spread), 1.0 + spread), "coverage": (max(0.0, NOMINAL - cover), min(1.0, NOMINAL + cover))}


def replicas(problem, result, sets: int, level: float = 1.0, seed: int = 1, planted=None, nuisances=None, model_of=None):
    """`sets` synthetic problems, each made when it is read, from its own seed: the 18th can be made without the 17 before
    it. `planted` is a discrepancy, one array a record (or None), added to the mean. `model_of(rng)` gives the model that
    one set is drawn from, for a quantity that is one number for a whole set and that the fit holds at its centre: the
    size of a term held at the centre of its theory, drawn from the law of that theory. The model of the problem is the
    one that every set is fitted with."""
    def one(rep: int):
        rng = np.random.default_rng([seed, int(round(1000 * level)), rep])
        model = problem.model if model_of is None else model_of(rng)            # drawn before the noise of the set, from the generator of the set
        recs = []
        for k, rec in enumerate(problem.records):
            design = rec.design if nuisances is None else nuisances(rng, rec)
            mean = np.asarray(model(rec.x, result.shared, result.own[k], design), dtype=float)
            if planted is not None and planted[k] is not None:
                mean = mean + planted[k]
            sd = level * rec.sd
            recs.append(Record(rec.x, mean + sd * rng.standard_normal(rec.x.size), sd, rec.design, rec.name))
        return problem.with_records(recs)
    return Lazy(one, sets)


def noiseless(problem, result, planted=None, max_iter: int = 80) -> dict:
    """The first rung of the ladder of noise: synthetic records of the fitted model with no noise drawn, fitted again with
    the weights of the stated noise. The search starts where the researcher's own fit starts, at the starts of the problem,
    and not at the values that the records were made at: started at the answer it could not fail (review 36). With no term
    planted the fit must end and return those values. A failure says that the search does not find them from those starts,
    or that the model cannot tell two of its parameters apart, and nothing further is worth running. What it cannot show: an
    error of the model function itself, since the records and the fit use the same one. With a term planted it returns
    what the term does to the estimates by itself, in their standard deviations."""
    recs = []
    for k, rec in enumerate(problem.records):
        mean = np.asarray(problem.model(rec.x, result.shared, result.own[k], rec.design), dtype=float)
        if planted is not None and planted[k] is not None:
            mean = mean + planted[k]
        recs.append(Record(rec.x, mean, rec.sd, rec.design, rec.name))
    got = fit(problem.with_records(recs), max_iter=max_iter)
    off = (got.shared - result.shared) / result.sd
    out = {"ended": bool(got.converged), "a_term_is_planted": planted is not None,
           "estimates_less_what_was_put_in_in_standard_deviations": dict(zip(problem.shared, off.tolist()))}
    if planted is None:
        out["returns_what_was_put_in"] = bool(got.converged and float(np.max(np.abs(off))) < 1e-3)
    return out


def _start(start: str) -> str:
    if start not in ("problem", "fit"):
        raise ValueError("a search starts at the starts of the problem or at the fitted values: start='problem' or start='fit'")
    return start


def _fit(twin, result, start: str, max_iter: int = 80):
    return fit(twin, max_iter=max_iter) if start == "problem" else fit(twin, result.shared, result.own, max_iter=max_iter)


def _admission_of_one_set(args):
    twin, result, start, max_iter = args
    try:
        got = _fit(twin, result, start, max_iter)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return None
    return (got.shared, got.sd) if got.converged else None


def measured(out: dict, b: dict) -> tuple:
    """What an admission measured, in words, and each statistic of it that falls outside its band, named with its value
    and the band (the cold tests of row M169: the page printed the bands alone, the same for every model, and a refusal
    named no statistic that failed)."""
    said, failed = [f"searches ended on {out['ended']:.2f} of the sets"], []
    if out["ended"] < ENDED:
        failed.append(f"searches ended on {out['ended']:.2f} of the sets, under {ENDED:g}")
    for name, c in out["parameters"].items():
        if c is None:
            failed.append(f"{name}: fewer than 30 searches ended, and nothing of it is read")
            continue
        parts = []
        if c["bias_over_spread"] is not None:
            parts.append(f"a bias of {c['bias_over_spread']:+.2f} of its spread")
            if not abs(c["bias_over_spread"]) < b["bias"]:
                failed.append(f"{name}: a bias of {c['bias_over_spread']:+.2f} of its spread, not under {b['bias']:.3g}")
        if c["spread_over_bound"] is not None:
            parts.append(f"a spread of {c['spread_over_bound']:.2f} of its bar of the fit")
            if not b["spread"][0] <= c["spread_over_bound"] <= b["spread"][1]:
                failed.append(f"{name}: a spread of {c['spread_over_bound']:.2f} of its bar of the fit, outside {b['spread'][0]:.3g} to {b['spread'][1]:.3g}")
        if c["coverage"] is not None:
            parts.append(f"an interval that holds the truth on {c['coverage']:.2f} of the sets")
            if not b["coverage"][0] <= c["coverage"] <= b["coverage"][1]:
                failed.append(f"{name}: an interval that holds the truth on {c['coverage']:.2f} of the sets, outside {b['coverage'][0]:.3g} to {b['coverage'][1]:.3g}")
        said.append(f"{name}: " + ", ".join(parts))
    return ". ".join(said), failed


def admission(problem, result, sets: int = 200, level: float = 1.0, seed: int = 1, nuisances=None, start: str = "problem", workers: int = 1,
              max_iter: int = 80, model_of=None) -> dict:
    """The fit on replicas of its own model: the rule that admits it at one level of noise. Each search starts at the starts
    of the problem, as the researcher's own fit does. `start="fit"` starts it at the fitted values, which is faster and says
    nothing of whether the search finds them. `workers`: the sets are read in that many forked processes, and the table
    is the same whatever their number (pool.py). `max_iter`: the limit of steps of each search, as in `fit`. A search
    that stops at the limit counts as one that did not end: where the records hold a curve with no visible peak the
    limit of 80 can be too low, and the table says how many searches ended."""
    if level <= 0:
        raise ValueError("the admission is a rule about spreads and needs noise: for the rung with no noise call noiseless(problem, result)")
    start = _start(start)
    truth = result.shared
    twins = replicas(problem, result, sets, level, seed, nuisances=nuisances, model_of=model_of)
    read = in_order(_admission_of_one_set, Lazy(lambda i: (twins[i], result, start, max_iter), sets), workers)
    est, bars = [r[0] for r in read if r is not None], [r[1] for r in read if r is not None]
    ended = len(est)
    est, bars = np.array(est), np.array(bars)
    out = {"sets": sets, "level": level, "ended": ended / sets, "nuisances_drawn": nuisances is not None, "searches_start_at": start, "parameters": {}}
    if model_of is not None:                                      # written only then: the table of a fit whose sets are those of its own model is what it was
        out["a_model_is_drawn_for_each_set"] = True
    ok = out["ended"] >= ENDED and len(est) >= 30
    bound = result.sd * level
    b = bands(len(est))
    out["bands"] = {"bias_over_spread_under": float(b["bias"]), "spread_over_bound": [float(v) for v in b["spread"]], "coverage": [float(v) for v in b["coverage"]]}
    out["bands_are_wider_than_stated"] = bool(b["bias"] > BIAS + 1e-12 or b["spread"][0] < SPREAD[0] - 1e-12 or b["coverage"][0] < COVERAGE[0] - 1e-12)
    for j, name in enumerate(problem.shared):
        if len(est) < 30:
            out["parameters"][name] = None
            continue
        e = est[:, j] - truth[j]
        spread = float(e.std(ddof=1))
        barred = np.isfinite(bars[:, j])                           # the sets on which the parameter is free: held at a bound it has no bar, and no interval to count
        held = float(1.0 - barred.mean())
        cell = {"bias_over_spread": float(e.mean() / spread) if spread > 0 else None,
                "spread_over_bound": float(spread / bound[j]) if np.isfinite(bound[j]) else None,
                "coverage": float(np.mean(np.abs(e[barred]) < Z * bars[barred, j])) if barred.any() else None,
                "spread_per_cent_of_the_value": float(100.0 * spread / abs(truth[j])) if truth[j] != 0 else None}
        if problem.bounded:                                       # a problem with no bounds holds nothing, and its table is what it was
            cell["held_at_a_bound_on"] = held
        if held > HELD or not np.isfinite(bound[j]):
            cell["passes"] = False                                # held on more sets than that, an estimate has not the law that the rule assumes
            out.setdefault("why_not", []).append(f"the parameter {name} is held at a bound " + ("in the fit of the records" if not np.isfinite(bound[j])
                                                 else f"on {held:.2f} of the sets whose search ended, more than {HELD:g}")
                                                 + ": there it has no bar, and the rule of the admission is about an estimate that is free. Move the bound away, "
                                                 "or read this parameter as one that the bound sets")
        else:                                                     # held on a set in twenty or fewer: the bands judge it, with the share said in the cell
            cell["passes"] = bool(abs(cell["bias_over_spread"]) < b["bias"] and b["spread"][0] <= cell["spread_over_bound"] <= b["spread"][1]
                                  and b["coverage"][0] <= cell["coverage"] <= b["coverage"][1])
        ok = ok and cell["passes"]
        out["parameters"][name] = cell
    out["admitted"] = bool(ok)
    out["measured"], out["failed"] = measured(out, b)
    out["rule"] = (f"searches ended on {ENDED:g} of the sets, a bias under {b['bias']:.3g} of the spread, a spread of {b['spread'][0]:.3g} to {b['spread'][1]:.3g} of the first order bound, "
                   f"an interval that holds the truth on {b['coverage'][0]:.3g} to {b['coverage'][1]:.3g} of the sets"
                   + (f". With {len(est)} sets each of these is known only within its standard error, so each band is widened to three of them. The bands as stated are a bias "
                      f"under {BIAS:g}, a spread of {SPREAD[0]:g} to {SPREAD[1]:g} and a coverage of {COVERAGE[0]:g} to {COVERAGE[1]:g}, and from 475 sets on they are held as stated"
                      if out["bands_are_wider_than_stated"] else ""))          # "made that wide" and "the three" named nothing that a reader could find (the cold tests of row M169, round 2)
    return out


def _is_its_own(cand, term: list) -> bool:
    """Whether a planted term lies in the span of the columns of one candidate: then the candidate comes first by construction."""
    used = [k for k, piece in enumerate(term) if piece is not None]
    if not used or np.all([np.all(term[k] == 0) for k in used]):
        return False
    y = np.concatenate([np.asarray(term[k], dtype=float) for k in used])
    cols = [np.concatenate([np.asarray(col[k], dtype=float) if col[k] is not None else np.zeros_like(term[k]) for k in used]) for col in cand.columns]
    A = np.column_stack(cols)
    left = y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
    return bool(np.linalg.norm(left) <= 1e-8 * np.linalg.norm(y))


def _planted_of_one_set(args):
    """One replica read as the records are: its fit, its list of candidates, its table. None where a search or a table failed."""
    twin, result, start, make, seed, symmetric_only, max_iter = args
    try:
        got = _fit(twin, result, start, max_iter)
        if not got.converged:
            return None
        lin = Linear(twin, got)
        cands = make(lin)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return None
    listed, keys = [c.name for c in cands], {c.key: c.name for c in cands}      # the name that each key has in this set: a name can hold a width that this fit computed
    try:
        table = propose(lin, cands, draws=1000, seed=seed, symmetric_only=symmetric_only)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return {"listed": listed, "keys": keys}
    return {"listed": listed, "keys": keys, "shared": got.shared, "leading": list(table["leading_set"])}


def planted(problem, result, columns: list, amplitude, sets: int = 100, level: float = 1.0, seed: int = 2, block: str | None = None,
            name: str | None = None, symmetric_only: bool = False, start: str = "problem", candidates=None, workers: int = 1, max_iter: int = 80) -> dict:
    """Replicas that hold a term which the working model lacks, read by the fit and by the table of candidates. `columns`
    and `amplitude` give the term as a candidate does. With an amplitude of zero the replicas are those of a right model,
    and the share of sets in which the table passes its family threshold is its rate of false alarms.

    `candidates` is a function of the linearised fit of a replica that returns the list of candidates of its table. It
    must be the list that is read on the records: a rate of false alarms holds for the list it was counted on, and a
    planted term is first or in the leading set only of a list that holds it. Left out, the list is the dictionary of the
    package alone. A `name` that the list does not hold is refused: counted against such a list a planted term is first
    on no set, and the table would say so without a word (the first trial of the package found it, trials/01).

    `name` is the name that the candidate has in the list read on the records. In the list of a replica the same
    candidate can have another name, since the name of a periodic term or of a slice holds a width that the fit of
    that replica computes: it is matched on its key, on every replica. A replica whose list does not hold the key at
    all is counted under `sets_whose_list_does_not_hold_it`."""
    if level <= 0:
        raise ValueError("replicas that are read by the table need noise: for a term with no noise call noiseless(problem, result, planted=...)")
    start = _start(start)
    amplitude = np.atleast_1d(np.asarray(amplitude, dtype=float))
    n_rec = len(problem.records)
    term = []
    for k in range(n_rec):
        parts = [a * np.asarray(col[k], dtype=float) for col, a in zip(columns, amplitude) if col[k] is not None]
        term.append(sum(parts) if parts else None)
    truth = result.shared
    est, seen, first, in_set, sizes = [], 0, 0, 0, []
    done = lacking = 0
    make = candidates if candidates is not None else (lambda lin: dictionary(lin, block=block))
    listed = None
    key = None
    its_own = None
    if name is not None:                                          # the planted candidate is known by its key: its name is the one that it has on the records
        own = make(Linear(problem, result))
        key = next((c.key for c in own if c.name == name), None)
        if key is None:
            raise ValueError(f"the planted candidate '{name}' is not in the list of {len(own)} candidates that the tables of the replicas are made from. "
                             "Pass candidates=, a function of the linearised fit that returns the list read on the records")
        its_own = _is_its_own(next(c for c in own if c.name == name), term)
    twins = replicas(problem, result, sets, level, seed, planted=term if np.any(amplitude != 0) else None)
    if workers <= 1 and not elsewhere():                          # one after another, each read when the loop below comes to it
        read = (_planted_of_one_set((twin, result, start, make, seed, symmetric_only, max_iter)) for twin in twins)
    else:
        read = in_order(_planted_of_one_set, Lazy(lambda i: (twins[i], result, start, make, seed, symmetric_only, max_iter), sets), workers)
    for r in read:
        if r is None:
            continue
        if listed is None:
            listed = r["listed"]
        if "leading" not in r:
            continue
        done += 1
        est.append((r["shared"] - truth) / (result.sd * level))
        leading = r["leading"]
        seen += int(bool(leading))
        sizes.append(len(leading))
        if name is not None:
            here = r["keys"].get(key)                             # looked for on every set, under the name that it has there
            if here is None:
                lacking += 1                                      # the list of this set does not hold it: said below, not counted in silence as not first
            elif leading:
                first += int(leading[0] == here)
                in_set += int(here in leading)
    est = np.array(est)
    out = {"sets": sets, "read": done, "level": level, "amplitude": amplitude.tolist(), "searches_start_at": start,
           "candidates_in_the_list": len(listed) if listed is not None else None, "the_list_is": "the one passed as candidates" if candidates is not None else "the dictionary of the package alone",
           "the_shares_are_over": "the sets that were read: those whose search ended and whose table was made",
           "share_of_sets_in_which_something_is_seen": seen / max(done, 1), "share_of_the_sets_asked_in_which_something_is_seen": seen / max(sets, 1),
           "median_size_of_the_leading_set": float(np.median(sizes)) if sizes else None,
           "estimates_less_truth_in_standard_deviations": {n: {"mean": float(est[:, j].mean()), "spread": float(est[:, j].std(ddof=1))} for j, n in enumerate(problem.shared)} if done > 2 else {}}
    if name is not None:
        out["the_planted_candidate"] = name
        out["share_of_sets_in_which_it_is_first"] = first / max(done, 1)
        out["share_of_sets_in_which_it_is_in_the_leading_set"] = in_set / max(done, 1)
        out["the_planted_term_is_the_named_candidate_s_own_column"] = bool(its_own)
        if its_own:                                               # the first trial kept in touch met it: a share of 1.00 read as evidence
            out["read_its_first_place_as"] = ("holding by construction: the term that was planted is this candidate's own column, so it comes first on nearly every set whatever the cause "
                                              "of the residuals of the records is. It says that the table can see a term of that shape and size, and not that this candidate is the cause")
        if lacking:                                               # written when it is not zero: such a set can show the candidate neither first nor in the leading set
            out["sets_whose_list_does_not_hold_it"] = lacking
    return out
