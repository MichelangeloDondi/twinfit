"""The candidates of a stage developed as fitted models, and the default choice among them.

At first order a candidate is ranked with no fit (propose.py). That is where a stage begins and not where it ends: a term
that is not small moves the estimates by more, or otherwise, than first order says. So the most plausible candidates,
three by default, are each made into a model with one parameter more (mechanisms.py), fitted from the current fit, and
set beside each other: what each gains, what its own parameter comes back as, what it does to every other parameter,
and what the criteria say (criteria.py).

    table = develop(problem, result, [mechanism, ...])            each candidate fitted, in the order given
    law = null_gains(problem, result, mechanisms, sets=30)       the gains that the same candidates reach on synthetic
                                                                 sets of the working model: the law of a gain when nothing is missing
    choice = decide(table, law)                                  the default choice of the stage
    held = carried(problem, result, mechanism)                   a term on theory that the records cannot fix: held at its
                                                                 theory size, with what the bar of the theory costs each estimate
    trial = on_the_twin(problem, result, table[0], sets=30)      the chosen candidate put into synthetic records: what the
                                                                 working model then shows, and what the whole procedure returns

`workers` above 1 runs the synthetic sets in several processes. The processes are forked, so a model written in a file
or in the session both work. Each process should use one thread of the linear algebra library: with many, a
fork can hang, and the sets would not run faster.

What a developed candidate does not settle. Two candidates with nearly the same gain cannot be told apart by these
records. A gain is the gain of the best form that was tried, and the true term may have another form. And an estimate
made after a candidate was chosen for its gain was chosen after looking: its bar keeps its level only once synthetic
records with the term have shown what the whole procedure returns.
"""
from __future__ import annotations

import numpy as np

from .criteria import choose, criteria, family_threshold, level_of
from .pool import Lazy
from .pool import in_order as _map
from .fit import fit
from .linear import Linear
from .twin import replicas


def working_fit_of(problem, result) -> str:
    """A short digest of a working fit: the abscissae, the values and the weights of its records, the names of its
    parameters and their fitted values. A law of the gains holds for the fit that it was simulated from, and `decide`
    reads this to refuse one that was simulated from other records, other weights or another working model. Two
    models that give the same fitted values to the last digit on the same records are not told apart: the model is a
    function, and a function has no digest. And the same model fitted again to another tolerance is another fit here,
    to its last digit: its law is simulated again, which costs a run and cannot mislead."""
    import hashlib
    h = hashlib.sha256()
    for rec in problem.records:
        for arr in (rec.x, rec.y, rec.sd):
            h.update(np.ascontiguousarray(arr, dtype=float).tobytes())
    h.update(("|".join(problem.shared) + "||" + "|".join(problem.own)).encode())
    h.update(np.ascontiguousarray(result.shared, dtype=float).tobytes())
    h.update(np.ascontiguousarray(result.own, dtype=float).tobytes())
    return h.hexdigest()[:16]


def develop_one(problem, result, mechanism, max_iter: int = 80) -> dict:
    """One candidate as a fitted model: the working model with the term, fitted from the current fit."""
    larger = mechanism.extend(problem, result)
    try:
        got = fit(larger, max_iter=max_iter)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError) as err:          # a model that refuses a point, as the toy's line does where a width overflows, fails this fit and not the loop (his message 156)
        return {"name": mechanism.name, "family": mechanism.family, "ground": mechanism.ground, "ended": False, "why_not": "its search failed: " + str(err),
                "theory_given": mechanism.theory is not None, "working_fit": working_fit_of(problem, result), "_mechanism": mechanism}
    n = len(problem.shared)
    value, bar = float(got.shared[n]), float(got.sd[n])
    at_its_bound = n in got.at_bound                        # the search ended where the term is absent, or as large as it may be
    out = {"name": mechanism.name, "family": mechanism.family, "ground": mechanism.ground, "ended": bool(got.converged), "steps": got.iterations,
           "at_its_bound": bool(at_its_bound), "theory_given": mechanism.theory is not None, "working_fit": working_fit_of(problem, result),
           "parameter": {"name": mechanism.parameter, "value": value, "bar": bar, "started_at": mechanism.start,
                         "from_its_start_in_bars": (value - mechanism.start) / bar if bar > 0 and not at_its_bound else float("nan")},
           "criteria": criteria(max(result.chi2 - got.chi2, 0.0) if at_its_bound else result.chi2 - got.chi2, 1, problem.samples, len(problem.records),
                                problem.samples - result.dof),   # the parameters that the working fit took out of the residuals
           "chi_square": got.chi2, "held": got.held,
           "what_it_does_to_each_parameter": {name: {"before": float(result.shared[j]), "after": float(got.shared[j]),
                                                     "moved_in_bars_of_before": float((got.shared[j] - result.shared[j]) / result.sd[j]),
                                                     "bar_before": float(result.sd[j]), "bar_after": float(got.sd[j])} for j, name in enumerate(problem.shared)}}
    if mechanism.theory is not None:
        centre, sd = mechanism.theory
        out["parameter"]["theory"] = {"centre": centre, "bar": sd, "from_theory_in_joint_bars": (value - centre) / float(np.hypot(bar, sd))}
    if not got.converged:                                     # a search that stopped before its limit is not helped by a higher limit: say which of the two it was
        out["why_not"] = (f"its search did not end in {max_iter} steps" if got.iterations >= max_iter else
                          f"its search stopped after {got.iterations} steps, before its limit of {max_iter}: no step lowered the chi square any more, "
                          "though the step that first order asks was not yet small. More steps do not help: look at the model of the term")
    out["_fit"], out["_problem"], out["_mechanism"] = got, larger, mechanism       # for the caller that goes on with this candidate. Not written to a table
    return out


def develop(problem, result, mechanisms: list, max_iter: int = 80) -> list:
    return [develop_one(problem, result, m, max_iter) for m in mechanisms]


def _gains_of_one_set(args) -> list:
    twin, mechanisms, max_iter = args
    try:
        base = fit(twin, max_iter=max_iter)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return [float("nan")] * len(mechanisms)
    row = []
    for m in mechanisms:
        try:
            got = fit(m.extend(twin, base), max_iter=max_iter)
            row.append(max(float(base.chi2 - got.chi2), 0.0) if base.converged and got.converged else float("nan"))     # a larger model fitted to its end never fits worse
        except (np.linalg.LinAlgError, FloatingPointError, ValueError):
            row.append(float("nan"))
    return row


def null_gains(problem, result, mechanisms: list, sets: int = 30, level: float = 1.0, seed: int = 5, max_iter: int = 80, workers: int = 1) -> dict:
    """The gain of each candidate on synthetic sets of the working model at its fit, where nothing is missing. Each set is
    fitted with the working model and then with each candidate."""
    twins = replicas(problem, result, sets, level=level, seed=seed)
    jobs = Lazy(lambda i: (twins[i], mechanisms, max_iter), sets)
    g = np.array(_map(_gains_of_one_set, jobs, workers), dtype=float)
    whole = g[np.all(np.isfinite(g), axis=1)]
    out = {"sets_asked": sets, "sets_that_ended": int(whole.shape[0]), "sets_dropped": int(sets - whole.shape[0]), "names": [m.name for m in mechanisms],
           "working_fit": working_fit_of(problem, result), "noise_level": float(level), "gains": whole.tolist(),
           "what_the_dropped_sets_cost": "a set is dropped when one of its searches did not end. Its other gains are lost with it, so the law is that of the sets that "
                                         "ended and can sit low: with a share of sets dropped, fit the candidate that fails apart or leave it out"}
    if whole.shape[0] >= 2:
        out.update(family_threshold(whole))
    return out


def _one_twin_of_the_larger(args) -> dict:
    working, mechanism, table_of, max_iter = args
    out = {"ended": False}
    try:
        base = fit(working, max_iter=max_iter)
        if not base.converged:
            return out
        tab = table_of(Linear(working, base)) if table_of is not None else None
        got = fit(mechanism.extend(working, base), max_iter=max_iter)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return out
    if not got.converged:
        return out
    out.update({"ended": True, "working": base.shared.tolist(), "larger": got.shared.tolist(), "bars": got.sd.tolist(), "gain": float(base.chi2 - got.chi2)})
    if tab is not None:
        ranked = [c for c in tab["candidates"] if c["status"] == "ranked"]
        out.update({"seen": bool(tab["leading_set"]), "first": ranked[0]["name"] if ranked else None, "family": ranked[0]["family"] if ranked else None})
    return out


def on_the_twin(problem, result, developed: dict, sets: int = 30, level: float = 1.0, seed: int = 6, table_of=None, max_iter: int = 80, workers: int = 1) -> dict:
    """A developed candidate put into synthetic records. The sets are drawn from the larger model at its fit, so the term
    is in them at the size that was fitted. Each is then read as the records were: the working model fitted, its
    table of candidates when `table_of(linear)` is given, the larger model fitted. What comes back:

      what the working model shows of the term: its estimates, to set beside those of the records, and how often its
        table sees that something is missing, with the family that it then puts first
      what the whole procedure returns: the larger model's estimates against the values that were put in, their
        scatter against their bars, and how often a bar holds the value. This is where a bar earns its level after
        a term was chosen for its gain
      the gains on these sets, and where the gain of the records stands among them

    A candidate that reproduces the readings of the records is not excluded. Nothing more can be said of it."""
    larger, at = developed["_problem"], developed["_fit"]
    mechanism = developed["_mechanism"]
    twins = replicas(larger, at, sets, level=level, seed=seed)
    jobs = Lazy(lambda i: (problem.with_records(twins[i].records), mechanism, table_of, max_iter), sets)
    rows = [r for r in _map(_one_twin_of_the_larger, jobs, workers) if r["ended"]]
    out = {"the_planted_candidate": developed["name"], "sets": sets, "read": len(rows), "level": level}
    if len(rows) < 2:
        return out
    names = list(larger.shared)
    est, bars = np.array([r["larger"] for r in rows]), np.array([r["bars"] for r in rows])
    work = np.array([r["working"] for r in rows])
    gains = np.array([r["gain"] for r in rows])
    scatter = est.std(axis=0, ddof=1)
    out["what_the_whole_procedure_returns"] = {
        n: {"put_in": float(at.shared[j]), "mean": float(est[:, j].mean()), "bias_in_scatters": float((est[:, j].mean() - at.shared[j]) / scatter[j]),
            "scatter": float(scatter[j]), "mean_bar_over_scatter": float(bars[:, j].mean() / scatter[j]),
            "inside_one_bar": float(np.mean(np.abs(est[:, j] - at.shared[j]) < bars[:, j]))} for j, n in enumerate(names)}
    out["what_the_working_model_shows"] = {
        n: {"on_the_records": float(result.shared[j]), "mean_on_these_sets": float(work[:, j].mean()), "scatter": float(work[:, j].std(ddof=1)),
            "the_records_in_scatters_from_that_mean": float((result.shared[j] - work[:, j].mean()) / work[:, j].std(ddof=1))} for j, n in enumerate(problem.shared)}
    recorded = developed["criteria"]["gain_in_chi_square"]
    out["gain"] = {"on_the_records": float(recorded), "mean": float(gains.mean()), "scatter": float(gains.std(ddof=1)),
                   "share_of_sets_with_a_smaller_gain": float(np.mean(gains < recorded))}
    if "seen" in rows[0]:
        fam = {}
        for r in rows:
            if r["seen"]:
                fam[r["family"]] = fam.get(r["family"], 0) + 1
        out["share_of_sets_in_which_something_is_seen"] = float(np.mean([r["seen"] for r in rows]))
        out["family_put_first_when_something_is_seen"] = dict(sorted(fam.items(), key=lambda kv: -kv[1]))
    return out


def _one_twin_of_the_held(args) -> dict:
    working, mechanism, centre, table_of, max_iter = args
    out = {"ended": False}
    try:
        base = fit(working, max_iter=max_iter)
        if not base.converged:
            return out
        tab = table_of(Linear(working, base)) if table_of is not None else None
        got = fit(mechanism.hold(working, centre), base.shared, base.own, max_iter=max_iter)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return out
    if not got.converged:
        return out
    out.update({"ended": True, "working": base.shared.tolist(), "larger": got.shared.tolist(), "bars": got.sd.tolist(), "gain": float(base.chi2 - got.chi2)})
    if tab is not None:
        ranked = [c for c in tab["candidates"] if c["status"] == "ranked"]
        out.update({"seen": bool(tab["leading_set"]), "first": ranked[0]["name"] if ranked else None, "family": ranked[0]["family"] if ranked else None})
    return out


def held_on_the_twin(problem, result, held: dict, sets: int = 30, level: float = 1.0, seed: int = 6, table_of=None, max_iter: int = 80, workers: int = 1,
                     nuisances=None) -> dict:
    """A term that is held at the size of its theory, put into synthetic records. `held` is what `carried` gives.

    A quantity that is given with a bar is drawn in the synthetic records and held in the fit. So each set holds the
    term at a size of its own, drawn from the normal law of the theory, its centre and its bar, inside the bounds of the
    term. Each set is then read as the records were: the working model fitted, its table of candidates when
    `table_of(linear)` is given, the model with the term held at the centre fitted. What comes back is what
    `on_the_twin` gives, with one thing more for each shared estimate: its scatter over these sets, which holds the
    noise and the bar of the theory together, over the bar of the fit, which holds the noise alone. Where that ratio is
    above one the bar of the fit is too short for a term that is held, and the scatter is the bar to quote.

    The sets are drawn at the fit of the model with the term held at its centre. What it cannot show: whether the
    centre of the theory is right. A term held at a wrong centre biases every set alike, and nothing here sees it.

    Two things to know of the draw. A size drawn outside the bounds of the term is put at the bound, so where the
    centre is within two or three bars of a bound the sizes pile up there and their scatter is under the bar of the
    theory: a size of centre 4 and bar 2 with a bound at 0 has 2.3 per cent of its draws at 0 and a scatter of 0.98 of
    its bar. And the size of a set is drawn from the generator of that set before its noise, so the noise of set i here
    is not the noise of set i of `on_the_twin` or of `admission` at the same seed: the runs are not paired."""
    mechanism, at, larger = held["_mechanism"], held["_fit"], held["_problem"]
    centre, bar = mechanism.theory

    def model_of(rng):
        return mechanism.hold(problem, float(np.clip(centre + bar * rng.standard_normal(), mechanism.lower, mechanism.upper))).model

    twins = replicas(larger, at, sets, level=level, seed=seed, nuisances=nuisances, model_of=model_of)
    jobs = Lazy(lambda i: (problem.with_records(twins[i].records), mechanism, centre, table_of, max_iter), sets)
    rows = [r for r in _map(_one_twin_of_the_held, jobs, workers) if r["ended"]]
    out = {"the_planted_candidate": held["name"], "held_at": float(centre), "its_size_is_drawn_with_a_bar_of": float(bar), "sets": sets, "read": len(rows), "level": level}
    if len(rows) < 2:
        return out
    names = list(problem.shared)
    est, bars = np.array([r["larger"] for r in rows]), np.array([r["bars"] for r in rows])
    work = np.array([r["working"] for r in rows])
    gains = np.array([r["gain"] for r in rows])
    scatter = est.std(axis=0, ddof=1)
    out["what_the_whole_procedure_returns"] = {
        n: {"put_in": float(at.shared[j]), "mean": float(est[:, j].mean()), "bias_in_scatters": float((est[:, j].mean() - at.shared[j]) / scatter[j]),
            "scatter": float(scatter[j]), "mean_bar_over_scatter": float(bars[:, j].mean() / scatter[j]),
            "scatter_over_the_bar_of_the_fit": float(scatter[j] / (at.sd[j] * level)),
            "inside_one_bar": float(np.mean(np.abs(est[:, j] - at.shared[j]) < bars[:, j]))} for j, n in enumerate(names)}
    out["what_the_working_model_shows"] = {
        n: {"on_the_records": float(result.shared[j]), "mean_on_these_sets": float(work[:, j].mean()), "scatter": float(work[:, j].std(ddof=1)),
            "the_records_in_scatters_from_that_mean": float((result.shared[j] - work[:, j].mean()) / work[:, j].std(ddof=1))} for j, n in enumerate(names)}
    recorded = held["gain_in_chi_square"]
    out["gain"] = {"on_the_records": float(recorded), "mean": float(gains.mean()), "scatter": float(gains.std(ddof=1)),
                   "share_of_sets_with_a_smaller_gain": float(np.mean(gains < recorded))}
    if "seen" in rows[0]:
        fam = {}
        for r in rows:
            if r["seen"]:
                fam[r["family"]] = fam.get(r["family"], 0) + 1
        out["share_of_sets_in_which_something_is_seen"] = float(np.mean([r["seen"] for r in rows]))
        out["family_put_first_when_something_is_seen"] = dict(sorted(fam.items(), key=lambda kv: -kv[1]))
    return out


def carried(problem, result, mechanism, max_iter: int = 80) -> dict:
    """A term that stands on theory and is held, not fitted: the working model with the term at the centre that the
    theory gives, and again one bar of the theory below and above it. The move of each shared estimate between those
    two is what the bar of the theory costs it: a systematic bar, half of that move, to be quoted beside the bar of the
    fit and never folded into it without a word. `mechanism.theory` gives the centre and the bar."""
    if mechanism.theory is None:
        raise ValueError(f"{mechanism.name}: a held term needs the centre and the bar that the theory gives, as `theory=(centre, bar)`")
    centre, bar = mechanism.theory
    low, high = max(centre - bar, mechanism.lower), min(centre + bar, mechanism.upper)     # a held value stays where the term is defined
    fits = {}
    for label, value in (("centre", centre), ("minus", low), ("plus", high)):
        held = mechanism.hold(problem, value)
        fits[label] = (held, fit(held, result.shared, result.own, max_iter=max_iter))
    at = fits["centre"][1]
    out = {"name": mechanism.name, "ground": mechanism.ground, "held_at": centre, "bar_of_the_theory": bar, "ended": all(f.converged for _, f in fits.values()),
           "held_between": (low, high), "interval_cut_at_a_bound": bool(low > centre - bar or high < centre + bar),
           "gain_in_chi_square": float(result.chi2 - at.chi2), "chi_square": at.chi2,
           "what_it_does_to_each_parameter": {
               name: {"before": float(result.shared[j]), "with_the_term": float(at.shared[j]), "bar_of_the_fit": float(at.sd[j]),
                      "moved_in_bars_of_before": float((at.shared[j] - result.shared[j]) / result.sd[j]),
                      "systematic_bar": float(0.5 * abs(fits["plus"][1].shared[j] - fits["minus"][1].shared[j]))} for j, name in enumerate(problem.shared)}}
    out["_fit"], out["_problem"], out["_mechanism"] = at, fits["centre"][0], mechanism
    return out


def decide(developed: list, law: dict | None = None, problem=None, result=None, max_iter: int = 80) -> dict:
    """The default choice, with what it rests on. The gain of each candidate is set against the law of the gains under a
    right model when one was simulated: the share of those sets on which the best candidate gained as much.

    A term on theory that could not be fitted is added held (criteria.choose). With `problem` and `result` given the
    held fit is made here, by `carried`: the choice then holds, under "held", the working model with the term at the
    size of the theory and what the bar of the theory costs each estimate. A held fit that fails is said, not hidden."""
    # The law is that of a list of candidates. It may be the list of every candidate here, or the list of the guesses alone: a term on theory is added
    # whatever a criterion says, so its gain is set against no law, and a stage simulates the law of its guesses
    every, guessed = sorted(d["name"] for d in developed), sorted(d["name"] for d in developed if not d.get("ground"))
    if law and law.get("names") is not None and sorted(law["names"]) not in (every, guessed):
        raise ValueError("REFUSED: the law of the gains was simulated for another list of candidates (" + ", ".join(law["names"]) + "), and the best gain of a list holds for "
                         "that list alone. Simulate it again for these, or for the guesses among them: " + ", ".join(d["name"] for d in developed))
    if law and law.get("working_fit") is not None and any(d.get("working_fit") not in (None, law["working_fit"]) for d in developed):
        raise ValueError("REFUSED: the law of the gains was simulated from other records, from these records with other weights, or from another working fit of them. "
                         "A law holds for the fit that it was simulated from: simulate it again from this one")
    if law and law.get("noise_level") not in (None, 1.0):
        raise ValueError(f"REFUSED: the law of the gains was simulated at {law['noise_level']:g} times the noise of the records, and the gains of the records are those of "
                         "their own noise. A choice is made against a law simulated at level 1")
    best = np.max(np.array(law["gains"], dtype=float), axis=1) if law and law.get("gains") else None
    if best is not None:
        for d in developed:
            if d.get("ended") and d["name"] in law.get("names", [d["name"]]):        # a level in the family is that of a member of the family
                d["level_in_the_family"] = level_of(d["criteria"]["gain_in_chi_square"], best)
    out = choose(developed, best)
    if out.get("how") == "held" and problem is not None and result is not None:
        chosen = next(d for d in developed if d["name"] == out["chosen"])
        try:
            out["held"] = carried(problem, result, chosen["_mechanism"], max_iter)
        except (np.linalg.LinAlgError, FloatingPointError, ValueError) as err:          # a model that refuses a point, as the toy's line does where a width overflows, fails this fit and not the loop (his message 156)
            out["held"] = {"name": out["chosen"], "ended": False, "why_not": "the fit with the term held failed: " + str(err)}
        if not out["held"].get("ended"):               # a term that could be neither fitted nor held is not added: the choice says so, and names it
            out["terms_on_theory_not_added_in_this_stage"] = [{"name": out["chosen"], "why": "it stands on theory, and the fit with it held at the size of the theory did not end: "
                                                              "nothing is added by default in this stage, and the researcher decides"}] + out["terms_on_theory_not_added_in_this_stage"]
            out.update(chosen=None, how=None, why="the term on theory could be neither fitted nor held, and no guess is put before it")
    if law:
        out["the_law_of_the_gains"] = {k: law.get(k) for k in ("sets_asked", "sets_that_ended", "sets_dropped", "noise_level", "threshold", "smallest_level_that_these_sets_allow",
                                                              "standard_error_of_the_level", "largest_gain_seen_under_a_right_model")}
    return out


def table(developed: list) -> list:
    """The developed candidates as a table that can be written to a file: without the fits and the models that each holds."""
    return [{k: v for k, v in d.items() if not k.startswith("_")} for d in developed]


def held_lines(held: dict | None) -> list:
    """What a term held at the size of its theory does to each estimate, with the bar that the theory's own bar adds."""
    if not held:
        return []
    if not held.get("ended"):
        return [f"The fit with {held['name']} held did not end: {held.get('why_not', 'one of its three searches stopped')}. The term is not added, and the stage says so."]
    change = -held["gain_in_chi_square"]
    out = [f"With {held['name']} held at {held['held_at']:g}, where the theory gives a bar of {held['bar_of_the_theory']:g}: the chi square changes by {change:+.1f} "
           f"against the working model, which is the model with the term absent. For one parameter that is {np.sqrt(abs(change)):.1f} standard deviations "
           + ("against the term at that size." if change > 0 else "for the term at that size.")
           + (f" The interval of the theory is cut at a bound of the term, to {held['held_between'][0]:g} and {held['held_between'][1]:g}: the systematic bars are those of the cut interval."
              if held.get("interval_cut_at_a_bound") else "")]
    for name, v in held["what_it_does_to_each_parameter"].items():
        out.append(f"  {name} = {v['with_the_term']:.5g}, with a bar of the fit of {v['bar_of_the_fit']:.2g} and a systematic bar of {v['systematic_bar']:.2g} from the bar of the theory. "
                   f"Before: {v['before']:.5g}.")
    return out


def lines(developed: list, choice: dict | None = None) -> str:
    """The table of the developed candidates, for a person."""
    rows = ["candidate | its parameter | gain | AIC | BIC (samples) | BIC (records) | level in the family | what it does to each parameter, in bars",
            "---|---|---|---|---|---|---|---"]
    for d in developed:
        if not d.get("ended"):
            rows.append(f"{d['name']} | not fitted: {d.get('why_not', 'the search did not end')} | | | | | |")
            continue
        p, c = d["parameter"], d["criteria"]
        moves = ", ".join(f"{k} {v['moved_in_bars_of_before']:+.1f}" for k, v in d["what_it_does_to_each_parameter"].items())
        level = f"{d['level_in_the_family']:.3f}" if "level_in_the_family" in d else "not simulated"
        ground = " (on theory)" if d.get("ground") else ""
        value = f"{p['name']} = {p['value']:.4g}, at its bound: the term is absent or as large as it may be" if d.get("at_its_bound") else f"{p['name']} = {p['value']:.4g} with a bar of {p['bar']:.2g}"
        by_records = f"{c['BIC_n_the_records']:+.1f}" if np.isfinite(c["BIC_n_the_records"]) else "one record: none"
        rows.append(f"{d['name']}{ground} | {value} | {c['gain_in_chi_square']:.1f} | {c['AIC']:+.1f} | {c['BIC_n_the_samples']:+.1f} | "
                    f"{by_records} | {level} | {moves}")
    if choice is not None:
        rows.append("")
        rows.append(f"By default: {choice['chosen'] or 'no term is added'}{', held at the size that the theory gives' if choice.get('how') == 'held' else ''}. {choice['why'][0].upper() + choice['why'][1:]}.")
        rows += held_lines(choice.get("held"))
        for other in choice.get("left_out", []):
            rows.append(f"Left out of the choice: {other['name']}. {other['why'][0].upper() + other['why'][1:]}.")
        for other in choice.get("terms_on_theory_not_added_in_this_stage", []):
            rows.append(f"On theory and not added in this stage: {other['name']}. {other['why'][0].upper() + other['why'][1:]}.")
    return "\n".join(rows)
