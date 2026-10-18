"""The output of the package, the owner's message 140: the best model after n iterations, its terms, its parameters with
proper care for their uncertainty, and what was tried and refused.

    out = best(iterations, sets=100)                         the output, from the models fitted in the iterations
    parts = uncertainty(problem, result, sets=100)           the five parts of the uncertainty of each shared parameter

An iteration is one model fitted on the records, a dict:

    "model"     its name
    "problem"   the model, and "result" its fit
    "added"     the term that it adds to the model that it extends, and for the first the terms that it starts with, or None
    "after"     the name of the model that it extends, by default that of the iteration before
    "tried"     the terms tried in the iteration, and "refused" each term that was refused, as {"term": ..., "why": ...}
    "held"      the values of theory that the model holds, {name: (centre, bar, at)}: `at(value)` is the model with
                that value held in place of the centre, as `refine.carried` holds a term

The best model is that of the smallest Bayesian criterion, the chi square plus the number of fitted parameters times
the logarithm of the number of samples. The models whose criterion is within `TIED` of it are those that the records
cannot tell apart from it. Every iteration must be fitted on the same records: a criterion compares models, not records.

The five parts of the uncertainty of each shared parameter, never folded into one number without a word:

    the bar of the fit                 the standard deviation of the fit, at first order
    the coverage on the twin           the share of synthetic sets of the fitted model on which the interval of 1.96
                                       bars, which holds the truth on 95 sets of 100 under a right model, holds the
                                       value that made them, with the bias of the estimate on those sets in bars of the
                                       fit, which the output gives beside the value and does not take out of it. Until
                                       0.0.4.dev9 this line said that it took it out, and no version did
    the bar of theory                  half the move of the parameter between the fits with each held value of theory one
                                       bar of the theory below and above its centre, for each held value and together
    the spread over the models         half the range of the parameter over the best model and those that the records
                                       cannot tell apart from it
    not told apart                     a parameter correlated with another beyond `TOLD_APART` gets no bar: the records
                                       fix a combination of the two and not each, and the output names the other

The package knows nothing of the field of the records. Which term of physical origin an iteration adds is the choice of
the researcher, or of the agent in the researcher's place.
"""
from __future__ import annotations

import math

import numpy as np

from .fit import fit
from .twin import admission

TIED = 2.0                  # two models whose Bayesian criteria differ by less than this are not told apart by the records
TOLD_APART = 0.999          # two parameters correlated beyond this get no bar


def theory_moves(result, held: dict | None, max_iter: int = 80) -> dict:
    """For each held value, {name: (centre, bar, at)}, half the move of each shared parameter between the fits made again
    at the value less and more one bar, from the fit `result`: the bar that the value puts on each. Each held value is
    put back at its centre after."""
    moved = {}
    for name, (centre, bar, at) in (held or {}).items():
        ends = []
        for value in (centre - bar, centre + bar):
            got = fit(at(value), result.shared, result.own, max_iter=max_iter)
            if not got.converged:
                at(centre)
                raise ValueError(f"the fit with {name} held at {value:g}, one bar of the theory from its centre, did not end: its bar cannot be carried")
            ends.append(got.shared)
        at(centre)
        moved[name] = 0.5 * np.abs(ends[1] - ends[0])
    return moved


def correlation(cov: np.ndarray, sd: np.ndarray) -> np.ndarray:
    """The correlations of the parameters whose standard deviation is a positive number, and none of the others: one
    parameter held at a bound, its deviation not a number, made the whole matrix the identity, so a pair at -0.9999 got
    two bars (the red team of the package, his message 181)."""
    sd = np.asarray(sd, dtype=float)
    out = np.eye(len(sd))
    free = np.isfinite(sd) & (sd > 0)
    if free.any():
        out[np.ix_(free, free)] = np.asarray(cov, dtype=float)[np.ix_(free, free)] / np.outer(sd[free], sd[free])
    return out


def uncertainty(problem, result, sets: int = 100, held: dict | None = None, alike=(), seed: int = 31, start: str = "problem", max_iter: int = 80, workers: int = 1) -> dict:
    """The five parts of the uncertainty of each shared parameter of a fit. `held`: the values of theory that the model
    holds, {name: (centre, bar, at)}. `alike`: the models that the records cannot tell apart from this one, each as
    (name, the names of its shared parameters, their values). `sets`, `seed`, `start`, `max_iter` and `workers` go to
    the admission of the twin, which needs 30 sets that end."""
    sd = result.sd
    corr = correlation(result.cov, sd)
    twin = admission(problem, result, sets=sets, seed=seed, start=start, max_iter=max_iter, workers=workers)
    moved = theory_moves(result, held, max_iter)
    out = {"admitted_on_the_twin": twin["admitted"], "rule_of_the_twin": twin["rule"], "sets_of_the_twin": twin["sets"],
           "searches_that_ended": int(round(twin["ended"] * twin["sets"])), "parameters": {}}          # 27 of 30, which best.md did not say (round 18 of the cold tests)
    for j, name in enumerate(problem.shared):
        partners = [problem.shared[k] for k in range(len(sd)) if k != j and abs(corr[j, k]) > TOLD_APART]
        cell = twin["parameters"].get(name) or {}
        bias = cell.get("bias_over_spread")
        bias = None if bias is None or cell.get("spread_over_bound") is None else float(bias * cell["spread_over_bound"])
        values = [float(result.shared[j])] + [float(dict(zip(names, vals))[name]) for _, names, vals in alike if name in names]
        theory = {k: float(v[j]) for k, v in moved.items()}
        out["parameters"][name] = {
            "value": float(result.shared[j]),
            "bar_of_the_fit": None if partners else float(sd[j]),
            "not_told_apart_from": partners,
            "coverage_on_the_twin": cell.get("coverage"),
            "bias_on_the_twin_in_bars": bias,
            "value_less_the_bias_of_the_twin": None if bias is None or partners else float(result.shared[j] - bias * sd[j]),
            "bar_of_theory": {"together": float(math.sqrt(sum(v * v for v in theory.values()))), "of_each": theory},
            "spread_over_the_models": float(0.5 * (max(values) - min(values))),
            "models_not_told_apart": [m for m, names, _ in alike if name in names],
            "absent_from": [m for m, names, _ in alike if name not in names]}          # a size of a term that a model not told apart lacks: the spread of 0 read as agreement (round 16 of the cold tests)
    return out


def own_bars(problem, result) -> list:
    """The own parameters of each record with their bars: for each record, {name: (value, bar)}. The bar is marginal over
    the shared parameters: the inverse of the record's own block, M, widened by what the shared ones leave uncertain,
    M^-1 + M^-1 B'A C A'B M^-1, with A and B the whitened columns of the shared and the own parameters and C the
    covariance of the shared ones. A record whose own block cannot be inverted gives no bar. A researcher of eclipses
    wanted the centre of each curve, the time of each eclipse, and no output held it (the cold tests of row M169)."""
    from .fit import derivatives
    out = []
    C = np.asarray(result.cov, dtype=float)
    for k in range(len(problem.records)):
        _, A, B = derivatives(problem, k, np.asarray(result.shared, dtype=float), np.asarray(result.own[k], dtype=float))
        try:
            Mi = np.linalg.inv(B.T @ B)
        except np.linalg.LinAlgError:
            out.append({name: (float(v), None) for name, v in zip(problem.own, result.own[k])})
            continue
        cov = Mi + (Mi @ B.T @ A @ C @ A.T @ B @ Mi if C.size else 0.0)
        out.append({name: (float(v), float(np.sqrt(cov[j, j])) if cov[j, j] > 0 else None) for j, (name, v) in enumerate(zip(problem.own, result.own[k]))})
    return out


def best(iterations: list, sets: int = 100, tied: float = TIED, seed: int = 31, start: str = "problem", max_iter: int = 80, workers: int = 1) -> dict:
    """The output after the iterations: the best model, its terms, its parameters with the five parts of their
    uncertainty, the models that the records cannot tell apart from it, and what each iteration tried and refused."""
    if not iterations:
        raise ValueError("no iteration was given: the output is that of at least one fitted model")
    first = iterations[0]["problem"]
    for it in iterations[1:]:
        p = it["problem"]
        if len(p.records) != len(first.records) or any(a.x.shape != b.x.shape or not (np.array_equal(a.x, b.x) and np.array_equal(a.y, b.y)) for a, b in zip(p.records, first.records)):
            raise ValueError(f"the model {it['model']} is fitted on other records than {iterations[0]['model']}: a criterion compares models on the same records")
    names = [it["model"] for it in iterations]
    if len(set(names)) != len(names):
        raise ValueError("two iterations have one name: the terms of a model are found through the names of the models that it extends")
    rows = []
    for i, it in enumerate(iterations):
        p, r = it["problem"], it["result"]
        fitted = p.samples - r.dof                  # the shared and the own parameters that the fit took out of the residuals
        rows.append({"model": it["model"], "added": it.get("added"), "after": it.get("after", names[i - 1] if i else None), "chi_square": float(r.chi2),
                     "fitted_parameters": int(fitted), "criterion": float(r.chi2 + fitted * math.log(p.samples)), "ended": bool(r.converged),
                     "tried": list(it.get("tried", ())), "refused": list(it.get("refused", ())), "aside": bool(it.get("aside"))})
    ended = [row for row in rows if row["ended"]]
    if not ended:
        raise ValueError("no fit of the iterations ended: there is no best model")
    top = min([row for row in ended if not row["aside"]] or ended, key=lambda row: row["criterion"])          # a model set aside, a refused term in the place of a kept one, is never the best by a hair: it is held in the spread, or named where it is better by more than tied
    for row in rows:
        row["criterion_less_that_of_the_best"] = row["criterion"] - top["criterion"]
    alike = [row["model"] for row in ended if row is not top and abs(row["criterion"] - top["criterion"]) < tied]
    better = [row["model"] for row in ended if row["aside"] and top["criterion"] - row["criterion"] >= tied]
    terms, at = [], top
    while at is not None:
        if at["added"]:
            terms.insert(0, at["added"])
        at = next((row for row in rows if row["model"] == at["after"]), None) if at["after"] else None
    chosen = iterations[rows.index(top)]
    others = [(it["model"], it["problem"].shared, it["result"].shared) for it in iterations if it["model"] in alike]
    parts = uncertainty(chosen["problem"], chosen["result"], sets=sets, held=chosen.get("held"), alike=others, seed=seed, start=start, max_iter=max_iter, workers=workers)
    return {"best": top["model"], "terms": terms, "models_not_told_apart": alike, "better_set_aside": better, "iterations": rows,
            "tried_and_refused": [dict(refused, iteration=row["model"]) for row in rows for refused in row["refused"]], **parts}
