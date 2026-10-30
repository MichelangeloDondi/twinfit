"""The analysis of the package, the same on the synthetic and on the recorded traces (the owner's message 119: "apply
always (for each tier of the model) the same analysis (of the generalized moments, chi functions, windowing/gridding,
bias map, etc) to the real traces and to the synthetic traces, and to compare them").

`of` reads one set of 8 traces under one model, through the joint fit of all 8 with the weights of the noise law:

    families    each family of statistics of the package: moments of the line under several kernels, windows in frequency
                and in signal level, the two sides of the line added and subtracted, the frequency at given levels,
                moments of the derivative, of the integral, of the logarithm and of the reciprocal. For each: what it
                reads in what the fit leaves, as a chi square with its degrees of freedom and its level, and its own
                estimate of every shared parameter with its first order uncertainty
    costs       each chi function: the fit under a cost of a residual that is not its square. Its estimate, its first
                order uncertainty, and what its residuals show on the plane of frequency and residual value
    windows     each pair of windows at equal distances from the centre of a line, the ranges in a geometric progression
                of the width: what the fit leaves there, both sides together and one side less the other
    usual       what an analyst runs with no method of the package: the neighbours, the runs of signs, equal bins
    seen        on synthetic traces, where the true line is known: the chance that each family reports what the model
                lacks, at first order, with no noise drawn
    combination for each term that the model lacks and the traces see, and for each shared parameter: the estimates of
                several families combined, with weights that add to one, so that the bias of that term cancels at first
                order (the owner's message 88: "generalised moments (and their bias free combinations)"). The model
                still lacks the term, and the estimate no longer carries its bias. `on_draws` checks it on many draws
    guess       what the package says that the model lacks, and how large it is (the owner's message 120: "to let a
                claude code agent to understand/guess blindly which model terms are missing and how big they are, going
                much beyond what the chi^2 can do"). Two lists, from one table of the package. With no hint: the generic
                terms of its dictionary, such as a width that changes along the temperature, a pedestal, a second
                component beside the centre, a tail on one side. With the names of the terms as hints: each term of the
                toy that the model does not free, at its size of config.toml, and the size that the residuals ask of it

`every` runs it for a subfolder: both sets of traces, the model with the new term and the model without it. The same
function reads both sets, so a difference between them is a difference between the traces and not between two programs.
The uncertainties are first order ones, and on the recorded traces they assume a right model.
"""
from __future__ import annotations

import numpy as np

import twinfit
from twinfit import families, losses, noise as N, propose, readings, standard, unsuspected
from twinfit.candidates import Candidate

from . import fits, model

SHORT = {"moments: powers 0 to 4": "moments, powers 0 to 4",
         "moments: powers 0 to 8": "moments, powers 0 to 8",
         "moments: fractional powers of the distance": "moments, fractional powers",
         "moments: inverse and logarithmic kernels": "moments, inverse and logarithm",
         "moments: cosines and sines": "moments, cosines and sines",
         "windows along the abscissa: pairs of intervals, the two sides added and subtracted": "windows in frequency",
         "the two sides added: pairs of intervals, the centre held": "windows, the two sides added",
         "the two sides subtracted: pairs of intervals": "windows, one side less the other",
         "windows along the value: bands between levels of the height, each side": "windows in signal level",
         "the abscissa against the value: the crossings of nine levels, each side": "frequency at nine levels",
         "moments of the derivative: powers 0 to 8": "moments of the derivative",
         "moments of the running integral: powers 0 to 8": "moments of the integral",
         "ratios of moments of the distance: orders 0.5 to 4 over order 0": "ratios of moments",
         "the logarithm of the value: power kernels 0 to 8": "moments of the logarithm",
         "the reciprocal of the value: power kernels 0 to 8": "moments of the reciprocal",
         "every sample: the fit": "chi square of the fit"}          # the families of the package, as a page names them
COSTS = {"square": "least squares", "1.5": "power 1.5", "huber": "Huber", "logarithm": "logarithm", "student": "Student"}          # the chi functions, as a page names them
CENTRE = "centre"          # the own parameter of a trace that places its line
LEVEL = 0.05          # a reading reports when its level is under this


def at(detail: int, rows: list, const: dict, config: dict, free: list | None = None) -> tuple:
    """The joint problem of some traces with the weights of the noise law, its fit, and the linear picture of the two."""
    problem = fits.problem_of(detail, rows, const, config, free)
    settled, result, _ = N.settle(problem, fits.variance_of(const))
    return settled, result, twinfit.Linear(settled, result)


def told_apart(detail: int, rows: list, const: dict, config: dict, free: list | None = None) -> tuple:
    """The free parameters of a model that some traces tell apart, and those that they do not. Where two free parameters
    do the same thing to every trace, as two widths over traces of one temperature, their correlation is beyond
    `fits.TIED` and no reading of that fit means anything: the one that enters the model later is then held at its value
    of config.toml, and the fit is made again, until no two are tied. Returns the free parameters that stay, the ones
    that were held, and for each held one the free parameter that it was tied with: that one then carries whatever the
    held value leaves, and its estimate is no estimate of its own term (review 8, moderate 2.1)."""
    free, held, tied = list(config["detail"][str(detail)]["free"] if free is None else free), [], {}
    while len(free) > 1:
        lin = at(detail, rows, const, config, free)[2]
        cov = np.linalg.pinv(lin.J.T @ lin.J)
        bars = np.sqrt(np.clip(np.diag(cov), 0.0, None))
        corr = np.abs(cov / np.where(np.outer(bars, bars) > 0.0, np.outer(bars, bars), np.inf))
        np.fill_diagonal(corr, 0.0)
        i, j = np.unravel_index(int(np.argmax(corr)), corr.shape)
        if not corr[i, j] > fits.TIED:
            break
        later = max((free[i], free[j]), key=lambda name: config["parameter"][name]["enters"])
        tied[later] = free[j] if later == free[i] else free[i]
        free.remove(later)
        held.append(later)
    return free, held, tied


def of_the_families(lin, lacks: list | None = None, names: list | None = None) -> list:
    """What each family of the package reads at a fit, and its own estimate. `lacks`: a discrepancy of the mean, one array
    a trace, for the chance that each family reports it. `names`: the families to read, all of them when None."""
    shared, out = list(lin.problem.shared), []
    for name in (names or list(SHORT)):
        fam = families.family(lin, name, CENTRE)
        reading, bounds = fam.reading(), fam.bounds()
        row = {"name": SHORT[name], "statistics": int(fam.size), "chi2": float(reading["chi_square"]), "dof": int(reading["degrees_of_freedom"]), "level": float(reading["level"]),
               "reports": bool(reading["level"] < LEVEL), "estimator": bool(bounds["estimator"]) and name != "every sample: the fit", "ended": False, "values": None, "bars": None,
               "seen": None if lacks is None else float(fam.power(lacks, LEVEL))}
        if row["estimator"]:
            got = fam.estimate()
            row.update(ended=bool(got.converged))
            if got.converged:
                row.update(values=dict(zip(shared, (float(v) for v in got.shared))), bars=dict(zip(shared, (float(v) for v in got.sd))))
        out.append(row)
    return out


def of_the_costs(problem, result, costs: list | None = None) -> list:
    """The fit under each chi function: its estimate, its first order uncertainty, which is its factor times that of
    the weighted fit, and the table of its residuals on the plane of frequency and residual value."""
    shared, out = list(problem.shared), []
    for cost in (costs or list(COSTS)):
        if not shared and cost != "square":          # a chi function gives another estimate of the shared parameters: a model with none has no other fit than the weighted one
            continue
        got = losses.chi_fit(problem, result, float(cost) if cost.replace(".", "").isdigit() else cost)
        row = {"name": COSTS[cost], "ended": bool(got.converged and got.fit.converged), "values": None, "bars": None, "factor": None, "plane": None}
        if row["ended"]:
            plane = losses.grid(problem, got.fit, got.scale)
            row.update(values=dict(zip(shared, (float(v) for v in got.fit.shared))), bars=dict(zip(shared, (float(v) for v in np.atleast_1d(got.sd)))), factor=float(got.factor),
                       plane={"chi2": float(plane["chi_square"]), "dof": int(plane["degrees_of_freedom"]), "level": float(plane["level"]), "reports": bool(plane["level"] < LEVEL),
                              "cells": [[float(v) for v in line] for line in plane["standardised_cells"]]})
        out.append(row)
    return out


def of_the_windows(lin) -> list:
    """What the fit leaves in each pair of windows about the centre of every line, the ranges in a geometric progression
    of the width of the line: both sides together, and one side less the other. Each is the root of a chi square with
    one degree of freedom, with the sign of the misfit: positive where the traces lie above the fitted line."""
    cuts = readings.slices(lin)          # in their order: for each range, both sides together, then one side less the other
    table = {c["name"]: c for c in propose.propose(lin, cuts)["candidates"]}          # the table ranks them: each is found again by its name
    edges, out = list(readings.GEOMETRIC), []
    for k, cut in enumerate(cuts):
        c = table[cut.name]
        scored = c.get("status") == "ranked"
        out.append({"from_widths": edges[k // 2], "to_widths": edges[k // 2 + 1], "sides": "together" if cut.parity == "even" else "one less the other",
                    "z": float(np.sign(c["amplitude"][0]) * np.sqrt(c["score"])) if scored else None, "level": float(c["level_on_its_own"]) if scored else None})
    return out


TOP = 8          # how many generic candidates a reading keeps, the best first
TERM = "a term of the model"          # the family of the candidates that are terms of the toy
HIDDEN = 0.05          # a term of which the fit leaves less than this share to be seen is taken by the model: no size is read from what is left of it, since a large misfit of another shape would be read as its size


CONVERGED = 1e-3          # the effect of a term at first order: its size is halved until the shape of what it adds to the lines changes by less than this share of its largest value
HALVINGS = 40          # and never more often than this: a term whose effect has no limit of small size is refused


def of_the_terms(problem, result, const: dict, config: dict, small: bool = False, detail: int | None = None) -> list:
    """The terms of the toy that the model does not hold, each as a candidate: what the fitted line of every trace gains
    when the term is added at its size of config.toml, the other parameters as fitted. Every term of the toy is absent
    at zero, so the line with one more term is the line of the last detail with the others at zero. A term that the
    model holds at a value is not in the list: a free parameter of the same kind takes it whole. Returns, for each
    term, its name, the size that was added and the candidate. With `small`, the candidate is the effect of the term at
    first order: the size is halved until the shape of what the term adds, over all the traces, changes by less than
    `CONVERGED` of its largest value from one halving to the next. That limit does not follow the size of config.toml,
    and no size is chosen: a threshold on the height of the change, which stood here for one candidate, decided the
    result (review 7, moderate 1.2). A combination of estimates is made from this one: its weights follow the shape of
    the effect, and with the size of config.toml they followed the size as well (review 6, moderate 1.1 and 1.2).
    `detail`: the detail of the model, where its free parameters do not say it: a model of which a sweep alone holds one
    parameter at its value still holds that term, and it is not offered as one that the model lacks."""
    last = max(int(d) for d in config["detail"])
    detail = max(config["parameter"][name]["enters"] for name in problem.shared) if detail is None else int(detail)
    inside = {name for name, e in config["parameter"].items() if e["enters"] <= detail}
    base = dict(model.truth(config), **{name: 0.0 for name in config["parameter"] if name not in inside})
    base.update(zip(problem.shared, (float(v) for v in result.shared)))

    def line(k: int, values: dict):
        rec, own = problem.records[k], result.own[k]
        u = rec.x - own[1]
        top = float(model.shape(last, 0.0, rec.design["T"], values, const))
        return own[0] * model.shape(last, u, rec.design["T"], values, const) / top + own[2] + own[3] * u
    out = []
    for name, e in config["parameter"].items():
        if name in inside:
            continue
        size = float(e["true"])

        def gained(size: float) -> list:
            more = dict(base, **{name: base[name] + size})
            return [line(k, more) - line(k, base) for k in range(len(problem.records))]
        column = gained(size)
        shape = lambda made: np.concatenate(made) / (float(np.max(np.abs(np.concatenate(made)))) or 1.0)          # noqa: E731
        for _ in range(HALVINGS if small else 0):
            half = gained(size / 2.0)
            settled = float(np.max(np.abs(shape(half) - shape(column)))) < CONVERGED
            size, column = size / 2.0, half
            if settled:
                break
        else:
            if small:
                raise SystemExit(f"REFUSED: what {name} adds to the lines has no limit of small size: its shape still changes after {HALVINGS} halvings")
        out.append((name, size, Candidate(f"{e['what'].split(':')[0]}", TERM, [column], "even")))
    return out


def of_the_guess(problem, result, lin, const: dict, config: dict, terms: list | None = None) -> dict:
    """What the package says that the model lacks, from the residuals of the fit: one table of `propose` over the generic
    shapes that the residuals hold and over the terms of the toy that the model does not free. A candidate has a score, a
    chi square with as many degrees of freedom as it has columns, whatever the fit absorbed, and its level in the
    family of all the candidates that were looked at. For a term of the toy, the amplitude of the candidate times the
    size that was added is the size that the residuals ask of it, at first order. A term that the model absorbs whole
    leaves no residual, and the table says so."""
    terms = of_the_terms(problem, result, const, config) if terms is None else terms
    bases = [c for c in unsuspected.held_by_the_residuals(lin) if c.family == "trend"]          # the one basis fixed before the residuals are read, the slow trend: a shape built from the residuals scores their own size there and put the collisions fifth behind four of them, and the best of a scan of periods led on the right model, its level paying nothing for the scan, which the twin of the report pays
    table = propose.propose(lin, [c for _, _, c in terms] + bases)          # beside the terms of the toy: the package holds no list of terms (the owner's message 340)
    rows = table["candidates"]

    def row(c: dict, rank: int) -> dict:
        seen = c["status"] == "ranked"
        return {"name": c["name"], "kind": c["family"], "sides": c["parity"], "rank": rank if seen else None, "seen": seen, "score": float(c["score"]) if seen else None,
                "dof": int(c["degrees_of_freedom"]) if seen else None, "level": float(c["family_level"]) if seen else None, "leads": bool(seen and c["name"] in table["leading_set"]),
                "amplitude": [float(v) for v in c["amplitude"]] if seen else None, "amplitude_bar": [float(v) for v in c["amplitude_standard_deviation"]] if seen else None,
                "alike": list(c.get("like", [])) if seen else [], "visible": float(c["visibility"]) if seen else 0.0}
    by_name = {c["name"]: (c, k + 1) for k, c in enumerate(rows)}
    named = []
    for name, size, cand in terms:
        got = row(*by_name[cand.name])
        told = got["seen"] and got["visible"] >= HIDDEN          # a term that a free parameter takes nearly whole has no size that these traces can tell, and does not lead
        a, bar = (got["amplitude"][0], got["amplitude_bar"][0]) if told else (None, None)
        if told and config["parameter"][name].get("effect") == "square":          # the line changes as the square of the size: the amplitude is the square of the size over the one that was added
            told = a > 0.0
            a, bar = (float(np.sqrt(a)), bar / (2.0 * float(np.sqrt(a)))) if told else (None, None)
        got.update(parameter=name, added=size, told=bool(told), leads=bool(got["leads"] and told), size=a * size if told else None, size_bar=bar * size if told else None,
                   absorbed=bool(not got["seen"] or got["visible"] < HIDDEN))          # no size for one of two reasons: the model takes the term nearly whole, or it comes out with the sign that the term cannot have
        named.append(got)
    generic = [row(c, k + 1) for k, c in enumerate(rows) if c["family"] != TERM][:TOP]
    return {"looked_at": int(table["candidates_looked_at"]), "leading": len(table["leading_set"]), "terms": named, "generic": generic}


COMBINED = ("every sample: the fit", "moments: powers 0 to 8", "moments: cosines and sines", "windows along the abscissa: pairs of intervals, the two sides added and subtracted",
            "windows along the value: bands between levels of the height, each side", "the abscissa against the value: the crossings of nine levels, each side",
            "the logarithm of the value: power kernels 0 to 8", "moments of the derivative: powers 0 to 8")          # the estimates that a combination is made of: the fit, and the families whose own search is fast
INFLATED = 10.0          # a combination whose uncertainty is more than this many times that of the fit is of no use: the families then share the bias of the term, and nothing of it cancels
DRAWS = 100          # the draws on which a combination is checked: as many as a model is checked on
BIAS = 0.25          # a combination passes with a bias under this share of its spread: the bound of a model in `fits.check`, which a test holds equal to this one


def combined(weights: dict, estimates: dict) -> float:
    """The weighted sum of the estimates, each weight with the estimate of its own name."""
    if set(weights) != set(estimates):
        raise ValueError(f"the weights and the estimates are not of the same families: {sorted(set(weights) ^ set(estimates))}")
    return float(sum(weights[name] * estimates[name] for name in weights))


def of_the_combination(lin, rows: list, terms: list, left_out: bool = False) -> list:  # noqa: C901
    """For each term given, as `of_the_terms` gives it, and for each shared parameter: the combination of the estimates of
    the families of `COMBINED` whose bias under that term is zero at first order, with the least variance among those.
    `rows` are the families as `of_the_families` read them: a family whose search did not end is left out. The weights
    and the uncertainty are the package's, with no draw. `useful` is False when the uncertainty is beyond `INFLATED`
    times that of the fit. With `left_out`, a useful combination is made again with each family but the fit left out
    in turn, and then with each further family of `rows` that gives an estimate added in turn, and the lowest and the
    highest of each kind are given: the set of families is a choice, and its uncertainty does not hold what the choice
    moves (review 6, moderate 1.3). One more family moves it more than one fewer (review 7, moderate 1.3)."""
    held = {SHORT[name]: name for name in COMBINED}
    mine = [f for f in rows if f["name"] in held and f["estimator"] and f["ended"]]
    names = [COMBINED[0]] + [held[f["name"]] for f in mine]
    shared = list(lin.problem.shared)
    bars = np.sqrt(np.diag(np.linalg.inv(lin.J.T @ lin.J)))
    out = []
    for parameter, _, cand in terms:
        for j, name in enumerate(shared):
            if len(names) < 3:
                continue          # with fewer than three estimates nothing is left to choose once the bias is cancelled
            free = families.free_of(lin, names, {parameter: cand.columns[0]}, CENTRE, j)
            estimates = dict([(COMBINED[0], float(lin.fit.shared[j]))] + [(held[f["name"]], float(f["values"][name])) for f in mine])
            over = float(free["standard_deviation_over_the_fit"])
            row = {"term": parameter, "parameter": name, "value": combined(free["weights"], estimates), "bar": over * float(bars[j]), "over_the_fit": over, "useful": bool(over <= INFLATED),
                   "fit": estimates[COMBINED[0]], "fit_bar": float(bars[j]), "estimates": len(names), "largest_weight": float(max(abs(w) for w in free["weights"].values()))}
            if left_out and row["useful"] and len(names) > 3:
                others = []
                for gone in names[1:]:
                    fewer = [n for n in names if n != gone]
                    others.append(combined(families.free_of(lin, fewer, {parameter: cand.columns[0]}, CENTRE, j)["weights"], {n: estimates[n] for n in fewer}))
                row.update(lowest_with_one_left_out=min(others), highest_with_one_left_out=max(others))
                more = []
                for f in rows:
                    one = next((n for n, short in SHORT.items() if short == f["name"]), None)
                    if one is None or one in names or not (f["estimator"] and f["ended"]):
                        continue
                    wider = names + [one]
                    more.append(combined(families.free_of(lin, wider, {parameter: cand.columns[0]}, CENTRE, j)["weights"], dict(estimates, **{one: float(f["values"][name])})))
                if more:
                    row.update(lowest_with_one_added=min(more), highest_with_one_added=max(more), added=len(more))
            out.append(row)
    return out


def one_draw(asked: tuple) -> list:
    """The combination on one draw of the synthetic traces of a detail, under the model that lacks its term: for a pool."""
    made, d, free, draw, const, config = asked
    problem, result, lin = at(d, fits.synthetic(made, const, config, draw), const, config, free)
    new = [t for t in of_the_terms(problem, result, const, config, small=True, detail=d) if config["parameter"][t[0]]["enters"] == made]
    return of_the_combination(lin, of_the_families(lin, None, list(COMBINED[1:])), new)


def of_the_draws(got: list, true: dict, draws: int) -> dict:
    """What many draws say of the fit and of the combination. `got`: for each draw that ended, the rows of
    `of_the_combination`. For each shared parameter and for each of the two estimates: its mean, its spread, its bias in
    spreads, the share of the draws whose interval of one uncertainty of its own holds the true value, the root of its
    mean square error, and whether it passes as a model has to before it goes on the recorded traces: every draw
    ended, a bias under `BIAS` of its spread, and a share of intervals that hold the truth inside the band of
    `fits.held`. A combination that is of no use does not pass, whatever its bias in its own large spread. The
    standard error of a bias in spreads is one over the root of the number of draws. The bound was
    three standard errors of the bias for one candidate, which is looser, and a check then passed on 48 draws and
    failed on 48 others (review 7, moderate 1.1). An estimate that is one and the same value in every draw has no
    spread: its bias in spreads is None and it does not pass. That is a fit on a bound of its parameter, and the
    page says so. Draws in which every estimate is one value are refused: they are one draw."""
    out = {"draws": draws, "ended": len(got), "parameters": {}}
    low, high = fits.held(len(got))
    varies = False
    for k, first in enumerate(got[0] if got else []):
        name = first["parameter"]
        said = {}
        for label, value, bar in (("fit", "fit", "fit_bar"), ("combination", "value", "bar")):
            v, b = np.array([rows[k][value] for rows in got]), np.array([rows[k][bar] for rows in got])
            spread = float(v.std(ddof=1)) if len(v) > 1 else 0.0
            varies = varies or spread > 0.0
            bias, cover = (float((v.mean() - true[name]) / spread) if spread > 0.0 else None), float(np.mean(np.abs(v - true[name]) <= b))
            said[label] = {"mean": float(v.mean()), "spread": spread, "bias_in_spreads": bias, "bias_standard_error": float(1.0 / np.sqrt(len(v))), "holds_the_truth": cover,
                           "root_mean_square_error": float(np.sqrt(np.mean((v - true[name]) ** 2))), "passes": bool(len(got) == draws and bias is not None and abs(bias) < BIAS and low <= cover <= high)}
        useful = bool(np.mean([rows[k]["useful"] for rows in got]) > 0.5)
        said["combination"]["passes"] = bool(said["combination"]["passes"] and useful)          # a combination of no use has a spread so large that any bias is small in it: it does not pass (review 8, minor 3.4)
        out["parameters"][name] = dict(said, true=true[name], term=first["term"], useful=useful)
    if got and got[0] and not varies:
        raise SystemExit(f"REFUSED: the {len(got)} draws give one and the same value of every estimate: they are not different draws")
    return out


def on_draws(made: int, d: int, free: list | None, const: dict, config: dict, draws: int = DRAWS, workers: int = 8, size: float | None = None) -> dict:
    """The combination against the fit on many draws of the 8 synthetic traces of the detail `made`, under the model of
    detail `d` with the free parameters `free`, which lacks the term that `made` brings: what `of_the_draws` says of
    them. `size`: the traces are drawn with the term of `made` at this size and not at the one of config.toml, as for
    the size that the recorded traces ask of it. The other true values stay."""
    import copy
    from concurrent.futures import ProcessPoolExecutor
    drawn = copy.deepcopy(config)
    terms = [name for name, e in drawn["parameter"].items() if e["enters"] == made]
    for name in terms if size is not None else []:
        drawn["parameter"][name]["true"] = float(size)
    with ProcessPoolExecutor(min(workers, draws)) as pool:
        got = [rows for rows in pool.map(one_draw, [(made, d, free, 1 + k, const, drawn) for k in range(draws)]) if rows]
    return dict(of_the_draws(got, model.truth(drawn), draws), sizes={name: float(drawn["parameter"][name]["true"]) for name in terms}, bias_bound=BIAS, band_of_the_share_that_holds_the_truth=list(fits.held(len(got))))


def of_the_usual(lin) -> dict:
    """The checks of the residuals that an analyst runs with no method of the package."""
    serial, runs, binned = standard.serial(lin), standard.runs(lin), standard.binned(lin)
    return {"neighbours": {"correlation": float(serial["correlation"]), "level": float(serial["level"])}, "runs": {"z": float(runs["z"]), "level": float(runs["level"])},
            "bins": {"chi2": float(binned["chi_square"]), "dof": int(binned["degrees_of_freedom"]), "level": float(binned["level"]), "reports": bool(binned["level"] < LEVEL)}}


def of(detail: int, rows: list, const: dict, config: dict, free: list | None = None, true: list | None = None, names: list | None = None, costs: list | None = None) -> dict:
    """The whole reading of one set of 8 traces under one model. `true`: the true line of each trace, where it is known,
    for the chance that each family reports what the model lacks. `names` and `costs`: the families and the chi
    functions to read, all of them when None."""
    problem, result, lin = at(detail, rows, const, config, free)
    lacks = None if true is None else [np.asarray(t, dtype=float) - problem.predict(k, result.shared, result.own[k]) for k, t in enumerate(true)]
    bars = np.sqrt(np.clip(np.diag(result.cov), 0.0, None))
    terms, read = of_the_terms(problem, result, const, config, detail=detail), of_the_families(lin, lacks, names)
    guess = of_the_guess(problem, result, lin, const, config, terms)
    told = {g["parameter"] for g in guess["terms"] if g["told"]}          # a combination is made for the terms that the traces see
    first_order = [t for t in of_the_terms(problem, result, const, config, small=True, detail=detail) if t[0] in told]
    return {"detail": detail, "free": list(problem.shared), "values": dict(zip(problem.shared, (float(v) for v in result.shared))), "bars": dict(zip(problem.shared, (float(v) for v in bars))),
            "chi2": float(result.chi2), "dof": int(result.dof), "families": read, "costs": of_the_costs(problem, result, costs),
            "windows": of_the_windows(lin), "usual": of_the_usual(lin), "guess": guess, "combination": of_the_combination(lin, read, first_order, left_out=True)}


def spread(one: dict, name: str | None = None) -> dict:
    """How far apart the estimates of one shared parameter lie, over the families and the chi functions of a reading:
    the lowest, the highest, and how many lie more than 3 of their own uncertainties from the fit by least squares."""
    name = name or one["free"][0]
    got = [(s["values"][name], s["bars"][name]) for s in one["families"] + [c for c in one["costs"] if c["name"] != COSTS["square"]] if s.get("values")]
    return {"parameter": name, "lowest": min(v for v, _ in got), "highest": max(v for v, _ in got), "far": int(sum(abs(v - one["values"][name]) > 3.0 * b for v, b in got)), "of": len(got)}


def largest(one: dict, sides: str) -> float:
    """The largest of what the fit leaves in a pair of windows, with its sign, for both sides together or for one less the other."""
    return max((w["z"] for w in one["windows"] if w["sides"] == sides and w["z"] is not None), key=abs)


def one(asked: tuple) -> tuple:
    """One reading, for a pool of processes."""
    key, detail, rows, const, config, free, true, names, costs = asked
    return key, of(detail, rows, const, config, free, true, names, costs)


def every(detail: int, which: dict, synthetic: list, recorded: list | None, const: dict, config: dict, keep=None, names: list | None = None, costs: list | None = None) -> dict:
    """The readings of a subfolder, side by side in processes: both sets of traces, under the model with the new term and
    under the model without it. `which`: the two models, each as (its detail, its free parameters). `keep(name, make)`
    reads a reading that was kept and makes it when it was not."""
    from concurrent.futures import ProcessPoolExecutor
    true = [model.line(detail, t["P"], t["T"], model.truth(config), const, config) for t in synthetic]          # the synthetic traces are those of the detail, under both models
    jobs = [((kind, word), d, rows, const, config, free, true if kind == "synthetic" else None, names, costs)
            for kind, rows in (("synthetic", synthetic), ("recorded", recorded)) if rows is not None for word, (d, free) in which.items()]
    out, todo = {"synthetic": {}, "recorded": {} if recorded is not None else None}, []
    for job in jobs:
        if keep is None:
            todo.append(job)
            continue
        got = keep(job, None)
        if got is None:
            todo.append(job)
        else:
            out[job[0][0]][job[0][1]] = got
    if todo:
        with ProcessPoolExecutor(len(todo)) as pool:
            for job, (key, got) in zip(todo, pool.map(one, todo)):
                out[key[0]][key[1]] = got if keep is None else keep(job, got)
    return out
