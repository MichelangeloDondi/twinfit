"""What can be said of the known parameters before the missing term is known.

The residuals of a fit may say that the model lacks something while nobody can name what. The orbit of Mercury did not
follow the model of the known planets, and for decades the term that was missing had no right name: another planet, or
another law. This file asks what a set of records still lets one say in that state, with no name of the term. It
answers four questions, each from the records, the model as it is, and synthetic sets of that model:

    level = level_for_all(problem, result, sets=200)        one level for all the readings together
    answer = missing(lin, level)                            is something missing: one answer, at that level
    kinds = classes(propose(lin, dictionary(lin)))          of which kind: the generic terms that pass, in classes
                                                            that these records cannot tell apart
    known = known_parameters(problem, result, lin, ...)     what can be said of each shared parameter: its estimate by
                                                            the fit, with a generic term added, fitted and at first
                                                            order, and by the combination of the families that is
                                                            free of that term, or of the leading class
    may = licence(problem, result, known, sets=100)         may first order be used on this problem: on its own
                                                            synthetic sets, first order beside the fit made again

    everything = read(problem, result, level)               the first four for one set of records

Nothing here knows the field of the problem. A generic term is a term of mechanisms.py, built from the fitted curves and
from the design: a wide copy under each curve, a kernel that widens it, a shared parameter that changes along a
covariate. None of them is the term that the records hold. Each stands for a class of terms that these records cannot
tell from it, and is added so that the known parameters stop carrying what the model lacks.

What it cannot do. A term that the fit absorbs whole leaves nothing to read: no level reports it, and every estimate
here is then off by what was absorbed, with no warning. `ledger.silent_zone` gives the size of that for a term that is
named. A parameter whose meaning needs the missing term, such as a value at zero of a covariate that the term follows,
is not freed by a generic term, which is centred on the records. And an estimate made after a term was chosen for its
gain was chosen after looking: `refine.on_the_twin` is where its bar earns its level.
"""
from __future__ import annotations

import numpy as np

from . import families
from .candidates import covariates_of, dictionary
from .fit import fit
from .linear import Linear
from .mechanisms import along, pedestal, widened
from .pool import Lazy, in_order
from .propose import propose
from .refine import develop_one
from .twin import replicas

EVERY = "every sample: the fit"          # the family of every sample: its estimate is the fit itself
INFLATED = 10.0          # a combination whose standard deviation is beyond this many times the fit's is of no use: the estimates then share the bias of the term, and nothing of it cancels
MOST = 2          # a combination is made free of at most this many members of a class, from the first: each member more costs variance, and the members of a class are nearly one curve
SAME = 0.5          # two estimates of a parameter closer than this, in the bar of the fitted one, are the same prediction
SHARE = 0.95          # first order has its licence when it gives the same prediction on at least this share of the synthetic sets


def readings(lin, names=None, centre=None) -> dict:
    """The level of each reading of the package at a fit: of the families of statistics that are named, or of the whole
    catalogue. A family that the package refuses on these records has no level: it is listed with the reason."""
    levels, refused = {}, {}
    for name in (list(names) if names is not None else list(families.catalogue(centre))):
        try:
            levels[name] = float(families.family(lin, name, centre).reading()["level"])
        except (ValueError, np.linalg.LinAlgError) as why:
            refused[name] = str(why)
    return {"levels": levels, "refused": refused}


def _levels_of_one_set(args):
    twin, result, names, centre, max_iter = args
    try:
        got = fit(twin, result.shared, result.own, max_iter=max_iter)
        return readings(Linear(twin, got), names, centre)["levels"] if got.converged else None
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        return None


def level_for_all(problem, result, names=None, centre=None, sets: int = 200, level: float = 0.05, seed: int = 11, max_iter: int = 80, workers: int = 1) -> dict:
    """One level for all the readings together. Each reading has a level of its own, and with many readings some reading
    reports on far more sets of a right model than its level says: with sixteen of them, on about a third. Here `sets`
    synthetic sets of the fitted model are each fitted and read as the records are, and the smallest level among the
    readings of each set is kept. The level for all is the k-th smallest of those, with k the whole part of
    level * (sets read + 1): a reading of the records under it is a finding at the level `level` for the whole of the
    readings, since on a further set of a right model the smallest level falls under it with that chance. Fewer sets
    than 1 / level give no such number, and are refused.

    What it assumes: the noise of the records is the stated one, independent from sample to sample, and the fitted
    model is near enough to the truth for its synthetic sets to stand for a right model. A noise that is correlated or
    of another size moves every level, this one too."""
    if int(level * (sets + 1)) < 1:
        raise ValueError(f"{sets} sets give no level of {level} for all the readings: at least {int(np.ceil(1.0 / level))} are needed")
    names = list(names) if names is not None else list(families.catalogue(centre))
    twins = replicas(problem, result, sets, 1.0, seed)
    read = [r for r in in_order(_levels_of_one_set, Lazy(lambda i: (twins[i], result, names, centre, max_iter), sets), workers) if r]
    k = int(level * (len(read) + 1))
    if k < 1:
        raise ValueError(f"the searches of {len(read)} of the {sets} synthetic sets ended: too few for a level of {level} for all the readings")
    smallest = np.sort([min(r.values()) for r in read])
    return {"level": float(level), "sets": int(sets), "read": len(read), "readings": names, "centre": centre, "seed": int(seed),
            "level_for_all": float(smallest[k - 1]),
            "any_reading_at_its_own_level": float(np.mean(smallest < level)),
            "each_reading_at_its_own_level": {name: float(np.mean([r[name] < level for r in read if name in r])) for name in names if any(name in r for r in read)}}


def missing(lin, calibrated: dict) -> dict:
    """Is something missing: one answer. The smallest level among the readings of the fit, against the level for all the
    readings together that `level_for_all` counted for this model. The readings and the centre are those of that count."""
    got = readings(lin, calibrated["readings"], calibrated["centre"])["levels"]
    if not got:
        raise ValueError("no reading of the package could be made on these records")
    name = min(got, key=got.get)
    return {"something_is_missing": bool(got[name] < calibrated["level_for_all"]), "smallest_level": float(got[name]), "of_the_reading": name,
            "level_for_all": float(calibrated["level_for_all"]), "level": float(calibrated["level"]), "levels": got}


def classes(table: dict) -> list:
    """Of which kind is what the model lacks. The table of `propose` ranks the candidates of a dictionary and marks, for
    each, the candidates whose visible parts are nearly the same curve: these records cannot tell those apart, nor a
    chain of such pairs. A class is such a chain. It leads when one of its members passes the threshold of the whole
    dictionary. For each shared parameter a class gives the range of the bias, in standard deviations of the fit, over
    its members: IF one of them is what the model lacks. A narrow range says that the correction is known though the
    member is not. A range that crosses zero says that the class does not even fix its sign. The classes come back
    with the one of the best member first."""
    rows = {row["name"]: row for row in table["candidates"] if row.get("status") == "ranked"}
    order = sorted(rows, key=lambda name: (rows[name]["family_level"], -rows[name]["score"] / max(rows[name]["degrees_of_freedom"], 1)))
    leading, seen, out = set(table["leading_set"]), set(), []
    for first in order:
        if first in seen:
            continue
        members, queue = [], [first]
        while queue:
            name = queue.pop()
            if name in seen or name not in rows:
                continue
            seen.add(name)
            members.append(name)
            queue += [other for other in rows[name].get("like", []) if other not in seen]
        members.sort(key=order.index)
        bias = {}
        for parameter in rows[first]["bias_if_this_is_what_is_missing"]:
            values = [rows[name]["bias_if_this_is_what_is_missing"][parameter]["in_its_standard_deviations"] for name in members]
            bias[parameter] = {"lowest": float(min(values)), "highest": float(max(values)), "of_the_first": float(values[0])}
        out.append({"members": members, "families": sorted({rows[name]["family"] for name in members}), "first": members[0], "leads": bool(leading & set(members)),
                    "level_among_all_the_candidates": float(rows[members[0]]["family_level"]), "bias_if_a_member_is_what_is_missing_in_standard_deviations": bias})
    return out


def generic_terms(problem, lin) -> list:
    """The generic terms that are fitted when none is given: a copy of each curve four times wider under it, a kernel
    with exponential wings of half the median width of the fitted curves, where the abscissa has equal steps, and each
    shared parameter along each numeric covariate of the design. `mechanisms.displaced` and `mechanisms.bent`, which
    lean a curve to one side, are not among them: they cost eight calls of the model for one, and are given by name
    when a lean is looked for."""
    width = float(np.median([lin.width_of(k) for k in range(len(problem.records))]))
    out = [pedestal(4.0), widened(size=lambda design, w=width: 0.5 * w)]
    for key in covariates_of(problem):
        out += [along(parameter, key) for parameter in problem.shared]
    return out


def at_first_order(problem, result, lin, mechanism, step: float = 1e-3) -> dict | None:
    """A generic term read at first order, with no fit: what it adds to the mean of each record for a small size is its
    column, the residuals give the size that they ask for, and the part of the column that the fit absorbed gives what
    each shared parameter is off by. Returns the size, the score of the term, a chi square of one degree of freedom
    under a right model, and for each shared parameter its value once that bias is taken off, with the standard
    deviation that it has when the term is free. None for a term that cannot be put on these records."""
    try:
        larger = mechanism.extend(problem, result)
        base = np.append(np.asarray(result.shared, dtype=float), 0.0)
        column = []
        for k, rec in enumerate(problem.records):
            up = np.asarray(larger.model(rec.x, base + np.eye(base.size)[-1] * step, result.own[k], rec.design), dtype=float)
            down = np.asarray(larger.model(rec.x, base - np.eye(base.size)[-1] * step, result.own[k], rec.design), dtype=float)
            column.append((up - down) / (2.0 * step))
        seen = lin.visible(column)
        norm2 = float(seen @ seen)
        if not np.isfinite(norm2) or norm2 <= 0.0:
            return None
        size = float(seen @ lin.r) / norm2
        moved = np.asarray(lin.absorbed(column), dtype=float) * size
        grows = lin.ledger.extended(lin.reduce(column))["sd_factor"] if len(problem.shared) else {}
    except (ValueError, np.linalg.LinAlgError, FloatingPointError):
        return None
    bars = [float(result.sd[j]) * float(grows[name]) for j, name in enumerate(problem.shared)]
    return {"term": mechanism.name, "family": mechanism.family, "size": size, "size_standard_deviation": float(1.0 / np.sqrt(norm2)), "score": float(size * size * norm2),
            "values": {name: float(result.shared[j] - moved[j]) for j, name in enumerate(problem.shared)}, "bars": {name: float(bars[j]) for j, name in enumerate(problem.shared)}, "_column": column}


def fitted(problem, result, mechanism, max_iter: int = 80) -> dict | None:
    """A generic term added to the model and fitted from the current fit: the gain in chi square, the size that comes
    back, and each shared parameter of the working model with its standard deviation in the larger fit. None when the
    term cannot be put on these records or its search does not end."""
    try:
        got = develop_one(problem, result, mechanism, max_iter=max_iter)
    except (ValueError, np.linalg.LinAlgError, FloatingPointError):
        return None
    if not got["ended"] or got.get("_fit") is None:
        return None
    larger = got["_fit"]
    return {"term": mechanism.name, "family": mechanism.family, "gain": float(got["criteria"]["gain_in_chi_square"]), "size": float(got["parameter"]["value"]), "size_standard_deviation": float(got["parameter"]["bar"]),
            "at_its_bound": bool(got["at_its_bound"]), "values": {name: float(larger.shared[j]) for j, name in enumerate(problem.shared)},
            "bars": {name: float(larger.sd[j]) for j, name in enumerate(problem.shared)}}


def estimates_of(lin, names=None, centre=None) -> dict:
    """The estimates of the shared parameters that are combined: the fit, and that of each family that estimates every
    shared parameter on these records and whose own search ends. `names` gives the families to take, all of the
    catalogue with none: the search of some of them is slow."""
    out = {EVERY: np.asarray(lin.fit.shared, dtype=float)}
    for name in (list(names) if names is not None else list(families.catalogue(centre))):
        if name == EVERY:
            continue
        try:
            fam = families.family(lin, name, centre)
            if not fam.bounds()["estimator"]:
                continue
            got = fam.estimate()
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            continue
        if got.converged:
            out[name] = np.asarray(got.shared, dtype=float)
    return out


def free_of(lin, estimates: dict, columns: dict, centre=None, inflated: float = INFLATED) -> dict:
    """The combination of several estimates of each shared parameter that is free, at first order, of the bias of each
    of some terms, with the least variance among those (families.free_of). `columns` gives each term as what it adds
    to the mean of every record, the first of them first. A combination whose standard deviation is beyond `inflated`
    times the fit's is of no use: it is made again free of one term fewer, from the last. Returns, for each shared
    parameter, the value and its standard deviation, or the reason why there is none."""
    out, terms = {}, list(columns)
    for j, parameter in enumerate(lin.problem.shared):
        out[parameter] = {"value": None, "bar": None, "why": "fewer than three estimates of this parameter: nothing is left to choose once a bias is cancelled"}
        for count in range(min(len(terms), len(estimates) - 2), 0, -1):
            try:
                free = families.free_of(lin, list(estimates), {name: columns[name] for name in terms[:count]}, centre, j)
            except (ValueError, np.linalg.LinAlgError):
                continue
            over = float(free["standard_deviation_over_the_fit"])
            if over <= inflated:
                out[parameter] = {"value": float(sum(free["weights"][name] * estimates[name][j] for name in estimates)), "bar": over * float(lin.fit.sd[j]), "over_the_fit": over, "free_of": count}
                break
            out[parameter]["why"] = f"its standard deviation is {over:.1f} times that of the fit, beyond {inflated:g}: it says too little to be of use"
    return out


def free_of_the_class(lin, kind: dict, candidates: list, estimates: dict, centre=None, inflated: float = INFLATED, most: int = MOST) -> dict:
    """The combination of the `estimates` that is free, at first order, of the bias of the first `most` members of a
    class: the members of a class are nearly one curve in what the residuals keep of them, and each member more costs
    variance. What the fit absorbed of each member is not one curve: a combination that is free of the first members
    may keep the bias of the term that the records hold, when that term is no member."""
    columns = {c.name: c.columns[0] for c in candidates}
    return free_of(lin, estimates, {name: columns[name] for name in [m for m in kind["members"] if m in columns][:most]}, centre, inflated)


def known_parameters(problem, result, lin, kind: dict | None = None, candidates: list | None = None, names=None, centre=None, terms: list | None = None, max_iter: int = 80) -> dict:
    """What can be said of each shared parameter of the model while the missing term has no name. Five estimates:

      the_fit                          the fit as it is: the right one when nothing is missing
      free_of_the_class                the combination of the families that is free of the leading class, at first order
      with_a_generic_term              the model with the generic term of the largest gain added and fitted
      with_that_term_at_first_order    the same term read with no fit: `licence` says whether this one may be used
      free_of_that_term                the combination of the families that is free of that term, at first order

    `kind` is the class that leads, from `classes`, with the `candidates` of its table: with none, the second estimate
    is not made. `names` are the families whose estimates are combined, with the fit. `terms` are the generic terms to
    try, those of `generic_terms` when none is given. The term of the last three estimates is the one whose fit gains
    the most chi square. With no shared parameter there is nothing to estimate, and the terms are still tried: their
    gains say of which kind the missing term is."""
    shared = list(problem.shared)
    kinds = ("free_of_the_class", "with_a_generic_term", "with_that_term_at_first_order", "free_of_that_term")
    out = {"parameters": {name: dict({"the_fit": {"value": float(result.shared[j]), "bar": float(result.sd[j])}}, **{k: None for k in kinds}) for j, name in enumerate(shared)},
           "class": None if kind is None else kind["members"], "terms_tried": [], "term": None}
    estimates = estimates_of(lin, names, centre) if shared else {}
    out["estimates_combined"] = len(estimates)
    if kind is not None and candidates is not None and shared:
        free = free_of_the_class(lin, kind, candidates, estimates, centre)
        for name in shared:
            out["parameters"][name]["free_of_the_class"] = free[name]
    tried = []
    for term in (generic_terms(problem, lin) if terms is None else terms):
        got = fitted(problem, result, term, max_iter)
        out["terms_tried"].append({"term": term.name, "family": term.family, "gain": None if got is None else got["gain"]})
        if got is not None:
            tried.append((got, term))
    if tried:
        best, term = max(tried, key=lambda pair: pair[0]["gain"])
        first = at_first_order(problem, result, lin, term)
        out["term"] = {"name": best["term"], "family": best["family"], "gain": best["gain"], "size": best["size"], "size_standard_deviation": best["size_standard_deviation"], "at_its_bound": best["at_its_bound"],
                       "size_at_first_order": None if first is None else first["size"], "score_at_first_order": None if first is None else first["score"]}
        out["_term"] = term
        of_it = free_of(lin, estimates, {term.name: first["_column"]}, centre) if first is not None and shared else {}
        for name in shared:
            out["parameters"][name]["with_a_generic_term"] = {"value": best["values"][name], "bar": best["bars"][name]}
            out["parameters"][name]["with_that_term_at_first_order"] = None if first is None else {"value": first["values"][name], "bar": first["bars"][name]}
            out["parameters"][name]["free_of_that_term"] = of_it.get(name)
    return out


def read(problem, result, calibrated: dict, candidates: list | None = None, names=None, terms: list | None = None, max_iter: int = 80) -> dict:
    """The first four answers for one set of records: is something missing, at the one level of `calibrated`. Of which
    kind, by the classes of the generic dictionary, or of the `candidates` that are given. And what can be said of the
    shared parameters. Where nothing is missing the fit stands and no term is tried."""
    lin = Linear(problem, result)
    answer = missing(lin, calibrated)
    candidates = dictionary(lin) if candidates is None else candidates
    kinds = classes(propose(lin, candidates))
    leading = next((kind for kind in kinds if kind["leads"]), None)
    if answer["something_is_missing"]:
        known = known_parameters(problem, result, lin, leading, candidates, names, calibrated["centre"], terms, max_iter)
    else:
        known = {"parameters": {name: {"the_fit": {"value": float(result.shared[j]), "bar": float(result.sd[j])}, "free_of_the_class": None, "with_a_generic_term": None, "with_that_term_at_first_order": None, "free_of_that_term": None}
                                for j, name in enumerate(problem.shared)}, "class": None, "terms_tried": [], "term": None}
    return {"missing": answer, "classes": kinds, "leading_class": leading, "known": known}


def _licence_of_one_set(args):
    twin, result, term, terms, max_iter = args
    try:
        got = fit(twin, result.shared, result.own, max_iter=max_iter)
        if not got.converged:
            return None
        lin = Linear(twin, got)
        first, again = at_first_order(twin, got, lin, term), fitted(twin, got, term, max_iter)
        if first is None or again is None:
            return None
        scores = {t.name: at_first_order(twin, got, lin, t) for t in terms}
        gains = {t.name: fitted(twin, got, t, max_iter) for t in terms}
    except (ValueError, np.linalg.LinAlgError, FloatingPointError):
        return None
    scores = {name: one["score"] for name, one in scores.items() if one is not None}
    gains = {name: one["gain"] for name, one in gains.items() if one is not None}
    return {"apart": {name: (first["values"][name] - again["values"][name]) / again["bars"][name] for name in again["values"]},
            "size_apart": (first["size"] - again["size"]) / again["size_standard_deviation"] if again["size_standard_deviation"] > 0 else float("inf"),
            "first_by_score": max(scores, key=scores.get) if scores else None, "first_by_gain": max(gains, key=gains.get) if gains else None}


def licence(problem, result, known: dict, terms: list | None = None, sets: int = 100, same: float = SAME, share: float = SHARE, seed: int = 13, max_iter: int = 80, workers: int = 1) -> dict:
    """May first order be used on this problem. First order is fast, and right for a term that is small. Whether the
    term of these records is small enough is not known beforehand, so it is counted. Synthetic sets are drawn from the
    model with the generic term of `known` at the size that was fitted: the nearest stand-in for the records that can be
    made with no name of the term. Each set is read twice, at first order and by fits made again, and two things are
    compared: the value of each shared parameter once the term is taken into account, and which of the generic terms
    comes first, by its score at first order and by its gain when fitted. First order has the licence for a number when
    the two readings agree on at least `share` of the sets: for a value, within `same` of the bar of the fitted one.

    The licence is of this problem at this size of the term. A larger term, another design or another model is another
    count."""
    term = known.get("_term")
    if term is None:
        raise ValueError("no generic term was fitted on these records: there is nothing whose first order could be set beside its fit")
    lin = Linear(problem, result)
    terms = generic_terms(problem, lin) if terms is None else list(terms)
    larger = term.extend(problem, result)
    grown = fit(larger, max_iter=max_iter)
    if not grown.converged:
        raise ValueError(f"the search of the model with {term.name} did not end: no synthetic set can be drawn from it")
    twins = replicas(larger, grown, sets, 1.0, seed)
    working = Lazy(lambda i: (problem.with_records(twins[i].records), result, term, terms, max_iter), sets)
    read_ = [r for r in in_order(_licence_of_one_set, working, workers) if r]
    if len(read_) < max(20, sets // 2):
        raise ValueError(f"only {len(read_)} of the {sets} synthetic sets could be read both ways: too few to give or refuse a licence")
    values = {name: float(np.mean([abs(r["apart"][name]) <= same for r in read_])) for name in problem.shared}
    first = float(np.mean([r["first_by_score"] == r["first_by_gain"] for r in read_]))
    return {"sets": int(sets), "read": len(read_), "term": term.name, "same_within_bars": float(same), "share_asked": float(share),
            "the_same_value": values, "the_same_term_comes_first": first,
            "largest_distance_in_bars": {name: float(np.max([abs(r["apart"][name]) for r in read_])) for name in problem.shared},
            "licence_for_the_values": {name: bool(v >= share) for name, v in values.items()}, "licence_for_the_term_that_comes_first": bool(first >= share)}
