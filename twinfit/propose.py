"""The candidates of a dictionary ranked on the residuals of a fit, with what each would mean for the estimates.

For a candidate with columns D, in units of the noise and freed of the own parameters, the residuals keep the visible part
V = (I - P) D. The score of the candidate is the part of the squared residual that V explains, a chi square with as many
degrees of freedom as V has independent columns when the working model is right. This is the canonical score test for an
added term, read here for a whole dictionary at once. For each candidate the table gives:

    score, degrees of freedom      the canonical test of the term, and its level on its own
    family level                   the same level once the whole dictionary was looked at: the share of draws of the
                                   residuals under a right model in which SOME candidate does as well. It is computed
                                   with no fit, from the correlations of the visible parts
    amplitude                      the size of the term that the residuals ask for, with its standard deviation
    bias                           by how much each shared parameter is off IF this term is what the model lacks: the
                                   absorbed part of the fitted term, in standard deviations and in per cent
    visibility                     the share of the fitted term that the residuals keep. The rest went into the estimates
    bias at detection              the bias of each shared parameter that a term of this shape causes at the size at which
                                   its score passes the 5 per cent point of the family with a stated probability: one
                                   half by default, and 0.8 and 0.95 beside it. It is not a bound: at the size of the
                                   first number the term is missed on half of the sets. It is given for every candidate,
                                   whatever its score, for the terms of this dictionary and for no others. It becomes
                                   part of an uncertainty only when the theory of the problem bounds the size of the term
    like                           the candidates whose visible parts are nearly the same curve: these records cannot tell
                                   them apart, and only theory or another design can
    leading set                    every candidate that passes the 5 per cent threshold of the family. The residuals say
                                   that something is missing and do not say which member of the set it is: their scores
                                   differ by less than the noise of a score. The biases that the members would mean are
                                   given as a range for each parameter, and they can differ in sign
    does the work of the first     a member of the set that, once taken out, leaves the first candidate under the 5 per
                                   cent point of its own test: the two explain the same part of the residuals

Three things the table cannot do, said here so that no reader takes it for more. A candidate that lies in the span of the
model is absorbed whole: it leaves no residual, its score is zero whatever its size, and it is listed as such. A term that
the dictionary does not hold is not named: the table then shows its nearest relatives. And a term chosen because it ranked
first is chosen after looking: an estimate made with it is not unbiased, and its interval does not keep its level, until
synthetic records with the term have shown what the whole procedure returns.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import chi2

from .candidates import ALONG, ALONG_SQUARE

from .ledger import size_at

ABSORBED = 1e-4          # a column whose visible share is under this is taken as inside the span of the model
RANK = 1e-9


def _combine(columns: list, weights) -> list:
    """The discrepancy sum of weights[j] * columns[j], record by record."""
    n_rec = len(columns[0])
    out = []
    for k in range(n_rec):
        parts = [w * np.asarray(col[k], dtype=float) for col, w in zip(columns, weights) if col[k] is not None]
        out.append(sum(parts) if parts else None)
    return out


def _one(lin, cand) -> dict:
    led = lin.ledger
    D = np.column_stack([lin.reduce(col) for col in cand.columns])
    V = np.column_stack([led.visible(D[:, j]) for j in range(D.shape[1])])
    total = np.array([lin.norm2(col) for col in cand.columns])
    share = np.sqrt(np.maximum((V * V).sum(axis=0), 0.0) / np.maximum(total, 1e-300))
    entry = {"name": cand.name, "family": cand.family, "parity": cand.parity, "columns": int(D.shape[1]),
             "visible_share_of_each_column": share.tolist()}
    if total.max() == 0.0 or share.max() < ABSORBED:
        entry.update({"status": "absorbed whole: no statistic of these records sees it", "W": None})
        return entry
    G = V.T @ V
    w, U = np.linalg.eigh(G)
    keep = w > RANK * w.max()
    W = V @ (U[:, keep] / np.sqrt(w[keep]))                      # an orthonormal basis of the visible part
    c = W.T @ lin.r
    amp = (U[:, keep] / w[keep]) @ (U[:, keep].T @ (V.T @ lin.r))
    amp_sd = np.sqrt(np.maximum(np.diag((U[:, keep] / w[keep]) @ U[:, keep].T), 0.0))
    fitted = _combine(cand.columns, amp)
    bias = lin.absorbed(fitted)
    values, sd = lin.fit.shared, led.sd
    visible = W @ c
    per_unit = np.column_stack([led.shift(D[:, j]) for j in range(D.shape[1])]) / sd[:, None]          # the move of each parameter, in its standard deviations, for one unit of each column
    Ginv = (U[:, keep] / w[keep]) @ U[:, keep].T
    entry.update({"status": "ranked", "degrees_of_freedom": int(keep.sum()), "score": float(c @ c),
                  "reach": np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", per_unit, Ginv, per_unit), 0.0)).tolist(),
                  "level_on_its_own": float(chi2.sf(float(c @ c), int(keep.sum()))),
                  "amplitude": amp.tolist(), "amplitude_standard_deviation": amp_sd.tolist(),
                  "visibility": float(np.sqrt(float(visible @ visible) / max(lin.norm2(fitted), 1e-300))),
                  "bias_if_this_is_what_is_missing": {n: {"in_its_standard_deviations": float(b / s), "per_cent_of_its_value": float(100.0 * b / abs(v)) if v != 0 else None}
                                                      for n, b, s, v in zip(lin.problem.shared, bias, sd, values)},
                  "harm_in_joint_standard_deviations": float(np.sqrt(bias @ np.linalg.solve(led.C, bias))),
                  "W": W, "unit": visible / max(np.linalg.norm(visible), 1e-300)})
    return entry


def _score_without(W: np.ndarray, other: np.ndarray, r: np.ndarray) -> float:
    """The score of the visible part W on the residuals once the visible part `other` has been taken out of both."""
    rest = r - other @ (other.T @ r)
    left = W - other @ (other.T @ W)
    keep = np.linalg.norm(left, axis=0) > 1e-6
    if not keep.any():
        return 0.0
    q, _ = np.linalg.qr(left[:, keep])
    return float(np.sum((q.T @ rest) ** 2))


def family_levels(entries: list, draws: int, seed: int) -> tuple:
    """For each ranked candidate, the share of draws of the residuals under a right model in which the best candidate of the
    whole dictionary reaches a level as small as this one's. First order: the scores are drawn from their joint normal law.
    Also the level on its own that the best candidate passes in 5 per cent of those draws: the threshold of the family."""
    ranked = [e for e in entries if e["status"] == "ranked"]
    if not ranked:
        return [], 0.05
    W = np.hstack([e["W"] for e in ranked])
    R = W.T @ W
    lam, E = np.linalg.eigh(R)
    root = E * np.sqrt(np.maximum(lam, 0.0))
    z = np.random.default_rng(seed).standard_normal((draws, R.shape[0])) @ root.T
    best = np.ones(draws)
    at = 0
    for e in ranked:
        m = e["degrees_of_freedom"]
        best = np.minimum(best, chi2.sf((z[:, at:at + m] ** 2).sum(axis=1), m))
        at += m
    return [float((np.sum(best <= e["level_on_its_own"]) + 1.0) / (draws + 1.0)) for e in ranked], float(np.quantile(best, 0.05))


POWERS = (0.5, 0.8, 0.95)       # the probabilities of being seen at which the bias of a candidate is given


def propose(lin, candidates: list, draws: int = 4000, seed: int = 7, similar: float = 0.9, symmetric_only: bool = False) -> dict:
    """The table of the candidates on the residuals of a fit. `symmetric_only`: the candidates that are odd about the centre
    of a curve are read as alarms and none of them is proposed as a term."""
    entries = [_one(lin, c) for c in candidates]
    levels, threshold = family_levels(entries, draws, seed)
    names = lin.problem.shared
    sizes = {}
    for e, level in zip([e for e in entries if e["status"] == "ranked"], levels):
        e["family_level"] = level
        reach = e.pop("reach")
        # the root of the noncentrality at which the score of this candidate passes the threshold of the family with each
        # probability. It depends on the degrees of freedom alone, so it is found once for each number of them
        m = e["degrees_of_freedom"]
        if m not in sizes:
            sizes[m] = {p: size_at(1.0, m, threshold, p) for p in POWERS}
        at = sizes[m]
        e["bias_at_which_it_is_seen_with_probability"] = {f"{p:g}": {n: float(r * at[p]) for n, r in zip(names, reach)} for p in POWERS}
    ranked = sorted([e for e in entries if e["status"] == "ranked"], key=lambda e: (e["level_on_its_own"], -e["score"]))
    for a in ranked:
        a["like"] = [b["name"] for b in ranked if b is not a and abs(float(a["unit"] @ b["unit"])) >= similar]
        a["proposed"] = not (symmetric_only and a["parity"] == "odd")
    leading = []
    if ranked:
        first = ranked[0]
        first["score_once_the_first_is_taken_out"] = None
        first["does_the_work_of_the_first"] = first["family_level"] <= 0.05
        for e in ranked[1:]:
            e["score_once_the_first_is_taken_out"] = _score_without(e["W"], first["W"], lin.r)
            e["does_the_work_of_the_first"] = False
            if first["family_level"] <= 0.05 and e["family_level"] <= 0.05:
                left = _score_without(first["W"], e["W"], lin.r)
                e["score_of_the_first_once_this_is_taken_out"] = left
                e["does_the_work_of_the_first"] = bool(left < chi2.isf(0.05, first["degrees_of_freedom"]))
        leading = [e for e in ranked if e["family_level"] <= 0.05 and e.get("proposed", True)]
    chi_sq = float(lin.r @ lin.r)
    table = []
    for e in ranked + [e for e in entries if e["status"] != "ranked"]:
        table.append({k: v for k, v in e.items() if k not in ("W", "unit")})
    allowance = {}
    for n in names:
        worst = max((e for e in ranked), key=lambda e: e["bias_at_which_it_is_seen_with_probability"][f"{POWERS[0]:g}"][n], default=None)
        if worst is not None:
            allowance[n] = {"in_its_standard_deviations": worst["bias_at_which_it_is_seen_with_probability"][f"{POWERS[0]:g}"][n], "from": worst["name"], "seen_with_probability": POWERS[0],
                            "the_same_seen_with_probability": {f"{p:g}": worst["bias_at_which_it_is_seen_with_probability"][f"{p:g}"][n] for p in POWERS[1:]}}
    spread = {}
    for n in names:
        if leading:
            vals = [e["bias_if_this_is_what_is_missing"][n]["in_its_standard_deviations"] for e in leading]
            spread[n] = {"from": float(min(vals)), "to": float(max(vals))}
    return {"model": lin.problem.name, "chi_square": chi_sq, "degrees_of_freedom": lin.dof,
            "threshold_of_the_family_as_a_level_on_its_own": threshold,
            "leading_set": [e["name"] for e in leading],
            "bias_over_the_leading_set_in_standard_deviations": spread,
            "bias_at_which_the_most_absorbed_candidate_is_seen": allowance,
            "chi_square_less_its_mean_in_its_standard_deviation": float((chi_sq - lin.dof) / np.sqrt(2.0 * lin.dof)),
            "candidates_looked_at": len(entries), "draws_for_the_family_level": draws, "similar_from": similar, "symmetric_only": symmetric_only,
            "candidates": table}


def _family_level(c: dict, floor: float) -> str:
    """A level in the family at the floor of its draws is a floor and is printed as one: fifteen candidates of the first
    trial printed the same 0.00025, which read as fifteen measured levels."""
    return f"under {floor:.2g}" if c["family_level"] <= floor * (1.0 + 1e-9) else f"{c['family_level']:.2g}"


def lines(table: dict, top: int = 8, show=()) -> str:
    """The table as text, for a person or an agent: the first candidates, then those named in `show` wherever they stand,
    then those that are absorbed whole. `show` is for the candidates that the researcher's own theory names: when the
    residuals do not see them they are not among the first rows, and that they are not seen is what must be said."""
    names = list(table["candidates"][0].get("bias_if_this_is_what_is_missing", {})) if table["candidates"] else []
    out = [f"{table['model']}: chi square {table['chi_square']:.1f} on {table['degrees_of_freedom']} degrees of freedom, "
           f"{table['chi_square_less_its_mean_in_its_standard_deviation']:+.1f} of its standard deviations from its mean. {table['candidates_looked_at']} candidates looked at.",
           "| candidate | score, degrees of freedom | level on its own | level in the family | bias if it is the missing term, in standard deviations: " + ", ".join(names) + " | kept by the residuals | cannot be told from |",
           "|---|---|---|---|---|---|---|"]
    read = []
    for c in table["candidates"]:                     # the keys of the design that the list reads as covariates, in the order of the list
        if c.get("family") == "grows with a covariate" and ALONG in c["name"]:
            key = c["name"].split(ALONG_SQUARE if ALONG_SQUARE in c["name"] else ALONG, 1)[1]
            if key not in read:
                read.append(key)
    out[0] += f" Read as covariates of the design: {', '.join(read) or 'none'}."
    ranked = [c for c in table["candidates"] if c["status"] == "ranked"]
    floor = 1.0 / (table.get("draws_for_the_family_level", 4000) + 1.0)
    asked = [c for c in ranked[top:] if c["name"] in set(show)]
    for c in ranked[:top] + asked:
        bias = ", ".join(f"{c['bias_if_this_is_what_is_missing'][n]['in_its_standard_deviations']:+.1f}" for n in names)
        mark = "" if c.get("proposed", True) else " (odd: an alarm, not proposed)"
        mark += " (shown because it was asked for)" if c in asked else ""
        out.append(f"| {c['name']}{mark} | {c['score']:.1f}, {c['degrees_of_freedom']} | {c['level_on_its_own']:.2g} | {_family_level(c, floor)} | {bias} | "
                   f"{c['visibility']:.2f} | {', '.join(c['like'][:3]) or 'none in the dictionary'} |")
    missing = [n for n in show if n not in {c["name"] for c in table["candidates"]}]
    if missing:
        out.append("Asked for and not in this table: " + ", ".join(missing) + ".")
    if table.get("leading_set"):
        rng = table["bias_over_the_leading_set_in_standard_deviations"]
        shown = table["leading_set"][:6]
        out.append(f"Something is missing: {len(table['leading_set'])} candidates pass the threshold of the family, and these records do not choose among them. The first: " + ", ".join(shown) + ".")
        out.append("The bias that they would mean, in standard deviations: " + ". ".join(f"{n} from {v['from']:+.1f} to {v['to']:+.1f}" for n, v in rng.items()) + ".")
    else:
        out.append("No candidate of the dictionary passes the threshold of the family: the residuals name nothing.")
    allowance = table.get("bias_at_which_the_most_absorbed_candidate_is_seen", {})
    if allowance:
        out.append("The bias at which the most absorbed term of this dictionary is seen on half of the sets, in standard deviations, and on 95 per cent of them: "
                   + ". ".join(f"{n} {a['in_its_standard_deviations']:.1f} and {a['the_same_seen_with_probability']['0.95']:.1f}, from {a['from']}" for n, a in allowance.items())
                   + ". It is not a bound: a term of the first size is missed on half of the sets.")
    gone = [c["name"] for c in table["candidates"] if c["status"] != "ranked"]
    if gone:
        out.append("Absorbed whole by the model, so that no statistic of these records sees them: " + ", ".join(gone) + ".")
    return "\n".join(out)
