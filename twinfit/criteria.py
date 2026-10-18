"""Information criteria for a model with terms more, and the threshold that a family of candidates must pass.

The errors are taken as normal with the standard deviations of the records, so minus twice the log likelihood is the chi
square up to a constant, and for a model with m parameters more whose chi square is lower by `gain`:

    AIC                 prefers the larger model when  gain > 2 m
    AICc                the same with the correction for a small number n of samples
    BIC                 when  gain > m ln n.  n is given twice: the number of samples, and the number of records, since
                        the samples of one record share its own parameters and the two counts give different answers
    Hannan and Quinn    when  gain > 2 m ln ln n

Each is given as a margin: the gain less what the criterion asks. A positive margin prefers the larger model.

None of these knows how many candidates were looked at. A term picked as the best of sixty passes AIC nearly always
under a right model. So the choice among candidates that rest on no theory is made in the family: the gain of the best
candidate is set against the best gains on synthetic sets of a right model, and the term is added when at most 5 per
cent of those sets do as well. That rule is this package's, from the plan of its author: the owner of the method asked
for the criteria to be shown and for a default choice.

And none knows the system. A term that stands on solid ground in the theory of the problem is added whatever these
numbers say: a criterion weighs a gain of fit against a count of parameters and nothing else. `choose` holds that rule.
"""
from __future__ import annotations

import math

import numpy as np


def criteria(gain: float, extra: int, samples: int, records: int, parameters: int) -> dict:
    """The margins of each criterion for a model with `extra` parameters more than one of `parameters`, whose chi square
    is lower by `gain`."""
    aicc = 2.0 * (parameters + extra) * (parameters + extra + 1) / (samples - parameters - extra - 1) - 2.0 * parameters * (parameters + 1) / (samples - parameters - 1) \
        if samples - parameters - extra - 1 > 0 else float("nan")
    return {"gain_in_chi_square": float(gain), "parameters_more": int(extra),
            "AIC": float(gain - 2.0 * extra),
            "AICc": float(gain - 2.0 * extra - aicc),
            "BIC_n_the_samples": float(gain - extra * math.log(samples)),
            "BIC_n_the_records": float(gain - extra * math.log(records)) if records > 1 else float("nan"),
            "Hannan_Quinn": float(gain - 2.0 * extra * math.log(math.log(samples))) if samples > 3 else float("nan")}


MIN_SETS = 19          # with n sets of a right model the smallest level that a gain can reach is 1 / (n + 1): 0.05 needs 19
GOOD_SETS = 211        # and a level of 0.05 is known to a standard error of 0.015 from this many


def level_of(gain: float, best) -> float:
    """The level of a gain in the family: the share of the sets of a right model on which the best candidate gained as
    much, counted so that it is never zero: (count + 1) / (sets + 1)."""
    best = np.asarray(best, dtype=float)
    return float((np.sum(best >= gain) + 1.0) / (best.size + 1.0))


def family_threshold(gains, level: float = 0.05) -> dict:
    """What the sets of a right model say of the best gain of a family of candidates. `gains` has one row a synthetic set
    of the working model and one column a candidate. The threshold is the gain that the best candidate passes on `level`
    of the sets: a number to print, with what it is worth. From 12 sets it lies anywhere between the gains that 1 and 25
    per cent of right models reach, and from 30 between those of 3 and 15 per cent (measured on the worked example). So the
    choice is not made on it: `choose` counts the sets, with `level_of`, and that count cannot fall under 1 / (sets + 1)."""
    g = np.atleast_2d(np.asarray(gains, dtype=float))
    best = g.max(axis=1)
    n = int(g.shape[0])
    return {"threshold": float(np.quantile(best, 1.0 - level)), "level_of_the_threshold": level, "sets": n, "candidates": int(g.shape[1]),
            "largest_gain_seen_under_a_right_model": float(best.max()), "best_gains": best.tolist(),
            "smallest_level_that_these_sets_allow": 1.0 / (n + 1.0), "enough_sets_for_a_level_of_0.05": bool(n >= MIN_SETS),
            "standard_error_of_the_level": float(math.sqrt(level * (1.0 - level) / n)),
            "sets_for_a_standard_error_of_0.015": GOOD_SETS}


def _unfitted(d: dict) -> str:
    """Why a term on theory has no fitted size: its search failed, or it ended at its bound, which is another thing."""
    if d.get("ended") and d.get("at_its_bound"):
        return "its own parameter ended at its bound on these records, where the term is absent or as large as the search lets it be"
    return "its own parameter could not be fitted on these records"


def _said(ground: str) -> str:
    """The ground of a term as the caller wrote it, with a capital: it follows a full stop or a colon in the sentence of the choice."""
    ground = str(ground).strip().rstrip(".")
    return ground[:1].upper() + ground[1:]


def choose(developed: list, best=None) -> dict:
    """The default choice of a stage. `best` holds, for each synthetic set of the working model, the largest gain that a
    candidate reached there: the law of a gain when nothing is missing.

    A term that stands on ground in the theory of the problem is added whatever a criterion says. Of several that were
    fitted, the one with the largest gain is added in this stage, and the others are named, since a stage adds one
    term: they are the next. One whose own parameter could not be fitted, because its search failed or ended at its
    bound, where the term is absent, is not dropped for a guess: when no other
    term on theory was fitted it is added held at the size that the theory gives, with the bar of the theory carried
    (refine.carried), and otherwise it is named for a later stage. One that could be neither fitted nor held, since
    no size was given for it, is named with what it needs, and the choice goes on among the others.

    With no term on theory left, the guess with the largest gain is added if its level in the family is 0.05 or under.
    With no law of the gains, or one of fewer than 19 sets, nothing is chosen among guesses: a gain with no law is not
    evidence, and on fewer sets no gain can reach that level."""
    out = {"chosen": None, "how": None, "by_default": True, "terms_on_theory_not_added_in_this_stage": [], "left_out": []}
    broken = [d for d in developed if d.get("ended") and not np.isfinite(d["criteria"]["gain_in_chi_square"])]      # a gain that is no number ranks nowhere, and is named
    developed = [d for d in developed if d not in broken]
    out["left_out"] = [{"name": d["name"], "why": "its gain is not a number"} for d in broken]
    grounded = [d for d in developed if d.get("ground")]
    fitted = [d for d in grounded if d.get("ended") and not d.get("at_its_bound")]
    unfitted = [d for d in grounded if d not in fitted]       # its search failed, or ended at its bound: the records put the term at nothing, and the theory says that it is there
    holdable = [d for d in unfitted if d.get("theory_given")]
    later = [{"name": d["name"], "why": "it stands on theory and " + _unfitted(d) + ": it is to be held at the size that the theory gives, with its bar carried"}
             for d in holdable]
    later += [{"name": d["name"], "why": "it stands on theory and its gain is not a number: fit it again before anything is added"} for d in broken if d.get("ground")]
    later += [{"name": d["name"], "why": "it stands on theory, its own parameter could not be fitted on these records, and no size was given at which to hold it: "
                                         "give the centre and the bar of the theory with the term, as theory=(centre, bar)"} for d in unfitted if not d.get("theory_given")]
    if fitted:
        first = max(fitted, key=lambda d: d["criteria"]["gain_in_chi_square"])
        out.update(chosen=first["name"], how="fitted", why="it stands on ground in the theory of the problem, so it is added whatever a criterion says: " + _said(first["ground"]),
                   rests_on={"ground": str(first["ground"])})
        out["terms_on_theory_not_added_in_this_stage"] = [{"name": d["name"], "why": "it also stands on theory: a stage adds one term, and this one is next"} for d in fitted if d is not first] + later
        return out
    if holdable:
        first = holdable[0]
        out.update(chosen=first["name"], how="held", why="it stands on ground in the theory of the problem, and " + _unfitted(first) + ": "
                   "it is added held at the size that the theory gives, with the bar of the theory carried. " + _said(first["ground"]),
                   rests_on={"ground": str(first["ground"])})
        out["terms_on_theory_not_added_in_this_stage"] = later[1:]
        return out
    out["terms_on_theory_not_added_in_this_stage"] = later            # those that could be neither fitted nor held: the choice goes on, and says so
    usable = [d for d in developed if d.get("ended") and not d.get("ground")]
    if not usable:
        out["why"] = "no candidate could be fitted"
        return out
    top = max(usable, key=lambda d: d["criteria"]["gain_in_chi_square"])
    gain = top["criteria"]["gain_in_chi_square"]
    if best is None or len(best) == 0:
        out["why"] = "the gains have no law yet: no synthetic set of the working model was fitted with these candidates, and nothing is chosen on a gain alone"
        return out
    n = len(best)
    if n < MIN_SETS:
        out["why"] = (f"the law of the gains rests on {n} sets, and on fewer than {MIN_SETS} no gain can reach a level of 0.05: the smallest that {n} sets allow is "
                      f"{1.0 / (n + 1.0):.3f}. Nothing is chosen")
        return out
    level = level_of(gain, best)
    count = int(np.sum(np.asarray(best, dtype=float) >= gain))
    if level <= 0.05:
        out.update(chosen=top["name"], how="fitted", why=f"its gain of {gain:.1f} is the largest, and the best candidate gained as much on {count} of {n} sets where nothing is missing: "
                   f"a level of {level:.3f} in the family", rests_on={"level_in_the_family": float(level), "sets": int(n)})
    else:
        out["why"] = (f"the largest gain, {gain:.1f}, is reached by the best candidate on {count} of {n} sets where nothing is missing: a level of {level:.3f} in the family, "
                      "above 0.05. The stage ends here")
    return out
