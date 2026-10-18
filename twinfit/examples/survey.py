"""The six hidden terms of two_scales.py, one line each: what the fit suffers, what the usual statistic says, what the table
of candidates says. And the twin of the right model: the admission of the fit and the false alarms of the table. Then the
same over many synthetic sets of each hidden term: how often the table sees something, how often a candidate of the hidden
kind is in the set that it cannot tell apart, how often it is first, and how often the range of bias over that set holds
the bias that the hidden term gives when no noise is drawn.

    python -m twinfit.examples.survey OUT.json [SETS]        everything: the single sets, the twin of the right model, the rates
    python -m twinfit.examples.survey --single OUT.json      the single sets and the twin again, the rates of OUT.json kept
    python -m twinfit.examples.survey --outside OUT.json [SETS]   adds the terms outside the dictionary to OUT.json
    python -m twinfit.examples.survey --derive OUT.json           adds the shares over the sets with a detection, from the counts that OUT.json holds

Every hidden term has the size of one per cent in the units of two_scales.py. The program that makes the records knows
the term. The fit and the table do not. No recorded data are read.
"""
from __future__ import annotations

import json
import sys

import numpy as np

from .. import Linear, fit
from ..candidates import dictionary
from ..propose import propose
from ..twin import admission, planted
from . import two_scales as E

RIGHT = {"second component": "a second component", "pedestal": "a pedestal", "periodic": "a periodic term",
         "scale grows with the square of the covariate": "k changes along the square of z", "one block wider": "the curves of the block z = 100 are wider"}


def main(sets: int = 120) -> dict:
    truth = np.array([E.TRUTH[n] for n in ("k", "g", "s")])
    out = {"what": "synthetic records of twinfit/examples/two_scales.py, 12 curves, noise of 0.004 a sample, seed 0, every hidden term at a size of 0.01 in its own unit",
           "what_the_size_means": {h: E.SIZES[h] for h in E.HIDDEN if h in E.SIZES},
           "hidden_terms": {}}
    for hidden in E.HIDDEN:
        problem = E.problem(E.records(hidden, size=0.01))
        result = fit(problem)
        lin = Linear(problem, result)
        table = propose(lin, dictionary(lin, block="block"))
        ranked = [c for c in table["candidates"] if c["status"] == "ranked"]
        off = (result.shared - truth) / result.sd
        cell = {"estimates_less_truth_in_standard_deviations": dict(zip(problem.shared, off.tolist())),
                "chi_square_less_its_mean_in_its_standard_deviation": table["chi_square_less_its_mean_in_its_standard_deviation"],
                "candidates_that_pass_the_family_threshold": len(table["leading_set"]),
                "bias_over_the_leading_set_in_standard_deviations": table["bias_over_the_leading_set_in_standard_deviations"],
                "bias_at_which_the_most_absorbed_candidate_is_seen_on_half_of_the_sets": {n: a["in_its_standard_deviations"] for n, a in table["bias_at_which_the_most_absorbed_candidate_is_seen"].items()}}
        if hidden in RIGHT:
            places = [i + 1 for i, c in enumerate(ranked) if c["name"].startswith(RIGHT[hidden])]
            cell["first_place_of_a_candidate_of_the_hidden_kind"] = places[0] if places else None
            cell["one_of_the_hidden_kind_is_in_the_leading_set"] = any(n.startswith(RIGHT[hidden]) for n in table["leading_set"])
        out["hidden_terms"][hidden] = cell
    problem = E.problem(E.records("none"))
    result = fit(problem)
    adm = admission(problem, result, sets=sets)
    alarms = planted(problem, result, [[None] * len(problem.records)], [0.0], sets=sets, block="block")
    out["twin_of_the_right_model"] = {"sets": sets, "the_fit_is_admitted": adm["admitted"], "parameters": adm["parameters"],
                                      "share_of_sets_in_which_the_table_passes_its_family_threshold": alarms["share_of_sets_in_which_something_is_seen"]}
    return out


def visible(hidden_terms=E.HIDDEN + E.OUTSIDE) -> dict:
    """For each hidden term at its size of the survey, with no noise drawn: the squared norm of what it adds to the records
    in units of the noise, and the part of that norm which stays in the residuals of the working model, which is the
    noncentrality of the statistic of fit. A term can go unseen because the model absorbs it, or because it is small:
    the two numbers tell which."""
    base = E.problem(E.records("none", noise=0.0))
    lin = Linear(base, fit(base))
    clean = [rec.y for rec in base.records]
    out = {}
    for hidden in hidden_terms:
        if hidden == "none":
            continue
        d = [rec.y - y0 for rec, y0 in zip(E.records(hidden, size=0.01, noise=0.0), clean)]
        total = float(lin.norm2(d))
        left = lin.ledger.visible(lin.reduce(d))
        seen = float(left @ left)
        out[hidden] = {"squared_norm_in_units_of_the_noise": total, "of_which_stays_in_the_residuals": seen, "share_of_the_norm_that_is_visible": float(np.sqrt(seen / total)) if total > 0 else None}
    return out


def of_the_kind(name: str, hidden: str) -> bool:
    return hidden in RIGHT and name.startswith(RIGHT[hidden])


def rates(sets: int = 60, hidden_terms=E.HIDDEN) -> dict:
    """Over `sets` synthetic sets of each hidden term, each with its own noise: what the table does. The bias that a hidden
    term gives is read from the fit of its records with no noise drawn, in the standard deviations of the fit. A term of
    E.OUTSIDE has no candidate of its kind in the dictionary: for it the table can see something and name nothing right."""
    truth = np.array([E.TRUTH[n] for n in ("k", "g", "s")])
    out = {"sets_of_each_hidden_term": sets, "size_of_each_hidden_term": 0.01, "what_the_size_means": {h: E.SIZES[h] for h in hidden_terms if h in E.SIZES}, "hidden_terms": {}}
    for hidden in hidden_terms:
        clean = fit(E.problem(E.records(hidden, size=0.01, noise=0.0)))
        bias = (clean.shared - truth) / clean.sd
        seen = first = inside = 0
        families = {}
        holds = np.zeros(3)
        sizes, off = [], []
        for seed in range(1, sets + 1):
            problem = E.problem(E.records(hidden, size=0.01, seed=seed))
            result = fit(problem)
            lin = Linear(problem, result)
            table = propose(lin, dictionary(lin, block="block"), seed=seed)
            off.append((result.shared - truth) / result.sd)
            lead = table["leading_set"]
            sizes.append(len(lead))
            if lead:
                seen += 1
                ranked = [c for c in table["candidates"] if c["status"] == "ranked"]
                first += int(of_the_kind(ranked[0]["name"], hidden))
                families[ranked[0]["family"]] = families.get(ranked[0]["family"], 0) + 1
                inside += int(any(of_the_kind(n, hidden) for n in lead))
                span = table["bias_over_the_leading_set_in_standard_deviations"]
                holds += np.array([span[n]["from"] <= bias[i] <= span[n]["to"] for i, n in enumerate(problem.shared)])
        off = np.array(off)
        out["hidden_terms"][hidden] = {
            "sets_in_which_something_is_seen": seen,
            "family_of_the_first_candidate_share_of_those": {f: n / seen for f, n in sorted(families.items())} if seen else None,
            "of_those_share_with_the_hidden_kind_in_the_set": inside / seen if seen and hidden in RIGHT else None,
            "of_those_share_with_the_hidden_kind_first": first / seen if seen and hidden in RIGHT else None,
            "bias_with_no_noise_in_standard_deviations": dict(zip(("k", "g", "s"), bias.tolist())),
            "mean_of_estimates_less_truth_in_standard_deviations": dict(zip(("k", "g", "s"), off.mean(axis=0).tolist())),
            "share_of_sets_in_which_the_table_sees_something": seen / sets,
            "median_number_of_candidates_in_the_set": float(np.median(sizes)),
            "share_of_sets_with_the_hidden_kind_in_the_set": inside / sets if hidden in RIGHT else None,
            "share_of_sets_with_the_hidden_kind_first": first / sets if hidden in RIGHT else None,
            "of_the_sets_in_which_something_is_seen_share_whose_range_of_bias_holds_the_bias": None if not seen else dict(zip(("k", "g", "s"), (holds / seen).tolist()))}
    return out


def derive(res: dict) -> dict:
    """The shares over the sets with a detection, for a block of rates that was written before rates() wrote them: they are
    counts of the same sets divided by another count, so no set is drawn again. What cannot be derived, the family of the
    first candidate, stays absent and the block says so."""
    block = res["over_many_sets"]
    sets = block["sets_of_each_hidden_term"]
    for hidden, c in block["hidden_terms"].items():
        if "sets_in_which_something_is_seen" in c:
            continue
        seen = int(round(c["share_of_sets_in_which_the_table_sees_something"] * sets))
        c["sets_in_which_something_is_seen"] = seen
        for key, share in (("of_those_share_with_the_hidden_kind_in_the_set", c["share_of_sets_with_the_hidden_kind_in_the_set"]), ("of_those_share_with_the_hidden_kind_first", c["share_of_sets_with_the_hidden_kind_first"])):
            c[key] = None if share is None or not seen else int(round(share * sets)) / seen
    block["the_shares_over_the_sets_with_a_detection"] = ("derived from the counts of the same sets, with no set drawn again: this block was written before rates() wrote them. "
                                                          "The family of the first candidate cannot be derived and is absent here: a new run of rates() writes it")
    return res


if __name__ == "__main__":
    if sys.argv[1:2] == ["--derive"]:                     # adds the shares over the sets with a detection to a block of rates that lacks them
        with open(sys.argv[2]) as fh:
            res = derive(json.load(fh))
        with open(sys.argv[2], "w") as fh:
            fh.write(json.dumps(res, indent=1) + "\n")
        for name, c in res["over_many_sets"]["hidden_terms"].items():
            print(f"{name:45s} seen on {c['sets_in_which_something_is_seen']} of {res['over_many_sets']['sets_of_each_hidden_term']} | of those, kind in the set {c['of_those_share_with_the_hidden_kind_in_the_set']} | kind first {c['of_those_share_with_the_hidden_kind_first']}")
        raise SystemExit(0)
    if sys.argv[1:2] == ["--outside"]:                    # adds the terms outside the dictionary to a file that holds the rest
        with open(sys.argv[2]) as fh:
            res = json.load(fh)
        res["outside_the_dictionary"] = rates(int(sys.argv[3]) if len(sys.argv) > 3 else 40, E.OUTSIDE)
        with open(sys.argv[2], "w") as fh:
            fh.write(json.dumps(res, indent=1) + "\n")
        for name, c in res["outside_the_dictionary"]["hidden_terms"].items():
            b, m = c["bias_with_no_noise_in_standard_deviations"], c["of_the_sets_in_which_something_is_seen_share_whose_range_of_bias_holds_the_bias"]
            print(f"{name:45s} bias " + ", ".join(f"{n} {v:+.1f}" for n, v in b.items()) + f" | seen on {c['share_of_sets_in_which_the_table_sees_something']:.2f} | median set {c['median_number_of_candidates_in_the_set']:.0f} | "
                  "the range holds the bias: " + ("nothing seen" if m is None else ", ".join(f"{n} {v:.2f}" for n, v in m.items())))
        raise SystemExit(0)
    if sys.argv[1:2] == ["--single"]:                     # the single sets and the twin again, the rates of the file kept as they are
        with open(sys.argv[2]) as fh:
            old = json.load(fh)
        res = main()
        res["what_each_hidden_term_leaves_visible_with_no_noise"] = visible()
        for part in ("over_many_sets", "outside_the_dictionary"):
            if part in old:
                res[part] = old[part]
        sys.argv = [sys.argv[0], sys.argv[2]]
    else:
        res = main()
        res["over_many_sets"] = rates(int(sys.argv[2]) if len(sys.argv) > 2 else 60)
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w") as fh:
            fh.write(json.dumps(res, indent=1) + "\n")
    for name, c in res["hidden_terms"].items():
        off = c["estimates_less_truth_in_standard_deviations"]
        print(f"{name:45s} off by " + ", ".join(f"{n} {v:+.1f}" for n, v in off.items()) + f" | chi square {c['chi_square_less_its_mean_in_its_standard_deviation']:+.1f} | "
              f"{c['candidates_that_pass_the_family_threshold']} pass | of its kind: place {c.get('first_place_of_a_candidate_of_the_hidden_kind')}, in the set {c.get('one_of_the_hidden_kind_is_in_the_leading_set')}")
    for name, c in res["over_many_sets"]["hidden_terms"].items():
        b, m = c["bias_with_no_noise_in_standard_deviations"], c["of_the_sets_in_which_something_is_seen_share_whose_range_of_bias_holds_the_bias"]
        print(f"{name:45s} bias " + ", ".join(f"{n} {v:+.1f}" for n, v in b.items()) + f" | seen on {c['share_of_sets_in_which_the_table_sees_something']:.2f} | "
              f"its kind in the set {c['share_of_sets_with_the_hidden_kind_in_the_set']}, first {c['share_of_sets_with_the_hidden_kind_first']} | median set {c['median_number_of_candidates_in_the_set']:.0f} | "
              "the range holds the bias: " + ("nothing seen" if m is None else ", ".join(f"{n} {v:.2f}" for n, v in m.items())))
    t = res["twin_of_the_right_model"]
    print(f"twin of the right model, {t['sets']} sets: admitted {t['the_fit_is_admitted']}, the table passes on {t['share_of_sets_in_which_the_table_passes_its_family_threshold']:.3f}")
