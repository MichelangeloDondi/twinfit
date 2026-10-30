"""One rung of the ladder, walked from its opening to its report, on the synthetic records of two_scales.py.

    python -m twinfit.examples.walkthrough "one block wider"            a hidden term of HIDDEN, or "none"
    python -m twinfit.examples.walkthrough "second component" DIR       the notebook is kept in DIR

It is what an agent does for one rung, in the order of agent/AGENT.md, with nothing of a field in it: the rung is opened
with a person's words, the design is read before the residuals, the fit and the table of candidates take their laws from
the twin, the records are read, the first candidate is tried on synthetic records before anyone adds it to the model, and
the report says what the rung does not show. The hidden term is known to the program that makes the records and to
nothing else: the walk does not use it, and prints it at the end beside what was found.
"""
from __future__ import annotations

import sys
import tempfile

import numpy as np

from .. import Linear, fit
from ..dictionary import candidates as dictionary          # the package's own list: this example walks the path that reads it
from ..identify import identify
from ..notebook import Notebook
from ..propose import lines, propose
from ..twin import admission, noiseless, planted
from . import two_scales as E

RUNG = "two scales"


def listed(lin) -> list:
    """The list of candidates of this walk. The table of the records and the tables of the replicas are made from the one
    list: a rate of false alarms holds for the list that it was counted on."""
    return dictionary(lin, block="block")


def walk(hidden: str = "one block wider", folder: str | None = None, sets: int = 60, say=print) -> dict:
    book = Notebook(folder or tempfile.mkdtemp(prefix="twinfit_"))
    problem = E.problem(E.records(hidden, size=0.01))
    book.open_rung(RUNG, by="the person who decides: start with the two scales and nothing else")

    result = fit(problem)
    book.step("intake", RUNG, {"records": len(problem.records), "samples": problem.samples, "ended": result.converged, "chi_square": result.chi2, "degrees_of_freedom": result.dof})
    say(f"1. Intake. {len(problem.records)} records, {problem.samples} samples. The fit ended: {result.converged}. Chi square {result.chi2:.0f} on {result.dof}.")

    lin = Linear(problem, result)
    table = identify(lin, block="block")
    book.step("identify", RUNG, table)
    say("2. What the design can fix, before any residual is read. " + ", ".join(f"{n} to {p['per_cent_of_its_value']:.1f} per cent" for n, p in table["parameters"].items())
        + f". The weakest is {table['weakest']}. Each estimate rests most on the block: " + ", ".join(f"{n} on {b}" for n, b in table["rests_on"].items()) + ".")

    book.stamp(f"twin of {RUNG}", "On replicas of the fitted model at the recorded noise the fit ends on 98 per cent of the sets or more, its bias is under a quarter of "
               "its spread, its spread is 0.85 to 1.15 of the bound, and its interval holds the truth on 92 to 98 per cent, each band made as wide as three standard errors of its statistic where the sets are few. The table of candidates passes its family "
               "threshold on at most 12 per cent of the sets. Could miss: no replica was drawn before this stamp.", rung=RUNG)
    clean = noiseless(problem, result)
    adm = admission(problem, result, sets=2 * sets)
    alarms = planted(problem, result, [[None] * len(problem.records)], [0.0], sets=sets, candidates=listed)
    rate = alarms["share_of_sets_in_which_something_is_seen"]
    holds = clean["returns_what_was_put_in"] and adm["admitted"] and rate <= 0.12          # a rate of false alarms has a ceiling and no floor
    quiet = {"sets": alarms["read"], "level": 1.0, "admitted": bool(rate <= 0.12), "false_alarm_rate": rate, "rule": "it passes its family threshold on at most 0.12 of the sets of a right model"}
    book.step("twin", RUNG, {"with_no_noise": clean, "admission": adm, "readings": {"the table of candidates": quiet}, "false_alarms": alarms}, stamp=f"twin of {RUNG}")
    book.grade(f"twin of {RUNG}", holds, f"with no noise the fit returns what was put in: {clean['returns_what_was_put_in']}, admitted: {adm['admitted']}, false alarms of the table: {rate:.3f}")
    book.admit("the fit", RUNG, 1.0)                          # read from the twin step that was just recorded, not handed over
    book.admit("the table of candidates", RUNG, 1.0)
    say(f"3. The twin. The fit is admitted: {adm['admitted']}. Under a right model the table of candidates raises a false alarm on {100 * rate:.0f} per cent of the sets.")

    book.read_records("the fit", RUNG)
    book.step("fit", RUNG, {"shared": dict(zip(problem.shared, result.shared.tolist())), "standard_deviations": dict(zip(problem.shared, result.sd.tolist()))})
    book.read_records("the table of candidates", RUNG)
    cands = listed(lin)
    proposal = propose(lin, cands)
    book.step("read", RUNG, {"chi_square": proposal["chi_square"], "degrees_of_freedom": proposal["degrees_of_freedom"]})
    book.step("propose", RUNG, proposal)
    say("4. The records, read. " + lines(proposal, top=5))

    tried = None
    if proposal["leading_set"]:
        first = next(c for c in proposal["candidates"] if c["name"] == proposal["leading_set"][0])
        cand = next(c for c in cands if c.name == first["name"])
        told = ", ".join(f"{n} {first['bias_if_this_is_what_is_missing'][n]['in_its_standard_deviations']:+.1f}" for n in problem.shared)
        book.stamp(f"candidate of {RUNG} on the twin", f"The candidate '{cand.name}', planted at the size that the residuals ask for, is seen on over half of the sets, "
                   f"and the estimates are then off by {told} standard deviations, each within one. Could miss: the table is first order and the size is a noisy estimate.", rung=RUNG)
        tried = planted(problem, result, cand.columns, first["amplitude"], sets=sets, name=cand.name, candidates=listed)
        book.step("candidate on the twin", RUNG, tried, note=cand.name, stamp=f"candidate of {RUNG} on the twin")
        came = tried["estimates_less_truth_in_standard_deviations"]
        near = all(abs(came[n]["mean"] - first["bias_if_this_is_what_is_missing"][n]["in_its_standard_deviations"]) < 1.0 for n in problem.shared)
        book.grade(f"candidate of {RUNG} on the twin", bool(tried["share_of_sets_in_which_something_is_seen"] > 0.5 and near),
                   f"seen on {tried['share_of_sets_in_which_something_is_seen']:.2f} of the sets, the estimates off by " + ", ".join(f"{n} {v['mean']:+.1f}" for n, v in came.items()))
        got = tried["estimates_less_truth_in_standard_deviations"]
        say(f"5. The first candidate, '{cand.name}', planted in synthetic records at the size that the residuals ask for. It is seen on "
            f"{100 * tried['share_of_sets_in_which_something_is_seen']:.0f} per cent of the sets and is in the leading set on "
            f"{100 * tried['share_of_sets_in_which_it_is_in_the_leading_set']:.0f}. The estimates are then off by "
            + ", ".join(f"{n} {v['mean']:+.1f}" for n, v in got.items()) + " standard deviations, where the table said "
            + ", ".join(f"{n} {first['bias_if_this_is_what_is_missing'][n]['in_its_standard_deviations']:+.1f}" for n in got) + ".")
    else:
        say("5. No candidate passes the threshold of the family, so none is tried on the twin. What cannot be excluded is in the table above.")
    book.step("report", RUNG, {"leading_set": proposal["leading_set"], "bias_range": proposal["bias_over_the_leading_set_in_standard_deviations"],
                               "bias_at_which_a_term_is_seen": proposal["bias_at_which_the_most_absorbed_candidate_is_seen"],
                               "what_it_does_not_show": "which member of a leading set is the cause, any term outside the dictionary, and any term that the model absorbs"})
    truth = np.array([E.TRUTH[n] for n in problem.shared])
    off = (result.shared - truth) / result.sd
    say(f"6. What the walk could not know: the hidden term was '{hidden}', and the estimates were off by " + ", ".join(f"{n} {v:+.1f}" for n, v in zip(problem.shared, off))
        + " standard deviations. What the rung does not show: which member of a leading set is the cause, and any term outside the dictionary.")
    say(f"The notebook of this walk is in {book.folder}.")
    return {"proposal": proposal, "tried": tried, "true_bias_in_standard_deviations": dict(zip(problem.shared, off.tolist())), "notebook": str(book.folder)}


if __name__ == "__main__":
    args = sys.argv[1:]
    walk(args[0] if args else "one block wider", args[1] if len(args) > 1 else None)
