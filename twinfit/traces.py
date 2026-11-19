"""The number of traces for a target bar, row M585 (the owner's message 5: the number of traces for a target
uncertainty).

    python -m twinfit traces DIR NAME=BAR...   the times each record of DIR is to be taken for each shared parameter
                                               NAME of the model that its term files give to reach BAR, read on the twin

The information of independent records adds, so under the working model the bar of the fit falls as one over the
square root of the times each record is taken, at the same conditions and each time with its own noise. The times are
the fit's own bar over the target, squared and rounded up: arithmetic, with no twin. The twin then draws sets of that
design from the working model at its fit and fits each again, and the bias and the scatter about the value that drew
them say whether the law holds where the fit is not linear, as on few samples or a parameter near a bound. It holds
under the working model only: a term that the records miss biases a parameter by as much whatever the times, so no
count of traces reaches a target below that bias, and the report of the loop says where one is missing.
"""
from __future__ import annotations

import math

import numpy as np

from .proposals import _designed, _on_design
from .seeds import seed_of


def traces_for(problem, result, targets: dict, sets: int = 50, seed: int | None = None, max_iter: int = 80, most: int = 64) -> dict:
    """For each shared parameter named in `targets`, {name: bar}: the times each record is to be taken, the traces that
    makes, the bar that the law gives there, and on `sets` sets of the twin of that design its bias, scatter and error,
    the root of the sum of their squares, with whether that error holds the target within one standard error of a
    scatter read on that many sets. Past `most` times a set holds too many records to fit, and the law is said alone."""
    seed = seed_of(problem, 47) if seed is None else seed          # traces_for: drawn from the records, as every Monte Carlo of the package (row M582)
    names = list(problem.shared)
    out = {}
    for name, want in targets.items():
        if name not in names:
            raise ValueError(f"{name} is no shared parameter of this model: its shared parameters are {', '.join(names)}")
        if not want > 0:
            raise ValueError(f"the target bar of {name} is {want}: a bar is larger than zero")
        j = names.index(name)
        bar = float(np.sqrt(result.cov[j, j]))
        if not (math.isfinite(bar) and bar > 0):
            out[name] = {"target": want, "bar_of_the_fit": bar, "times": None, "why_not": "the fit gives it no bar"}
            continue
        times = max(1, math.ceil((bar / want) ** 2 - 1e-9))          # the information adds: the bar of the fit over the square root of the times
        row = {"target": want, "bar_of_the_fit": bar, "times": times, "traces": times * len(problem.records), "by_the_law": bar / math.sqrt(times)}
        if times > most:
            row["why_not"] = f"a set of the twin would hold {times * len(problem.records)} records, past {most} times each"
        else:
            got = _on_design(problem, result, _designed(problem, result, {"kind": "repeats", "factor": times}), None, sets, seed, max_iter)
            if got is None:
                row["why_not"] = "fewer than half the fits of the twin ended"
            else:
                read = got["with_the_working_model"][name]
                band = 1.0 / math.sqrt(2.0 * (got["sets"] - 1))          # the relative standard error of a scatter read on that many sets: two of them let an error of 1.2 times the target pass on 50 sets (the red team of his 531)
                row.update({"sets": got["sets"], "on_the_twin": read, "standard_error": read["error"] * band, "tolerance": want * band, "holds": bool(read["error"] <= want * (1.0 + band))})          # the standard error of the error read, and apart the tolerance of the check at the target (the second reading of 0.1.3rc2)
        out[name] = row
    return out


def page(out: dict, records: int, used: list) -> str:
    """The lines that the command prints: each parameter, its times and traces, and what the twin read."""
    from .loop import _pm
    lines = [f"The model that the term files give: {', '.join(used) or 'the first model alone'}, on {records} records."]
    for name, r in out.items():
        if r["times"] is None:
            lines.append(f"{name}: {r['why_not']}.")
            continue
        head = (f"{name}: its bar is {r['bar_of_the_fit']:.3g} on these records, at or under the target of {r['target']:.3g}: they reach it already."
                if r["times"] == 1 else
                f"{name}: its bar is {r['bar_of_the_fit']:.3g} on these records. To reach {r['target']:.3g}, each record is taken {r['times']} times, "
                f"{r['traces']} traces, where the law gives {r['by_the_law']:.3g}.")
        if "on_the_twin" in r:
            t = r["on_the_twin"]
            error = " ± ".join(_pm(t["error"], r["standard_error"]))          # the error at the decimals of its standard error, two significant digits (review 233, N4)
            head += (f" On {r['sets']} sets of the twin of that design: bias {t['bias']:+.2g}, scatter {t['scatter']:.2g}, error {error}, "
                     + ("which holds the target." if r["holds"] else "which does not hold the target: the fit is not linear enough there for the law, so take more."))
        else:
            head += f" The twin is not run: {r['why_not']}."
        lines.append(head)
    lines.append("This holds under the model alone: a term that the records miss biases a parameter by as much whatever the number of traces. "
                 "Read the report of the loop first.")
    return "\n".join(lines) + "\n"
