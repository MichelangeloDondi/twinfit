"""The report of a stage: what a researcher needs in order to keep control of an analysis that an agent runs.

After each stage, and before any term is added, the report says five things, in this order:

    1  the results so far, each with two bars: the bar of the fit, checked on synthetic sets where they were run, and the
       bias that a term of the dictionary could cause before the residuals show it
    2  the structure of the residuals: the chi square against its law, the residuals by position about the centre of
       the curves, even and odd, and the candidates that the residuals rank first
    3  the candidates that were developed as fitted models, three by default: what each gains, what its own parameter
       comes back as, what the criteria say, and what it does to every other parameter
    4  the default choice, and the rule that made it
    5  what the stage does not show

    text = stage_report(problem, result, table, developed, choice, law=law, admission=admission, title="stage 1")
    write(folder, "stage 1", text, problem, result)        the page, and a picture of the residuals where matplotlib is there

The researcher can stop at any report, or go back: a default is a default. The report gives numbers and names no cause:
which mechanism of the system a term stands for is the researcher's to say.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .linear import Linear
from .propose import lines as table_lines
from .refine import held_lines
from .refine import lines as developed_lines

BANDS = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)        # distances from the centre of a curve, in widths of the curve
LETTERS = 13          # the size of the words of the figure of a report, in points: 9 and 10 until his message 191


def by_position(problem, result) -> list:
    """The whitened residuals in pairs of bands about the centre of each curve, over all the records: the sum of both
    sides and the difference of the two. Each is read through the linear picture of the fit: the part of the band that
    the parameters can take is taken out of it, and what the residuals keep is divided by its own standard deviation
    after the fit. So each is a standard normal number at first order when nothing is missing, whatever the fit took.

    Until candidate 41 each sum was divided by the root of its samples, the standard deviation with no fit. Under a
    right model the twelve numbers of the worked example then had standard deviations well under one, since the own
    parameters of every record are fitted to those same residuals: 0.41 to 0.94 on the 60 sets of the test of this
    function, 0.46 to 1.02 on the 120 sets of review 40, 0.38 to 0.92 on the 80 sets of another reader of the export.
    `share_kept` says how much of the standard deviation of a band the fit leaves."""
    lin = Linear(problem, result)
    unit = float(np.median([lin.width_of(k) for k in range(len(problem.records))]))
    centre = [lin.centre_of(k) for k in range(len(problem.records))]
    rows = []
    for lo, hi in zip(BANDS, BANDS[1:]):
        even, lean, n = [], [], 0
        for k, rec in enumerate(problem.records):
            u = (rec.x - centre[k]) / unit
            inside = (np.abs(u) >= lo) & (np.abs(u) < hi)
            even.append(rec.sd * inside)
            lean.append(rec.sd * inside * np.sign(u))
            n += int(inside.sum())
        if n:
            (a, ka), (b, kb) = _after_the_fit(lin, even, n), _after_the_fit(lin, lean, n)
            rows.append({"from": lo, "to": hi, "samples": n, "both_sides": a, "one_side_less_the_other": b, "share_kept": (ka, kb)})
    return rows


def by_position_and_record(problem, result) -> list:
    """The readings of `by_position`, record by record: for each pair of bands about the centre and each record, the sum
    of both sides and the difference of the two over the samples of that record alone, each read through the linear
    picture of the fit of every record and divided by its own standard deviation after that fit. A standard normal
    number at first order when nothing is missing, as `by_position` gives one over all the records. How the structure of
    what the fit leaves changes from one record to the next, along the settings of the design, is read on these
    (the owner's message 170). A record whose share of a band the fit takes whole reads nan there."""
    lin = Linear(problem, result)
    unit = float(np.median([lin.width_of(k) for k in range(len(problem.records))]))
    centre = [lin.centre_of(k) for k in range(len(problem.records))]
    rows = []
    for lo, hi in zip(BANDS, BANDS[1:]):
        both, lean, n = [], [], 0
        for k, rec in enumerate(problem.records):
            u = (rec.x - centre[k]) / unit
            inside = (np.abs(u) >= lo) & (np.abs(u) < hi)
            if not inside.any():
                both.append(float("nan"))
                lean.append(float("nan"))
                continue
            alone = lambda w: [w if j == k else None for j in range(len(problem.records))]          # noqa: E731
            both.append(_after_the_fit(lin, alone(rec.sd * inside), int(inside.sum()))[0])
            lean.append(_after_the_fit(lin, alone(rec.sd * inside * np.sign(u)), int(inside.sum()))[0])
            n += int(inside.sum())
        if n:
            rows.append({"from": lo, "to": hi, "both_sides": both, "one_side_less_the_other": lean})
    return rows


def _margin(v: float) -> str:
    """A margin of a criterion as the report prints it. A criterion that these records do not define, such as the one that
    counts the records when there is one record, is said to be none."""
    return f"{v:+.1f}" if np.isfinite(v) else "none here"


def _after_the_fit(lin, per_record, samples: int) -> tuple:
    """A sum of whitened residuals over a set of samples, over its standard deviation after the fit, and the share of its
    standard deviation with no fit that the fit leaves. A set that the fit takes whole has no reading: nan."""
    v = lin.visible(per_record)
    norm = float(np.sqrt(v @ v))
    kept = norm / np.sqrt(samples)
    return (float(v @ lin.r) / norm if kept > 1e-6 else float("nan")), kept


def by_block(problem, result, key: str) -> list:
    """The chi square of the records of each value of one covariate of the design, less the mean that it has after the
    fit when nothing is missing: its samples, less the own parameters that the fit took from those records, less the
    share of the shared parameters that falls to them by their samples. The last is an approximation, of at most the
    number of shared parameters. Until candidate 41 the samples alone were taken off, and the number sat under zero."""
    lin = Linear(problem, result)
    shared, total = len(problem.shared), sum(rec.x.size for rec in problem.records)
    rows = []
    for value, members in problem.blocks(key).items():
        c = sum(float(np.sum(problem.whitened(k, result.shared, result.own[k]) ** 2)) for k in members)
        n = sum(problem.records[k].x.size for k in members)
        mean = n - sum(lin.parts[k][1].shape[1] for k in members) - shared * n / total
        rows.append({"value": value, "records": len(members), "chi_square_less_its_samples": c - n, "chi_square_less_its_mean": c - mean,
                     "in_its_standard_deviations": (c - mean) / np.sqrt(2.0 * mean) if mean > 0 else float("nan")})
    return rows


def stage_report(problem, result, table: dict, developed: list | None = None, choice: dict | None = None, law: dict | None = None,
                 admission: dict | None = None, title: str = "A stage", blocks: tuple = (), allowance: dict | None = None) -> str:
    out = [f"# {title}: {problem.name}", ""]
    out += ["## 1. The results so far", ""]
    allow = table.get("bias_at_which_the_most_absorbed_candidate_is_seen", {})
    out += ["| parameter | value | bar of the fit | its scatter on synthetic sets, over that bar | bias that a term of the dictionary could cause before half of the sets show it, in bars |", "|---|---|---|---|---|"]
    for j, name in enumerate(problem.shared):
        cell = (admission or {}).get("parameters", {}).get(name) if admission else None
        seen = ("none: held at a bound" if cell and cell.get("spread_over_bound") is None else f"{cell['spread_over_bound']:.2f}" if cell
                else "too few searches ended" if admission is not None else "not simulated")
        a = allow.get(name)
        bar = f"{result.sd[j]:.2g}" if np.isfinite(result.sd[j]) else "none: held at its bound"
        out.append(f"| {name} | {result.shared[j]:.5g} | {bar} | {seen} | {a['in_its_standard_deviations']:.1f}, from {a['from']} |" if a
                   else f"| {name} | {result.shared[j]:.5g} | {bar} | {seen} | none ranked |")
    second = " The search with the bounds of the problem had stopped before its end: this is the end of a second search with no bounds, which lies inside them." if getattr(result, "finished_without_bounds", False) else ""
    out += ["", f"The search {'ended' if result.converged else 'DID NOT END: nothing below is to be trusted'} in {result.iterations} steps." + second
            + (f" {result.held} own parameters sit at a bound and are held." if result.held else "")
            + (f" {result.undetermined} own parameters are not fixed by their record: they mean nothing." if result.undetermined else "")
            + (" HELD AT A BOUND, with no bar: " + ", ".join(problem.shared[j] for j in result.at_bound) + ". The bars of the others are those of a fit with these held."
               if result.at_bound else "")]
    if admission is not None:
        out.append(f"On {admission['sets']} synthetic sets at level {admission['level']:g} the fit is {'admitted' if admission['admitted'] else 'NOT admitted'}: {admission['rule']}."
                   + "".join(f" {w[0].upper() + w[1:]}." for w in admission.get("why_not", [])))
    out += ["", "## 2. The structure of the residuals", ""]
    said = ""
    if allowance:
        each = ", ".join(f"one bar of {name} moves that mean by {v['in_standard_deviations_of_the_chi_square']:+.1f}" for name, v in allowance["parts"].items())
        said = (f" The weights come from a law of the noise that has bars: {each}, {allowance['together_in_standard_deviations']:.1f} together, in the same units. "
                f"With its own scatter that is {np.sqrt(1.0 + allowance['together_in_standard_deviations'] ** 2):.1f}, and this chi square sits at "
                f"{allowance['in_standard_deviations_of_its_scatter_and_the_bars_together']:+.1f} of that."
                + (" The bars of the law, more than the noise, set how far the chi square of these records can be read: a missing term shows there only if it moves "
                   "the chi square by more than they do. Read the structure of the residuals." if allowance["together_in_standard_deviations"] > 1.0 else ""))
    out += [f"The chi square is {result.chi2:.1f} on {result.dof} degrees of freedom: {(result.chi2 - result.dof) / np.sqrt(2.0 * result.dof):+.1f} of its standard deviations from its mean." + said,
            "", "By position about the centre of each curve, over all the records. Each number is the sum of the whitened residuals of a band, over its standard deviation "
            "after the fit: at first order a standard normal number when nothing is missing. The last column is the share of the standard deviation of the band that the fit "
            "leaves, for the two sums: where it is small the parameters take most of what such a band could show, and a misfit there is seen late. "
            "A table with no number beyond 2 does not show that nothing is missing.", "",
            "| distance from the centre, in widths | samples | both sides together | one side less the other | share that the fit leaves |", "|---|---|---|---|---|"]
    for r in by_position(problem, result):
        out.append(f"| {r['from']:g} to {r['to']:g} | {r['samples']} | {r['both_sides']:+.1f} | {r['one_side_less_the_other']:+.1f} | {r['share_kept'][0]:.2f} and {r['share_kept'][1]:.2f} |")
    for key in blocks:
        out += ["", f"By {key}: the chi square of its records less the mean that it has after the fit when nothing is missing, in standard deviations of that difference.", "",
                f"| {key} | records | in standard deviations |", "|---|---|---|"]
        for r in by_block(problem, result, key):
            shown = f"{r['value']:.6g}" if isinstance(r["value"], float) else r["value"]          # a value of the design is shown to six digits: 403.15 as it is, 7 / 9 as 0.777778
            out.append(f"| {shown} | {r['records']} | {r['in_its_standard_deviations']:+.1f} |")
    out += ["", "What the dictionary ranks first on these residuals, at first order and with no fit:", "", table_lines(table, top=6)]
    out += ["", "## 3. The candidates developed as fitted models", ""]
    if developed:
        out += [developed_lines(developed), ""]
        for d in developed:
            if not d.get("ended"):
                continue
            p, c = d["parameter"], d["criteria"]
            th = p.get("theory")
            came = (f"Its parameter ended at its bound, {p['value']:.4g}: the records put the term at nothing, or as large as the search lets it be, and it has no bar"
                    + (f". The theory gives {th['centre']:g} with {th['bar']:g}. " if th else ". ")) if d.get("at_its_bound") else (
                    f"Its parameter came back at {p['value']:.4g} with a bar of {p['bar']:.2g}"
                    + (f", where the theory gives {th['centre']:g} with {th['bar']:g}: {th['from_theory_in_joint_bars']:+.1f} joint bars from it. " if th else ". "))
            out.append(f"- **{d['name']}.** " + (f"On theory: {d['ground']}. " if d.get("ground") else "An educated guess, with no ground in the theory given. ") + came
                       + f"Margins, positive for the larger model: AIC {_margin(c['AIC'])}, AICc {_margin(c['AICc'])}, BIC with the samples {_margin(c['BIC_n_the_samples'])}, "
                       f"BIC with the records {_margin(c['BIC_n_the_records'])}, Hannan and Quinn {_margin(c['Hannan_Quinn'])}. "
                       + "With it: " + ", ".join(f"{k} = {v['after']:.5g} " + (f"with a bar of {v['bar_after']:.2g}" if np.isfinite(v["bar_after"]) else "held at its bound")
                                                 for k, v in d["what_it_does_to_each_parameter"].items()) + ".")
        if law and "threshold" in law:
            n = law["sets_that_ended"]
            out += ["", f"The law of a gain when nothing is missing: on {n} synthetic sets of the working model"
                    + (f", at {law['noise_level']:g} times the noise of the records," if law.get("noise_level") not in (None, 1.0) else "")
                    + " the best of these candidates gained more than "
                    f"{law['threshold']:.1f} on 5 per cent of the sets, and {law['largest_gain_seen_under_a_right_model']:.1f} at most. "
                    f"What these {n} sets are worth: the smallest level that a gain can reach on them is {law['smallest_level_that_these_sets_allow']:.3f}, and a level "
                    f"of 0.05 is known to a standard error of {law['standard_error_of_the_level']:.3f}. "
                    + ("" if law.get("enough_sets_for_a_level_of_0.05") else "That is too few for a choice among guesses: none is made. ")
                    + (f"{law['sets_dropped']} of the {law['sets_asked']} sets asked were dropped, since one of their searches did not end: the law is that of the sets that ended. "
                       if law.get("sets_dropped") else "")
                    + "The level in the family of the table above is the share of these sets on which the best candidate gained as much, and the default choice is made on that count."]
        else:
            out += ["", "The gains have no law yet: no synthetic set of the working model was fitted with these candidates."]
    else:
        out.append("None was developed in this stage.")
    out += ["", "## 4. The default choice", ""]
    if choice:
        head = "No term is added" if choice["chosen"] is None else "Next: " + choice["chosen"] + (", held at the size that the theory gives" if choice.get("how") == "held" else "")
        out.append(f"**{head}.** {choice['why'][0].upper() + choice['why'][1:]}. "
                   "This is a default: the researcher can stop here, choose another candidate, or go back a stage.")
        out += held_lines(choice.get("held"))
        for other in choice.get("left_out", []):
            out.append(f"- Left out of the choice: {other['name']}. {other['why'][0].upper() + other['why'][1:]}.")
        for other in choice.get("terms_on_theory_not_added_in_this_stage", []):
            out.append(f"- On theory and not added in this stage: {other['name']}. {other['why'][0].upper() + other['why'][1:]}.")
        if choice["chosen"] is None and not choice.get("the_law_of_the_gains"):
            out.append("To go on by default the gains need their law: `refine.null_gains` fits these candidates on synthetic sets of the working model, 30 by default.")
    else:
        out.append("No choice was made in this stage.")
    out += ["", "## 5. What this stage does not show", ""]
    gone = [c["name"] for c in table["candidates"] if c["status"] != "ranked"]
    out += ["- A term that the working model absorbs whole leaves no residual. No statistic of these records sees it, whatever its size"
            + (": " + ", ".join(gone[:4]) + (" and others" if len(gone) > 4 else "") + (" are such terms here." if len(gone) > 1 else " is such a term here.") if gone else "."),
            "- A term that is in neither the dictionary nor the list of developed candidates is not named: the tables show its nearest relatives.",
            "- Two candidates whose gains are close cannot be told apart by these records. Theory, or another design, can.",
            "- A candidate chosen for its gain was chosen after looking: an estimate made with it keeps the level of its bar only once synthetic records with the term have shown what the whole procedure returns.",
            "- The bars are those of the working model. If a term is missing, the bias that it causes is not in them."]
    if admission is None:
        out.append("- The bar of the fit was not checked on synthetic sets in this stage.")
    return "\n".join(out) + "\n"


def write(folder, name: str, text: str, problem=None, result=None) -> list:
    """Writes the page of a stage, and beside it a picture of the mean residual by position where matplotlib is there."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    stem = "".join(ch if ch.isalnum() else "_" for ch in name).strip("_")
    page = folder / f"report_{stem}.md"
    made = [page]
    if problem is not None and result is not None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            text += "\nNo picture: matplotlib is not installed where this ran.\n"
        else:
            lin = Linear(problem, result)
            unit = float(np.median([lin.width_of(k) for k in range(len(problem.records))]))
            edges = np.linspace(-8.0, 8.0, 65)
            where = [np.digitize((rec.x - lin.centre_of(k)) / unit, edges) - 1 for k, rec in enumerate(problem.records)]
            value = np.full(64, np.nan)
            for b in range(64):
                count = int(sum(int((idx == b).sum()) for idx in where))
                if count:
                    value[b] = _after_the_fit(lin, [rec.sd * (idx == b) for rec, idx in zip(problem.records, where)], count)[0]
            fig, ax = plt.subplots(figsize=(9.0, 4.6))          # larger, and its letters larger, his message 191: "make the fonts bigger to make the figures and the slides more clear"
            mid = 0.5 * (edges[:-1] + edges[1:])
            use = np.isfinite(value)
            ax.axhspan(-2.0, 2.0, color="0.85", lw=0)
            ax.axhline(0.0, color="0.4", lw=0.8)
            ax.plot(mid[use], value[use], color="#B5482A", lw=1.8)
            ax.set_xlabel("distance from the centre of the curve, in widths of the curve", fontsize=LETTERS)
            ax.set_ylabel("sum of the residuals of a band\nover its standard deviation after the fit", fontsize=LETTERS)
            ax.set_title("The residuals of all the records by position.\nThe grey band is two standard deviations.", fontsize=LETTERS + 1)
            ax.tick_params(labelsize=LETTERS - 1)
            fig.tight_layout()
            fig.savefig(folder / f"report_{stem}.png", dpi=150)
            plt.close(fig)
            made.append(folder / f"report_{stem}.png")
            text += (f"\n![The residuals by position](report_{stem}.png)\n\n**How to read it.** Along the bottom is the distance from the centre of a curve, in widths of the curve, "
                     "left of the centre to the left. Up is the sum of the whitened residuals of all the records in a band of that distance, over its standard deviation "
                     "after the fit: at first order it stays inside the grey band on 95 bands in 100 when nothing is missing. "
                     "A shape that leaves the band is what the model lacks, summed over the records. A flat line does not show that nothing is missing.\n")
    page.write_text(text)
    return made
