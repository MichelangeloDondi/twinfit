"""The visual log of a run of the loop: for each model that the loop reached, each of the eight phases of its turn first
as one plot, row M600, then beside four of them the older slides, the equation of the model and the effect of each of its
terms, every curve in a row, the answers of the report and each model's criterion, thirteen kinds of slide a model in
`SLIDES`, a kind on more than one page where its curves, its terms or its proposals ask
it, and the best model at the end with the
parts of the bar of each of its first ten shared parameters that have a size. The owner's message 309: "show clearly also the equations used and the
effects of the terms of the models, in each turn, with at least a visual log slide for each phase".

    python -m twinfit loop DIR
    python -m twinfit log DIR

It reads what the run left in DIR/loop: equations.json and effects_K.npz, the notebook of each phase, residuals_K.csv,
model_K.json and best.json, beside DIR/description.toml for the settings and DIR/terms/ for the files of the terms. It fits nothing, so it takes seconds. It writes DIR/loop/log/, one picture a slide,
DIR/loop/log.pdf with every slide in order, and DIR/loop/log.md, which lists the slides with the one sentence that each
shows. A run made before 0.1.0rc8 left no equations.json: its slides say so where the equation and the effects stand.

The phases of a turn are those of the notebook, in its order:

    1 the records            what was read: the curves, how many, how many samples
    2 the model              its shared parameters with their bars, then its equation, term by term, and what each term does
    3 the twin               whether synthetic sets of the fitted model give back what made them, so that the fit may be read
    4 the fit                the joint fit with what it leaves above it, then every curve in a row with what it leaves
    5 what the fit leaves    the mean square of each curve's residuals over its noise, against the setting that its sweep changes
    6 where to search        where the residuals leave the twin's band, then the answers of the report and the six proposals
    7 the next term          the term tried after this model, on synthetic sets that hold it
    8 the decision           what the next term changes, then each model's criterion and whether the term was kept, and why
"""
from __future__ import annotations

import csv
import json
import math
import re
import textwrap
import tomllib
from pathlib import Path

import numpy as np

from .best import TOLD_APART
from .notebook import ORDER

TITLES = {"intake": "The records", "identify": "The model and its equation", "twin": "The twin: may this fit be read",
          "fit": "The joint fit", "read": "What the fit leaves", "propose": "Where to search",
          "candidate on the twin": "The next term, tried on the twin", "report": "The decision"}
SHORT = {"intake": "records", "identify": "model", "twin": "twin", "fit": "fit", "read": "residuals", "propose": "where",
         "candidate on the twin": "next term", "report": "decision"}
COLOURS = ("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#000000", "#999933", "#882255", "#44AA99")
ACCENT, GOOD, BAD, GREY = "#0072B2", "#009E73", "#D55E00", "#8a8f98"
SIZE = (13.33, 7.5)          # a slide of 16 by 9
TEXT, SMALL, TITLE = 17, 15, 26          # no letter under 15 points, row M583 and the toy's rule
LETTERS = ((r"\\ell", "ℓ the straight line under the curve"), (r"\\Theta", "Θ one side of the centre"), (r"w\(z\)", "w(z) the width at the setting"), (r"G_", "G a Gaussian of standard deviation σ"), (r"L_", "L a Lorentzian of full width Γ"),
           (r"(?<![a-z\\])v(?![a-z])", "v the share of the displacement"), (r"(?<![a-z\\])p ", "p the parameter that it changes"),
           (r"(?<![a-z\\])z", "z the setting, in standard deviations from its mean"), (r"(?<![a-z\\])c(?![a-z])", "c the centre of the curve"), (r"(?<![a-z\\])s(?![a-z])", "s its stretch"),
           (r"(?<![a-z\\])D(?![a-z])", "D its distance"), (r"(?<![a-z\\])m(?![a-z])", "m the middle of the abscissas"), (r"(?<![a-z\\])h(?![a-z])", "h their half span"),
           (r"(?<![a-z\\])P(?![a-z])", "P the period"), (r"(?<![a-z\\])b(?![a-z])", "b the scale of the kernel"), (r"(?<![a-z\\])r(?![a-z])", "r its reach"),
           (r"(?<![a-z\\])u(?![a-z])", "u the distance within the kernel"), (r"(?<![a-z\\])d(?![a-z])", "d the largest displacement"), (r"(?<![a-z\\])q(?![a-z])", "q the setting over its mean over the records"))          # the letters of the formulas of the package's terms, each said on the slide where it stands (review 84, minor 3)
MOST_PANELS = 6              # the effects drawn on one slide, the largest first: a slide of twelve panels reads as none
TOP = 0.74                   # the top of the figures of a slide, under its title and its two lines of words
WIDTH = 104                  # the characters of a line of the words of a slide


def _rc() -> dict:
    return {"font.size": TEXT, "axes.titlesize": TEXT, "axes.labelsize": TEXT, "xtick.labelsize": SMALL, "ytick.labelsize": SMALL,
            "legend.fontsize": SMALL, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white",
            "mathtext.fontset": "cm"}          # the formulas in Computer Modern, the letters of LaTeX (his message 314: "use a proper latex writing for the formulas")


def read_run(folder, description=None) -> dict:
    """What a run of the loop left in `folder`/loop, read once: the equations and effects of each model, the notebook of
    each phase by model, the residuals of each model by curve, the table of each model, best.json, and the settings of
    each curve from the description, `folder`/description.toml or else the file `description`, as for a loop kept in a
    folder of its own beside the project. Raises where the folder holds no run."""
    folder = Path(folder)
    described = folder / "description.toml" if (folder / "description.toml").is_file() else Path(description) if description else None
    loop = folder / "loop"
    if not (loop / "notebook").is_dir():
        raise ValueError(f"{loop} holds no run of the loop: run `python -m twinfit loop {folder}` first")
    phases: dict = {}
    for path in sorted((loop / "notebook").glob("*_model_*_*.json")):
        m = re.match(r"\d+_model_(\d+)_(.+)\.json$", path.name)
        if m:
            phases.setdefault(int(m.group(1)), {})[m.group(2).replace("_", " ")] = json.loads(path.read_text())
    equations = json.loads((loop / "equations.json").read_text()) if (loop / "equations.json").is_file() else []
    residuals, tables, effects = {}, {}, {}
    for k in phases:
        path = loop / f"residuals_{k}.csv"
        if path.is_file():
            by: dict = {}
            with path.open() as fh:
                for row in csv.DictReader(fh):
                    by.setdefault(row["curve"], []).append([float(row[c]) for c in ("abscissa", "value", "fitted", "residual", "noise")])
            residuals[k] = {name: np.array(rows) for name, rows in by.items()}
        if (loop / f"model_{k}.json").is_file():
            tables[k] = json.loads((loop / f"model_{k}.json").read_text())
        if (loop / f"effects_{k}.npz").is_file():
            with np.load(loop / f"effects_{k}.npz") as held:
                effects[k] = {key: held[key] for key in held.files}
    settings, units = {}, {}
    if described is not None and described.is_file():
        d = tomllib.loads(described.read_text())
        units = {key: str(c.get("unit", "")) for key, c in d.get("covariates", {}).items() if isinstance(c, dict)}          # in the order of the description
        for t in d.get("trace", []):
            own = {key: v for key, v in t.items() if key not in ("file", "columns") and isinstance(v, (int, float, str))}
            settings[t.get("file", "")] = {key: own[key] for key in [k for k in units if k in own] + [k for k in own if k not in units]}
        problem = d.get("problem", {})
    else:
        problem = {}
    best = json.loads((loop / "best.json").read_text()) if (loop / "best.json").is_file() else {}
    return {"folder": folder, "loop": loop, "phases": phases, "equations": equations, "residuals": residuals, "tables": tables,
            "effects": effects, "settings": settings, "units": units, "best": best, "problem": problem}


def label_of(run: dict, curve: str, room: int = 32) -> str:
    """A curve named by its settings where the description gives them, in its order, each number at three significant
    figures with its unit, as `290 K`, or with its name where it has no unit, as `pressure 10`, as many of them as
    `room` characters hold, or else by its file. A line of eight curves of four settings each, one of them at
    eight figures, was named past the edge of its slide and its panels' titles lay on one another (0.1.0rc8)."""
    units = run.get("units", {})
    parts = [(f"{v:.3g} {units[k]}" if units.get(k) else f"{k} {v:.3g}") if isinstance(v, (int, float)) else f"{k} {v}" for k, v in run["settings"].get(curve, {}).items()]
    out = ""
    for part in parts:
        if out and len(out) + 2 + len(part) > room:
            break
        out = part if not out else out + ", " + part
    return (out if len(out) <= room else out[:room - 3] + "...") or Path(curve).stem


def chain_of(run: dict, k: int) -> list:
    """The rows of equations.json from the first model to model k, through the model that each extends."""
    rows = {row["model"]: row for row in run["equations"]}
    at, out = rows.get(f"model {k}"), []
    while at:
        out.insert(0, at)
        at = rows.get(at.get("after")) if at.get("after") else None
    return out


def pm(value: float, bar: float) -> str:
    """A value with its bar, the bar at two significant digits and the value at its decimals."""
    if not (math.isfinite(value) and math.isfinite(bar)) or bar <= 0:
        return f"{value:.4g}"
    digits = max(0, 1 - int(math.floor(math.log10(bar))))
    return f"{value:.{digits}f} ± {bar:.{digits}f}"


def _frame(title: str, k: int, phase: str | None, run: dict, sentence: str = ""):
    """A slide: its title, the model and the phase it belongs to, the one sentence it shows, and at its foot the eight
    phases of a turn with the current one filled."""
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=SIZE)
    fig.text(0.035, 0.94, title, fontsize=TITLE, fontweight="bold", va="top")
    where = f"model {k}" if k >= 0 else "the output"
    if phase is not None:
        where += f"  ·  phase {ORDER.index(phase) + 1} of {len(ORDER)}"
    fig.text(0.965, 0.94, where, fontsize=SMALL, color=GREY, ha="right", va="top")
    if sentence:
        fig.text(0.035, 0.865, two_lines(sentence, WIDTH), fontsize=TEXT, va="top", linespacing=1.3)          # cut with its sign: two lines cut without one dropped the verdict of the best model (the slides read against their words, 0.1.3rc4)
    for j, name in enumerate(ORDER):
        x = 0.035 + j * 0.1175
        on = phase == name
        fig.patches.append(plt.Rectangle((x, 0.018), 0.108, 0.045, transform=fig.transFigure, facecolor=ACCENT if on else "white",
                                         edgecolor=ACCENT, linewidth=1.2))
        fig.text(x + 0.054, 0.0405, f"{j + 1} {SHORT[name]}", fontsize=SMALL, ha="center", va="center", color="white" if on else ACCENT)
    fig.text(0.965, 0.075, run["folder"].name, fontsize=SMALL, color=GREY, ha="right")
    return fig


def two_lines(text: str, width: int) -> str:
    """A text on two lines of `width` characters at most, cut at a word with its sign where it is longer, shortened until
    its two lines hold it: wrapped and then cut to two lines, a third line and its sign were dropped with no word (the
    fourteenth reading of 0.1.3rc4)."""
    lines = textwrap.wrap(text, width)
    room = 2 * width
    while len(lines) > 2 and room > 10:
        room -= 4
        lines = textwrap.wrap(textwrap.shorten(text, room, placeholder="..."), width)
    return "\n".join(lines[:2])


def _reason(why, room: int) -> str:
    """A reason as a slide says it: cut at a word with its sign where it is long, closed by one full stop where it is not.
    Cut at 140 characters with no mark and a full stop after it, a refused term's reason read as whole (the twelfth reading
    of 0.1.3rc4)."""
    said = textwrap.shorten(str(why), room, placeholder="...")
    return said if said.endswith("...") else said.rstrip(".") + "."


def _no_data(fig, words: str) -> None:
    fig.text(0.5, 0.45, "\n".join(textwrap.wrap(words, 80)), fontsize=TEXT, ha="center", va="center", color=GREY)          # wrapped by its characters: a long sentence ran past both edges (the blind red team of 0.1.0rc12)


def slide_records(run: dict, k: int):
    """The slide of phase 1 before row M600: every curve that the model reads, each in its colour, with the count of
    records and samples. The log no longer draws it: its phase 1 is `slide_curves`, `SLIDES`. It stays for a reader who
    draws it alone, and for its tests."""
    intake = run["phases"].get(k, {}).get("intake", {})
    res = run["residuals"].get(k, {})
    fig = _frame(TITLES["intake"], k, "intake", run, f"{intake.get('records', len(res))} records, {intake.get('samples', '?')} samples, read as the description gives them."
                 + (" The same records as model 0." if k else ""))
    if not res:
        _no_data(fig, "The run kept no residuals of this model.")
        return fig
    ax = fig.add_axes([0.08, 0.17, 0.86, TOP - 0.17])
    for j, (curve, a) in enumerate(res.items()):
        ax.plot(a[:, 0], a[:, 1], ".", ms=2.5, color=COLOURS[j % len(COLOURS)], label=label_of(run, curve))
    ax.set_xlabel(axis_names(run)[0])
    ax.set_ylabel(axis_names(run)[1])
    ax.legend(loc="upper left", frameon=True, markerscale=5, fontsize=SMALL)          # inside the plot, where its corner is empty: beside it, the names of eight curves left the slide (0.1.0rc8)
    return fig


def parameter_rows(params: dict, room: int = 30) -> list:
    """The rows of the table of parameters of a slide: each name on lines of `room` characters at most, its value with
    its bar, and the bar as a share of the value. A name of 42 characters ran past its cell in a blind run of 0.1.0rc10
    (a time scale of the model over its theory, squared)."""
    return [["\n".join(textwrap.wrap(name, room)) or name, pm(v["value"], v["standard_deviation"]), f"{v['per_cent_of_its_value']:.2g} %"] for name, v in params.items()]


LETTER_LINES = 3          # the lines of letters under a formula, the last saying how many more where they do not fit: two lines left five of the first CO2 model's ten letters out with no word (the slides read against their words, 0.1.3rc4)
LETTER_STEP = 0.03          # the height of each line of letters after the first
TALL = ("\\frac", "\\sum", "\\int", "\\prod")          # what stands taller than a line in a formula: under a sum of HITRAN's lines its letters lay on the sum's lower limit (the slides read against their words, 0.1.3rc4)


def _tall(formula) -> bool:
    return any(x in str(formula or "") for x in TALL)


def _rise(formula) -> float:
    """The room a tall formula takes under its line: a sum's lower limit reaches lower than a fraction's denominator."""
    return 0.045 if "\\sum" in str(formula or "") or "\\prod" in str(formula or "") else 0.025 if _tall(formula) else 0.0


def name_labels(names: list, room: int) -> list:
    """The names of parameters on a narrow axis, each on two lines, the first of `room` characters and the second of four
    more at most, cut at a word, and
    numbered by their place where two still read alike: cut on one line, "temperature exponent" read "..." and the
    amplitudes of the fringe and of its fundamental both "amplitude..." (the slides read against their words, 0.1.3rc4)."""
    out = []
    for n in names:
        lines = textwrap.wrap(n, room) or [n]
        rest = " ".join(lines[1:])
        if len(rest) > room + 4:          # cut: its leading article left out, where "the..." alone told the fundamental from nothing
            rest = " ".join(rest.split()[1:]) if rest.split()[0] in ("the", "a", "an") and len(rest.split()) > 1 else rest
        out.append(lines[0] if len(lines) == 1 else lines[0] + "\n" + textwrap.shorten(rest, room + 4, placeholder="..."))
    return [x + (f" ({i + 1})" if out.count(x) > 1 else "") for i, x in enumerate(out)]


def name_size(labels: list) -> int:
    """The size of names on an axis: 15 points, the least of his rule, where names of two lines stand on more than five
    rows, else TICK. At 20 points the two lines of neighbouring rows touched with no gap (the slides read against their
    words, 0.1.3rc4)."""
    return 15 if len(labels) > 5 and any("\n" in x for x in labels) else TICK


def term_letters(run: dict | None, j: int, row: dict) -> tuple:
    """The definitions of the letters of line j of the equation, and the file that defines them where the slide cannot:
    the model's own `LETTERS` for the first model, a term's own letters where it gives them, the letters of the package's
    mechanisms where the term is one of them, and none for a term that writes its own formula, whose letters only its
    file defines. Read by patterns, the CO2 terms of HITRAN's lines called c_2 the centre of the curve and P/P_0 the
    period (the slides read against their words, 0.1.3rc4)."""
    formula = str(row.get("formula") or "")
    if not formula.startswith("$"):
        return [], None
    if not j or row.get("letters"):
        return [f"{a} {b}" for a, b in (row.get("letters") or {}).items()], None
    path = Path(run["folder"]) / "terms" / str(row.get("file")) if run and run.get("folder") and row.get("file") else None
    try:
        own = path is not None and "formula=" in path.read_text()
    except OSError:
        own = False
    if own:
        return [], f"terms/{row['file']}"
    return letters_list(formula, row.get("parameter")), None


def _letter_lines(said: list, where: str | None = None) -> list:
    """The lines of letters a slide shows, `LETTER_LINES` at most, the last naming how many letters it leaves out."""
    lines = letters_lines(said)
    if len(lines) <= LETTER_LINES:
        return lines
    m = len(said)
    while m and len(letters_lines(said[:m], width=80)) > LETTER_LINES:
        m -= 1
    shown = letters_lines(said[:m], width=80)
    shown[-1] += f", and {len(said) - m} more in {where or 'the LETTERS of its file'}"
    return shown


def _height_of(j: int, row: dict, run: dict | None = None) -> float:
    """The height on a slide of line j of the equation: its head, its formula, and its letters where the formula is
    written in mathematics, a line of letters more for each beyond the first."""
    said, where = term_letters(run, j, row)
    extra = LETTER_STEP * (max(len(_letter_lines([x for x in said if x.split(" ")[0] not in ("f", "a")])), 1) - 1)
    lettered = (said or where) and str(row.get("formula") or "").startswith("$")
    head = f"$f_{{{j}}}$: " + ("the first model, " + str(row.get("added", "")) if j == 0 else str(row.get("added", "")))
    return 0.055 + 0.03 * two_lines(head, 64).count("\n") + (0.045 + extra if lettered else 0.0) + 0.07 + _rise(row.get("formula"))          # a fraction stands taller


def slides_equation(run: dict, k: int) -> list:
    """Phase 2, the equation of model k on as many slides as its terms ask, every term shown and the parameters of the
    fit on the last: a model of three terms left its third off its slide, "and 1 more terms, in equations.json", over
    the table (the report of the blind run of his message 310)."""
    out, first = [], 0
    while True:
        fig = slide_equation(run, k, first)
        out.append(fig)
        if fig.rest is None or fig.rest <= first:
            return out
        first = fig.rest


WORDS = 60                   # the words of a slide at most, his message 149 of the repository it is made in: a slide of the equation filled to its foot held 103


def _lettered(text: str) -> int:
    """The words of a text as a slide counts them: each token that holds a letter, so that a number is no word."""
    return sum(1 for w in str(text).split() if any(c.isalpha() for c in w))


def _words_of(j: int, row: dict, run: dict | None = None) -> int:
    """The words that line `j` of an equation puts on its slide: its heading, its value, the letters under it, and a
    formula's code and its note where it has no EQUATION."""
    letters, where = term_letters(run, j, row)
    letters = letters + ([f"its letters as {where} defines them"] if where else [])
    head = f"f: {'the first model, ' if j == 0 else ''}{row.get('added', '')}"
    letters = [x for x in letters if x.split(" ")[0] not in ("f", "a")]          # said once by the sentence at the top of the slide
    formula = str(row.get("formula") or "")
    shown = (formula + " no EQUATION in its model file: its code is shown" if formula and not formula.startswith("$")
             else "" if formula else "its formula is not written: a term of one's own takes formula= in its Mechanism")          # the code and its note count, as the slide shows them: uncounted, a slide of code reached 96 words (the slides read against their words, 0.1.3rc4)
    return _lettered(head) + (1 if j else 0) + (3 if row.get("held_at") else 0) + _lettered(" ".join(letters)) + _lettered(shown)


def slide_equation(run: dict, k: int, first: int = 0):
    """Phase 2, the slides of the equation that `slides_equation` draws after `slide_fixed`: the equation of the model, the
    first model and then each term on the model before it, with
    the value of its parameter, and the parameters of the fit with their bars. From line `first` of the equation, as many
    lines as the slide holds: the table where the equation ends on this slide, else none, and `fig.rest` the line that
    the next slide starts from, or None."""
    ident = run["phases"].get(k, {}).get("identify", {})
    chain = chain_of(run, k)
    alone = bool(first) and first >= len(chain)          # its caption in log.md the parameters, the twelfth reading of 0.1.3rc4. a slide of the table alone: "continued" promised an equation it did not hold (the slides read against their words, 0.1.3rc4)
    fig = _frame("The model: its shared parameters" if alone else TITLES["identify"] + (", continued" if first else ""), k, "identify", run, r"Each line turns $f$, the model above it, into the next, with $a$ its parameter." if not first else "")
    fig.rest = None
    fig.caption = "the shared parameters of the fit" if alone else None
    params = ident.get("parameters", {})
    if not chain:
        _no_data(fig, "This run left no equations.json: it was made before the package wrote the equation of each model.")
    y = TOP
    table = 3 + sum(_lettered(n) for n in list(params)[:6])          # the words of the table of the parameters: its heads and the names, its numbers no words
    box = 0.045 * (sum(1 + r[0].count("\n") for r in parameter_rows(params)[:6]) + 1) + 0.03          # by its lines: a name on two lines doubles its row, and six such rows ran into the phases at the foot          # the height of the table, by its rows: a fixed box of six rows sent a table of two to a slide of its own
    last = (y - sum(_height_of(j, row, run) for j, row in enumerate(chain) if j >= first) >= 0.14 + box
            and sum(_words_of(j, row, run) for j, row in enumerate(chain) if j >= first) + table <= WORDS)          # the rest of the equation ends on this slide, above the table, in its room and its words
    said_letters = set() if first else {"f", "a"}          # the sentence at the top of the first slide says f and a
    spent = 0
    for j, row in enumerate(chain):
        if j < first:
            continue
        if not last and y - _height_of(j, row, run) < 0.17 and j > first:          # a line that does not reach the foot goes on this slide, and the table on the next: a slide of two terms left its lower half blank (the cold round 22)
            fig.rest = j
            return fig
        if j > first and spent + _words_of(j, row, run) > WORDS:          # nor past the words of a slide
            fig.rest = j
            return fig
        spent += _words_of(j, row, run)
        head = f"$f_{{{j}}}$: " + ("the first model, " + row["added"] if j == 0 else row["added"])
        said = fig.text(0.04, y, two_lines(head, 64), fontsize=TEXT, fontweight="bold", va="top")          # on two lines with the sign of a cut: on one, no term was named (the slides read against their words, 0.1.3rc4)          # cut at a word: cut at 61 characters, "its Lorentzian wid..." named no term
        if j and (row.get("parameter") in params or row.get("held_at")):
            v = params.get(row["parameter"]) or {"value": row["held_at"][0], "standard_deviation": row["held_at"][1]}
            end = said.get_window_extent(renderer=fig.canvas.get_renderer()).x1 / fig.bbox.width          # the value beside the words, where they end
            fig.text(end + 0.015, y - 0.03 * two_lines(head, 64).count("\n"), f"$a$ = {pm(v['value'], v['standard_deviation'])}" + (", held from theory" if row.get("parameter") not in params else ""), fontsize=TEXT, va="top", color=ACCENT)          # a held term's value is its theory (the cold round 22)
        y -= 0.055 + 0.03 * two_lines(head, 64).count("\n")          # a second line of the head, counted by `_height_of`
        formula = row.get("formula") or ""
        if formula.startswith("$"):
            fig.text(0.07, y, formula, fontsize=TEXT + 5, va="top")          # the formula larger than its heading, as a displayed equation
            listed, where = term_letters(run, j, row)
            if listed or where:          # the letters of the first model from LETTERS of its file (the blind red team of 0.1.0rc12: h, x_0, w and b were defined nowhere)
                y -= 0.045 + _rise(formula)          # under a fraction, which stands taller: its letters lay on its denominator
                new = [x for x in listed if x.split(" ")[0] not in said_letters]          # each letter defined once on a slide (his message 311: fewer words a slide)
                said_letters.update(x.split(" ")[0] for x in new)
                shown = _letter_lines(new, where) if new else [f"its letters as {where} defines them"] if where else []          # nothing where every letter was said above it on the slide
                fig.text(0.07, y, "\n".join(shown), fontsize=SMALL, color=GREY, va="top")
                y -= LETTER_STEP * max(len(shown) - 1, 0)          # a line of letters more, its room counted by `_height_of`          # the letters of the formula on its slide, each typeset as mathematics
        elif formula:          # a model file with no EQUATION: its Python line, and how to have it typeset (his message 314)
            code = textwrap.wrap(formula, 64, break_on_hyphens=False) or [formula]
            shown = code[:2] if len(code) <= 2 else [code[0], textwrap.shorten(" ".join(code[1:]), 64, placeholder="...")]          # cut at 49 characters, ** 2 and + own[2] were left out under "its code is shown" (the slides read against their words, 0.1.3rc4)
            fig.text(0.07, y, "\n".join(shown), fontsize=TEXT - 2, family="monospace", va="top")
            fig.text(0.07, y - 0.04 - 0.035 * (len(shown) - 1), "no EQUATION in its model file: its code is shown" + (", cut" if len(code) > 2 else ""), fontsize=SMALL, color=GREY, va="top")
            y -= 0.035 * (len(shown) - 1)          # read as a promise that it was typeset (the blind red team of 0.1.0rc12, stumble 7)
        else:
            fig.text(0.07, y, "its formula is not written: a term of one's own takes formula= in its Mechanism", fontsize=SMALL, color=GREY, va="top")
        y -= 0.07 + (_rise(formula) if not any(term_letters(run, j, row)) else 0.0)          # a fraction stands taller, as `_height_of` counts it, its room taken above its letters where it has them
    if not last:          # the equation ended low on a slide with no room left for the table: the table on a slide of its own
        fig.rest = len(chain) if first < len(chain) else None
        if fig.rest is not None:
            return fig
    ax = fig.add_axes([0.04, 0.12, 0.62, box])
    ax.axis("off")
    rows = parameter_rows(params)
    if rows:
        tab = ax.table(cellText=rows[:6], colLabels=["parameter", "value of this turn", "bar"], loc="upper left", cellLoc="left", colLoc="left", colWidths=[0.46, 0.39, 0.15])          # the values of a turn read as results, which best.md alone quotes (the blind red team of 0.1.0rc16, mislead 11), and then a value printed under a head that said it was not quoted (the blind red team of 0.1.0rc19, mislead 1)
        tab.auto_set_font_size(False)
        tab.set_fontsize(SMALL)
        tab.scale(1.0, 1.9)
        for (r, c), cell in tab.get_celld().items():
            if r and "\n" in rows[r - 1][0]:          # a name on two lines, its row twice as high
                cell.set_height(1.8 * cell.get_height())
            cell.set_edgecolor("#d0d4da")
            if r == 0:
                cell.set_text_props(fontweight="bold")
        if len(rows) > 6:          # the seventh parameter of the CO2 model with HITRAN's lines left the table with no word (the slides read against their words, 0.1.3rc4)
            fig.text(0.68, 0.14, f"and {len(rows) - 6} more in loop/model_{k}.md", fontsize=SMALL, color=GREY, va="bottom")
    weak = ident.get("weakest")
    pairs = _correlated(ident.get("correlations", {}))
    told = [f"{a} and {b}" for a, b, r in pairs if abs(r) >= TOLD_APART]          # not told apart where best.md gives no bar, and correlated below: the slide said "not told apart" of two values that best.md quotes with bars (the cold round 22)
    close = [f"{a} and {b}, r = {r:.3f}" for a, b, r in pairs if abs(r) < TOLD_APART]
    note = (f"Least fixed: {weak}." if weak and len(params) > 1 else "") + (f" Not told apart: {'; '.join(told)}." if told else "") + (f" Correlated: {'; '.join(close)}." if close else "")          # the least fixed of one parameter says nothing (his message 311: fewer words)
    if note:
        fig.text(0.68, 0.12 + box, "\n".join(textwrap.wrap(note, 40)[:4]), fontsize=SMALL, va="top")          # beside the top of the table, whatever its rows
    return fig


def letters_of(formula: str, parameter: str | None = None) -> str:
    """The letters of a term's formula in words, those that stand in it, after f, the model before the term, and a, its
    parameter. A change along a setting names its parameter and its setting, read from the name of its parameter, as
    `width along concentration`."""
    return ", ".join(letters_list(formula, parameter))


MATH = {"ℓ": r"\ell", "Θ": r"\Theta", "σ": r"\sigma", "Γ": r"\Gamma", "G": r"G_\sigma", "L": r"L_\Gamma"}          # the symbols of the definitions of letters, as LaTeX writes them


def math_letter(said: str) -> str:
    """A definition of `letters_list` as the slide sets it, its symbol and every Greek letter in it typeset as
    mathematics: `ℓ the straight line` as `$\\ell$ the straight line` (his message 314: "use a proper latex writing for the
    formulas")."""
    head, _, rest = said.partition(" ")
    words = [f"${MATH[w.rstrip(',')]}${',' if w.endswith(',') else ''}" if w.rstrip(",") in MATH and w.rstrip(",") not in ("G", "L") else w for w in rest.split(" ")]
    return f"${MATH.get(head, head)}$ " + " ".join(words)


def letters_lines(said: list, width: int = 100) -> list:
    """The definitions of letters typeset by `math_letter` on lines of about `width` characters as they read, a line never
    cut inside a definition: a mathematical symbol cut in two does not typeset, and one line of six definitions ran off
    the slide (his message 314)."""
    import re as re_
    lines, plain = [], 0
    for x in said:
        typeset, length = math_letter(x), len(re_.sub(r"\$[^$]*\$", "xx", math_letter(x))) + 2          # a symbol reads about two characters wide
        if lines and plain + length > width:
            lines.append(typeset)
            plain = length
        else:
            lines[-1:] = [lines[-1] + ", " + typeset] if lines else [typeset]
            plain += length
    return lines


def letters_list(formula: str, parameter: str | None = None) -> list:
    """The letters of `letters_of`, each its own definition beginning with its symbol, so that a slide defines each
    letter once (his message 311)."""
    import re as re_
    said = ["f the model above it", "a the parameter of the line"]
    body = formula.replace("$", "")
    for pattern, words in LETTERS:
        if re_.search(pattern, body + " "):
            said.append(words)
    if "b^2 = a" in body:          # the kernel's s is the size that it was given, not a stretch
        said = [w.replace("s its stretch", "s the size that it was given") for w in said]
    if parameter and " along " in parameter:
        changed, setting = parameter.split(" along ", 1)
        said = [w.replace("p the parameter that it changes", f"p the {changed}").replace("z the setting,", f"z the {setting.rstrip('^2')},") for w in said]
    return list(said)          # each definition whole, with a comma inside it


def _correlated(correlations, beyond: float = 0.9) -> list:
    out = []
    if isinstance(correlations, dict):
        for key, r in correlations.items():
            if isinstance(r, (int, float)) and abs(r) >= beyond:
                a, _, b = str(key).partition(" and ") if " and " in str(key) else str(key).partition(", ")          # a pair named either way: the loop writes "a, b" (the run of his message 310 read "and ?")
                out.append((a, b or "?", float(r)))
    return out


def axis_names(run: dict) -> tuple:
    """The names of the abscissa and of the value of the curves, each with its unit, as the description gives them."""
    p = run.get("problem") or {}          # a run read without its description names its axes "abscissa" and "value"
    return (f"{p.get('abscissa', 'abscissa')}" + (f", {p['abscissa_unit']}" if p.get("abscissa_unit") else ""),
            f"{p.get('value', 'value')}" + (f", {p['value_unit']}" if p.get("value_unit") else ""))


def slide_effects(run: dict, k: int):
    """Phase 2, its last slide, after `slide_fixed` and the equation: what each term of the model does to the curves at
    the fit, the fitted model less the same
    model with that term at zero, on every record, beside one noise of a sample. For the first model, which has no term,
    each shared parameter moved by one bar of its fit."""
    rows = next((row for row in run["equations"] if row["model"] == f"model {k}"), None)
    arrays = run["effects"].get(k)
    first = k == 0 or not rows or not rows.get("effects") or not rows.get("after")          # a term moved by one bar where at zero it leaves no curve is no first model
    moved = not first and any(e.get("how") not in ("at zero", "held") for e in (rows or {}).get("effects", []))          # a term whose parameter is moved by one bar and not set at zero: the speed dependence of the public CO2 lines (the slides read against their words, 0.1.3rc4)
    sentence = ("The first model has no term. Each panel: the curves when one of its shared parameters moves by one bar of its fit, less the fitted curves, on every record. Grey: one noise, or its size in the title." if first else
                "Each panel: the fitted model less the same model with that term at zero, or one bar off where it says one bar, all else held, on every record. Grey: one noise, or its size in the title." if moved else
                "Each panel: the fitted model less the same model with that term at zero, all else held, on every record. Grey: one noise, or its size in the title.")
    fig = _frame("What each term does to the curves", k, "identify", run, sentence)
    if not rows or arrays is None or not rows.get("effects"):
        _no_data(fig, "This run left no effects: it was made before the package kept the effect of each term.")
        return fig
    effects = rows["effects"]
    size = {e["key"]: max(float(np.max(np.abs(arrays[f"e{e['key']}_{r}"]))) for r in range(len(rows["records"])) if f"e{e['key']}_{r}" in arrays) for e in effects}
    shown = sorted(effects, key=lambda e: -size.get(e["key"], 0.0))[:MOST_PANELS]
    n = len(shown)
    cols = 1 if n == 1 else 2 if n <= 4 else 3
    lines = math.ceil(n / cols)
    noise = float(np.median(rows["noise"])) if rows.get("noise") else 0.0
    w, h, left, top = 0.86 / cols, (TOP - 0.12) / lines, 0.1, TOP - 0.02
    room = {1: 60, 2: 40, 3: 26}[cols]          # the characters of a panel's title: on three columns titles of 44 ran into each other (the slides read against their words, 0.1.3rc4)
    keyed = len(rows["records"]) <= 6          # a legend of 25 records ran off the top of the slide over its title
    said_noise = any(noise and size.get(e["key"], 0.0) and not size.get(e["key"], 0.0) / 40 <= noise <= 5 * size.get(e["key"], 0.0) for e in shown)          # a title of four lines, the noise's size the last
    for i, e in enumerate(shown):
        c, r = i % cols, i // cols
        ax = fig.add_axes([left + c * w + 0.03, top - (r + 1) * h + 0.05, w - (0.11 if cols > 1 else 0.06), h - (0.17 if lines > 1 else 0.09) - 0.04 * bool(said_noise)])          # a line lower where a title says the noise: its fourth line met the sentence above it (the slides read against their words, 0.1.3rc4)          # room above a lower panel for a name on two lines          # room between panels for the ticks of the next: they lay on the panel before
        for j in range(len(rows["records"])):
            key = f"e{e['key']}_{j}"
            if key in arrays:
                ax.plot(arrays[f"x{j}"], arrays[key], color=COLOURS[j % len(COLOURS)], lw=1.6, label=label_of(run, rows["records"][j], room=22))
        biggest = size.get(e["key"], 0.0)
        far = ""
        if noise and biggest / 40 <= noise <= 5 * biggest:          # a band of a hundredth of the change was a hairline that no reader saw (the slides read against their words, 0.1.3rc4)
            ax.axhspan(-noise, noise, color=GREY, alpha=0.2, lw=0)
        elif noise and biggest:
            far = f"\nnoise: {noise / biggest:.2g} {'times' if noise > biggest else 'of'} the change"          # as narrow as the value above it: "one noise of a sample: 0.0105 times the largest change" ran into the next panel's title (the slides read against their words, 0.1.3rc4)          # under the title: inside the panel at 15 points it lay on the curves (row M583)
        ax.axhline(0.0, color=GREY, lw=0.8)
        moved = e["how"] not in ("at zero", "held")
        tail = f"+ one bar, {pm(e['value'], e['bar'])}"          # whole: cut at the panel's room, "+ one bar" and its value were the words lost (the slides read against their words, 0.1.3rc4)
        said = {"at zero": f"$a$ = {pm(e['value'], e['bar'])}", "held": f"$a$ = {pm(e['value'], e['bar'])}, held from theory"}.get(
            e["how"], tail)
        ax.set_title("\n".join((textwrap.wrap(e["parameter"] if moved else e["name"], room) or [""])[:1] + [textwrap.shorten(" ".join(textwrap.wrap(e["parameter"] if moved else e["name"], room)[1:]), room, placeholder="...")] * (len(textwrap.wrap(e["parameter"] if moved else e["name"], room)) > 1)) + "\n" + said + far, fontsize=SMALL, loc="left")
        if i == 0 and keyed:
            ax.legend(fontsize=SMALL, frameon=True, loc="lower right")          # which curve is which record (review 84, minor 3)
        if r == lines - 1 or i + cols >= n:          # the names of the axes, on the lowest panel of each column and the first of each line (review 85, minor 1, and the red team of message 310, M7)
            ax.set_xlabel(axis_names(run)[0], fontsize=SMALL)
        ax.tick_params(labelsize=SMALL)
    fig.text(0.012, top - lines * h / 2, "change of the " + axis_names(run)[1], fontsize=SMALL, rotation=90, va="center", ha="left")          # one name for every row: on two rows it stood past the slide's left edge
    if not keyed:
        fig.text(0.965, 0.905, f"a colour a record, again after {len(COLOURS)}", fontsize=SMALL, color=GREY, ha="right", va="top")          # under the phase, clear of the names of the axes
    if len(effects) > n:
        fig.text(0.965, 0.1, f"the {n} largest of {len(effects)}", fontsize=SMALL, color=GREY, ha="right")
    return fig


def slide_twin(run: dict, k: int):
    """The slide of phase 3 before row M600: the admission of the fit by its twin, each of its first twelve shared parameters'
    bias over its spread, spread over its bar, and the share of sets whose interval holds the truth, against the bands of
    the rule. The log no longer draws it: its phase 3 is `slide_admitted`, `SLIDES`. It stays for a reader who draws it
    alone, and for its tests."""
    adm = run["phases"].get(k, {}).get("twin", {}).get("admission", {})
    verdict = adm.get("admitted")
    fig = _frame(TITLES["twin"], k, "twin", run, f"Synthetic sets drawn from the fitted model, {adm.get('sets', '?')} of them, fitted again. The dots must stand in the white bands.")
    if not adm.get("parameters"):
        _no_data(fig, "The run kept no twin of this model.")
        return fig
    fig.text(0.965, 0.80, "admitted" if verdict else "not admitted", fontsize=TITLE - 2, fontweight="bold", ha="right", va="top", color=GOOD if verdict else BAD)
    bands = adm.get("bands", {})
    names = list(adm["parameters"])[:12]
    under = float(bands.get("bias_over_spread_under", 0.5))
    what = [("bias over spread", "bias_over_spread", (-under, under)),
            ("spread over bar", "spread_over_bound", tuple(bands.get("spread_over_bound", (0.6, 1.4)))),
            ("share of sets held", "coverage", tuple(bands.get("coverage", (0.83, 1.0))))]
    for i, (label, key, (lo, hi)) in enumerate(what):
        ax = fig.add_axes([0.14 + i * 0.29, 0.25, 0.23, TOP - 0.29])
        ys = np.arange(len(names))
        vals = [float((adm["parameters"][n] or {}).get(key, np.nan) or np.nan) for n in names]          # a parameter that the twin could not read is None: the log ended there under a twin that did not admit the fit (the cold round 22)
        ax.axvspan(lo, hi, color="white", zorder=0)
        ax.set_facecolor("#f2e3dc")
        ax.scatter(vals, ys, s=70, color=[GOOD if lo <= v <= hi else BAD for v in vals], zorder=3)
        ax.set_yticks(ys)
        ax.set_yticklabels(["\n".join(textwrap.wrap(n, 14)) for n in names] if i == 0 else [""] * len(names), fontsize=SMALL)          # wrapped: a long name was cut at the left edge (the cold round 22)
        ax.set_title(label, fontsize=TEXT)
        span = max(hi - lo, 1e-9)
        finite = [v for v in vals if math.isfinite(v)] or [lo]
        ax.set_xlim(min(lo - 0.3 * span, min(finite) - 0.1 * span), max(hi + 0.3 * span, max(finite) + 0.1 * span))
        ax.invert_yaxis()
    fig.text(0.04, 0.195, "\n".join(textwrap.wrap("Bias over spread: the mean fit of the sets less the value that drew them, in their scatter. "
                                                  "Spread over bar: their scatter over the bar of one fit. Share held: the sets whose interval holds that value. "
                                                  "Admitted: the fit returns what drew the sets, so its bar is honest if the model is right, which this does not test.", 116)),          # three lines above the phases: four of 104 characters reached them at 15 points (row M583)
             fontsize=SMALL, color=GREY, va="top", linespacing=1.3)          # the three measures were named and not said (the cold round 22)
    return fig


def chi_said(read: dict) -> str:
    """The chi square for each degree of freedom of a fit, with the bar of a right model, as the brief of the model gives
    it, `loop.chi_bar`: that of independent samples, and where neighbouring residuals correlate, that bar times the root
    of the samples over which they correlate. The bar of independent samples alone, 0.015, read 3.219 as 150 bars from 1,
    where the guide gave 0.096 (the blind red teams of 0.1.0rc19, mislead 2 and mislead 1)."""
    chi, dof = read.get("chi_square"), read.get("degrees_of_freedom")
    if not chi or not dof:
        return ""
    from .loop import chi_bar
    said = chi_bar(chi, dof, read.get("correlation_of_neighbours"), read.get("samples_of_the_correlation"))
    return said[:1].upper() + said[1:]


LEVEL = 0.03          # a level of a sweep spans no more than three per cent of the larger of its first value and the value that joins it, its gaps under a tenth of the largest gap between the setting's values over the run, `levels_of`: a measured temperature of 330.97 K and one of 331.02 K are one setting, his message 604, and the public CO2 line's 49.51 to 50.59 Torr, 2.1 per cent of the larger, is one pressure: at two per cent its 246 K curve stood alone, the grid was not whole and the row split its 50 Torr sweep (the slides read against their words, 0.1.3rc4)


def correlated(run: dict, k: int) -> float:
    """The samples over which the neighbouring residuals of model k correlate, as `loop.chi_bar` reads them from the phase
    read, where their correlation passes 0.3, else 1: the band of a right model of independent samples, times its root,
    is that of a right model whose noise correlates as these residuals do. The slides of the fit and of what it leaves
    drew the band of independent samples alone, 1.96 standard deviations of 0.020, 0.04, where `chi_said` gave one
    standard deviation of 0.28 for samples that correlate over 192, so 0.55 at 1.96 (the eighth reading of 0.1.3rc4)."""
    read = (run.get("phases") or {}).get(k, {}).get("read", {})          # a run read without its notebook has no phases
    rho, span = read.get("correlation_of_neighbours"), read.get("samples_of_the_correlation")
    return float(span) if isinstance(rho, (int, float)) and abs(rho) > 0.3 and isinstance(span, (int, float)) and span > 1 else 1.0


def levels_of(values) -> dict:
    """Each value of a setting mapped to its level. Two values are one level where their gap is under a tenth of the
    largest gap between the values, so small against the steps of the sweep itself, and the level spans no more than
    LEVEL of the larger of its first value and the value that joins it. So the jitter of a measured setting, 330.97 K and 331.02 K, is one level, and a fine
    sweep in steps of 3 K is thirty (review 265, moderate 1: a gap judged against LEVEL alone, value by value, chained
    them into one). The level is named by the median of its values."""
    ordered = sorted(set(values))
    widest = max((b - a for a, b in zip(ordered, ordered[1:])), default=0.0)
    groups = []
    for v in ordered:
        if groups and v - groups[-1][-1] <= 0.1 * widest and abs(v - groups[-1][0]) <= LEVEL * max(abs(v), abs(groups[-1][0]), 1e-12):
            groups[-1].append(v)
        else:
            groups.append([v])
    return {v: float(np.median(g)) for g in groups for v in g}


def crowding(numeric: dict, level: dict, keys: list) -> dict:
    """How many curves stand in each cell of the grid of settings, each setting read at its level."""
    out = {}
    for c, s in numeric.items():
        cell = tuple(level[k].get(s.get(k)) for k in keys)
        out[cell] = out.get(cell, 0) + 1
    return out


def sweeps_of(run: dict, curves: list) -> list:
    """The curves in sweeps, row M583 (the owner's messages 118, 135, 142 and 149: "each sweep in one plot"): a sweep is
    the curves that share every numeric setting of the description but one, each setting read within its level, `levels_of`:
    values one level where their gap is under a tenth of the largest gap of the setting and the level spans no more than
    three per cent of the larger of its first value and the value that joins it, ordered along that one, as (the setting, the
    settings that are fixed, the curves). A curve may stand in two sweeps, as the corner of an L of settings does. Where no
    setting varies, or one alone, every curve is one sweep, and a curve that shares no sweep stands alone in its own."""
    numeric = {c: {k: v for k, v in run["settings"].get(c, {}).items() if isinstance(v, (int, float))} for c in curves}
    keys = list(dict.fromkeys(k for s in numeric.values() for k in s))
    level = {k: levels_of([s[k] for s in numeric.values() if k in s]) for k in keys}
    exact = {k: {v: v for v in level[k]} for k in keys}
    if max(crowding(numeric, level, keys).values(), default=1) > max(crowding(numeric, exact, keys).values(), default=1):
        level = exact          # a merge that puts two curves of other settings in one cell of the grid joined a sweep, not the jitter of a measure (review 266, moderate 1)
    varying = [k for k in keys if len(set(level[k].values())) > 1]
    if len(varying) <= 1:
        along = varying[0] if varying else None
        return [(along, {}, sorted(curves, key=lambda c: numeric[c].get(along, 0.0)) if along else list(curves))]
    out, placed = [], set()
    for along in varying:
        groups = {}
        for c in curves:
            if along in numeric[c]:
                groups.setdefault(tuple((k, level[k].get(numeric[c].get(k))) for k in varying if k != along), []).append(c)
        for fixed, members in groups.items():
            if len(members) > 1:
                out.append((along, dict(fixed), sorted(members, key=lambda c: numeric[c][along])))
                placed.update(members)
    return out + [(None, numeric[c], [c]) for c in curves if c not in placed]


def in_a_row(run: dict, curves: list) -> list:
    """The curves in the order of the owner's message 118, for any sweeps: "start from left with 70C-225mW, then increase
    temperature up to 130C-225mW, then decrease power down to 130C-25mW". Along each sweep from its far end to the curve
    that it shares with the next, then along the next away from it, each curve once and the curves of no sweep after.
    Where sweeps cross, a curve stands in the first sweep that reaches it, and a later sweep shows only its curves not yet
    shown. In a full grid of two settings the row runs along the sweeps of one setting, each whole (the slides read
    against their words, 0.1.3rc4)."""
    sweeps = sweeps_of(run, curves)
    out = []
    for along, _, members in reversed(sweeps):
        if not along or len(members) < 2:
            continue
        if out and members[-1] in out and members[0] not in out:          # the next sweep starts where the last ended, at the curve they share
            members = members[::-1]
        out += [c for c in members if c not in out]
    return out + [c for c in curves if c not in out]


def _set_on(a) -> tuple:
    """A curve set on the others as the toy sets the traces of a sweep (the owner's message 118, row M583): its base off,
    over its height, about the abscissa of its extremum. The base is the mean of the fitted curve over a twentieth of
    its samples at each end, and its height the fitted value farthest from that base less the base: read from the fitted
    curve alone, for a model of any field. Returns the abscissa, the values, the fitted values and the residuals over
    their noise, in the order of the abscissa."""
    order = np.argsort(a[:, 0])
    x, y, f = a[order, 0], a[order, 1], a[order, 2]
    z = a[order, 3] / np.where(a[order, 4] > 0, a[order, 4], 1.0)
    edge = max(1, x.size // 20)
    base = 0.5 * (f[:edge].mean() + f[-edge:].mean())
    j = int(np.argmax(np.abs(f - base)))
    height = f[j] - base if f[j] != base else 1.0
    return x - x[j], (y - base) / height, (f - base) / height, z


def _cells(x, values, n: int = 60) -> tuple:
    """`values` in `n` cells of one width along the abscissa: the middle of each cell that holds samples, the mean of its
    values and their count."""
    edges = np.linspace(x.min(), x.max(), n + 1)
    k = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, n - 1)
    held = np.unique(k)
    return 0.5 * (edges[held] + edges[held + 1]), np.array([values[k == i].mean() for i in held]), np.array([(k == i).sum() for i in held])


def _strip(ax, x, z, colour, reach: float, n: int = 60) -> None:
    """What the fit leaves on one curve over its noise, summed in `n` cells of the abscissa and divided by the root of
    their count, so that each bar has a scatter of one noise on a right model of independent noise and stands inside the
    grey band of 1.96 on 19 cells of 20 (the toy's strips, the owner's message 118)."""
    middle, mean, count = _cells(x, z, n)
    ax.axhspan(-1.96, 1.96, color=GREY, alpha=0.18, lw=0)
    ax.bar(middle, mean * np.sqrt(count), width=(x.max() - x.min()) / n, color=colour, lw=0)
    ax.axhline(0.0, color=GREY, lw=0.8)
    _seen(ax, reach)
    ax.tick_params(labelbottom=False, labelsize=SMALL)


SEEN = 20.0          # past this many noises an axis turns to symlog, linear within 1.96 and logarithms past it, as the slide of where does: on axes of 75 to 14 000 noises the band of a right model, 1.96, was a hairline (the slides read against their words, 0.1.3rc4)


def _seen(ax, reach: float, low: float | None = None, axis: str = "y", linthresh: float = 1.96) -> None:
    """An axis on which the band of a right model can be seen: where its reach passes `SEEN` noises, symlog, linear
    within `linthresh`, 1.96 for noises, and a scale of logarithms past it, which the slides say as "log past 2". Else
    linear."""
    top = max(abs(reach), abs(low if low is not None else reach))
    if top > SEEN:
        from matplotlib.ticker import FixedLocator, FuncFormatter
        (ax.set_yscale if axis == "y" else ax.set_xscale)("symlog", linthresh=linthresh, linscale=0.6)
        mark = 10.0 ** math.floor(math.log10(top))
        line = ax.yaxis if axis == "y" else ax.xaxis
        box = ax.get_position()
        short = (box.height * ax.figure.get_figheight() if axis == "y" else box.width * ax.figure.get_figwidth()) * 72 < 110          # an axis under 110 points: 1.96 touched 10 on the strips of the row of curves
        line.set_major_locator(FixedLocator([-mark, 0.0, linthresh] + [mark] * (not bool(short))))          # four ticks, three on a short axis: a symlog axis gave six, and 1 and -1 lay on one another
        line.set_major_formatter(FuncFormatter(lambda v, _: "" if v == 0 else f"{v:g}"))          # the edge of the band named: between 0 and 100 it read as 25 and not 1.96 (the slides read against their words, 0.1.3rc4)
    (ax.set_ylim if axis == "y" else ax.set_xlim)(-reach if low is None else low, reach)


def _logged(label: str, reach: float) -> str:
    """The name of an axis of `_seen`, with its scale where it turns to logarithms: the band of 1.96 read as 5 to 25 against
    ticks of 10 and 100 where nothing said so (the slides read against their words, 0.1.3rc4)."""
    return label + ", log past 2" if abs(reach) > SEEN else label


def _reach(arrays) -> float:
    """The half height of the strips of a slide, one for all of them: 4.5 noises, or more where a bar reaches past."""
    bars = [np.abs(m * np.sqrt(c)).max() for x, z in arrays for _, m, c in [_cells(x, z)]]
    return max(4.5, 1.1 * max(bars, default=0.0))


def slide_row(run: dict, k: int):
    """Phase 4, the first kind of slide of the owner's message 118: every curve in a row, in the order of its sweeps
    (`in_a_row`), each set on one scale with its fit and what the fit leaves under it. Ten curves a slide at most, at letters of
    15 points, on as many slides as the curves need (his message 604)."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator          # matplotlib is asked only by the log
    res = run["residuals"].get(k, {})
    fig = _frame(TITLES["fit"] + ": the curves in a row", k, "fit", run,
                 "Each curve with its fit, once, along its sweeps, a curve of two sweeps in the first, base off, over its height. Under it what the fit leaves in noises, summed in cells: " + ("grey, a right model in white noise." if correlated(run, k) > 1 else "inside the grey band on a right model."))
    if not res:
        _no_data(fig, "The run kept no residuals of this model.")
        return fig
    curves = in_a_row(run, list(res))
    lengths = {len(members) for along, _, members in sweeps_of(run, curves) if along and len(members) > 1}
    row = next(iter(lengths)) if len(lengths) == 1 and next(iter(lengths)) <= 10 and len(curves) > 10 else 0          # the sweeps of a grid, all of one length: a page of whole sweeps, where pages of nine split a sweep of five across two (the slides read against their words, 0.1.3rc4)
    n_pages = max(math.ceil(len(curves) / 10), 1)          # never a division by nought (review 265)
    size = math.ceil(len(curves) / n_pages)          # pages of equal size, ten panels at most: the rest were left to the report, fifteen of the twenty five public CO2 curves, his message 604
    if row:
        size = row * (10 // row)
        n_pages = math.ceil(len(curves) / size)
    figs = []
    for page_n, start_at in enumerate(range(0, len(curves), size)):
        shown = curves[start_at:start_at + size]
        if n_pages > 1:          # which curves, in the title: under the plots the words met the name of the project
            plt.close(fig)
            fig = _frame(TITLES["fit"] + f": the curves in a row, {start_at + 1} to {start_at + len(shown)} of {len(curves)}", k, "fit", run,
                         "Each curve with its fit, once, along its sweeps, a curve of two sweeps in the first, base off, over its height. Under it what the fit leaves in noises, summed in cells: " + ("grey, a right model in white noise." if correlated(run, k) > 1 else "inside the grey band on a right model."))
        set_on = {c: _set_on(res[c]) for c in shown}
        reach = _reach([(x, z) for x, _, _, z in set_on.values()])
        w = 0.9 / len(shown)
        for i, curve in enumerate(shown):
            x, y, f, z = set_on[curve]
            shades = [c for c in COLOURS if c != "#000000"]          # not black under the black fit (the slides read against their words, 0.1.3rc4)
            colour = shades[i % len(shades)]
            top = fig.add_axes([0.07 + i * w, 0.36, w - 0.03, TOP - 0.42])          # apart by 0.03: the last tick of a panel met the first of the next
            low = fig.add_axes([0.07 + i * w, 0.15, w - 0.03, 0.17], sharex=top)
            middle, mean, _ = _cells(x, y)
            top.plot(middle, mean, ".", ms=4, color=colour)
            top.plot(x, f, "-", color="black", lw=1.0)
            top.set_ylim(-0.25, 1.25)
            top.set_title(label_of(run, curve, room=24).replace(", ", "\n"), fontsize=SMALL, loc="left")          # every setting, one a line: "750 Torr" alone named five curves of a grid
            top.tick_params(labelbottom=False, labelleft=i == 0, labelsize=SMALL)
            _strip(low, x, z, colour, reach)
            low.tick_params(labelbottom=True, labelleft=i == 0)
            low.xaxis.set_major_locator(MaxNLocator(2))
            if len(shown) > 6:          # two numbers a curve, 0 and its far end: on ten narrow plots three ran together, and none left the unit of the name with no scale (the slides read against their words, 0.1.3rc4)
                lo, hi = low.get_xlim()
                ends = [v for v in MaxNLocator(3).tick_values(lo, hi) if lo <= v <= hi]
                far = max(ends, key=abs, default=None)
                low.set_xticks([0.0, far] if far and lo <= 0.0 <= hi else [far] if far else [])
                low.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: "" if v == 0 else f"{v:g}"))          # the far end alone in numbers, 0 the extremum that the name of the axis says: "-0.15 0" read "-0.150" (the eighth reading of the slides of 0.1.3rc4)
            if i == 0:
                top.set_ylabel("over its height", fontsize=SMALL)
                low.set_ylabel(_logged("in noise", reach), fontsize=SMALL)
        fig.text(0.5, 0.075, f"{axis_names(run)[0]}, from each curve's extremum", fontsize=SMALL, ha="center", va="bottom")          # the abscissa named: on pages of narrow plots it had no name and no numbers (the slides read against their words, 0.1.3rc4)
        figs.append(fig)
    return figs[0] if len(figs) == 1 else figs

def slide_fit(run: dict, k: int):
    """The slides of phase 4 before row M600, the other kinds of slide of the owner's message 118: each sweep alone on a
    slide of its own (`sweeps_of`), its curves together on one axis, each set on the others (`_set_on`), and above them a
    strip for each curve of what the fit leaves. The log no longer draws them: its phase 4 is `slide_fitted` and
    `slide_row`, `SLIDES`. They stay for a reader who draws one sweep alone, and for their tests. Returns one slide, or a
    list of one a sweep."""
    res = run["residuals"].get(k, {})
    read = run["phases"].get(k, {}).get("read", {})
    said = chi_said(read)
    whole = "One fit of every curve together: its shared parameters are one value for all. " + said
    sentence = whole if len(textwrap.wrap(whole, WIDTH)) <= 2 else said          # the frame shows two lines, and the bar of correlated neighbours was cut at "if its noise" (the blind red team of 0.1.0rc20)
    if not res:
        fig = _frame(TITLES["fit"], k, "fit", run, sentence)
        _no_data(fig, "The run kept no residuals of this model.")
        return fig
    units = run.get("units", {})
    figs = []
    for along, fixed, members in sweeps_of(run, list(res)):
        head = (f"along {along}" if along else "one curve" if len(members) == 1 else f"{len(members)} curves") + ("" if not fixed else ", " + ", ".join(f"{key} {v:.3g}" for key, v in list(fixed.items())[:1]))
        fig = _frame(TITLES["fit"] + (f": {head}" if along else ""), k, "fit", run, sentence)
        set_on = {c: _set_on(res[c]) for c in members}
        reach = _reach([(x, z) for x, _, _, z in set_on.values()])
        n = len(members)
        h = min(0.07, 0.30 / n)          # a strip for each curve, above them: at most 0.30 of the slide for all
        low_top = TOP - n * h - 0.035
        main = fig.add_axes([0.08, 0.17, 0.76, low_top - 0.17])          # its label above the phases at the foot
        for j, curve in enumerate(members):
            x, y, f, z = set_on[curve]
            colour = COLOURS[j % len(COLOURS)]
            value = run["settings"].get(curve, {}).get(along) if along else None
            name = (f"{value:.3g} {units.get(along, '')}".strip() if isinstance(value, (int, float)) else label_of(run, curve, room=16))
            middle, mean, _ = _cells(x, y)
            main.plot(middle, mean, "o", ms=3.5, color=colour, alpha=0.8, label=name)
            main.plot(x, f, "-", color=colour, lw=1.6)
            strip = fig.add_axes([0.08, TOP - (j + 1) * h, 0.76, h - 0.006], sharex=main)
            _strip(strip, x, z, colour, reach)
            strip.set_yticks([])
            fig.text(0.85, TOP - (j + 0.5) * h, name, fontsize=SMALL, color=colour, va="center")
            if j == 0:
                strip.set_title(head[:40], fontsize=SMALL, loc="left")
        fig.text(0.045, TOP - 0.5 * n * h, "in noise", fontsize=SMALL, rotation=90, ha="center", va="center")          # the unit of the strips, whose band is 1.96 either side
        main.tick_params(labelsize=SMALL)
        main.set_ylabel("over its height", fontsize=SMALL)
        main.set_xlabel("from the extremum of each curve", fontsize=SMALL)
        if n > 1:
            main.legend(fontsize=SMALL, frameon=True, loc="upper right", markerscale=1.5, handletextpad=0.2, borderpad=0.3)
        figs.append(fig)
    return figs[0] if len(figs) == 1 else figs


def slide_read(run: dict, k: int):
    """The slide of phase 5 before row M600: the residuals over the noise of each sample against the normal law, and the
    mean of their square in bands of the height of the curve, which is 1 where the noise is the one stated. The log no
    longer draws it: its phase 5 is `slide_left`, `SLIDES`. It stays for a reader who draws it alone, and for its tests."""
    res = run["residuals"].get(k, {})
    table = run["tables"].get(k, {})
    fig = _frame(TITLES["read"], k, "read", run, "Left: what the fit leaves over the noise, against the normal law. Right: its mean square along the height of the curve, 1 where the model and the noise are right, above 1 where either is wrong.")          # above 1 read as a wrong noise alone (the blind red team of 0.1.0rc12)
    if not res:
        _no_data(fig, "The run kept no residuals of this model.")
        return fig
    z = np.concatenate([a[:, 3] / np.where(a[:, 4] > 0, a[:, 4], 1.0) for a in res.values()])
    ax = fig.add_axes([0.07, 0.16, 0.38, TOP - 0.18])
    ax.hist(z, bins=np.linspace(-5, 5, 61), density=True, color=ACCENT, alpha=0.7)
    u = np.linspace(-5, 5, 400)
    ax.plot(u, np.exp(-0.5 * u * u) / math.sqrt(2 * math.pi), color="black", lw=1.6)
    ax.set_xlabel("residual over noise")
    ax.set_ylabel("share of the samples, a density")          # no label on the axis of the histogram (the cold round 22)
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1)
    bands = (table.get("noise") or {}).get("bands") if isinstance(table.get("noise"), dict) else None
    ax2 = fig.add_axes([0.57, 0.16, 0.38, TOP - 0.18])
    if bands:
        mids = [0.5 * (b["height_from"] + b["to"]) for b in bands]
        ax2.errorbar(mids, [b["mean_squared_residual"] for b in bands], yerr=[b.get("standard_error", 0) for b in bands], fmt="o", ms=8, color=ACCENT, capsize=4)
        ax2.axhline(1.0, color="black", lw=1.2)
        ax2.set_xlabel("height of the curve")
        ax2.set_ylabel("mean square over the noise")
    else:
        _no_data_in(ax2, "no bands of the noise in this run")
    return fig


def _no_data_in(ax, words: str) -> None:
    ax.axis("off")
    ax.text(0.5, 0.5, words, ha="center", va="center", color=GREY, fontsize=SMALL, transform=ax.transAxes)


def answer_label(a: dict, room: int = 40, lines: int = 2) -> str:
    """An answer of the report as the slide of where to search names its bar: its question and its answer on two lines
    of `room` characters at most, cut only past them. 0.1.0rc10 cut each at 36 characters, so that "a Gaussian width of
    0.03 widths" and "of 0.3 widths" read alike (review 85 of the repository it is made in, minor 1). Where the whole does
    not fit and the answer names a family before a colon, the family alone, whole: its statistic after the colon was the
    words cut away (the slides read against their words, 0.1.3rc4)."""
    said = textwrap.wrap(f"{a['question']}: {a['answer']}", room)
    if len(said) > lines and ": " in str(a["answer"]):
        said = textwrap.wrap(f"{a['question']}: {str(a['answer']).split(': ', 1)[0]}", room)
    return "\n".join(said[:lines]) + ("..." if len(said) > lines else "")


def slide_where(run: dict, k: int):
    """Phase 6: the answers of the report, each with its excess, the line that a right model passes on 5 sets of 100,
    and the three proposals for the model and the three for the experiment, those for the experiment that the slide has
    no room for whole on a page of their own after it: under the proposals for the model and their term files, no slide
    of eight read had room for one (the slides read against their words, 0.1.3rc4). Returns one slide, or two."""
    table = run["tables"].get(k, {})
    where = table.get("where", {})
    every = sorted([a for a in where.get("answers", []) if isinstance(a.get("excess"), (int, float)) and math.isfinite(a["excess"])], key=lambda a: -a["excess"])
    places, setting, _, _, _ = place_read(table)
    key = (("changes along", f"{places[0]} along {setting}") if setting else ("where", places[0])) if places else None
    pick = next((a for a in every if key and (a.get("question"), a.get("answer")) == key), None)          # the answer whose excess slide 6a writes: it stood on none of these bars (the slides read against their words, 0.1.3rc4)
    answers = every[:4] if pick is None or pick in every[:4] else every[:3] + [pick]          # four of three lines at 15 points, the rest in the report: eight of two lines of 40 characters left the slide at its left edge, and five of three lines touched (row M583)
    line = where.get("largest_excess_on_the_twin")
    fig = _frame(TITLES["propose"], k, "propose", run, (f"The {len(answers)} largest of the report's {len(every)} answers, the largest first. " if answers == every[:len(answers)] else f"The 3 largest of the report's {len(every)} answers, then the place of the slide before, rank {every.index(answers[-1]) + 1}. ") + "The line: the excess a right model passes on 5 sets of 100, all answers together. Blue beyond it, grey within.")          # the colours were said nowhere (the cold round 22)
    if not answers:
        _no_data(fig, "The run kept no report of this model.")
        return fig
    ax = fig.add_axes([0.3, 0.22, 0.24, TOP - 0.22])          # its label in two lines under it, above the phases
    ys = np.arange(len(answers))
    values = [max(a["excess"], 1e-2) for a in answers]
    ax.barh(ys, values, color=[ACCENT if a.get("beyond_the_twin") else GREY for a in answers])
    ax.set_yticks(ys)
    ax.set_yticklabels([answer_label(a, 30, 3) for a in answers], fontsize=SMALL)
    ax.invert_yaxis()
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter          # matplotlib is asked only by the log
    ax.set_xscale("log")
    decades = math.log10(max(max(values), line or 1.0) / max(min(min(values), line or math.inf), 1e-2))          # the span of the axis, the red line among it: bars that all stood near 10 000 beside a line at 2 read as half a decade, and their ticks ran together again (the blind red team of 0.1.0rc13)
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)) if decades <= 2.0 else LogLocator(base=10, subs=(1.0,), numticks=4))          # four decades at most: "10 100 1000100001000000" ran together over five          # 1, 2 and 5 of each decade over two decades at most: over four, "5 1020 50 100200" ran together (the blind red team of 0.1.0rc12)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}" if decades <= 2.0 else f"$10^{{{round(math.log10(v))}}}$"))          # one form on an axis: "1000100000 1e+071e+09" mixed two
    ax.xaxis.set_minor_formatter(NullFormatter())          # the labels of the minor ticks lay on one another at 30 sets (the cold round 22)
    if line:
        ax.axvline(line, color=BAD, lw=2)
    ax.set_xlabel("each answer's excess over its\ntwin, 1 at its 95th percentile", fontsize=SMALL)          # two lines as wide as its axis: two of 46 characters reached the words on its right (row M583)          # in two lines under its axis: aligned on the right of a narrow axis, its first words left the slide (the blind red team of 0.1.0rc13)          # the line is the largest over every answer, the unit is each answer's own (the blind red team of 0.1.0rc12)
    props = table.get("proposals", {})
    renderer = fig.canvas.get_renderer()

    def put(x, y, text, **style):          # a text at y, and the y under it by its measured height: three proposals of three lines each left the slide by fixed steps (review 85, moderate 1)
        said = fig.text(x, y, text, va="top", **style)
        return y - said.get_window_extent(renderer=renderer).height / fig.bbox.height

    beyond = any(a.get("beyond_the_twin") for a in answers)
    heading = ("For the model, by each term's fall of chi square, not by the bars" if beyond else "For the model, though none is asked")          # the bars put the core first and this list the flanks, with no word that they rank two things (the fresh red team of 0.1.0rc20, mislead 1)
    y = put(0.58, TOP, "\n".join(textwrap.wrap(heading, 38)), fontsize=SMALL + 1, fontweight="bold") - 0.015          # three proposals under grey bars read as asked (the cold round 22)
    for m, p in enumerate(props.get("for_the_model", [])[:3]):          # not k, the model's number, which the page after it names (the slides read against their words, 0.1.3rc4)
        y = put(0.59, y, "\n".join(textwrap.wrap(f"{m + 1}. {p.get('term', '')}", 44)[:2]), fontsize=SMALL)          # each proposal under its own words, its term file beneath where it has one (review 84, minor 3)
        if p.get("as_a_term_file"):
            y = put(0.61, y - 0.004, "\n".join(textwrap.wrap(p["as_a_term_file"], 40, break_on_hyphens=False)[:2]), fontsize=SMALL, family="monospace", color=GREY)          # whole, where 0.1.0rc10 cut it in a word
        y -= 0.012
    if any("width" in str(p.get("term", "")) for p in props.get("for_the_model", [])[:3]):          # three widths on one page, a parameter, the unit of the names and a deviation, and the names said none (the blind red team of 0.1.0rc14 on a fresh machine)
        y = put(0.59, y, "\n".join(textwrap.wrap("A width in these names: the full width of the fitted curves at half their height, not a parameter of model.py.", 46)), fontsize=SMALL, color=GREY) - 0.012
    y = put(0.58, y - 0.02, "For the experiment", fontsize=SMALL + 1, fontweight="bold") - 0.015
    experiment = props.get("for_the_experiment", [])[:3]
    for i, p in enumerate(experiment):
        said = fig.text(0.59, y, "\n".join(textwrap.wrap(experiment_line(p), 46)), fontsize=SMALL, va="top")
        if y - said.get_window_extent(renderer=renderer).height / fig.bbox.height < 0.13:          # by its measured height, above the name of the run and the phases at the foot
            said.remove()
            fig.text(0.59, y, f"{'all' if i == 0 else 'the rest'} on the next slide", fontsize=SMALL, color=GREY, va="top")
            break
        y -= said.get_window_extent(renderer=renderer).height / fig.bbox.height + 0.012          # whole: four lines at most cut the ratio that the twin measured (the cold round 22)
    else:
        return fig
    more = _frame(TITLES["propose"] + ": for the experiment", k, "propose", run,
                  "The proposals of the report for the experiment that the slide before had no room for, whole, with what the twin measured of each.")
    y = TOP
    for j, p in enumerate(experiment[i:], start=i + 1):
        said = more.text(0.06, y, "\n".join(textwrap.wrap(f"{j}. {experiment_line(p)}", 90)), fontsize=TEXT, va="top")
        y -= said.get_window_extent(renderer=more.canvas.get_renderer()).height / more.bbox.height + 0.03
    return [fig, more]


PLACE_ROWS = 6          # the curves of one slide of where, his message 604
BIG, LABEL, TICK = 28, 22, 20          # the plotted slides of row M600 (his messages 586 to 588): one line of words at 28 points, axes at 22, ticks at 20


def width_of(x, f) -> float:
    """The full width at half height of a curve set on one scale by `_set_on`, its top at 1 and its base at 0 at the
    abscissa 0: between its two crossings of half its height about the top, each found linearly between the samples on
    either side of it. From the samples outside the crossings, a Lorentzian sampled every 0.05 read 4.10 for 4."""
    j = int(np.argmin(np.abs(x)))
    low = np.flatnonzero(f[:j] < 0.5)
    high = np.flatnonzero(f[j:] < 0.5)
    if low.size:
        i = low[-1]
        left = x[i] + (0.5 - f[i]) * (x[i + 1] - x[i]) / (f[i + 1] - f[i])
    else:
        left = x[0]
    if high.size:
        i = j + high[0]
        right = x[i - 1] + (f[i - 1] - 0.5) * (x[i] - x[i - 1]) / (f[i - 1] - f[i])
    else:
        right = x[-1]
    return float(right - left) if right > left else float(x[-1] - x[0])


def place_ticks(reach: float, tall: bool) -> list:
    """The ticks of the top strip of slide 6a, symlog within 1.96 noises: the edge of the band always, and past 10 noises
    the decade under the reach on both sides, its upper one only on a strip of 110 points or more. Under 10 noises the
    two edges of the band alone: a decade of 2 printed its 2 on the 1.96 of the band's edge (the twentieth reading of
    0.1.3rc4)."""
    if reach < 10:
        return [-1.96, 0.0, 1.96]
    mark = 10.0 ** math.floor(math.log10(reach))
    return [-mark, 0.0, 1.96] + [mark] * bool(tall)


def place_read(table: dict) -> tuple:
    """The place where the report reads the residuals beyond the twin, as `(places, setting, excess, line, seen_by)`: the
    places that lead, tied ones among them, the setting along which the place changes where it changes along one, the
    excess of the first place, the excess that a right model passes on 5 sets of 100, the 95th percentile on the twin, and the family of
    statistics that sees most. No place where none stands beyond the twin."""
    rep = table.get("where") or {}
    first = rep.get("first_of_each_question") or {}
    answers = rep.get("answers", [])
    along = rep.get("where_read_along")
    places, setting = [], None
    if along and " along " in along:
        place, setting = along.split(" along ", 1)
        places = [place]
        excess = next((a["excess"] for a in answers if a["question"] == "changes along" and a["answer"] == along), None)
    else:
        if first.get("where"):
            places = [first["where"]] + [q for q in ((rep.get("tied_with_the_first") or {}).get("where") or []) if q != first["where"]]
        excess = next((a["excess"] for a in answers if a["question"] == "where" and places and a["answer"] == places[0]), None)
    return places, setting, excess, rep.get("largest_excess_on_the_twin"), first.get("seen by")


def slide_place(run: dict, k: int):
    """Phase 6 as one plot, row M600: every curve's residuals over its noise, summed in cells of the distance from its
    centre in widths of the curve and divided by the root of their count, inside the grey band that a right model of
    independent noise fills on 19 cells of 20, with the place that the report reads beyond the twin shaded where it is the core, the flanks or
    the wings, a place of another kind, as the asymmetry, named in the title and left unshaded, the cells beyond the
    band in colour, and the excess of that place over the twin written on the plot with the family of statistics that
    sees the residuals most, the report's first answer to its question seen by, as long as the room beside the note of
    the band holds it. Where more than six curves make a grid of more than three sweeps along one setting, `_place_grid` draws them
    instead, one colour map a level, which shades no place and draws no grey band, the core and the flanks dotted, its
    colours pale inside the band of a right model, and writes the excess and its line but not the family that saw it.
    Where no place stands beyond the twin alone and the report reads one along a setting, `where_read_along`, the title
    says that the place changes along it, the curves stand in its order, and the excess is that of the answer of the
    question changes along (the eighth reading of 0.1.3rc4). With no place, the title says whether another answer of
    the report stands beyond the twin. The strips are on a scale of logarithms past 2 noises, and their unit says so."""
    import matplotlib.pyplot as plt
    from .where import CORE, WINGS          # the places of the report, in widths of the curve
    res = run["residuals"].get(k, {})
    places, setting, excess, line, seen = place_read(run["tables"].get(k, {}))
    title = (f"{ORDER.index('propose') + 1}. Beyond the twin: {places[0]}, changing along {setting}" if places and setting
             else f"{ORDER.index('propose') + 1}. Beyond the twin's band: {' and '.join(places)}") if places else (
        f"{ORDER.index('propose') + 1}. No place beyond the twin, other answers beyond it" if any(a.get("beyond_the_twin") for a in ((run["tables"].get(k, {}).get("where") or {}).get("answers") or []))
        else f"{ORDER.index('propose') + 1}. No answer of the report beyond its twin")          # from every answer: read from the place alone, the title said the residuals stayed inside where the shape and the noise stood beyond (the eleventh reading of 0.1.3rc4)          # short: "where the residuals leave the twin's band: the asymmetry" ran off the slide on the public CO2 lines, his 604
    fig = plt.figure(figsize=SIZE)
    fig.text(0.035, 0.95, title, fontsize=BIG, fontweight="bold", va="top")
    curves = in_a_row(run, list(res))
    if setting:
        curves = sorted(curves, key=lambda c: run["settings"].get(c, {}).get(setting, 0.0))
    if not curves:
        _no_data(fig, "The run kept no residuals of this model.")
        return fig
    sets = []
    for c in curves:
        x, _, f, z = _set_on(res[c])
        sets.append((c, x / width_of(x, f), z))
    span = max(abs(u).max() for _, u, _ in sets)
    reach = _reach([(u, z) for _, u, z in sets])
    bounds = {"the core": [(-CORE, CORE)], "the flanks": [(-WINGS, -CORE), (CORE, WINGS)], "the wings": [(-span, -WINGS), (WINGS, span)]}
    grid = plotted_sweeps(run, curves) if len(sets) > PLACE_ROWS else []
    if len(grid) > 3 and len({along for along, _, _ in grid}) == 1 and sum(len(m) for _, _, m in grid) == len(sets):
        return _place_grid(run, k, fig, dict((c, (u, z)) for c, u, z in sets), grid, places, excess, line, seen)
    n_pages = max(math.ceil(len(sets) / PLACE_ROWS), 1)          # never a division by nought (review 265)
    size = math.ceil(len(sets) / n_pages)          # pages of equal size, at most six curves: twenty five on one slide were rows too thin to read, the public CO2 lines of his 604, and pages of six left a last page of one
    pages = [sets[i:i + size] for i in range(0, len(sets), size)]
    figs = []
    for n, page in enumerate(pages):
        fig = fig if n == 0 else plt.figure(figsize=SIZE)
        if n:
            fig.text(0.035, 0.95, title, fontsize=BIG, fontweight="bold", va="top")
        top, foot = 0.83, 0.23          # room under the name of the abscissa for the family on two lines (the slides read against their words, 0.1.3rc4)
        each = (top - foot) / size          # every slide of a run at one height of row
        for i, (c, u, z) in enumerate(page):
            ax = fig.add_axes([0.2, top - (i + 1) * each, 0.73, each * 0.92])          # the scale of the top strip on its right, inside the slide
            for place in places:
                for a, b in bounds.get(place, []):
                    ax.axvspan(a, b, color=ACCENT, alpha=0.16, lw=0)
            middle, mean, count = _cells(u, z)
            bars = mean * np.sqrt(count)
            ax.axhspan(-1.96, 1.96, color=GREY, alpha=0.22, lw=0)
            ax.bar(middle, bars, width=2 * span / 60, color=[BAD if abs(v) > 1.96 else GREY for v in bars], lw=0)
            ax.axhline(0.0, color=GREY, lw=0.8)
            ax.set_xlim(-span, span)
            ax.set_yscale("symlog", linthresh=1.96, linscale=0.6)          # the band of a right model as wide as a decade beyond it: on a linear axis a misfit of 50 noises left the band a line
            ax.set_ylim(-reach, reach)
            ax.set_yticks(place_ticks(reach, bool(each * fig.get_figheight() * 72 >= 110)) if i == 0 else [])          # the band's edge named: between 0 and 10 it read as 5 (the slides read against their words, 0.1.3rc4)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: "" if v == 0 else f"{v:g}"))
            ax.yaxis.set_minor_locator(plt.NullLocator())
            ax.yaxis.tick_right()
            ax.tick_params(axis="y", labelsize=TICK)
            for side in ("top", "left") + (("right",) if i else ()):
                ax.spines[side].set_visible(False)
            ax.tick_params(labelbottom=i == len(page) - 1, labelsize=TICK)
            ax.text(-0.01, 0.5, label_of(run, c, 16), transform=ax.transAxes, ha="right", va="center", fontsize=TICK)
        first = n * size + 1
        ax.set_xlabel("distance from the centre, in widths of the curve" + (f", curves {first} to {first + len(page) - 1} of {len(sets)}" if len(pages) > 1 else ""), fontsize=LABEL)          # which curves, under the plot: beside a long title it touched
        if places and excess is not None and line:
            said = f"excess {excess:.0f} over the twin, a right model under {line:.1f}"          # short of the unit of the scale on its right: on the public CO2 spectra the longer line ran into it
            fig.text(0.035, 0.875, said, fontsize=LABEL, color=ACCENT, va="top")
        note = fig.text(0.965, 0.03, f"grey: ±1.96 noises, orange past, model {k}" if correlated(run, k) > 1 else f"grey: a right model, orange past it, model {k}", fontsize=TICK, color=GREY, ha="right", va="bottom")          # the band of independent noise, where the slide of the fit widens its own, said by its size and not as a right model's (the ninth reading of 0.1.3rc4)
        renderer = fig.canvas.get_renderer()
        limit = note.get_window_extent(renderer=renderer).x0 / fig.bbox.width - 0.02
        family = f"all answers seen most by {seen}" if seen and len(f"all answers seen most by {seen}") <= 64 else f"all answers seen most by {seen.split(': ', 1)[0]}" if seen else ""          # the family of all the answers, not of this place's (the slides read against their words, 0.1.3rc4)          # the family alone, without what follows its colon, where it is long
        for room in range(len(family), 12, -2):          # as much of it as the room beside the note holds: cut at 38 characters, two families both read "the chi functions in..." (the tenth reading of 0.1.3rc4)
            lines = textwrap.wrap(family, room)
            said_family = fig.text(0.035, 0.03, "\n".join(lines[:1] + [textwrap.shorten(" ".join(lines[1:]), room, placeholder="...")] * (len(lines) > 1)), fontsize=TICK, va="bottom", linespacing=1.1)          # two lines: on one, two families both read "the chi functions..." (the slides read against their words, 0.1.3rc4)
            if said_family.get_window_extent(renderer=renderer).x1 / fig.bbox.width <= limit:
                break
            said_family.remove()          # the band of independent noise, where the slide of the fit widens its own (the ninth reading of 0.1.3rc4)
        fig.text(0.965, 0.875, "residual in noises, log past 2", fontsize=TICK, color=GREY, ha="right", va="top")          # its strips are symlog at every reach          # the unit of the scale on the right of the top strip
        figs.append(fig)
    return figs[0] if len(figs) == 1 else figs


PANELS = 6          # the plots of one slide of sweeps, the levels of a grid of two settings


def panel_title(run: dict, along, fixed: dict, low: float, high: float) -> str:
    """The title of one plot of a sweep: the settings that it holds fixed where it has them, one level of a grid, as
    "245 K", or the setting that it changes and its range. Three titles of a range, "along pressure, 49.7 to 750 Torr",
    ran into each other on the public CO2 lines (his message 604)."""
    units = run.get("units", {})
    if fixed:
        return ", ".join(f"{v:.3g} {units.get(key, '')}".strip() for key, v in list(fixed.items())[:2])
    return f"along {along}, {low:.3g} to {high:.3g} {units.get(along, '')}".strip() if along else "the curves"


def _place_grid(run: dict, k: int, fig, sets: dict, grid: list, places: list, excess, line, seen):
    """Phase 6 for a grid of two settings, his message 605: one colour map a level of the held setting, one row a curve
    along the swept one, one column a cell of the distance from the centre in widths of the curve, coloured by what the fit
    leaves there over the noise, red above the fit and blue below, white at nought and pale inside the band of a right model, 1.96 either side, the colour growing on a scale of logarithms beyond it. Twenty five
    strips on five slides were still not readable at the size that a reader sees them."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from .where import CORE, WINGS
    reach = min(max(np.percentile(np.abs(u), 95) for u, _ in sets.values()), 3.0 * WINGS)          # the core, the flanks and the near wings: on fifteen widths a side the curves of high pressure were thin stripes, his message 605
    edges = np.linspace(-reach, reach, 41)
    maps, top = [], 2.0
    for along, fixed, members in grid:
        rows = []
        for c in members:
            u, z = sets[c]
            idx = np.digitize(u, edges) - 1
            row = np.full(len(edges) - 1, np.nan)
            for j in range(len(edges) - 1):
                held = z[idx == j]
                if held.size:
                    row[j] = held.sum() / np.sqrt(held.size)          # one noise of scatter a cell on a right model
            rows.append(row)
        maps.append(np.array(rows))
        top = max(top, float(np.nanmax(np.abs(maps[-1]))))
    norm = SymLogNorm(linthresh=1.96, linscale=0.5, vmin=-top, vmax=top)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#eeeeee")
    n = len(grid)
    w = 0.72 / n
    units = run.get("units", {})
    for i, ((along, fixed, members), m) in enumerate(zip(grid, maps)):
        ax = fig.add_axes([0.13 + i * w, 0.17, w - 0.03, 0.6])          # apart, so that the numbers under one map do not meet the next
        image = ax.imshow(m, aspect="auto", cmap=cmap, norm=norm, extent=(-reach, reach, len(members) - 0.5, -0.5), interpolation="nearest")
        for edge in (-WINGS, -CORE, CORE, WINGS):
            if abs(edge) < reach:
                ax.axvline(edge, color="black", lw=1.2, ls=":")
        ax.set_title(panel_title(run, along, fixed, 0.0, 0.0), fontsize=LABEL)
        ax.set_yticks(range(len(members)))
        if i == 0:
            ax.set_yticklabels([f"{run['settings'][c][along]:.3g}" for c in members], fontsize=TICK)
            ax.set_ylabel(f"{along}, {units.get(along, '')}".strip(", "), fontsize=LABEL)
        else:
            ax.set_yticklabels([])
        ax.set_xticks([-round(reach / 2), 0, round(reach / 2)])          # three numbers a map
        ax.tick_params(labelsize=TICK)
    bar = fig.add_axes([0.875, 0.17, 0.015, 0.6])
    marks = [v for v in (-100, -10, -2, 0, 2, 10, 100) if abs(v) <= top]
    fig.colorbar(image, cax=bar, ticks=marks)
    bar.set_yticklabels([f"{v:g}" for v in marks])          # plain numbers: the powers of ten of a symmetric log scale were a crowd
    bar.tick_params(labelsize=TICK)
    bar.set_ylabel("residual in noises, log past 2", fontsize=TICK)          # the colours of SymLogNorm
    fig.text(0.13 + 0.36, 0.03, f"distance from the centre, in widths of the curve: dotted, the core and the flanks\ndrawn to {reach:.0f} widths a side, light grey: no sample", fontsize=LABEL, ha="center", va="bottom")          # the reach of the maps and their light grey said (the slides read against their words, 0.1.3rc4)
    if places and excess is not None and line:
        fig.text(0.035, 0.875, f"excess {excess:.0f} over the twin, a right model under {line:.1f}", fontsize=LABEL, color=ACCENT, va="top")
    return fig


def plotted_sweeps(run: dict, curves: list) -> list:
    """The sweeps of the plotted slides of row M600, three at most, or up to `PANELS`, six, where a grid of two settings
    is drawn a level a plot, the curves of further levels, or of no sweep beside three, left out of those plots, and
    drawn on `slide_row`, as (the setting, the settings fixed, the curves in
    its order): those of `sweeps_of`, and the curves of no sweep, whose settings change together, as a density with its
    temperature, along the first setting that changes among them, with the curve of another sweep that shares their
    fixed settings, the corner of an L of settings (his message 118)."""
    every = [s for s in sweeps_of(run, curves) if len(s[2]) > 1]
    by_along = {}
    for s in every:
        by_along.setdefault(s[0], []).append(s)
    grid = next((v for v in by_along.values() if len(v) >= 3 and len({c for s in v for c in s[2]}) == len(curves)), None)
    if grid:          # a grid of two settings: every sweep along one of them, each plot one level of the other, each curve once (the public CO2 lines, 5 temperatures by 5 pressures, his message 604)
        return sorted(grid, key=lambda s: tuple(s[1].values()))[:PANELS]
    sweeps = every[:3]
    placed = {c for s in sweeps for c in s[2]}
    rest = [c for c in curves if c not in placed]
    if not rest:
        return sweeps or [(None, {}, list(curves))]
    numeric = {c: {key: v for key, v in run["settings"].get(c, {}).items() if isinstance(v, (int, float))} for c in curves}
    keys = list(numeric[rest[0]])
    along = next((key for key in keys if len({numeric[c].get(key) for c in rest}) > 1), None)
    fixed = {key: numeric[rest[0]][key] for key in keys if len({numeric[c].get(key) for c in rest}) == 1}
    corner = [c for c in placed if fixed and all(numeric[c].get(key) == v for key, v in fixed.items())]
    members = sorted(rest + corner, key=lambda c: numeric[c].get(along, 0.0)) if along else rest
    return (sweeps + [(along, fixed, members)])[:3]

def slide_curves(run: dict, k: int):
    """Phase 1 as plots, row M600: the curves as they were recorded, one plot a sweep of `plotted_sweeps`, which draws six
    levels of a grid or three sweeps at most and leaves the curves beyond them out, `slide_row` drawing every curve, each
    curve coloured by the setting that the sweep changes, the number of curves and of samples written on the plot."""
    import matplotlib.pyplot as plt
    from matplotlib import cm
    res = run["residuals"].get(k, {})
    curves = list(res)
    fig = plt.figure(figsize=SIZE)
    fig.text(0.035, 0.95, f"{ORDER.index('intake') + 1}. The curves as recorded, one plot a sweep", fontsize=BIG, fontweight="bold", va="top")
    if not curves:
        _no_data(fig, "The run kept no curves of this model.")
        return fig
    sweeps = plotted_sweeps(run, curves)
    samples = sum(len(res[c]) for c in curves)
    drawn = len({c for _, _, members in sweeps for c in members})
    fig.text(0.035, 0.875, f"{len(curves)} curves, {samples} samples" + (f", {drawn} drawn here, every one in a row in phase 4" if drawn < len(curves) else ""), fontsize=LABEL, color=ACCENT, va="top")          # a grid of 64 curves said 64 and drew 48 (the tenth reading of 0.1.3rc4)
    names = axis_names(run)
    width = 0.9 / len(sweeps)
    for i, (along, fixed, members) in enumerate(sweeps):
        ax = fig.add_axes([0.1 + i * width, 0.16, width - 0.08, 0.62])          # at 0.08 "fluorescence, a.u." lost its first letter at the slide's left edge
        values = [run["settings"].get(c, {}).get(along) for c in members] if along else list(range(len(members)))
        low, high = min(values), max(values)
        for c, v in zip(members, values):
            a = res[c][np.argsort(res[c][:, 0])]
            ax.plot(a[:, 0], a[:, 1], color=cm.viridis(0.85 * (v - low) / (high - low) if high > low else 0.0), lw=1.6)
        ax.set_title(panel_title(run, along, fixed, low, high), fontsize=LABEL)
        if i == 0:
            ax.set_ylabel(names[1], fontsize=LABEL)
        ax.tick_params(labelsize=TICK)
        if len(sweeps) > 2:
            from matplotlib.ticker import MaxNLocator
            lo, hi = ax.get_xlim()
            ends = [t for t in MaxNLocator(4).tick_values(lo, hi) if lo <= t <= hi]
            ax.set_xticks([ends[0], ends[-1]] if len(ends) > 1 else ends)          # two ticks: three, "-0.4 0.0 0.4", still touched on five plots
            ax.tick_params(axis="x", labelsize=15)          # at the 15 points of his rule, the least: at 20 three ticks did not fit a plot of a fifth of the slide
    fig.text(0.5, 0.03, names[0], fontsize=LABEL, ha="center", va="bottom")          # one label under the plots: two of them ran into each other
    alongs = {along for along, _, _ in sweeps if along}
    fig.text(0.965, 0.875, f"dark to light: low to high {next(iter(alongs)) if len(alongs) == 1 else 'setting'}, model {k}", fontsize=TICK, color=GREY, ha="right", va="top")
    return fig

def slide_fixed(run: dict, k: int):
    """Phase 2 as one plot, row M600: the model by its name, or by its first model and its terms where it has terms, and each of its first ten shared parameters by the bar of the joint
    fit in per cent of its value, on a scale of powers of ten, the parameters of a bar under the value read and those above it not.
    Its number is the largest correlation of two parameters, which says which the curves cannot tell apart."""
    import matplotlib.pyplot as plt
    ident = run["phases"].get(k, {}).get("identify", {})
    fig = plt.figure(figsize=SIZE)
    fig.text(0.035, 0.95, f"{ORDER.index('identify') + 1}. The model and the bar of its shared parameters", fontsize=BIG, fontweight="bold", va="top")
    params = {n: v for n, v in (ident.get("parameters") or {}).items() if isinstance(v, dict) and v.get("per_cent_of_its_value") is not None}
    if not params:
        _no_data(fig, "The run kept no parameters of this model.")
        return fig
    chain = chain_of(run, k) if run.get("equations") else []
    named = (f"the first model and {len(chain) - 1} term{'s' * (len(chain) > 2)}, the last {chain[-1].get('added', '')}" if len(chain) > 1
             else ident.get("model", f"model {k}"))          # the terms named: the first model's long name cut at 80 characters named none of the five terms of model 5 (the slides read against their words, 0.1.3rc4)
    fig.text(0.035, 0.875, textwrap.shorten(named, 80, placeholder="..."), fontsize=LABEL, color=ACCENT, va="top")
    names = list(params)[:10]
    shares = [max(abs(float(params[n]["per_cent_of_its_value"])), 1e-6) for n in names]
    ax = fig.add_axes([0.3, 0.2, 0.62, 0.56])
    ax.axvspan(1e-6, 100.0, color=GREY, alpha=0.15, lw=0)
    ax.scatter(shares, range(len(names)), s=180, color=[GOOD if v < 100.0 else BAD for v in shares], zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(min(min(shares) / 3.0, 0.01), max(max(shares) * 3.0, 300.0))
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(name_labels(names, 22), fontsize=name_size(name_labels(names, 22)))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))          # 0.01 to 100, not powers of ten
    ax.set_xlabel("its bar, in per cent of its value", fontsize=LABEL)
    ax.set_title("in the grey, under 100 per cent, green, and past it orange", fontsize=TICK, loc="left")          # the shade and the colours said above the plot, where under it they ran into the correlation (the slides read against their words, 0.1.3rc4)          # the shade and the colours said (the slides read against their words, 0.1.3rc4)
    ax.tick_params(axis="x", labelsize=TICK)          # the x axis alone: on both, it set the names back to 20 points
    ax.invert_yaxis()
    pairs = {pair: float(r) for pair, r in (ident.get("correlations") or {}).items() if isinstance(r, (int, float))}
    if pairs:
        pair, r = max(pairs.items(), key=lambda kv: abs(kv[1]))
        fig.text(0.035, 0.03, textwrap.shorten(f"largest correlation {r:.2f}, {pair}", 90, placeholder="..."), fontsize=TICK, va="bottom")
    return fig

def slide_admitted(run: dict, k: int):
    """Phase 3 as plots, row M600: synthetic sets drawn from the fitted model and fitted again, and for each of its first
    eight parameters its bias over its spread, its spread over its bar and the share of sets whose interval holds the
    value that drew them, each against the white band of the rule. Its number is how many of those stand inside every
    band."""
    import matplotlib.pyplot as plt
    adm = run["phases"].get(k, {}).get("twin", {}).get("admission", {})
    fig = plt.figure(figsize=SIZE)
    fig.text(0.035, 0.95, f"{ORDER.index('twin') + 1}. The twin: does the fit return what drew it?", fontsize=BIG, fontweight="bold", va="top")
    if not adm.get("parameters"):
        _no_data(fig, "The run kept no twin of this model.")
        return fig
    bands = adm.get("bands", {})
    names = list(adm["parameters"])[:8]
    under = float(bands.get("bias_over_spread_under", 0.5))
    what = [("bias over spread", "bias_over_spread", (-under, under)),
            ("spread over bar", "spread_over_bound", tuple(bands.get("spread_over_bound", (0.6, 1.4)))),
            ("share of sets held", "coverage", tuple(bands.get("coverage", (0.83, 1.0))))]
    inside = {n: True for n in names}
    for i, (label, key, (lo, hi)) in enumerate(what):
        ax = fig.add_axes([0.25 + i * 0.255, 0.16, 0.2, 0.6])          # the names on two lines inside the slide: at 0.2 the last lost its first letters
        ys = np.arange(len(names))
        vals = [float((adm["parameters"][n] or {}).get(key, np.nan) or np.nan) for n in names]
        for n, v in zip(names, vals):
            inside[n] = inside[n] and lo <= v <= hi
        ax.axvspan(lo, hi, color="white", zorder=0)
        ax.set_facecolor("#f2e3dc")
        ax.scatter(vals, ys, s=160, color=[GOOD if lo <= v <= hi else BAD for v in vals], zorder=3)
        ax.set_yticks(ys)
        ax.set_yticklabels(name_labels(names, 14) if i == 0 else [""] * len(names), fontsize=name_size(name_labels(names, 14)))
        ax.set_title(label, fontsize=LABEL)
        span = max(hi - lo, 1e-9)
        finite = [v for v in vals if math.isfinite(v)] or [lo]
        ax.set_xlim(min(lo - 0.3 * span, min(finite) - 0.1 * span), max(hi + 0.3 * span, max(finite) + 0.1 * span))
        ax.set_xticks([lo, hi])          # the edges of the band alone: three ticks of four digits ran together at 20 points
        ax.set_xticklabels([f"{lo:.2g}", f"{hi:.2g}"], fontsize=TICK)
        ax.invert_yaxis()
    held = sum(inside.values())
    fig.text(0.035, 0.875, f"{held} of {len(names)} parameters inside every white band, on {adm.get('sets', '?')} synthetic sets", fontsize=LABEL, color=GOOD if held == len(names) else BAD, va="top")
    fig.text(0.965, 0.015, "green: inside the white band, orange: outside it", fontsize=TICK, color=GREY, ha="right", va="bottom")          # the colours said (the fifteenth reading of 0.1.3rc4)
    return fig

def slide_fitted(run: dict, k: int):
    """Phase 4 as plots, row M600: the curves of each sweep of `plotted_sweeps`, six levels of a grid or three sweeps at
    most, with the joint fit drawn through them, and above them what the fit leaves on each over its noise, summed in cells and divided by the root of their count, inside the grey band
    that a right model of independent noise fills on 19 cells of 20. Its number is the chi square over its degrees of
    freedom, beside the band of a right model, 1.96 times the root of 2 over the degrees, times the root of the samples
    over which neighbouring residuals correlate where they do, `correlated`: 1.96 of the standard deviations that
    `chi_said` gives."""
    import matplotlib.pyplot as plt
    from matplotlib import cm
    res = run["residuals"].get(k, {})
    curves = list(res)
    table = run["tables"].get(k, {})
    fig = plt.figure(figsize=SIZE)
    fig.text(0.035, 0.95, f"{ORDER.index('fit') + 1}. The joint fit and what it leaves", fontsize=BIG, fontweight="bold", va="top")
    if not curves:
        _no_data(fig, "The run kept no residuals of this model.")
        return fig
    chi, dof = table.get("chi_square"), table.get("degrees_of_freedom")
    if isinstance(chi, (int, float)) and isinstance(dof, (int, float)) and dof > 0:
        span = correlated(run, k)
        fig.text(0.035, 0.875, (f"chi square a degree {chi / dof:.2f}, a right model 1 within {1.96 * math.sqrt(2.0 / dof * span):.2f}, {span:.0f} samples correlated" if span > 1
                                else f"chi square over its degrees {chi / dof:.2f}, a right model 1 within {1.96 * math.sqrt(2.0 / dof):.2f}"), fontsize=LABEL, color=ACCENT, va="top")          # the line of correlated noise ran past the slide's edge at its first form
    sweeps = plotted_sweeps(run, curves)
    alongs = {along for along, _, _ in sweeps if along}
    drawn = len({c for _, _, members in sweeps for c in members})
    fig.text(0.965, 0.975, f"dark to light:\nlow to high {next(iter(alongs)) if len(alongs) == 1 else 'setting'}" + (f"\n{drawn} of {len(curves)} curves" if drawn < len(curves) else ""), fontsize=TICK, color=GREY, ha="right", va="top")          # what the colours say, as on the slide of the curves: here nothing said it (the slides read against their words, 0.1.3rc4)
    names = axis_names(run)
    width = 0.85 / len(sweeps)
    reach = 4.5
    for along, _, members in sweeps:
        for c in members:
            x, _, _, z = _set_on(res[c])
            _, mean, count = _cells(x, z)
            reach = max(reach, 1.1 * float(np.abs(mean * np.sqrt(count)).max()))
    for i, (along, fixed, members) in enumerate(sweeps):
        values = [run["settings"].get(c, {}).get(along, j) if along else j for j, c in enumerate(members)]
        low, high = min(values), max(values)
        top = fig.add_axes([0.12 + i * width, 0.58, width - 0.05, 0.2])          # room on the left for the names of the scales, cut on the public CO2 lines
        bottom = fig.add_axes([0.12 + i * width, 0.14, width - 0.05, 0.42], sharex=top)
        top.axhspan(-1.96, 1.96, color=GREY, alpha=0.25, lw=0)
        for c, v in zip(members, values):
            colour = cm.viridis(0.85 * (v - low) / (high - low) if high > low else 0.0)
            a = res[c][np.argsort(res[c][:, 0])]
            z = a[:, 3] / np.where(a[:, 4] > 0, a[:, 4], 1.0)
            middle, mean, count = _cells(a[:, 0], z)
            top.step(middle, mean * np.sqrt(count), where="mid", color=colour, lw=1.6)
            bottom.plot(a[:, 0], a[:, 1], ".", ms=2, color=colour, alpha=0.5)
            bottom.plot(a[:, 0], a[:, 2], "-", lw=1.4, color="black")
        _seen(top, reach)
        top.axhline(0.0, color=GREY, lw=0.8)
        top.tick_params(labelbottom=False, labelleft=i == 0 or len(sweeps) <= 3, labelsize=TICK)          # one scale of noises for every strip: on five plots the numbers of each ran into the next, his message 604
        bottom.tick_params(labelsize=TICK)
        top.set_title(panel_title(run, along, fixed, low, high), fontsize=LABEL)
        if i == 0:
            top.set_ylabel("noises", fontsize=TICK)          # the unit alone: the strip is too short for more, and "noises, log past 2" ran onto the name below it
            bottom.set_ylabel(names[1], fontsize=TICK)
    fig.text(0.5, 0.03, names[0] + ", the fit in black, grey: 1.96 noises" + (", log past 2" if reach > SEEN else ""), fontsize=LABEL, ha="center", va="bottom")
    return fig

def slide_left(run: dict, k: int):
    """Phase 5 as one plot, row M600: what the fit leaves on each curve, the mean square of its residuals over their
    noise, against the setting that each sweep changes, the sweeps of `plotted_sweeps`. A right model reads 1 within the
    grey band, 1.96 standard errors of a mean square of n normal residuals, the root of 2 over n, times the root of the
    samples over which neighbouring residuals correlate where they do, `correlated`, and the largest over
    every curve is written on the plot, its curve drawn or not. Where more than
    three sweeps run along one setting, a grid of two settings, the slide is one plot of lines, a line a level, on a scale
    of powers of ten, with the line of a right model at 1 and no band, his message 605."""
    import matplotlib.pyplot as plt
    res = run["residuals"].get(k, {})
    curves = list(res)
    fig = plt.figure(figsize=SIZE)
    fig.text(0.035, 0.95, f"{ORDER.index('read') + 1}. What the fit leaves on each curve, along its setting", fontsize=BIG, fontweight="bold", va="top")
    if not curves:
        _no_data(fig, "The run kept no residuals of this model.")
        return fig
    square = {}
    for c in curves:
        a = res[c]
        z = a[:, 3] / np.where(a[:, 4] > 0, a[:, 4], 1.0)
        square[c] = (float(np.mean(z * z)), z.size)
    sweeps = plotted_sweeps(run, curves)
    span = correlated(run, k)
    top = max(max(m for m, _ in square.values()), 1.0 + 1.96 * math.sqrt(2.0 * span / min(n for _, n in square.values())))
    lines_drawn = False
    if len(sweeps) > 3 and len({along for along, _, _ in sweeps}) == 1:          # a grid of two settings as one plot, a line a level: five small plots of dots were not readable, his message 605
        from matplotlib import cm
        ax = fig.add_axes([0.12, 0.16, 0.6, 0.62])
        along = sweeps[0][0]
        for j, (_, fixed, members) in enumerate(sweeps):
            xs = [run["settings"][c][along] for c in members]
            ax.plot(xs, [square[c][0] for c in members], "o-", lw=2.4, ms=10, color=cm.viridis(0.85 * j / max(len(sweeps) - 1, 1)), label=panel_title(run, along, fixed, 0.0, 0.0))
        ax.set_yscale("log")
        ax.axhline(1.0, color="black", lw=1.2)
        ax.text(ax.get_xlim()[0], 1.15, "  a right model", fontsize=TICK, va="bottom", ha="left")          # inside the plot: at its edge it stood outside the axes (review 265)
        unit = run.get("units", {}).get(along, "")
        ax.set_xlabel(f"{along}, {unit}" if unit else along, fontsize=LABEL)
        ax.set_ylabel("mean square over the noise", fontsize=LABEL)
        ax.tick_params(labelsize=TICK)
        ax.legend(fontsize=TICK, frameon=True, loc="upper left", bbox_to_anchor=(1.02, 1.0))
        sweeps, lines_drawn = [], True
    width = 0.84 / max(len(sweeps), 1)
    for i, (along, fixed, members) in enumerate(sweeps):
        ax = fig.add_axes([0.12 + i * width, 0.16, width - 0.04, 0.6])          # room on the left for the scale's name, which ran off the slide on the public CO2 lines
        xs = [run["settings"].get(c, {}).get(along, j) if along else j for j, c in enumerate(members)]
        for x, c in zip(xs, members):
            m, n = square[c]
            band = 1.96 * math.sqrt(2.0 * span / n)
            ax.add_patch(plt.Rectangle((x - 0.02 * (max(xs) - min(xs) or 1.0), 1.0 - band), 0.04 * (max(xs) - min(xs) or 1.0), 2 * band, color=GREY, alpha=0.3, lw=0))
            ax.plot([x], [m], "o", ms=12, color=BAD if abs(m - 1.0) > band else ACCENT)
        ax.axhline(1.0, color="black", lw=1.2)
        ax.set_ylim(0.0, 1.08 * top)
        unit = run.get("units", {}).get(along, "")
        named = f"{along}, {unit}" if along and unit else (along or "curve")
        if len(sweeps) <= 3:
            ax.set_xlabel(named, fontsize=LABEL)
        elif i == 0:
            fig.text(0.12 + 0.42, 0.03, named, fontsize=LABEL, ha="center", va="bottom")          # one name under the plots: five of them ran into each other, his message 604
        if fixed:
            ax.set_title(panel_title(run, along, fixed, 0.0, 0.0), fontsize=LABEL)          # which level of the grid, his message 604
        if i == 0:
            ax.set_ylabel("mean square over the noise", fontsize=LABEL)
        else:
            ax.tick_params(labelleft=False)          # one scale for all, on the first plot
        ax.tick_params(labelsize=TICK)
        if len(sweeps) > 3:
            ax.xaxis.set_major_locator(plt.MaxNLocator(3))          # three ticks on a narrow plot
    worst = max(curves, key=lambda c: square[c][0])
    drawn = len({c for _, _, members in plotted_sweeps(run, curves) for c in members})          # from the sweeps that the slide draws, whichever form: counted over `sweeps`, which the grid's plot of lines empties, it read 0 of 25 (the sixteenth reading of 0.1.3rc4)
    if drawn < len(curves):          # the curves left out said, as slides 1 and 4a say them (the fifteenth reading of 0.1.3rc4)
        fig.text(0.965, 0.015, f"{drawn} of {len(curves)} curves drawn, every one in a row in phase 4", fontsize=TICK, color=GREY, ha="right", va="bottom")
    fig.text(0.035, 0.875, f"largest {square[worst][0]:.2f} on {label_of(run, worst, 24)}, a right model " + ("at 1" if lines_drawn else "1 within the grey band, orange past it"), fontsize=LABEL, color=ACCENT, va="top")          # the plot of lines draws no grey band (review 265)
    return fig

def slide_tried(run: dict, k: int):
    """Phase 7 as plots, row M600: the next term tried on synthetic sets that hold it. Left, the normal law of the gain in
    chi square that the term brings on the sets, drawn from their mean and scatter, its band of 1.96 scatters about the
    mean, and the gain on the records as a line. Right, each of the first six shared parameters that the sets return with a scatter, or the first five and the last, the
    term's own, where there are more, against the value that drew them, in their scatter, orange past two. Its number is
    the gain on the records beside that of the sets."""
    import matplotlib.pyplot as plt
    tried = run["phases"].get(k, {}).get("candidate on the twin", {})
    fig = plt.figure(figsize=SIZE)
    fig.text(0.035, 0.95, f"{ORDER.index('candidate on the twin') + 1}. The next term, tried on its own twin", fontsize=BIG, fontweight="bold", va="top")
    gain = tried.get("gain") or {}
    if not isinstance(gain.get("mean"), (int, float)) or not isinstance(gain.get("on_the_records"), (int, float)):
        this = next((r for r in (run.get("best") or {}).get("iterations", []) if r.get("model") == f"model {k}"), None) or {}
        gone = this.get("refused") or []
        _no_data(fig, f"The term tried after this model, {gone[0]['term']}, was refused before it could be planted on its twin: {_reason(gone[0]['why'], 140)}" if gone
                 else "No next term was tried on the twin of this model.")          # a refused term said as none tried, beside a decision that named it refused (the tenth reading of 0.1.3rc4, as the blind red team of 0.1.0rc13 found it on slide_next)
        return fig
    fig.text(0.035, 0.875, two_lines(f"gain {gain['on_the_records']:.0f} on the records, {gain['mean']:.0f} on {tried.get('sets', '?')} sets: {tried.get('the_planted_candidate', 'the next term')}", 75), fontsize=LABEL, color=ACCENT, va="top")          # two lines with the sign of a cut: cut at 80 the term was never named whole
    ax = fig.add_axes([0.08, 0.16, 0.34, 0.6])          # its right edge clear of the names of the plot beside it
    mean, scatter = float(gain["mean"]), float(gain.get("scatter") or 0.0)
    u = np.linspace(mean - 5 * scatter, mean + 5 * scatter, 400) if scatter > 0 else np.array([mean])
    if scatter > 0:
        ax.plot(u, np.exp(-0.5 * ((u - mean) / scatter) ** 2), color=GREY, lw=2)
        ax.axvspan(mean - 1.96 * scatter, mean + 1.96 * scatter, color=GREY, alpha=0.2, lw=0)
    ax.axvline(float(gain["on_the_records"]), color=ACCENT, lw=3)
    ax.set_yticks([])
    ax.set_xlabel("gain in chi square, grey: the sets' law\nand its 1.96 scatters, blue: the records", fontsize=TICK)          # lines that fit the plot: "of the sets'" lost its o at the slide's edge          # a law drawn from two numbers, not the sets' own gains (the twelfth reading of 0.1.3rc4)
    ax.tick_params(labelsize=TICK)
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(3))          # three ticks: "445000447500450000452500455000" ran together (the slides read against their words, 0.1.3rc4)
    back = tried.get("what_the_whole_procedure_returns") or {}
    every = [n for n, v in back.items() if isinstance(v, dict) and v.get("scatter")]
    names = every if len(every) <= 6 else every[:5] + every[-1:]          # the first five and the last, the term's own: the first six left out the parameter of the term tried (the slides read against their words, 0.1.3rc4)
    side = fig.add_axes([0.62, 0.16, 0.34, 0.6])
    if names:
        moves = [(back[n]["mean"] - back[n]["put_in"]) / back[n]["scatter"] for n in names]
        side.axvspan(-1.0, 1.0, color=GREY, alpha=0.25, lw=0)
        side.scatter(moves, range(len(names)), s=160, color=[BAD if abs(m) > 2 else GOOD for m in moves], zorder=3)
        side.set_yticks(range(len(names)))
        side.set_yticklabels(name_labels(names, 12), fontsize=name_size(name_labels(names, 12)))
        side.set_xlabel("returned less put in, scatters", fontsize=TICK)
        if len(every) > 6:
            side.set_title(f"the first five of {len(every)} and the last", fontsize=TICK, loc="left")          # above the plot: under it, it ran into the key of the colours
        fig.text(0.965, 0.015, "grey: one scatter, orange: past two", fontsize=TICK, color=GREY, ha="right", va="bottom")          # the colours said
        side.set_xlim(-max(3.0, 1.2 * max(abs(m) for m in moves)), max(3.0, 1.2 * max(abs(m) for m in moves)))
        side.tick_params(axis="x", labelsize=TICK)          # the x axis alone: on both, it set the names back to 20 points
        side.invert_yaxis()
    else:
        _no_data_in(side, "no parameter returned")
    return fig

def pooled_strip(res: dict, n: int = 60) -> tuple:
    """What a fit leaves on all its curves together, whatever its sign on each: each curve's residuals over their noise
    placed in widths of the curve from its centre, and in each of `n` cells of one width the mean of their squares less
    one, over its standard error, the root of 2 over the samples of the cell, so that a right model of independent noise
    stands in the band of 1.96 on 19 cells of 20. Summed with their signs, a width that changes along the temperature cancelled between the
    hot and the cold curves, and the slide of the decision read the first model of the tutorial inside the band (row
    M600). Returns the middles of the cells, their values and the half span."""
    places = []
    for a in res.values():
        x, _, f, z = _set_on(a)
        places.append((x / width_of(x, f), z))
    u = np.concatenate([q for q, _ in places])
    z = np.concatenate([w for _, w in places])
    span = float(np.abs(u).max())
    edges = np.linspace(-span, span, n + 1)
    k = np.clip(np.searchsorted(edges, u, side="right") - 1, 0, n - 1)
    held = np.unique(k)
    values = np.array([(np.mean(z[k == i] ** 2) - 1.0) / math.sqrt(2.0 / (k == i).sum()) for i in held])
    return 0.5 * (edges[held] + edges[held + 1]), values, span


def slide_chosen(run: dict, k: int):
    """Phase 8 as plots, row M600: the term that the next model adds, what the fit leaves on all the curves before it in
    grey and after it in colour, inside the band of a right model of independent noise, and the move of each of the first six shared parameters that both share with a bar, from one model to the next,
    in its own bar. Its number is the chi square over its degrees before and after."""
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=SIZE)
    fig.text(0.035, 0.95, f"{ORDER.index('report') + 1}. The decision: what the next term changes", fontsize=BIG, fontweight="bold", va="top")
    before, after = run["residuals"].get(k, {}), run["residuals"].get(k + 1, {})
    if not before or not after:
        _no_data(fig, "No model was made after this one: it is the last of the run.")
        return fig
    tables, phases = run.get("tables") or {}, run.get("phases") or {}
    name0 = (phases.get(k, {}).get("identify") or {}).get("model") or tables.get(k, {}).get("model", "")
    name1 = (phases.get(k + 1, {}).get("identify") or {}).get("model") or tables.get(k + 1, {}).get("model", "")
    added = name1[len(name0):].strip(" ,") if name0 and name1.startswith(name0) else name1
    added = added[5:] if added.startswith("with ") else added
    said = []
    for j in (k, k + 1):
        t = tables.get(j, {})
        if isinstance(t.get("chi_square"), (int, float)) and t.get("degrees_of_freedom"):
            said.append(f"{t['chi_square'] / t['degrees_of_freedom']:.2f}")
    fig.text(0.035, 0.875, two_lines(f"the next model adds {added}", 72), fontsize=LABEL, color=ACCENT, va="top", linespacing=1.1)          # the term on two lines with the sign of a cut
    ax = fig.add_axes([0.15, 0.22, 0.4, 0.5])          # its right edge clear of the names of the plot beside it, its own name inside the slide
    if len(said) == 2:
        ax.set_title(f"chi square over its degrees {said[0]} then {said[1]}", fontsize=TICK, loc="left")          # above the plot: inside it the text lay on the curves
    ax.axhspan(-1.96, 1.96, color=GREY, alpha=0.25, lw=0)
    reach, low = 4.5, -4.5
    for res, colour, label in ((before, GREY, "before"), (after, ACCENT, "after")):
        middle, values, span = pooled_strip(res)
        ax.step(middle, values, where="mid", color=colour, lw=2.2, label=label)
        reach, low = max(reach, 1.1 * float(values.max())), min(low, 1.1 * float(values.min()))
    ax.axhline(0.0, color="black", lw=0.8)
    _seen(ax, reach, low)          # from the lowest cell to the highest: a mean square falls below 1 by little, and a symmetric axis left its lower half empty
    ax.set_xlabel("distance from the centre, in widths", fontsize=LABEL)
    ax.set_ylabel(_logged("mean square less 1,\nin errors", max(reach, -low)), fontsize=TICK)          # short enough for its axis with the scale's words
    ax.tick_params(labelsize=TICK)
    ax.legend(fontsize=TICK, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)          # under the plot: inside it, it lay on the curves wherever it stood
    p0 = (phases.get(k, {}).get("identify") or {}).get("parameters") or {}
    p1 = (phases.get(k + 1, {}).get("identify") or {}).get("parameters") or {}
    shared = [n for n in p0 if n in p1 and p0[n].get("standard_deviation")][:6]
    side = fig.add_axes([0.76, 0.22, 0.2, 0.5])          # level with the plot beside it
    if shared:
        moves = [(p1[n]["value"] - p0[n]["value"]) / p0[n]["standard_deviation"] for n in shared]
        side.axvspan(-1.0, 1.0, color=GREY, alpha=0.25, lw=0)
        side.scatter(moves, range(len(shared)), s=160, color=[BAD if abs(m) > 2 else ACCENT for m in moves], zorder=3, clip_on=False)          # whole at the axis's edge
        side.set_yticks(range(len(shared)))
        side.set_yticklabels(name_labels(shared, 12), fontsize=name_size(name_labels(shared, 12)))
        if max(abs(m) for m in moves) > SEEN:
            _seen(side, 1.15 * max(abs(m) for m in moves), axis="x", linthresh=1.0)          # the band of one bar seen beside moves of 200 bars, four ticks, three on a short axis
        side.set_xlabel("move in bars" + (",\nlog past 1" if max(abs(m) for m in moves) > SEEN else ""), fontsize=TICK)
        fig.text(0.965, 0.015, "grey: 1.96 errors left, one bar right, blue: under two bars, orange: past two", fontsize=TICK, color=GREY, ha="right", va="bottom")          # both bands named: the left one read as one bar (the fifteenth reading of 0.1.3rc4)          # the colours said (the slides read against their words, 0.1.3rc4)          # inside the slide: at 22 points "log past 1" lost its 1
        side.tick_params(axis="x", labelsize=TICK)          # the x axis alone: on both, it set the names back to 20 points
        side.invert_yaxis()
    else:
        _no_data_in(side, "no shared parameter")
    return fig

def experiment_line(p: dict) -> str:
    """A proposal for the experiment as its slide says it: the change, and what the twin measured of it, the error of the
    known values with the change over that under the records, where it was measured (message 142 of the repository it is
    made in: each measured on the twin)."""
    said = f"{p.get('lever', '')}: {p.get('words', '')}"
    ratio = p.get("error_over_that_under_the_records")
    names = ", ".join(p.get("of_each") or {})          # "the known values" were named nowhere (the blind red team of 0.1.0rc16, mislead 2)
    return said + (f". On the twin, the error of the shared parameters{f' ({names})' if names else ''} times {ratio:.2f}" if isinstance(ratio, (int, float)) else "")


def slide_next(run: dict, k: int):
    """The slide of phase 7 before row M600, which the log no longer draws, its phase 7 being `slide_tried`, `SLIDES`, and
    which stays for a reader who draws it alone and for its tests: the term tried after this model, planted at its fitted
    size on synthetic sets and found again by the whole procedure: the bias of each of the first ten values in scatters, the share of sets in which one bar holds it, and the fall of
    chi square that the term gives on the records beside the fall on the sets."""
    cand = run["phases"].get(k, {}).get("candidate on the twin")
    fig = _frame(TITLES["candidate on the twin"], k, "candidate on the twin", run,
                 f"The term {cand['the_planted_candidate']}, planted at its fitted size on {cand.get('read', '?')} synthetic sets, and found again." if cand else "")
    if not cand:
        this = next((r for r in run.get("best", {}).get("iterations", []) if r.get("model") == f"model {k}"), None) or {}
        gone = this.get("refused") or []
        _no_data(fig, (f"The term tried after this model, {gone[0]['term']}, was refused before it could be planted on the twin: {_reason(gone[0]['why'], 140)} " if gone
                       else "No term was tried after this model: the loop ended here. ")
                 + "To go on, write the next file of terms/ from a proposal of phase 6 and run the loop again.")          # a term refused was said to be no term tried, beside a slide that said it was refused (the blind red team of 0.1.0rc13)
        return fig
    returns = cand.get("what_the_whole_procedure_returns", {})
    names = list(returns)[:10]
    ax = fig.add_axes([0.2, 0.18, 0.32, TOP - 0.2])
    ys = np.arange(len(names))
    bias = [returns[n].get("bias_in_scatters", np.nan) for n in names]
    ax.axvspan(-1, 1, color="#e8f1f8")
    ax.scatter(bias, ys, s=80, color=ACCENT, zorder=3)
    ax.axvline(0, color=GREY, lw=0.8)
    ax.set_yticks(ys)
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.set_xlabel("bias, in scatters of the sets")
    for y, n in zip(ys, names):
        share = returns[n].get("inside_one_bar")
        if share is not None:
            ax.text(1.02, y, f"{share:.0%} of sets\nwithin one bar", transform=ax.get_yaxis_transform(), va="center", fontsize=SMALL)          # "63% held" said nothing of what held (the cold round 22)
    gain = cand.get("gain", {})
    ax2 = fig.add_axes([0.66, 0.18, 0.3, TOP - 0.2])
    if gain.get("scatter"):
        u = np.linspace(gain["mean"] - 4 * gain["scatter"], gain["mean"] + 4 * gain["scatter"], 300)
        ax2.fill_between(u, np.exp(-0.5 * ((u - gain["mean"]) / gain["scatter"]) ** 2), color=GREY, alpha=0.4, lw=0, label="on the sets")
        ax2.axvline(gain["on_the_records"], color=BAD, lw=2.5, label="on the records")
        ax2.set_yticks([])
        ax2.set_xlabel("fall of chi square with the term")
        ax2.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=True)          # above the axes: it lay on the curve and the line (the cold round 22)
    else:
        _no_data_in(ax2, "no gain kept")
    return fig


def kept_said(this: dict | None, kept: dict) -> str:
    """What the decision slide says of a term kept after model `this`: its fall of the criterion, or its rise where its
    file's keep = True kept it, and no number where best.json holds either criterion not: 0 read in its place printed
    1697.9 for a rise of 7.2 (the second reading of 0.1.3rc2)."""
    if this is None or "criterion" not in this or "criterion" not in kept:
        return f"Kept: {kept['added']}."
    fall = this["criterion"] - kept["criterion"]
    return f"Kept: {kept['added']}. Its criterion is {fall:.1f} lower." if fall >= 0 else f"Kept: {kept['added']}, by its file's keep = True. Its criterion is {-fall:.1f} higher."


def decision_rows(run: dict, k: int) -> tuple:
    """The rows of the decision slide of model `k`: its own, the one of the model made after it, and every model that
    ended, each with its criterion. best.json compares the models from a held or kept term on alone, so the earlier
    slides read their criteria from each model's own table and what was kept after them from equations.json: they said
    "one model so far" and "no term was tried" on a run of three that kept one (the third reading of 0.1.3rc2)."""
    rows = run["best"].get("iterations", [])
    crit = {f"model {j}": float(t["criterion"]) for j, t in (run.get("tables") or {}).items() if isinstance(t, dict) and t.get("criterion") is not None}
    this = next((r for r in rows if r["model"] == f"model {k}"), None)
    if this is None and f"model {k}" in crit:
        this = {"model": f"model {k}", "criterion": crit[f"model {k}"]}
    kept = next((r for r in rows if r.get("after") == f"model {k}" and not r.get("aside")), None)
    if kept is None:
        made = next((e for e in run.get("equations") or [] if e.get("after") == f"model {k}"), None)
        if made is not None:
            kept = {"model": made["model"], "added": made.get("added"), "after": made["after"]} | ({"criterion": crit[made["model"]]} if made["model"] in crit else {})
    ended = [r for r in rows if r.get("ended") and r.get("criterion") is not None and re.fullmatch(r"model \d+", str(r.get("model", "")))]
    if len(crit) > len(ended):          # every model of the run, where best.json holds fewer
        ended = [{"model": m, "criterion": c, "ended": True} for m, c in sorted(crit.items(), key=lambda mc: int(mc[0].split()[-1]))]
    return this, kept, ended


def slide_decision(run: dict, k: int):
    """Phase 8: the criterion of each model so far, this one marked, and what became of the term tried after it."""
    this, kept, ended = decision_rows(run, k)
    tried = (this or {}).get("tried", [])
    refused = (this or {}).get("refused", [])
    if kept:
        said, colour = kept_said(this, kept), GOOD
    elif refused:
        said, colour = f"Refused: {refused[0]['term']}. " + _reason(refused[0]["why"], 160), BAD
    elif tried:
        said, colour = f"Tried: {', '.join(tried)}.", GREY
    else:
        said, colour = "No term was tried after this model.", GREY
    fig = _frame(TITLES["report"], k, "report", run, "The Bayesian criterion of each model, lower is better. A term is kept where its fit ends inside its bounds, its twin admits it and the criterion falls by more than 2, the criterion waived by keep = True.")          # the three conditions of the tutorial, where one sentence named two and a term that lowered the criterion by 0.53 read as kept (the blind red team of 0.1.0rc13)
    fig.text(0.04, 0.24, "\n".join(textwrap.wrap(said, WIDTH - 12)[:3]), fontsize=TEXT, color=colour, fontweight="bold", va="top", linespacing=1.3)          # bold letters are wider: the line ran past the edge (the cold round 22)
    # the models that ended, from `decision_rows`, where best.json held fewer than the run made
    if len(ended) == 1:          # one model is no comparison: a plot of one point was a graph of nothing (the blind red team of 0.1.0rc12)
        fig.text(0.08, 0.5, f"One model so far, its criterion {ended[0]['criterion']:.1f}" + (": the next term is judged against it." if kept or refused or tried else "."), fontsize=TEXT, color=GREY)          # with no term tried, nothing is judged against it: the slide said both (the thirteenth reading of 0.1.3rc4)
    elif ended:
        ax = fig.add_axes([0.11, 0.32, 0.57, TOP - 0.34])          # its name of the axis on the slide: at 0.08 "criterion" stood half past the left edge
        idx = [int(r["model"].split()[-1]) for r in ended]
        ax.plot(idx, [r["criterion"] for r in ended], "o-", color=ACCENT, ms=9)
        if k in idx:
            ax.plot([k], [ended[idx.index(k)]["criterion"]], "o", ms=18, mfc="none", mec=BAD, mew=2.5, label=f"ringed: model {k}, this slide's")
            ax.legend(fontsize=SMALL, frameon=False, loc="upper right")          # the ring said: it marked one model while the sentence spoke of the next (the slides read against their words, 0.1.3rc4)
        ax.set_xticks(idx)
        ax.set_xlabel("model")
        ax.set_ylabel("criterion")
    return fig


def terms_said(terms: list) -> str:
    """The terms of the best model in the sentence of its slide, their count first and where they stand whole, since the
    slide cuts the sentence where it is long: "Its terms:" held two of six, cut, and never said how many (the eighth
    reading of the slides of 0.1.3rc4)."""
    return f"Its {len(terms)} term{'' if len(terms) == 1 else 's'}, whole in best.md: " + "; ".join(terms) + "."


def slide_best(run: dict):
    """The output: the best model, its terms, and each of its first ten shared parameters with the parts of its uncertainty that have a size side by side, of the four bars, the fit's, the bias on the twin, theory and the spread over the models, the coverage in words,
    in units of the bar of the fit."""
    best = run["best"]
    quotable = best.get("quotable") or {}
    verdict = ("Every value may be quoted with its bar of the fit." if quotable and all(quotable.values())
               else "No value is to be quoted yet: best.md says why." if quotable and not any(quotable.values())
               else "Only the values marked quoted may be quoted: best.md says why." if quotable else "")
    fig = _frame(f"The best model: {best.get('best', '?')}", -1, None, run, (verdict + " " if verdict else "") + terms_said(best.get("terms", [])))          # the verdict first, where a cut cannot take it          # the slide dropped the warning of best.md, and a reader quoted a value that it forbids (the blind red team of 0.1.0rc14 on a fresh machine)
    params = best.get("parameters", {})
    if not params:
        _no_data(fig, "The run left no best.json.")
        return fig
    names = list(params)[:10]
    ax = fig.add_axes([0.3, 0.16, 0.28, TOP - 0.18])
    ys = np.arange(len(names))
    parts = [("fit", lambda p: 1.0), ("bias on the twin, its size", lambda p: abs(p.get("bias_on_the_twin_in_bars") or 0.0)),          # its sign left out, and said so (the slides read against their words, 0.1.3rc4)
             ("theory", lambda p: (p.get("bar_of_theory", {}).get("together") or 0.0) / (p.get("bar_of_the_fit") or 1.0)),
             ("over the models", lambda p: (p.get("spread_over_the_models") or 0.0) / (p.get("bar_of_the_fit") or 1.0))]
    drawn = [(i, label, [f(params[n]) for n in names]) for i, (label, f) in enumerate(parts)]
    drawn = [(i, label, sizes) for i, label, sizes in drawn if not i or any(sizes)]          # a part of no size is no bar and no line of the legend: "theory" and "over the models" stood empty and unexplained (the blind red team of 0.1.0rc12)
    for place, (i, label, sizes) in enumerate(drawn):
        ax.barh(ys + (place - (len(drawn) - 1) / 2) * 0.2, sizes, height=0.2, color=COLOURS[i], label=label)          # centred on its name by the parts drawn: offsets for four left a pair of two above its name (the slides read against their words, 0.1.3rc4)
    ax.set_yticks(ys)
    ax.set_yticklabels(name_labels(names, 20))          # on two lines inside the slide: seven names lost their first letters at its left edge
    ax.invert_yaxis()
    ax.set_xlabel("in bars of the fit")
    key = ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=True)
    read = best.get("searches_that_ended") or best.get("sets_of_the_twin")          # the sets that the run read: "held on 87 of 100 sets" was a coverage times 100 on a run of 30 (the cold round 22)
    paragraphs = []
    for n in names:
        cov = params[n].get("coverage_on_the_twin")
        said = (f"{n}: {pm(params[n]['value'], params[n].get('bar_of_the_fit') or float('nan'))}"
                + (f", its interval holds the truth on {round(cov * read)} of {read} sets of the twin of best.md, drawn apart from those of the slides of the turns" if isinstance(cov, (int, float)) and read else "")
                + ((", quoted" if quotable[n] else ", not yet quoted") if n in quotable else "")
                + (f", not told apart from {', '.join(params[n]['not_told_apart_from'])}" if params[n].get("not_told_apart_from") else ""))
        paragraphs.append(textwrap.wrap(said, 40, subsequent_indent="    "))          # wrapped inside the slide: its last words were cut at the right edge (the cold round 22)
    under = key.get_window_extent(renderer=fig.canvas.get_renderer()).y0 / fig.bbox.height - 0.03          # down from under the legend: up from the foot, the words of two parameters reached it at 15 points (row M583)
    renderer = fig.canvas.get_renderer()
    for i, para in enumerate(paragraphs):          # a paragraph a parameter, while the slide has room: seven ran past its foot and over the phases (the slides read against their words, 0.1.3rc4)
        said = fig.text(0.62, under, "\n".join(para), fontsize=SMALL, va="top", ha="left", linespacing=1.3)
        low = under - said.get_window_extent(renderer=renderer).height / fig.bbox.height
        if low < 0.13:
            said.remove()
            fig.text(0.62, under, f"and {len(paragraphs) - i} more in best.md", fontsize=SMALL, color=GREY, va="top")
            break
        under = low - 0.012
    return fig


SLIDES = {"intake": (slide_curves,), "identify": (slide_fixed, slides_equation, slide_effects), "twin": (slide_admitted,), "fit": (slide_fitted, slide_row),
          "read": (slide_left,), "propose": (slide_place, slide_where), "candidate on the twin": (slide_tried,), "report": (slide_chosen, slide_decision)}          # each phase first as one plot, row M600, his messages 586 to 588: the equation, the effects, every curve in a row, the six proposals and the decision's reasons stay beside it
PAGE_CAPTIONS = {slide_where: "the proposals for the experiment that the slide before had no room for"}          # the pages after the first, where they hold other things
CAPTIONS = {slide_curves: "the curves as recorded, one plot a sweep, three sweeps or six levels of a grid at most", slide_fixed: "the model and the bar of its first ten shared parameters", slides_equation: "the equation, term by term",
            slide_effects: "what each term does", slide_admitted: "the twin: does the fit return what drew it", slide_fitted: "the joint fit, and what it leaves above it",
            slide_row: "every curve in a row, what the fit leaves under each", slide_left: "what the fit leaves on each curve, along its setting",
            slide_place: "where the residuals leave the twin's band", slide_where: "the answers of the report and its proposals",
            slide_tried: "the next term, tried on its own twin", slide_chosen: "the decision: what the next term changes", slide_decision: "the decision: each model's criterion and why"}          # the line of log.md for each slide


def log(folder, dpi: int = 110) -> list:
    """The visual log of the last run of the loop in `folder`: DIR/loop/log/, one picture a slide at `dpi`, DIR/loop/log.pdf,
    and DIR/loop/log.md. Returns the paths of the pictures in order."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    run = read_run(folder)
    out = run["loop"] / "log"
    out.mkdir(exist_ok=True)
    for old in out.glob("*.png"):
        old.unlink()
    made, lines = [], [f"# The visual log of {run['folder'].name}", "",
                       "Each phase of each turn of the loop first as one plot, then beside four of them the equation and the effect of each term, every curve in a row, the answers of the report and each model's criterion. `log.pdf` holds them in order.", ""]
    with plt.rc_context(_rc()), PdfPages(run["loop"] / "log.pdf") as pdf:
        for k in sorted(run["phases"]):
            lines += [f"## Model {k}", ""]
            for p, phase in enumerate(ORDER):
                for s, draw in enumerate(SLIDES[phase]):
                    drawn = draw(run, k)
                    for i, fig in enumerate(drawn if isinstance(drawn, list) else [drawn]):          # a phase may need several slides: the equation of a model of many terms (his message 311)
                        path = out / f"{k:02d}_{p + 1}{'abc'[s] if len(SLIDES[phase]) > 1 else ''}{i + 1 if i else ''}_{SHORT[phase].replace(' ', '_')}.png"
                        fig.savefig(path, dpi=dpi)
                        pdf.savefig(fig)
                        plt.close(fig)
                        made.append(path)
                        said = getattr(fig, "caption", None) or (PAGE_CAPTIONS.get(draw, CAPTIONS.get(draw, TITLES[phase])) if i else CAPTIONS.get(draw, TITLES[phase]))          # a page after the first says what it holds: the page of the proposals for the experiment was listed as the answers and six proposals (the tenth reading of 0.1.3rc4)
                        lines.append(f"- phase {p + 1}, {said}: ![]({path.relative_to(run['loop']).as_posix()})")
            lines.append("")
        if run["best"]:
            fig = slide_best(run)
            path = out / "99_best.png"
            fig.savefig(path, dpi=dpi)
            pdf.savefig(fig)
            plt.close(fig)
            made.append(path)
            lines += ["## The best model", "", f"![]({path.relative_to(run['loop']).as_posix()})", ""]
    (run["loop"] / "log.md").write_text("\n".join(lines))
    return made
