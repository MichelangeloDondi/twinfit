"""The visual log of a run of the loop: for each model that the loop reached, one slide for each of the eight phases of
its turn, with the equation of the model and the effect of each of its terms drawn on the curves, and the best model
at the end with the five parts of each bar. The owner's message 309: "show clearly also the equations used and the
effects of the terms of the models, in each turn, with at least a visual log slide for each phase".

    python -m twinfit loop DIR
    python -m twinfit log DIR

It reads only what the run left in DIR/loop: equations.json and effects_K.npz, the notebook of each phase, residuals_K.csv,
model_K.json and best.json. It fits nothing, so it takes seconds. It writes DIR/loop/log/, one picture a slide,
DIR/loop/log.pdf with every slide in order, and DIR/loop/log.md, which lists the slides with the one sentence that each
shows. A run made before 0.1.0rc8 left no equations.json: its slides say so where the equation and the effects stand.

The phases of a turn are those of the notebook, in its order:

    1 the records            what was read: the curves, how many, how many samples
    2 the model              its equation, term by term, its parameters with their bars, and what each term does to the curves
    3 the twin               whether synthetic sets of the fitted model give back what made them, so that the fit may be read
    4 the fit                each curve with the model, and what the model leaves in units of the noise
    5 what the fit leaves    the residuals against the noise: their spread, and their spread along the height of the curve
    6 where to search        the answers of the report beyond the twin, and the three proposals for the model and the experiment
    7 the next term          the term tried after this model, on synthetic sets that hold it
    8 the decision           the criterion of each model, and whether the next term was kept or refused, and why
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
TEXT, SMALL, TITLE = 17, 14, 26
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
        fig.text(0.035, 0.865, "\n".join(textwrap.wrap(sentence, WIDTH)[:2]), fontsize=TEXT, va="top", linespacing=1.3)
    for j, name in enumerate(ORDER):
        x = 0.035 + j * 0.1175
        on = phase == name
        fig.patches.append(plt.Rectangle((x, 0.018), 0.108, 0.045, transform=fig.transFigure, facecolor=ACCENT if on else "white",
                                         edgecolor=ACCENT, linewidth=1.2))
        fig.text(x + 0.054, 0.0405, f"{j + 1} {SHORT[name]}", fontsize=12, ha="center", va="center", color="white" if on else ACCENT)
    fig.text(0.965, 0.075, run["folder"].name, fontsize=11, color=GREY, ha="right")
    return fig


def _no_data(fig, words: str) -> None:
    fig.text(0.5, 0.45, "\n".join(textwrap.wrap(words, 80)), fontsize=TEXT, ha="center", va="center", color=GREY)          # wrapped by its characters: a long sentence ran past both edges (the blind red team of 0.1.0rc12)


def slide_records(run: dict, k: int):
    """Phase 1: every curve that the model reads, each in its colour, with the count of records and samples."""
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
    ax.legend(loc="upper left", frameon=True, markerscale=5, fontsize=11)          # inside the plot, where its corner is empty: beside it, the names of eight curves left the slide (0.1.0rc8)
    return fig


def parameter_rows(params: dict, room: int = 30) -> list:
    """The rows of the table of parameters of a slide: each name on lines of `room` characters at most, its value with
    its bar, and the bar as a share of the value. A name of 42 characters ran past its cell in a blind run of 0.1.0rc10
    (a time scale of the model over its theory, squared)."""
    return [["\n".join(textwrap.wrap(name, room)) or name, pm(v["value"], v["standard_deviation"]), f"{v['per_cent_of_its_value']:.2g} %"] for name, v in params.items()]


def _height_of(j: int, row: dict) -> float:
    """The height on a slide of line j of the equation: its head, its formula, and its letters where the formula is
    written in mathematics."""
    return 0.055 + (0.045 if (j or row.get("letters")) and str(row.get("formula") or "").startswith("$") else 0.0) + 0.07 + (0.025 if "\\frac" in str(row.get("formula") or "") else 0.0)          # a fraction stands taller


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


def _words_of(j: int, row: dict) -> int:
    """The words that line `j` of an equation puts on its slide: its heading, its value and the letters under it."""
    formula = row.get("formula") or ""
    letters = (letters_list(formula, row.get("parameter")) if j else [f"{k} {v}" for k, v in (row.get("letters") or {}).items()]) if formula.startswith("$") else []
    head = f"f: {'the first model, ' if j == 0 else ''}{row.get('added', '')}"
    letters = [x for x in letters if x.split(" ")[0] not in ("f", "a")]          # said once by the sentence at the top of the slide
    return _lettered(head) + (1 if j else 0) + (3 if row.get("held_at") else 0) + _lettered(" ".join(letters))


def slide_equation(run: dict, k: int, first: int = 0):
    """Phase 2, first slide: the equation of the model, the first model and then each term on the model before it, with
    the value of its parameter, and the parameters of the fit with their bars. From line `first` of the equation, as many
    lines as the slide holds: the table where the equation ends on this slide, else none, and `fig.rest` the line that
    the next slide starts from, or None."""
    ident = run["phases"].get(k, {}).get("identify", {})
    chain = chain_of(run, k)
    fig = _frame(TITLES["identify"] + (", continued" if first else ""), k, "identify", run, r"Each line turns $f$, the model above it, into the next, with $a$ its parameter." if not first else "")
    fig.rest = None
    params = ident.get("parameters", {})
    if not chain:
        _no_data(fig, "This run left no equations.json: it was made before the package wrote the equation of each model.")
    y = TOP
    table = 3 + sum(_lettered(n) for n in list(params)[:6])          # the words of the table of the parameters: its heads and the names, its numbers no words
    box = 0.045 * (min(len(params), 6) + 1) + 0.03          # the height of the table, by its rows: a fixed box of six rows sent a table of two to a slide of its own
    last = (y - sum(_height_of(j, row) for j, row in enumerate(chain) if j >= first) >= 0.14 + box
            and sum(_words_of(j, row) for j, row in enumerate(chain) if j >= first) + table <= WORDS)          # the rest of the equation ends on this slide, above the table, in its room and its words
    said_letters = set() if first else {"f", "a"}          # the sentence at the top of the first slide says f and a
    spent = 0
    for j, row in enumerate(chain):
        if j < first:
            continue
        if not last and y - _height_of(j, row) < 0.17 and j > first:          # a line that does not reach the foot goes on this slide, and the table on the next: a slide of two terms left its lower half blank (the cold round 22)
            fig.rest = j
            return fig
        if j > first and spent + _words_of(j, row) > WORDS:          # nor past the words of a slide
            fig.rest = j
            return fig
        spent += _words_of(j, row)
        head = f"$f_{{{j}}}$: " + ("the first model, " + row["added"] if j == 0 else row["added"])
        said = fig.text(0.04, y, head if len(head) < 64 else head[:61] + "...", fontsize=TEXT, fontweight="bold", va="top")
        if j and (row.get("parameter") in params or row.get("held_at")):
            v = params.get(row["parameter"]) or {"value": row["held_at"][0], "standard_deviation": row["held_at"][1]}
            end = said.get_window_extent(renderer=fig.canvas.get_renderer()).x1 / fig.bbox.width          # the value beside the words, where they end
            fig.text(end + 0.015, y, f"$a$ = {pm(v['value'], v['standard_deviation'])}" + (", held from theory" if row.get("parameter") not in params else ""), fontsize=TEXT, va="top", color=ACCENT)          # a held term's value is its theory (the cold round 22)
        y -= 0.055
        formula = row.get("formula") or ""
        if formula.startswith("$"):
            fig.text(0.07, y, formula, fontsize=TEXT + 5, va="top")          # the formula larger than its heading, as a displayed equation
            if j or row.get("letters"):          # the letters of the first model from LETTERS of its file (the blind red team of 0.1.0rc12: h, x_0, w and b were defined nowhere)
                y -= 0.045 + (0.025 if "\\frac" in formula else 0.0)          # under a fraction, which stands taller: its letters lay on its denominator
                new = [x for x in (letters_list(formula, row.get("parameter")) if j else [f"{a} {b}" for a, b in row["letters"].items()]) if x.split(" ")[0] not in said_letters]          # each letter defined once on a slide (his message 311: fewer words a slide)
                said_letters.update(x.split(" ")[0] for x in new)
                fig.text(0.07, y, "\n".join(letters_lines(new)[:2]), fontsize=12, color=GREY, va="top")          # the letters of the formula on its slide, each typeset as mathematics
        elif formula:          # a model file with no EQUATION: its Python line, and how to have it typeset (his message 314)
            fig.text(0.07, y, formula if len(formula) < 52 else formula[:49] + "...", fontsize=TEXT - 2, family="monospace", va="top")
            fig.text(0.07, y - 0.04, "no EQUATION in its model file: its code is shown", fontsize=11, color=GREY, va="top")          # read as a promise that it was typeset (the blind red team of 0.1.0rc12, stumble 7)
        else:
            fig.text(0.07, y, "its formula is not written: a term of one's own takes formula= in its Mechanism", fontsize=SMALL, color=GREY, va="top")
        y -= 0.07 + (0.025 if "\\frac" in formula and not (j or row.get("letters")) else 0.0)          # a fraction stands taller, as `_height_of` counts it, its room taken above its letters where it has them
    if not last:          # the equation ended low on a slide with no room left for the table: the table on a slide of its own
        fig.rest = len(chain) if first < len(chain) else None
        if fig.rest is not None:
            return fig
    ax = fig.add_axes([0.04, 0.12, 0.6, box])
    ax.axis("off")
    rows = parameter_rows(params)
    if rows:
        tab = ax.table(cellText=rows[:6], colLabels=["parameter", "value of this turn", "bar"], loc="upper left", cellLoc="left", colLoc="left", colWidths=[0.48, 0.34, 0.18])          # the values of a turn read as results, which best.md alone quotes (the blind red team of 0.1.0rc16, mislead 11), and then a value printed under a head that said it was not quoted (the blind red team of 0.1.0rc19, mislead 1)
        tab.auto_set_font_size(False)
        tab.set_fontsize(SMALL)
        tab.scale(1.0, 1.9)
        for (r, c), cell in tab.get_celld().items():
            if r and "\n" in rows[r - 1][0]:          # a name on two lines, its row twice as high
                cell.set_height(1.8 * cell.get_height())
            cell.set_edgecolor("#d0d4da")
            if r == 0:
                cell.set_text_props(fontweight="bold")
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
    p = run["problem"]
    return (f"{p.get('abscissa', 'abscissa')}" + (f", {p['abscissa_unit']}" if p.get("abscissa_unit") else ""),
            f"{p.get('value', 'value')}" + (f", {p['value_unit']}" if p.get("value_unit") else ""))


def slide_effects(run: dict, k: int):
    """Phase 2, second slide: what each term of the model does to the curves at the fit, the fitted model less the same
    model with that term at zero, on every record, beside one noise of a sample. For the first model, which has no term,
    each shared parameter moved by one bar of its fit."""
    rows = next((row for row in run["equations"] if row["model"] == f"model {k}"), None)
    arrays = run["effects"].get(k)
    first = k == 0 or not rows or not rows.get("effects") or not rows.get("after")          # a term moved by one bar where at zero it leaves no curve is no first model
    sentence = ("The first model has no term. Each panel: the curves when one of its parameters moves by one bar of its fit, less the fitted curves, on every record, beside the noise of a sample." if first else
                "Each panel: the fitted model less the same model with that term at zero, all else held, on every record. Grey: one noise of a sample.")
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
    w, h, left, top = 0.9 / cols, (TOP - 0.12) / lines, 0.06, TOP - 0.02
    for i, e in enumerate(shown):
        c, r = i % cols, i // cols
        ax = fig.add_axes([left + c * w + 0.03, top - (r + 1) * h + 0.05, w - 0.06, h - 0.09])
        for j in range(len(rows["records"])):
            key = f"e{e['key']}_{j}"
            if key in arrays:
                ax.plot(arrays[f"x{j}"], arrays[key], color=COLOURS[j % len(COLOURS)], lw=1.6, label=label_of(run, rows["records"][j], room=22))
        biggest = size.get(e["key"], 0.0)
        if noise and noise <= 5 * biggest:
            ax.axhspan(-noise, noise, color=GREY, alpha=0.2, lw=0)
        elif noise and biggest:
            ax.text(0.98, 0.96, f"one noise of a sample: {noise / biggest:.3g} times the largest change", transform=ax.transAxes, ha="right", va="top", fontsize=12, color=GREY)
        ax.axhline(0.0, color=GREY, lw=0.8)
        said = {"at zero": f"$a$ = {pm(e['value'], e['bar'])}", "held": f"$a$ = {pm(e['value'], e['bar'])}, held from theory"}.get(e["how"], f"{e['parameter']} + one bar, {pm(e['value'], e['bar'])}")
        ax.set_title((e["name"] if len(e["name"]) < 44 else e["name"][:41] + "...") + "\n" + said, fontsize=SMALL, loc="left")
        if i == 0:
            ax.legend(fontsize=9, frameon=True, loc="lower right")          # which curve is which record (review 84, minor 3)
        if r == lines - 1 or i + cols >= n:          # the names of the axes, on the lowest panel of each column and the first of each line (review 85, minor 1, and the red team of message 310, M7)
            ax.set_xlabel(axis_names(run)[0], fontsize=12)
        if c == 0:
            ax.set_ylabel("change of the " + axis_names(run)[1], fontsize=12)
        ax.tick_params(labelsize=12)
    if len(effects) > n:
        fig.text(0.965, 0.1, f"the {n} largest of {len(effects)}", fontsize=SMALL, color=GREY, ha="right")
    return fig


def slide_twin(run: dict, k: int):
    """Phase 3: the admission of the fit by its twin, each parameter's bias over its spread, spread over its bar, and the
    share of sets whose interval holds the truth, against the bands of the rule."""
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
        ax.set_yticklabels(["\n".join(textwrap.wrap(n, 16)) for n in names] if i == 0 else [""] * len(names), fontsize=SMALL - 2)          # wrapped: a long name was cut at the left edge (the cold round 22)
        ax.set_title(label, fontsize=TEXT)
        span = max(hi - lo, 1e-9)
        finite = [v for v in vals if math.isfinite(v)] or [lo]
        ax.set_xlim(min(lo - 0.3 * span, min(finite) - 0.1 * span), max(hi + 0.3 * span, max(finite) + 0.1 * span))
        ax.invert_yaxis()
    fig.text(0.04, 0.185, "Bias over spread: how far the mean of the fits of the sets is from the value that drew them, in their scatter.\n"
             "Spread over bar: their scatter over the bar that one fit quotes. Share held: the sets whose interval holds the value that drew them.\n"
             "Admitted: on sets drawn from this model, the fit returns what drew them, so its bar is honest if the model is right. It does not test the model.",
             fontsize=SMALL - 4, color=GREY, va="top", linespacing=1.4)          # the three measures were named and not said (the cold round 22)
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


def slide_fit(run: dict, k: int):
    """Phase 4: each curve with the model fitted to all of them together, and under it what the model leaves in units
    of the noise of each sample."""
    res = run["residuals"].get(k, {})
    read = run["phases"].get(k, {}).get("read", {})
    said = chi_said(read)
    whole = "One fit of every curve together: its shared parameters are one value for all. " + said
    fig = _frame(TITLES["fit"], k, "fit", run, whole if len(textwrap.wrap(whole, WIDTH)) <= 2 else said)          # the frame shows two lines, and the bar of correlated neighbours was cut at "if its noise" (the blind red team of 0.1.0rc20, mislead 2)
    if not res:
        _no_data(fig, "The run kept no residuals of this model.")
        return fig
    curves = list(res.items())[:8]
    cols = min(4, len(curves))
    lines = math.ceil(len(curves) / cols)
    w, h = 0.92 / cols, (TOP - 0.1) / lines
    for i, (curve, a) in enumerate(curves):
        c, r = i % cols, i // cols
        x0, y0 = 0.06 + c * w, TOP + 0.02 - (r + 1) * h
        top = fig.add_axes([x0, y0 + 0.38 * h, w - 0.03, 0.48 * h])
        low = fig.add_axes([x0, y0 + 0.1 * h, w - 0.03, 0.24 * h], sharex=top)
        colour = COLOURS[i % len(COLOURS)]
        top.plot(a[:, 0], a[:, 1], ".", ms=2, color=colour, alpha=0.5)
        top.plot(a[:, 0], a[:, 2], "-", color="black", lw=1.2)
        top.set_title(label_of(run, curve, room=22), fontsize=12, loc="left")          # a quarter of the slide holds 22 characters: 32 lay on the next title (0.1.0rc8)
        top.tick_params(labelbottom=False, labelsize=11)
        low.plot(a[:, 0], a[:, 3] / np.where(a[:, 4] > 0, a[:, 4], 1.0), ".", ms=1.5, color=colour)
        low.axhline(0.0, color=GREY, lw=0.8)
        low.set_ylim(-5, 5)
        low.tick_params(labelsize=11)
        if c == 0:
            low.set_ylabel("in noise", fontsize=12)
    return fig


def slide_read(run: dict, k: int):
    """Phase 5: the residuals over the noise of each sample against the normal law, and the mean of their square in
    bands of the height of the curve, which is 1 where the noise is the one stated."""
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


def answer_label(a: dict, room: int = 40) -> str:
    """An answer of the report as the slide of where to search names its bar: its question and its answer on two lines
    of `room` characters at most, cut only past them. 0.1.0rc10 cut each at 36 characters, so that "a Gaussian width of
    0.03 widths" and "of 0.3 widths" read alike (review 85 of the repository it is made in, minor 1)."""
    lines = textwrap.wrap(f"{a['question']}: {a['answer']}", room)
    return "\n".join(lines[:2]) + ("..." if len(lines) > 2 else "")


def slide_where(run: dict, k: int):
    """Phase 6: the answers of the report, each with its excess, the line that a right model passes on 5 sets of 100,
    and the three proposals for the model and the three for the experiment."""
    table = run["tables"].get(k, {})
    where = table.get("where", {})
    answers = sorted([a for a in where.get("answers", []) if isinstance(a.get("excess"), (int, float)) and math.isfinite(a["excess"])], key=lambda a: -a["excess"])[:8]
    line = where.get("largest_excess_on_the_twin")
    fig = _frame(TITLES["propose"], k, "propose", run, "The answers of the report, the largest first. The vertical line: the largest excess of all the answers together that a right model passes on 5 sets of 100. Blue: beyond it. Grey: within it.")          # the colours were said nowhere (the cold round 22)
    if not answers:
        _no_data(fig, "The run kept no report of this model.")
        return fig
    ax = fig.add_axes([0.3, 0.2, 0.24, TOP - 0.2])
    ys = np.arange(len(answers))
    values = [max(a["excess"], 1e-2) for a in answers]
    ax.barh(ys, values, color=[ACCENT if a.get("beyond_the_twin") else GREY for a in answers])
    ax.set_yticks(ys)
    ax.set_yticklabels([answer_label(a) for a in answers], fontsize=10)
    ax.invert_yaxis()
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter          # matplotlib is asked only by the log
    ax.set_xscale("log")
    decades = math.log10(max(max(values), line or 1.0) / max(min(min(values), line or math.inf), 1e-2))          # the span of the axis, the red line among it: bars that all stood near 10 000 beside a line at 2 read as half a decade, and their ticks ran together again (the blind red team of 0.1.0rc13)
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0) if decades <= 2.0 else (1.0,)))          # 1, 2 and 5 of each decade over two decades at most: over four, "5 1020 50 100200" ran together (the blind red team of 0.1.0rc12)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())          # the labels of the minor ticks lay on one another at 30 sets (the cold round 22)
    if line:
        ax.axvline(line, color=BAD, lw=2)
    ax.set_xlabel("the excess of each answer over its own twin,\n1 at the 95th percentile of that answer alone", fontsize=10)          # in two lines under its axis: aligned on the right of a narrow axis, its first words left the slide (the blind red team of 0.1.0rc13)          # the line is the largest over every answer, the unit is each answer's own (the blind red team of 0.1.0rc12)
    props = table.get("proposals", {})
    renderer = fig.canvas.get_renderer()

    def put(x, y, text, **style):          # a text at y, and the y under it by its measured height: three proposals of three lines each left the slide by fixed steps (review 85, moderate 1)
        said = fig.text(x, y, text, va="top", **style)
        return y - said.get_window_extent(renderer=renderer).height / fig.bbox.height

    beyond = any(a.get("beyond_the_twin") for a in answers)
    heading = ("For the model, by each term's fall of chi square, not by the bars" if beyond else "For the model, though none is asked")          # the bars put the core first and this list the flanks, with no word that they rank two things (the fresh red team of 0.1.0rc20, mislead 1)
    y = put(0.58, TOP, "\n".join(textwrap.wrap(heading, 38)), fontsize=SMALL + 1, fontweight="bold") - 0.015          # three proposals under grey bars read as asked (the cold round 22)
    for k, p in enumerate(props.get("for_the_model", [])[:3]):
        y = put(0.59, y, "\n".join(textwrap.wrap(f"{k + 1}. {p.get('term', '')}", 50)[:2]), fontsize=SMALL - 2)          # each proposal under its own words, its term file beneath where it has one (review 84, minor 3)
        if p.get("as_a_term_file"):
            y = put(0.61, y - 0.004, "\n".join(textwrap.wrap(p["as_a_term_file"], 50, break_on_hyphens=False)[:2]), fontsize=SMALL - 4, family="monospace", color=GREY)          # whole, where 0.1.0rc10 cut it in a word
        y -= 0.012
    if any("width" in str(p.get("term", "")) for p in props.get("for_the_model", [])[:3]):          # three widths on one page, a parameter, the unit of the names and a deviation, and the names said none (the blind red team of 0.1.0rc14 on a fresh machine)
        y = put(0.59, y, "A width in these names: the full width of the fitted\ncurves at half their height, not a parameter of model.py.", fontsize=SMALL - 4, color=GREY) - 0.012
    y = put(0.58, y - 0.02, "For the experiment", fontsize=SMALL + 1, fontweight="bold") - 0.015
    for p in props.get("for_the_experiment", [])[:3]:
        if y < 0.16:          # above the phases at the foot: the rest are in model_K.md
            fig.text(0.59, y, "and more in the report of this model", fontsize=SMALL - 3, color=GREY, va="top")
            break
        y = put(0.59, y, "\n".join(textwrap.wrap(experiment_line(p), 54)), fontsize=SMALL - 2) - 0.012          # whole: four lines at most cut the ratio that the twin measured (the cold round 22)
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
    """Phase 7: the term tried after this model, planted at its fitted size on synthetic sets and found again by the
    whole procedure: the bias of each value in scatters, the share of sets in which one bar holds it, and the fall of
    chi square that the term gives on the records beside the fall on the sets."""
    cand = run["phases"].get(k, {}).get("candidate on the twin")
    fig = _frame(TITLES["candidate on the twin"], k, "candidate on the twin", run,
                 f"The term {cand['the_planted_candidate']}, planted at its fitted size on {cand.get('read', '?')} synthetic sets, and found again." if cand else "")
    if not cand:
        this = next((r for r in run.get("best", {}).get("iterations", []) if r.get("model") == f"model {k}"), None) or {}
        gone = this.get("refused") or []
        _no_data(fig, (f"The term tried after this model, {gone[0]['term']}, was refused before it could be planted on the twin: {gone[0]['why'][:140]}. " if gone
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
            ax.text(1.02, y, f"{share:.0%} of sets\nwithin one bar", transform=ax.get_yaxis_transform(), va="center", fontsize=11)          # "63% held" said nothing of what held (the cold round 22)
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


def slide_decision(run: dict, k: int):
    """Phase 8: the criterion of each model so far, this one marked, and what became of the term tried after it."""
    rows = run["best"].get("iterations", [])
    this = next((r for r in rows if r["model"] == f"model {k}"), None)
    kept = next((r for r in rows if r.get("after") == f"model {k}" and not r.get("aside")), None)
    tried = (this or {}).get("tried", [])
    refused = (this or {}).get("refused", [])
    if kept:
        said, colour = f"Kept: {kept['added']}. Its criterion is {(this or {}).get('criterion', 0) - kept.get('criterion', 0):.1f} lower.", GOOD
    elif refused:
        said, colour = f"Refused: {refused[0]['term']}. " + refused[0]["why"][:160], BAD
    elif tried:
        said, colour = f"Tried: {', '.join(tried)}.", GREY
    else:
        said, colour = "No term was tried after this model.", GREY
    fig = _frame(TITLES["report"], k, "report", run, "The Bayesian criterion of each model, lower is better. A term is kept where its fit ends inside its bounds, the criterion falls by more than 2 and its twin admits it.")          # the three conditions of the tutorial, where one sentence named two and a term that lowered the criterion by 0.53 read as kept (the blind red team of 0.1.0rc13)
    fig.text(0.04, 0.24, "\n".join(textwrap.wrap(said, WIDTH - 12)[:3]), fontsize=TEXT, color=colour, fontweight="bold", va="top", linespacing=1.3)          # bold letters are wider: the line ran past the edge (the cold round 22)
    ended = [r for r in rows if r.get("ended") and r.get("criterion") is not None and re.fullmatch(r"model \d+", str(r.get("model", "")))]          # a term tried alone in the place of another is no model of the run: "model 1 with 02.py in the place of 01.py" ended the log (the cold round 22)
    if len(ended) == 1:          # one model is no comparison: a plot of one point was a graph of nothing (the blind red team of 0.1.0rc12)
        fig.text(0.08, 0.5, f"One model so far, its criterion {ended[0]['criterion']:.1f}: the next term is judged against it.", fontsize=TEXT, color=GREY)
    elif ended:
        ax = fig.add_axes([0.08, 0.32, 0.6, TOP - 0.34])
        idx = [int(r["model"].split()[-1]) for r in ended]
        ax.plot(idx, [r["criterion"] for r in ended], "o-", color=ACCENT, ms=9)
        if k in idx:
            ax.plot([k], [ended[idx.index(k)]["criterion"]], "o", ms=18, mfc="none", mec=BAD, mew=2.5)
        ax.set_xticks(idx)
        ax.set_xlabel("model")
        ax.set_ylabel("criterion")
    return fig


def slide_best(run: dict):
    """The output: the best model, its terms, and each parameter with the five parts of its uncertainty side by side,
    in units of the bar of the fit."""
    best = run["best"]
    quotable = best.get("quotable") or {}
    verdict = ("Every value may be quoted with its bar of the fit." if quotable and all(quotable.values())
               else "No value is to be quoted yet: best.md says why." if quotable and not any(quotable.values())
               else "Only the values marked quoted may be quoted: best.md says why." if quotable else "")
    fig = _frame(f"The best model: {best.get('best', '?')}", -1, None, run, "Its terms: " + "; ".join(best.get("terms", [])) + ". " + verdict)          # the slide dropped the warning of best.md, and a reader quoted a value that it forbids (the blind red team of 0.1.0rc14 on a fresh machine)
    params = best.get("parameters", {})
    if not params:
        _no_data(fig, "The run left no best.json.")
        return fig
    names = list(params)[:10]
    ax = fig.add_axes([0.25, 0.16, 0.33, TOP - 0.18])
    ys = np.arange(len(names))
    parts = [("fit", lambda p: 1.0), ("bias on the twin", lambda p: abs(p.get("bias_on_the_twin_in_bars") or 0.0)),
             ("theory", lambda p: (p.get("bar_of_theory", {}).get("together") or 0.0) / (p.get("bar_of_the_fit") or 1.0)),
             ("over the models", lambda p: (p.get("spread_over_the_models") or 0.0) / (p.get("bar_of_the_fit") or 1.0))]
    for i, (label, f) in enumerate(parts):
        sizes = [f(params[n]) for n in names]
        if i and not any(sizes):          # a part of no size is no bar and no line of the legend: "theory" and "over the models" stood empty and unexplained (the blind red team of 0.1.0rc12)
            continue
        ax.barh(ys + (i - 1.5) * 0.2, sizes, height=0.2, color=COLOURS[i], label=label)
    ax.set_yticks(ys)
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.set_xlabel("in bars of the fit")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=True)
    read = best.get("searches_that_ended") or best.get("sets_of_the_twin")          # the sets that the run read: "held on 87 of 100 sets" was a coverage times 100 on a run of 30 (the cold round 22)
    lines = []
    for n in names:
        cov = params[n].get("coverage_on_the_twin")
        said = (f"{n}: {pm(params[n]['value'], params[n].get('bar_of_the_fit') or float('nan'))}"
                + (f", its interval holds the truth on {round(cov * read)} of {read} sets of the twin of best.md, drawn apart from those of the slides of the turns" if isinstance(cov, (int, float)) and read else "")
                + ((", quoted" if quotable[n] else ", not yet quoted") if n in quotable else "")
                + (f", not told apart from {', '.join(params[n]['not_told_apart_from'])}" if params[n].get("not_told_apart_from") else ""))
        lines += textwrap.wrap(said, 34, subsequent_indent="    ")          # wrapped inside the slide: its last words were cut at the right edge (the cold round 22)
    fig.text(0.62, 0.16, "\n".join(lines), fontsize=12, va="bottom", ha="left", linespacing=1.4)
    return fig


SLIDES = {"intake": (slide_records,), "identify": (slides_equation, slide_effects), "twin": (slide_twin,), "fit": (slide_fit,),
          "read": (slide_read,), "propose": (slide_where,), "candidate on the twin": (slide_next,), "report": (slide_decision,)}


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
                       "One slide for each phase of each turn of the loop, with the equation of each model and the effect of each term. `log.pdf` holds them in order.", ""]
    with plt.rc_context(_rc()), PdfPages(run["loop"] / "log.pdf") as pdf:
        for k in sorted(run["phases"]):
            lines += [f"## Model {k}", ""]
            for p, phase in enumerate(ORDER):
                for s, draw in enumerate(SLIDES[phase]):
                    drawn = draw(run, k)
                    for i, fig in enumerate(drawn if isinstance(drawn, list) else [drawn]):          # a phase may need several slides: the equation of a model of many terms (his message 311)
                        path = out / f"{k:02d}_{p + 1}{'ab'[s] if len(SLIDES[phase]) > 1 else ''}{i + 1 if i else ''}_{SHORT[phase].replace(' ', '_')}.png"
                        fig.savefig(path, dpi=dpi)
                        pdf.savefig(fig)
                        plt.close(fig)
                        made.append(path)
                        lines.append(f"- phase {p + 1}, {TITLES[phase] if s == 0 else 'what each term does'}: ![]({path.relative_to(run['loop']).as_posix()})")
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
