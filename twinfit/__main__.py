"""The command line of the package.

    python -m twinfit guide                 prints the protocol for an agent (agent/AGENT.md)
    python -m twinfit guide manual          prints the older manual path, for a model that one file model.py cannot hold
    python -m twinfit guide tutorial        prints the tutorial, twinfit/TUTORIAL.md, and `guide readme` the README, with
                                            the path of the file where the package is installed
    python -m twinfit guide rubidium        prints a real case, recorded traces of a line of rubidium as plots, and where they are
    python -m twinfit guide examples        prints the gallery, twinfit/EXAMPLES.md: projects of other fields, what to do
                                            and what not to do on each, and what the loop gave beside their truth
    python -m twinfit walk [HIDDEN] [DIR]   walks one rung of the older manual path on synthetic records with a hidden term
    python -m twinfit init DIR              writes into DIR what a project needs to start: description.toml and model.py
                                            to fill, traces/ for the curves, terms/ for the terms, the rules of the
                                            folder and the skill for an agent
    python -m twinfit example NAME DIR [--noise L] [--kind K]
                                            writes into DIR a worked project of another field: chromatography or nmr, or one
                                            of the gallery that `guide examples` prints, as clock or gel. Each
                                            holds its curves, the description of their collection and the first model. A
                                            project of the gallery takes its noise at L times its own, 0 for none, and of
                                            the kind K: white, correlated over neighbouring samples, or that of a counter
    python -m twinfit terms                 prints the terms that a term file can name, and the form of a term file
    python -m twinfit scan DIR TERM X...    fits TERM, a term with x in the place of its distance, stretch or period, as
                                            "second(x)", at each value X beside the model that the term files of DIR give
    python -m twinfit counts DIR [--floor F] writes a copy of DIR beside it whose curves carry the noise of a counter from
                                            the fitted values of the best model of the last run, the column sd
    python -m twinfit thin DIR N            writes three copies of DIR beside it, each keeping one sample in every N of each
                                            curve, from three first samples, for a noise correlated over N samples
    python -m twinfit loop DIR [--n N] [--sets S] [--max-seconds T]
                                            the loop without an agent on the description in DIR and its term files, up
                                            to N of them, each twin on S synthetic sets (loop.py). With T, a run that
                                            passes T seconds of the clock stops at the end of the step in course and is
                                            refused, and the last whole run stays in loop/
    python -m twinfit.tutorial --check      runs every command of the tutorial again and says whether an output changed
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MARK = "{{WHERE_THE_PACKAGE_IS}}"      # in the three starting files: the place of the package, which `init` writes in
PYTHON = "{{THE_PYTHON}}"               # in the rules and the skill of an agent: the Python that ran `init`, which has numpy and scipy, where the plain python of a machine may not (round 12 of the cold tests)


def thin(folder: Path, every: int) -> list:
    """Three copies of a folder beside it, `DIR_thin_N_S`, each keeping one sample in every `every` of each curve from
    the first sample S, for S at nought, a third and two thirds of `every`, made without `loop/`, `loop.before/` and
    `agent/`. Where the report leads with a noise correlated over N samples, the samples so kept are nearly apart, and a
    term that one copy keeps and another refuses is not to be believed. A blind agent wrote the copies with numpy of its
    own (round 15 of the cold tests). Returns the copies made."""
    import shutil
    import tomllib
    from .loop import THIN_KEPT as KEPT
    folder, every = Path(folder), int(every)
    if every < 2:
        raise ValueError(f"thin keeps one sample in every N, with N of 2 or more: {every} was asked")
    if not (folder / "description.toml").is_file():
        raise ValueError(f"{folder} holds no description.toml")
    files = [t["file"] for t in tomllib.loads((folder / "description.toml").read_text()).get("trace", [])]
    shortest = min((len((folder / name).read_text().splitlines()) - 1 for name in files if (folder / name).is_file()), default=0)
    if shortest and shortest // every < KEPT:          # thin 400 on curves of 201 samples wrote copies of one sample each, with no word (round 18 of the cold tests)
        raise ValueError(f"thin {every} keeps {shortest // every} of the {shortest} samples of the shortest curve, under {KEPT}, too few to fit: "
                         f"the largest N that keeps {KEPT} is {max(shortest // KEPT, 1)}")
    made = []
    for first in sorted({0, every // 3, (2 * every) // 3}):
        copy = folder.parent / f"{folder.name}_thin_{every}_{first}"
        if copy.exists():
            raise ValueError(f"{copy} is there already: thin writes new copies only")
        shutil.copytree(folder, copy, ignore=shutil.ignore_patterns("loop", "loop.before", "agent", "__pycache__"))
        (copy / "agent").mkdir()
        for name in files:
            lines = (copy / name).read_text().splitlines()
            (copy / name).write_text("\n".join([lines[0]] + lines[1:][first::every]) + "\n")
        made.append(copy)
    return made


def counts(folder: Path, floor: float = 1.0) -> Path:
    """A copy of a folder of counted records, `DIR_counts`, whose curves carry a column `sd`, the noise of each sample
    as a counter gives it, the root of the fitted value of the best model of the last run of the loop, never under the
    root of `floor`, and whose description reads the noise from that column. The root of the counts of the records
    themselves gives the lowest counts too small a noise and moves the values, as section 9 of the tutorial says, and
    all three testers of round 18 wrote this column with a script of their own. Made without `loop/`, `loop.before/` and
    `agent/`. Refused with no run of the loop, or where a sample of a curve has no fitted value."""
    import csv
    import math
    import shutil
    import tomllib
    folder = Path(folder)
    if not (folder / "loop" / "best.json").is_file():
        raise ValueError(f"{folder} holds no run of the loop: `python -m twinfit loop {folder}` first, with a model whose chi square is near 1 for each degree of freedom")
    best = json.loads((folder / "loop" / "best.json").read_text())["best"]
    table = folder / "loop" / f"residuals_{str(best).split()[-1]}.csv"
    if not table.is_file():
        raise ValueError(f"{table} is not there: the loop writes the residuals of each model it fits")
    data = tomllib.loads((folder / "description.toml").read_text())
    offset = float(data.get("problem", {}).get("abscissa_offset", 0.0))
    fitted = {}
    with table.open() as fh:
        for row in csv.DictReader(fh):
            fitted[(row["curve"], round(float(row["abscissa"]) + offset, 9))] = float(row["fitted"])
    copy = folder.parent / f"{folder.name}_counts"
    if copy.exists():
        raise ValueError(f"{copy} is there already: counts writes a new copy only")
    shutil.copytree(folder, copy, ignore=shutil.ignore_patterns("loop", "loop.before", "agent", "__pycache__", ".twinfit.lock"))
    (copy / "agent").mkdir()
    for trace in data.get("trace", []):
        name, (xcol, _) = trace["file"], trace["columns"][:2]
        with (copy / name).open() as fh:
            rows = list(csv.DictReader(fh))
        for row in rows:
            key = (trace.get("name", name), round(float(row[xcol]), 9))
            key = key if key in fitted else (name, key[1])
            if key not in fitted:
                shutil.rmtree(copy)
                raise ValueError(f"the sample at {row[xcol]} of {name} has no fitted value in {table.name}: the run of the loop read other curves than these")
            row["sd"] = f"{math.sqrt(max(fitted[key], floor)):.6g}"
        with (copy / name).open("w", newline="") as fh:
            out = csv.DictWriter(fh, fieldnames=[k for k in rows[0] if k != "sd"] + ["sd"])
            out.writeheader()
            out.writerows(rows)
    text = (copy / "description.toml").read_text()
    noise = '[noise]\nfrom = "column"\ncolumn = "sd"\n'
    import re as _re
    text = _re.sub(r"(?ms)^\[noise\].*?(?=^\[|\Z)", noise + "\n", text) if _re.search(r"(?m)^\[noise\]", text) else text.rstrip("\n") + "\n\n" + noise
    (copy / "description.toml").write_text(text)
    return copy


def init(folder: Path) -> list:
    """Writes the starting files of a project, a folder that `loop` reads once it is filled: the description of the
    data collection, the first model, a folder for the curves and one for the terms, the rules of the folder and the
    skill of an agent. A file that is already there is left as it is and is said to be. Until version 0.0.4.dev9 it
    wrote a file of the older manual ladder, which `loop` could not read (the red team of his message 169)."""
    folder.mkdir(parents=True, exist_ok=True)
    for name in ("traces", "terms", "agent"):          # agent/: the record of an agent's iterations, which the guide names
        (folder / name).mkdir(exist_ok=True)
    plan = [(HERE / "agent" / "description_template.toml", folder / "description.toml"), (HERE / "agent" / "model_template.py", folder / "model.py"),
            (HERE / "agent" / "CLAUDE_template.md", folder / "CLAUDE.md"), (HERE / "agent" / "skill" / "SKILL.md", folder / ".claude" / "skills" / "twinfit" / "SKILL.md")]
    out = []
    where = str(HERE.parent)                                   # the folder that holds the package: its path goes into each file
    held = folder / "description.toml"
    if held.exists() and not re.search(r"^\s*model\s*=", held.read_text(), re.M):          # a description that names no model starts from the generic one: a template written beside it made the loop refuse that start (round 17 of the cold tests)
        plan = [(src, folder / "model_template.py" if dst.name == "model.py" else dst) for src, dst in plan]          # the template beside, under its own name: with none an agent had nothing to start from, and as model.py it made the loop refuse the generic start (rounds 17 and 18 of the cold tests)
        out.append(f"not written: {folder / 'model.py'}, since the description names no first model and the loop then starts from the generic one of the package: "
                   f"to start from a model of your own, write model.py from {folder / 'model_template.py'}, and name it under [addenda], model = \"model.py\"")
    for src, dst in plan:
        if dst.exists():
            out.append(f"kept as it is: {dst}")
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(src.read_text().replace(MARK, where).replace(PYTHON, sys.executable))
        out.append(f"written: {dst}")
    ready = all(line.startswith(("kept as it is", "not written")) for line in out)          # "Fill description.toml" said of a project whose files were all filled (round 12 of the cold tests)
    out.append(f"made: {', '.join(name + '/' for name in ('traces', 'terms', 'agent'))}, for the curves, the term files and the record of an agent")
    out.append(f"the package is imported from {where}. It needs numpy and scipy, which {sys.executable} has. "
               + (f"Every starting file was there already and is kept: run python -m twinfit loop {folder}" if ready else f"Fill description.toml and model.py, put the curves in traces/, then: python -m twinfit loop {folder}"))
    return out


TERM_FILE = '''A term file is DIR/terms/01.py, then 02.py, and so on, one term each, read by the loop in that order:

    from twinfit import mechanisms
    term = mechanisms.along("width", "temperature")
    why = "the report: the core changes along temperature, and a width changes with the temperature of the sample"

`term` is one of the terms above, or a mechanism of your own, and `why` says why you add it, in six words or more. A
file that is kept may be written again later: the loop fits every file again at each run. The loop fits the
model with the term and keeps it when three things hold: its fit ends inside its bounds, its Bayesian criterion falls
by more than {TIED}, and the twin admits the larger model. The Bayesian criterion is the chi square plus the number of
fitted parameters times the logarithm of the number of samples, so a fall of more than {TIED} is a gain that one more
parameter does not buy alone. Otherwise it refuses the term, with the reason.

A term that the physics asks for and the records cannot show is held at its value of theory and not fitted: give the
term its theory, as mechanisms.pedestal(4.0, theory=(0.05, 0.01)), its size and the bar of that size, and write a line
hold = True in its file. The loop keeps it at that size whatever the criterion says, and carries its bar into the bar of
theory of every other value.'''


NAMED = (("widened(0.5)", (lambda design: 0.5,)), ("displaced(0.5)", (lambda design: 0.5,)), ("pedestal(8)", (8,)), ("bent()", ()), ("second(1.5)", (1.5,)),
         ("tail(2, 'upper')", (2, "upper")), ("level(1)", (1,)), ("periodic(3.0)", (3.0,)), ("periodic(3.0, 'sine')", (3.0, "sine")), ("along('width', 'temperature')", ("width", "temperature")))          # each term by the name that the reports print: a tester mapped "a kernel with exponential wings" to widened by hand (round 18 of the cold tests)


def terms() -> str:
    """The terms that a term file can name, read from the docstring of mechanisms.py, and the form of a term file
    (row M168: the output of the loop sends the researcher here to read what the parameter of a term means)."""
    from . import mechanisms
    blocks = mechanisms.__doc__.split("\n\n")
    table = next(b for b in blocks if "along(parameter, key" in b)
    at = next(k for k, b in enumerate(blocks) if b.startswith("Which term answers"))
    from .best import TIED
    named = "The name that a report and best.md give a term, and the call of a term file that writes it:\n\n" + "\n".join(f"    {getattr(mechanisms, call.split('(')[0])(*args).name:<46}{call}" for call, args in NAMED)
    return "\n\n".join(["The terms of the package, each with the one parameter a that it adds:", table] + blocks[at:at + 3] + [named, TERM_FILE.replace("{TIED}", f"{TIED:g}")])


def _span(lo: float, hi: float) -> tuple:
    """Two ends of a range, each to the digits that its half width gives at two significant digits: a width from
    0.988371 to 1.00289 was printed where every other bar of the package has two (round 8 of the cold tests)."""
    import math
    half = abs(hi - lo) / 2.0
    if not half > 0.0 or not math.isfinite(half):
        return f"{lo:g}", f"{hi:g}"
    d = max(0, 1 - int(math.floor(math.log10(half))))
    return f"{lo:.{d}f}", f"{hi:.{d}f}"


GUIDES = {"": ("agent", "AGENT.md"), "manual": ("agent", "MANUAL.md"), "examples": ("EXAMPLES.md",), "tutorial": ("TUTORIAL.md",), "readme": ("README.md",), "rubidium": ("RUBIDIUM.md",)}          # the pages that `guide` prints


def guide_page(args: list) -> Path | None:
    """The page that `python -m twinfit guide` prints for its arguments, or None for a topic it does not hold: until
    round 21 of the tests of the tutorial an unknown topic, as `tutorial`, printed the protocol of an agent."""
    topic = " ".join(args)
    return HERE.joinpath(*GUIDES[topic]) if topic in GUIDES else None


def main(argv: list) -> int:
    if set(argv) & {"--help", "-h", "help"} or not argv:          # asked for help: the usage, and no command run on a folder named --help (the cold tests of row M169)
        print(__doc__)
        return 0
    if argv == ["terms"]:
        print(terms())
        return 0
    if argv[:1] == ["counts"]:
        floor = 1.0
        if len(argv) == 4 and argv[2] == "--floor":
            floor = float(argv[3])
        elif len(argv) != 2:
            print("REFUSED: counts takes the folder, as `python -m twinfit counts DIR`, and with `--floor F` the least count whose root is a noise")
            return 1
        try:
            made = counts(Path(argv[1]), floor)
        except ValueError as why:
            print(f"REFUSED: {why}")
            return 1
        print(f"written: {made}")
        print("Each curve of the copy has a column sd, the root of the fitted value of the best model, and its description reads the noise from it. Run the loop on the copy.")
        return 0
    if argv[:1] == ["thin"]:
        if len(argv) != 3 or not argv[2].isdigit():
            print("REFUSED: thin takes the folder and N, as `python -m twinfit thin DIR 10`, and keeps one sample in every N")
            return 1
        try:
            made = thin(Path(argv[1]), int(argv[2]))
        except ValueError as why:
            print(f"REFUSED: {why}")
            return 1
        print("\n".join(f"written: {m}" for m in made))
        print(f"Run the loop on each copy. A term that one keeps and another refuses is not to be believed. Quote the values of one copy with the half range of the {len(made)} copies beside them: "
              f"each copy holds one sample in {argv[2]}, so its bars are about the root of {argv[2]} times those of the whole folder, which were too small.")
        return 0
    if argv[:1] == ["scan"] and len(argv) in (2, 3):          # a scan with no values printed the whole usage (round 12 of the cold tests)
        print("REFUSED: scan takes the folder, the term with x in the place of the number scanned, and the values of x, as `python -m twinfit scan DIR \"second(x)\" 1 1.5 2`")
        return 1
    if argv[:1] == ["scan"] and len(argv) >= 4:
        import warnings
        from .loop import _pm, scan
        options = [v for v in argv[3:] if v.startswith("--")]
        if options:          # REFUSED: could not convert string to float: '--sets' (round 8 of the cold tests)
            print(f"REFUSED: scan takes the folder, the term and its values alone, and no option: {', '.join(options)}. It runs no twin, so it reads no sets")
            return 1
        try:
            with warnings.catch_warnings():          # a fit that cannot give a bar printed a warning of numpy into the table (round 3 of the cold tests)
                warnings.simplefilter("ignore", RuntimeWarning)
                out = scan(Path(argv[1]), argv[2], [float(v) for v in argv[3:]])
        except ValueError as why:
            print(f"REFUSED: {why}")
            return 1
        except TypeError as why:          # a term written with too few numbers ended in a traceback (round 12 of the cold tests)
            print(f"REFUSED: {argv[2]} is not a term of the package with x in the place of one of its numbers: {why}. `python -m twinfit terms` lists each term with its numbers, as along(\"width\", \"temperature\", x), whose x is the power")
            return 1
        print(f"{out['expression']} beside the model of {argv[1]}" + (f" with {', '.join(out['beside'])}" if out["beside"] else ", its first model") + ", at each x, in the unit that the term takes: "
              "a distance or a period in the unit of the abscissa, a stretch or a power as a number. The fall of the Bayesian criterion is from the model without the term, "
              "the size is that of the term and the bar its bar of the fit. No twin is run here: the loop decides.")
        print("| x | its fit | fall of the Bayesian criterion, negative where it rises | size | bar |")
        print("|---|---|---|---|---|")
        for r in out["rows"]:
            fall = "none" if r["fall_of_the_criterion"] is None else f"{r['fall_of_the_criterion']:.1f}"
            size = "none" if r["size"] is None else " | ".join(_pm(r["size"], r["bar"]))
            print(f"| {r['x']:g} | {r['how']} | {fall} | {size} |" if r["size"] is not None else f"| {r['x']:g} | {r['how']} | {fall} | none | none |")
        ended = [r for r in out["rows"] if r.get("ended") and r.get("fall_of_the_criterion") is not None]
        if len(ended) >= 2 and max(r["fall_of_the_criterion"] for r in ended) - min(r["fall_of_the_criterion"] for r in ended) <= 1e-6 * (1.0 + max(abs(r["fall_of_the_criterion"]) for r in ended)):
            print("Every value of x gave the same fit: the term does not read x. A term of your own is scanned through term_at(x), which writes x where the number goes, "
                  "and a number that the term reads from HELD is not moved by the scan: scan first, then hold the number that the scan gives.")          # four rows alike and "scan further" for a term that read HELD (round 14 of the cold tests)
            return 0
        b = out["best"]
        from .best import TIED
        if b is None:
            print("No fit ended inside its bounds.")
        elif b["fall_of_the_criterion"] <= TIED:
            print(f"No value of x lowers the Bayesian criterion by more than {TIED:g}, the least that the loop keeps: these records do not ask for this term at any x scanned.")
        else:
            lo, hi = out["within_one"]
            at = {row["x"]: row.get("shared") or {} for row in out["rows"]}
            moved = {name: [at[x][name] for x in (lo, hi) if name in at.get(x, {})] for name in (at.get(lo) or {})} if lo != hi else {}
            xs = sorted(r["x"] for r in out["rows"] if r.get("ended"))
            near = b.get("scanned", b["x"])
            step = min((abs(x - near) for x in xs if x != near), default=float("inf"))
            fine = out.get("from_a_parabola") and step <= hi - lo          # the same advice at a step of 0.002 as at 0.05, and no stop (round 12 of the cold tests)
            gaps = [b2 - a2 for a2, b2 in zip(xs, xs[1:])]
            uneven = len(gaps) >= 2 and max(gaps) > 1.5 * min(gaps)          # 1.53 from -2 -1.4 -1 1 1.4 2, where an even finer scan gave 1.4 (round 15 of the cold tests)
            if out.get("from_a_parabola") and fine:
                print(f"Only x = {near:g} lies within 1 of the largest fall, and the values scanned about it are spaced within the range that a parabola through the three falls about the largest gives: "
                      f"its top, {b['x']:g}, is the value to write." + (" The values scanned are spaced unevenly, which moves the top of a parabola: scan values evenly spaced about it before it is written." if uneven else ""))
            elif out.get("from_a_parabola"):
                print(f"Only x = {near:g} lies within 1 of the largest fall: the values scanned are farther apart than the range, which a parabola through the three falls about the largest gives, "
                      f"and its top, {b['x']:g}. A parabola through three coarse falls gives the top roughly: scan values finer than that range about the top, and write the top of that finer scan.")          # "no finer scan is needed" on a range that missed the truth (round 11)
            elif lo == hi:
                print(f"Only x = {lo:g} lies within 1 of the largest fall: the values scanned are too far apart to give a range, so scan a finer set about it.")
            print(f"The largest fall: at x = {b['x']:g}. The values of x whose fall is within 1 of it, which is one standard deviation of x as the profile gives it, the criterion of the best fit at each x, if the model is otherwise right: {' to '.join(_span(lo, hi))}."
                  + (" That range reaches an end of the values scanned: scan further." if out["at_the_edge"] else "")
                  + " Where the range is wide, these records fix the size and x together and not each. "
                  + ("After the finer scan, write its top in the place of this one. " if out.get("from_a_parabola") and not fine else "")
                  + (f"To use it, write `{b['line']}` in {argv[2]}, then run the loop." if argv[2].endswith(".py") else f"To try it, write `{b['line']}` in the next term file, then run the loop."))
            if out.get("next_to_no_end"):          # the best beside values whose fits did not end, written as a top (round 17 of the cold tests)
                print(f"The largest fall stands beside x = {', '.join(f'{x:g}' for x in out['next_to_no_end'])}, where the fit did not end: the best may lie where these records give no fit, "
                      "so read it as the edge of the values that the fits reach, and not as a top, before it is written.")
            if any(len(v) == 2 for v in moved.values()):          # the other values across the range: a width quoted at one distance moved 1.2 bars across it (round 5 of the cold tests)
                print("Across that range the other values of the model move: " + ", ".join(f"{name} from {' to '.join(_span(v[0], v[1]))}" for name, v in moved.items() if len(v) == 2)
                      + ". Half of each move is the part of its uncertainty that the chosen x carries: quote it beside the five parts of best.md.")
        return 0
    if argv[:1] == ["guide"]:
        page = guide_page(argv[1:])
        if page is None:
            print(f"REFUSED: no guide {' '.join(argv[1:])!r}: `python -m twinfit guide` takes nothing, or one of {', '.join(t for t in GUIDES if t)}", file=sys.stderr)
            return 2
        print(page.read_text())
        if argv[1:] in (["tutorial"], ["readme"], ["examples"], ["rubidium"]):
            print(f"This page is the file {page}.")          # where an installed package keeps it: a tester of round 21 found no twinfit/ folder beside the one that pip wrote
        if argv[1:] == ["rubidium"]:
            print(f"Its figures are the files of the folder {page.parent / 'examples' / page.stem.lower()}.")          # what the page and the README promise: until review 30 the command printed no folder of the figures
        return 0
    if argv[:1] == ["walk"]:
        from .examples.walkthrough import walk
        try:
            walk(argv[1] if len(argv) > 1 else "one block wider", argv[2] if len(argv) > 2 else None)
        except KeyError as why:          # KeyError: 'nonsense', the only trace of the code in 25 mistakes of a newcomer (round 9)
            from .examples.two_scales import HIDDEN
            print(f"REFUSED: walk hides no term {why}: it hides one of {', '.join(repr(h) for h in HIDDEN)}")
            return 1
        return 0
    if argv[:1] == ["loop"] and len(argv) >= 2 and argv[1].startswith("--"):
        print(f"REFUSED: the folder comes first, then the options: python -m twinfit loop DIR {' '.join(argv[1:])}")
        return 1
    if argv[:1] == ["loop"] and len(argv) >= 2:
        from . import loop as the_loop
        from .loop import run
        the_loop.TALK = lambda text: print(text, file=sys.stderr, flush=True)          # each step as it starts, on the error stream, so that the output stays the summary
        opts = dict(zip(argv[2::2], argv[3::2]))
        if len(argv[2:]) % 2 or set(opts) - {"--n", "--sets", "--max-seconds"}:
            print(__doc__)
            return 1
        import warnings
        try:
            with warnings.catch_warnings():          # a start that a search leaves printed a warning of numpy before the refusal in words (round 5 of the cold tests)
                warnings.simplefilter("ignore", RuntimeWarning)
                out = run(Path(argv[1]), n=int(opts.get("--n", 10)), sets=int(opts.get("--sets", 30)), max_seconds=float(opts["--max-seconds"]) if "--max-seconds" in opts else None)          # a budget of seconds, his message 181
        except KeyError as why:          # REFUSED: 'pressure', where a term of the researcher read a setting that the curves do not carry (round 7 of the cold tests)
            print(f"REFUSED: the model or a term reads {why} where it is not, most often from `design`, the settings of a curve: the settings are those that the description gives each trace")
            return 1
        except (ValueError, FileNotFoundError, ArithmeticError, RuntimeError, __import__('numpy').linalg.LinAlgError) as why:          # a refusal is said in its words, with no trace of the code (review 21), and a fault that the model raises too (round 3)
            print(f"REFUSED: {why}")
            return 1
        kept = ", ".join(k["term"] + (", held at its theory" if "held_at" in k else "") + (", which no answer of the report before it asked for" if k.get("asked") is False and "held_at" not in k else "")
                         + (", of a negative size, which takes its shape away from the curve: best.md says what that asks" if (out.get("output") or {}).get("parameters", {}).get(k.get("parameter"), {}).get("value", 0.0) < 0.0 and k.get("parameter", "").startswith(("tail ", "pedestal ", "second component at ")) else "") for k in out["kept"]) or "nothing"          # outside the f-string: a line break inside one is the grammar of Python 3.12, and the package runs from 3.11
        print(f"kept: {kept}. Refused: {', '.join(r['term'] for r in out['refused']) or 'nothing'}. The best model: {out['best']}."
              + (" best.md gives the reason of each refusal." if out["refused"] else "")
              + ((" The check of the noise of the best model fails with every band under 1: the noise is stated too large"
                  + (", as the differences of neighbouring samples give it on a sharp curve" if out.get("noise_from", "neighbours") == "neighbours" else f", as the description gives it from its {out.get('noise_from')}")          # a noise given as a value was said to come from neighbouring samples (round 19)
                  + ", so the bars are too large and the twin sees less than the records hold. "
                     "Estimate it again from the residuals of this model, as section 9 of the tutorial says, before a refusal is believed."
                     if all(r["mean_squared_residual"] < 1.0 for r in out["noise"].get("bands", [])) and out["noise"].get("bands") else
                     " The check of the noise of the best model fails, which a missing term causes as well as a noise stated wrong, so what was kept is not to be believed before that is read: best.md says how.")
                    if out["kept"] and not out["noise"]["right"] else "")          # too large said apart from too small, where neighbours stated it 1.2 times the truth (round 13 of the cold tests)
              + (f" A refusal by the twin on {out['sets']} sets can turn on more: run with --sets 100 before setting that term aside." if out.get("sets", 100) < 100 and any("the twin does not admit" in r["why"] or "synthetic sets, and a term is read on" in r["why"] for r in out["refused"]) else ""))          # the noise steered to, where a term was missing, and a refusal at 30 sets that 100 turned (round 8)          # a cosine kept on 16 spikes with no word (round 7)
        if out.get("own_trends"):          # the report reads the shared parameters alone, and said the records ask for no more beside a place that moved 400 bars (round 9)
            nxt_own = " An own value moves from curve to curve, which no question of the report reads: best.md says which, and how to write it."
        else:
            nxt_own = ""
        if out.get("noise_leads") and out["write_next"]:          # "write terms/01.py" where the first kind beyond the twin was glitches (round 7 of the cold tests)
            nxt = (f" The first kind beyond the twin in the last report is {out['noise_leads']}, which a term file cannot hold: {out['report']} says what to do, in a copy of this folder,"
                   " and the loop is run on the copy.")
        elif out.get("nothing_beyond_the_twin"):
            nxt = " No answer of where, what or along in the last report stands beyond its twin: the records ask for no more, and a further term rests on theory alone." + (
                f" Families of statistics stand beyond it all the same, {', '.join(out['seen_by_beyond'][:3])}, which name no term: a noise taken from the records themselves, as the root of their counts, "
                "moves with the noise of the records and is read so. Where the noise is so given, give it from the fitted values of the model, as section 9 of the tutorial says, before this is believed." if out.get("seen_by_beyond") else "")
        elif out.get("one_answer_just_past"):          # its words were the guide's, which a researcher of the loop had not met (round 3 of the cold tests)
            nxt = (" Only one answer of where, what, along or changes along, the questions that name a term, stands beyond its twin in the last report, and just past its line: a model that lacks nothing shows one such answer on some of its synthetic sets, "
                   "so the records ask for no more here, and the loop stops." + (f" A term that the physics asks for and the records cannot show is held at its value of theory: "
                   f"{out['write_next']} would hold it, with hold = True and the term's theory, as `python -m twinfit terms` says." if out["write_next"] else ""))
        else:
            nxt = f" Write {Path(argv[1]) / out['write_next']} with the next term, then run the loop again." if out["write_next"] else ""          # a path of the project, read from the folder above it (round 16 of the cold tests)
        if out.get("seen_by_beyond") and out.get("noise_from_a_column") and not out.get("nothing_beyond_the_twin"):          # a column of the root of the observed counts, and the command said write the next term (round 17 of the cold tests)
            nxt += (f" Families of statistics read beyond the twin, {', '.join(out['seen_by_beyond'][:3])}, with a noise given as a column: where it comes from the values of the records, "
                    "as the root of their counts, those families read it, so give it from the fitted values, as section 9 of the tutorial says, before writing a term for them.")
        if out.get("ties"):          # the records ask for no more, said beside two causes that tie (round 16 of the cold tests)
            nxt += f" Two causes fit alike, {', '.join(out['ties'])}: the records do not choose between them, best.md says which, and the physics or a design that moves one of them alone chooses."
        if out.get("better_set_aside"):
            nxt += f" A refused term is better than a kept one, {', '.join(out['better_set_aside'])}: best.md says how to take the kept one out."
        hair = [r["term"] for r in out["refused"] if r.get("by_a_hair")]
        if hair:          # "write terms/03.py" after a right term refused by 0.005 of coverage (round 13 of the cold tests)
            nxt = (f" The twin refused {', '.join(hair)} by less than the standard error of the statistics that failed: run the loop again with --sets 500 before that term is set aside"
                   + (f" and before {out['write_next']} is written." if out["write_next"] else "."))
        print(f"The report of each model and the output are in {Path(argv[1]) / 'loop'}." + nxt + nxt_own)
        if out.get("misnamed"):
            print(f"Not read: {', '.join('terms/' + n for n in out['misnamed'])}, since a term file is named with two digits, as terms/02.py.")
        if out.get("not_read_past_n"):
            print(f"Not read: {', '.join('terms/' + n for n in out['not_read_past_n'])}, since --n {out['n']} reads {out['n']} term files.")
        if out.get("not_read"):
            print(f"Not read: {', '.join(out['not_read'])}, since {Path(out['write_next']).name} is missing and the loop reads the term files in their order up to the first that is.")
        return 0
    if argv[:1] == ["example"] and len(argv) in (3, 5, 7):
        from .examples import fields, gallery
        options = dict(zip(argv[3::2], argv[4::2]))          # --noise L and --kind K, his message 194
        if set(options) - {"--noise", "--kind"} or (options and argv[1] not in gallery.EXAMPLES):
            print(f"REFUSED: an example of the gallery takes --noise, a multiple of its own noise, 0 for none, and --kind, one of {', '.join(gallery.KINDS)}: python -m twinfit example clock DIR --noise 3 --kind correlated")
            return 1
        if argv[1] in gallery.EXAMPLES:          # an example of the gallery, his message 188
            try:
                level = float(options.get("--noise", 1.0))
                truth = gallery.write(argv[1], Path(argv[2]), level, options.get("--kind"))
            except ValueError as why:
                print(f"REFUSED: {why}")
                return 1
            noise = ", with no noise drawn" if level == 0 else f", its noise {gallery.KIND_WORDS[options.get('--kind') or gallery.kinds_of(argv[1])[0]]}" + ("" if level == 1 else f", {level:g} times as large as its own noise")
            print(f"written into {argv[2]}{noise}: description.toml, model.py, an empty terms/ and the curves in traces/. Its truth and what the loop gives are in twinfit/EXAMPLES.md. Next: python -m twinfit loop {argv[2]}")
            return 0
        if argv[1] not in fields.FIELDS:          # the refusal named the two examples of fields.py and not the nine of the gallery (round 20 of the tests of the tutorial)
            print(f"REFUSED: no example {argv[1]!r}: the examples are {', '.join(list(fields.FIELDS) + list(gallery.EXAMPLES))}, and `python -m twinfit guide examples` says what each holds")
            return 1
        try:
            written = fields.write(argv[1], Path(argv[2]))
        except ValueError as why:
            print(f"REFUSED: {why}")
            return 1
        curves = sum(1 for p in written if p.parent.name == "traces")
        print(f"written into {argv[2]}: description.toml, model.py, an empty terms/, and {curves} curves in traces/. Next: python -m twinfit loop {argv[2]}")
        return 0
    if argv[:1] == ["init"] and len(argv) == 2:
        print("\n".join(init(Path(argv[1]))))
        return 0
    print(__doc__)
    return 1


def entry() -> None:
    """The command `twinfit` of an installed package, the same as `python -m twinfit`."""
    sys.exit(main(sys.argv[1:]))          # the words after the command


if __name__ == "__main__":
    entry()
