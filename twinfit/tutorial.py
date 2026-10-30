"""The tutorial of the package, written from runs (row M168, his message 169: "prepare a tutorial for the researchers help
them to understand how to use it. I'll be the beta tester of the tutorial, and I'm so severe, so make sure it should be
flawless and straightforward before showing it to me").

    python -m twinfit.tutorial --write    runs every command of the tutorial in an empty folder and writes TUTORIAL.md
    python -m twinfit.tutorial --check    runs them again and says whether any output differs from TUTORIAL.md

With no flag it prints these two lines and runs nothing: a newcomer who looked for a help text rewrote the tutorial
that he was reading (round 7 of the cold tests).

Every output in the tutorial is what its command printed, and every number in its words is read from a file that a
command wrote: none is typed. A sentence that reads a result, as "the report puts the misfit along the temperature",
is a claim, and the builder refuses to write the page when the files do not hold it. A command is shown as the
researcher types it, `python -m twinfit ...` in the folder of the tutorial, and is run with the interpreter that runs
this file, the package on its path. The seconds that the commands took are said once, at the head, and the check
leaves them out, since they change from run to run.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAGE = HERE / "TUTORIAL.md"

CAPTIONS = {          # what the reader takes from each slide of the visual log, in italics under it, the style of the public repository of the analysis (his message 432). No number: a slide drawn again leaves its caption true
    "00_1_records.png": "The four spectra of the project, one for each temperature, as its description reads them.",
    "00_2a_model.png": "The equation of the first model, one Lorentzian line of one width for all four spectra.",
    "00_2b_model.png": "What each term of the first model does to the curves.",
    "00_3_twin.png": "On synthetic sets drawn from its fit, the first model returns what drew them: the twin admits it.",
    "00_4_fit.png": "The joint fit of the four spectra by the first model.",
    "00_5_residuals.png": "What the first model leaves exceeds the noise at the top of the lines.",
    "00_6_where.png": "The answers beyond the twin, the core and the flanks changing along the temperature first.",
    "00_7_next_term.png": "The width along the temperature, planted on the twin and found again: its fall of chi square on the records lies within the twin's.",
    "00_8_decision.png": "The width along the temperature lowers the criterion and is kept.",
    "01_1_records.png": "The same four spectra, read by model 1.",
    "01_2a_model.png": "The equation of model 1, the width of the first model times one plus a times z.",
    "01_2b_model.png": "What the new term does: the curves wider on the hot records and narrower on the cold ones.",
    "01_3_twin.png": "The twin admits model 1, as a term is kept only where it does.",
    "01_4_fit.png": "The joint fit of the four spectra by model 1.",
    "01_5_residuals.png": "What model 1 leaves, against the noise.",
    "01_6_where.png": "No answer stands beyond the twin after model 1, so no term is asked.",
    "01_7_next_term.png": "No term is tried after model 1, and the loop ends.",
    "01_8_decision.png": "Model 1 holds the lowest criterion of the run.",
    "99_best.png": "The best model, model 1, each value with its bar of the fit and its bias on the twin.",
}

SECONDS = re.compile(r"ran in [0-9.]+ seconds")


def _run(folder: Path, command: str) -> str:
    """The output of `command`, typed as the researcher types it, run in `folder` with this interpreter and the package
    on its path. Its standard output and its standard error, as the terminal shows them, with the place of the package
    and that of the Python written as words, since they differ from machine to machine."""
    words = _split(command)
    if words[0] != "python":
        raise ValueError(f"a command of the tutorial starts with python: {command!r}")
    argv = [sys.executable] + words[1:]
    env = dict(os.environ, PYTHONPATH=str(HERE.parent), OMP_NUM_THREADS="1")
    done = subprocess.run(argv, cwd=folder, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)          # one stream, in the order that a terminal shows: the steps of the loop came after its summary (round 14 of the cold tests)
    out = done.stdout.replace(str(HERE.parent), "/the/folder/that/holds/twinfit").replace(sys.executable, "/the/python/that/runs/it").replace(str(folder) + "/", "").replace(str(folder), ".")
    if done.returncode != 0:
        raise RuntimeError(f"the command {command!r} of the tutorial failed:\n{out}")
    return out.rstrip("\n")


def _split(command: str) -> list:
    """A command split as a shell splits it, for the quoting that the tutorial uses: double quotes around one word."""
    return [a.strip('"') for a in re.findall(r'"[^"]*"|\S+', command)]


def _tree(folder: Path, top: str) -> str:
    """The files and folders under `top`, one a line, sorted, as `find` would list them."""
    base = folder / top
    return "\n".join([top + "/"] + [str(p.relative_to(folder)) + ("/" if p.is_dir() else "") for p in sorted(base.rglob("*")) if "__pycache__" not in p.parts])


def _part(text: str, start: str, end: str | None) -> str:
    """The lines of `text` from the line that starts with `start` to the line before the one that starts with `end`."""
    lines = text.splitlines()
    i = next(k for k, s in enumerate(lines) if s.startswith(start))
    j = next((k for k in range(i + 1, len(lines)) if end is not None and lines[k].startswith(end)), len(lines))
    return "\n".join(lines[i:j]).rstrip()


def _json(folder: Path, path: str) -> dict:
    return json.loads((folder / path).read_text())


def _claim(holds: bool, what: str) -> None:
    if not holds:
        raise AssertionError(f"the tutorial says that {what}, and the files of its run do not hold it: the words are to be written again")


def build(folder: Path, figures: Path | None = None) -> str:
    """Runs the tutorial in `folder`, which is empty, and returns its page. With `figures`, the slides of the visual log
    that the page shows are copied there, as `--write` does into examples/figures of the package."""
    from .examples import fields
    from .best import TIED
    from .where import MISFIT_SETS, SMOOTH, TIE as TIE_PLACES
    from .loop import _pm
    began = time.time()
    out = []

    def say(text):
        out.append(text.strip("\n"))

    def run(command):
        said = _run(folder, command)
        out.append(f"```bash\n{command}\n```\n\n```\n{said}\n```")
        return said

    def show(path, text=None, kind=""):
        out.append(f"`{path}`:\n\n```{kind}\n{(folder / path).read_text().rstrip() if text is None else text}\n```")

    def write(path, text):
        (folder / path).write_text(text)
        out.append(f"`{path}`:\n\n```python\n{text.rstrip()}\n```")

    say("""# The tutorial of twinfit

This tutorial takes you from nothing to the best model of a set of curves, with the uncertainty of each of its parameters and six proposals, three for the model and three for the experiment. It runs on two worked projects that the package writes, one of nuclear magnetic resonance and one of chromatography, so nothing of yours is needed to follow it. Their curves were drawn from a known truth, so you can see whether the package finds it.

It builds on the five minutes of the README and sets out the whole path of the package, its words, its reports and its choices. Not covered here: the protocol of an agent, which `agent/AGENT.md` gives, and the projects of other fields, which `EXAMPLES.md` gives.

> Section 1 below defines every word that the package and its other pages use.

Every command below was run as it is written, in an empty folder, and every output is what it printed. `python -m twinfit.tutorial --check`, with matplotlib installed by `pip install "twinfit[log]"`, runs them all again, in a temporary folder, and says whether any output changed. Here the whole of it {SECONDS}. A machine busy with other work takes longer.""")
    say(f"""## 1. Purpose and vocabulary of the package

You give it curves recorded across a design, such as spectra at several temperatures, and the model that should describe them. It fits the model to every curve at once. It reads the structure of what the fit leaves, the residuals, with families of statistics: generalized moments of the residuals and of their transformations, and windows along the abscissa and along the value. It reads how that structure changes from one setting of the design to the next.

It sets every reading against a twin: synthetic curves drawn from the fitted model with the noise of the records, fitted and read the same way. The twin draws the noise of each sample from a normal law of the deviation that the description gives it, apart from the others: where a counter gives few counts, as under about 20, its noise is skewed, and the twin does not hold that.

The loop, `python -m twinfit loop`, does all of it on a folder: it fits the first model, reads it against its twin, writes a report of where to search, then fits each term that you write in turn, keeps or refuses it, and writes the best model with the uncertainty of each parameter. Where the records stand beyond the twin, the model lacks something, and the report says where in the curves, of which kind, and along which setting. You write a term in answer. The package fits it, keeps or refuses it, and says what it does to the other parameters.

The package names a place, a kind and a setting. It never names a cause: which term of physical origin fills them is your choice, or the choice of an agent in your place.

The words that the rest of this page uses, each once:

- **A curve** is one recording, a spectrum or a chromatogram, as a column of abscissas and one of values. The description calls its file a trace, and the package calls a curve with its noise a record: the three words name one thing.
- **A setting** is a quantity under which a curve was taken and that the description gives for it, as the temperature of the sample. The package calls a setting that changes from curve to curve a covariate. The settings of one curve, by name, are its `design`, the dictionary that the model receives.
- **A bar** is the standard deviation of a parameter that a fit gives from the curvature of its chi square. **First order** is the model taken as straight in its parameters about the fit, so that bar is a first order bar.
- **A shared parameter** is one value for every curve, as a width. **An own parameter** is one value for each curve, as its height. The known parameters are the shared ones: those that the analysis is after, known in the sense that the model names them, and not that their values are known. The value of an own parameter in one curve is an own value, and the pages say an own value moves where it changes from curve to curve beyond its bars.
- **An ansatz** is a value held before the fit, with its uncertainty, from outside the records: section 9 writes one. **A lever** is a change of the experiment that the proposals weigh on the twin, as a quieter detector: section 8 reads them.
- **The coverage and the bias on the twin** of a parameter: on synthetic sets drawn from the fitted model, the share on which its interval of 1.96 bars of the fit holds the value that made them, near 0.95 where the bar is right, and the mean of its fitted values less that value, in bars. Section 6 reads them.
- **A scatter** is the standard deviation of a parameter over the synthetic sets of the twin, the spread of the estimate. A report gives a bias on the twin in scatters, and `best.md` in bars of the fit.
- **The level** of a curve is its baseline.
- **A generalized moment** is a sum over the samples of a residual times a weight of the abscissa: a power of the distance from the center of the curve, a fractional power, an inverse or a logarithmic weight, a cosine or a sine. It is also taken of the derivative of the curve, of its running integral, of its logarithm, and as a ratio of two moments.
- **A window** keeps the samples of a pair of intervals at equal distances from the center, a window along the abscissa, or the samples whose fitted value lies between two levels of its height, a window along the value.
- **A kernel** is the weight of a generalized moment. The report names each family of statistics by its kernels and its windows, as "moments: powers 0 to 8", "inverse and logarithmic kernels", "the two sides subtracted: pairs of intervals" or "the crossings of nine levels", the abscissa where the curve crosses each of nine levels of its height.
- **The excess** of an answer is its statistic on the records less the median of the same statistic on the synthetic sets, over the span from that median to their 95th percentile. An answer is beyond the twin where its excess passes a line, the largest excess that the synthetic sets of a right model reach on 95 of 100 of them: each report prints it, read on the twin of its model.
- **Two kinds tie** where their excesses lie within a factor of 1.5 of each other, and **two places tie** where their excesses differ by less than {TIE_PLACES:g} standard deviations of their difference, read on {MISFIT_SETS} sets of a twin that carries the misfit, the fit plus what it leaves averaged over {SMOOTH} samples. The records then do not tell them apart, the report names both, and a term of either is as asked.
- **The five parts of the uncertainty** of a value: its bar of the fit, its coverage and its bias on the twin, its bar of theory where values are held, and the spread over the models that the records cannot tell apart. Section 6 reads them.
- **The Bayesian criterion** of a model is its chi square plus the number of its fitted parameters times the logarithm of the number of samples. A model with one term more is kept when its criterion falls by more than {TIED:g}.
- **The gain in chi square** of a proposal is how far the chi square falls when the shape that the residuals hold there is added to the model with one free size and every parameter fitted again. It ranks the directions beyond the twin, and it is not the fall for a term of your own, which the next run of the loop gives.
- **A profile** of a number that is chosen and not fitted, as the distance of a second component, is the best criterion at each value of it: the values within 1 of the best are its standard deviation. `python -m twinfit scan` draws it, section 10.
- **A Mechanism** is a term written in Python, section 11. **A rung** is one model of the loop: model 0, the first, and one rung more for each term kept on it.
- **Glitches** are a few samples far from the curve in one record. **A noise correlated from sample to sample** is a noise that moves neighboring samples together. The report names both among the kinds of a misfit, and no term of the curve fills them: they are read in the noise.
- **A chi function** is the cost that a fit gives each residual over its noise: its square in a plain fit, another function of it in a fit that weighs large residuals less. Where the residuals are those of the noise, the two fits agree.""")
    say("""## 2. Requirements and installation

Python 3.12 or later, with numpy 2.5 or later and scipy 1.16 or later, and matplotlib for the visual log of a run: numpy 2.5 convolves about ten times faster than numpy 1.26 on a Mac with an arm64 processor, and about as fast on Linux with an x86_64 processor, and the oldest versions that the package has run on are Python 3.12, numpy 2.5.0 and scipy 1.16.0.

Where the Python of your machine is called `python3`, or lives at a path of its own, type that in the place of `python` in every command of this page. Install the package from PyPI, `pip install twinfit`, or `pip install "twinfit[log]"` for matplotlib with it, or from a folder that holds it, `pip install /the/folder/that/holds/twinfit`, or put that folder on the path of Python, `export PYTHONPATH=/the/folder/that/holds/twinfit`. Then this prints the version of the package:""")
    run('python -c "import twinfit; print(twinfit.__version__)"')
    say("""## 3. A worked project

This writes a project of nuclear magnetic resonance into a new folder, `my_nmr`:""")
    run("python -m twinfit example nmr my_nmr")
    say("It writes what a researcher writes before any fit:")
    out.append(f"```\n{_tree(folder, 'my_nmr')}\n```")
    say("**The curves**, one file a curve in `traces/`, with a line of names and then two columns, the abscissa and the value. The first lines of one:")
    show("my_nmr/traces/temperature_290.csv", "\n".join((folder / "my_nmr/traces/temperature_290.csv").read_text().splitlines()[:4]))
    say("""**The description of the data collection**, `description.toml`. It says what was measured, the abscissa and the value with their units, and each setting under which a curve was taken, here the temperature of the sample. Then one table a curve, with its file, its columns and its setting, and the noise of each sample, here one standard deviation for every sample. Then the addenda: the first model, which every description names, and what is optional, here three changes of the experiment that the proposals for the experiment weigh against each other.

A key that the package does not read is refused, so that a misspelled one is not dropped in silence. The order in which the curves were taken goes under `[problem]` as `order`, since a setting taken in the order of the time cannot be told from a drift:""")
    show("my_nmr/description.toml", kind="toml")
    say("""**The first model**, `model.py`: the curve as a function of the abscissa `x`. Its parameters are `SHARED` by every curve or `OWN` to each. Here it is one Lorentzian line of one width for all four spectra, each with its own height, place and baseline. `CENTRE` names the own parameter that places a curve, and `centre` returns it, so that the report can tell the core of a curve from its wings. `START_SHARED` and `start_own` give the fit its starts, and `BOUNDS_SHARED` the bounds of the shared parameters, a list of the lower bounds and one of the upper.

`model`, `SHARED`, `OWN`, `START_SHARED` and `start_own` are required, `CENTRE` and `centre` come together or not at all, and the rest is optional. `NAME` says the model in a few words, as the reports name it, and `HELD` holds the values that it reads and does not fit, section 9. `record` holds a curve, `record.x` its abscissa and `record.y` its values, and `design` its settings by name, as the description gives them:""")
    show("my_nmr/model.py", kind="python")
    say("## 4. The first run")
    say("The loop runs here in the foreground, since each run of this section takes a minute or two. A run of your own records may take an hour: start it in the background with `nohup`, `nohup python -m twinfit loop DIR > DIR/run.log 2>&1 &`, and read its end with `python -m twinfit wait DIR`, as the card of the guide for agents runs it in the background, so that a closed window or the limit of a tool does not stop it.")          # the card ran it in the background and the tutorial in the foreground, and neither said why (the fresh red team of 0.1.0rc19, mislead 4)
    run("python -m twinfit loop my_nmr")
    m0 = _json(folder, "my_nmr/loop/model_0.json")
    first = m0["where"]["first_of_each_question"]
    _claim(first["changes along"] == "the core along temperature" and first["where"] == "the core" and m0["where"]["where_read_along"] == "the core along temperature" and first["what"] is None and first["shape"] is None,
           "the report puts the misfit in the core changing along the temperature, names the core as its place, read along that setting, and has no answer of shape or of the noise beyond the twin")
    p0 = m0["proposals"]["for_the_model"][0]
    _claim(p0["stands_for"].endswith("change along temperature") and p0.get("as_a_term_file") is None, "the first proposal for the model is a direction along the temperature, and names no term")
    w = p0["known_parameters"]["width"]
    _claim(abs(w["with_the_term"]["bias"]) < abs(w["without_the_term"]["bias"]),
           "on the twin of the first proposal the width is biased less with its shape than without it")
    page0 = (folder / "my_nmr/loop/model_0.md").read_text()
    say("""The loop fitted the first model, `model 0`, and wrote its report to `my_nmr/loop/model_0.md`. The report opens with the words that it uses:""")
    show("my_nmr/loop/model_0.md", _part(page0, "## The words", "## The fit"), "markdown")
    _claim(not m0["noise"]["right"] and m0["chi_square"] / m0["degrees_of_freedom"] > 1.0 + 3.0 * (2.0 / m0["degrees_of_freedom"]) ** 0.5,
           "the chi square of the first model is far above 1 and its check of the noise fails")
    say(f"""Then the fit. Its chi square is {m0['chi_square']:.1f} on {m0['degrees_of_freedom']} degrees of freedom, {m0['chi_square'] / m0['degrees_of_freedom']:.2f} for each, where a right model with a right noise gives 1 with a standard deviation of {(2.0 / m0['degrees_of_freedom']) ** 0.5:.2g}. So something is missing. The page says it in its first line:""")
    show("my_nmr/loop/model_0.md", _part(page0, "## The fit", "## Where to search"), "markdown")
    off_bands = [k for k, r in enumerate(m0["noise"]["bands"]) if abs(r["mean_squared_residual"] - 1.0) > 3.0 * r["standard_error"]]
    _claim(m0["admitted_on_the_twin"] and off_bands == [len(m0["noise"]["bands"]) - 1], "the twin admits the first model and only the highest band of the check of the noise fails")
    say("""The twin admits the fit: it gives back the width that made it. But the check of the noise fails in its highest band, where the line is tallest. That check cannot tell a noise stated wrong from a model that misses. While the report finds an answer beyond the twin, as here, read it as the model that misses, and read the check again once the term is in, as section 6 does. Until then `best.md` marks the coverage and the bias on the twin as not to be read.""")
    say("""Then the first answer of each question. An answer stands beyond the twin past a line. An excess of 1 is the 95th percentile of one answer on the synthetic sets, but the report reads many answers at once, so its line is higher: the largest excess of all the answers together that a set of a right model reaches on 95 of 100 sets. Some answer then passes the line on 5 of 100 sets of a right model, and not each answer on 5 of 100:""")
    show("my_nmr/loop/model_0.md", _part(page0, "## Where to search", "Every answer"), "markdown")
    say("""Read it row by row. No place stands beyond the twin on its own, no shape and no law of the noise, but what the fit leaves changes along the temperature, in the core of the lines and in their flanks, so the row of where names the core, read along that setting. The package says no more than that: it names no term. The term is yours, from the physics of the curves.

Where two places stand beyond the twin and the records scatter their excesses alike, the row names both, tied: the report reads the scatter of each excess on sets drawn from the fit plus what the fit leaves, and the records do not tell which of the two leads.

In a resonance, what changes the core and the flanks of a line together with the temperature of the sample is its width, which the motion of the molecules sets. The proposals for the model are the directions of the report, each measured with the shape that the residuals hold there, fitted beside the model. On their twin that shape tells what a term of that place and that setting would cost the width if it were real:""")
    show("my_nmr/loop/model_0.md", _part(page0, "## Three proposals for the model", "## Three proposals for the experiment"), "markdown")
    say("""A width in the page of the proposals is a unit of their own: the full width of the fitted curves at half their height above their median value, as the page says above its list. It is not the parameter `width` of `model.py`, which here is a full width too, and in section 10 a standard deviation.

The first direction is the one whose shape lowers the chi square most, and a term that takes part of another's misfit can lead: on records that lack a slope of the baseline, a second component written first can be kept, and fall to nothing once the slope is written. A term kept first is not looked at again when a later one arrives, so `best.md` names a kept term that has fallen within two bars of none, and the model without it is tried in a copy, as section 7 says.

""")
    from .loop import matrix as read_by_families
    mx = read_by_families(folder / "my_nmr")
    readings = {n: r for n, r in mx["readings"].items() if n != "every sample: the fit"}
    top = max(readings, key=lambda n: readings[n]["chi_square"] - readings[n]["degrees_of_freedom"])
    every = mx["readings"]["every sample: the fit"]
    _claim(readings[top]["chi_square"] / readings[top]["degrees_of_freedom"] > 5.0 * every["chi_square"] / every["degrees_of_freedom"], "a family of few directions reads the misfit many times over what the chi square of every sample reads")
    say("""The row seen by names the family of statistics that sees the misfit most. A family turns the residuals into a few numbers: the moments of the abscissa, and of the distance from the center to whole, fractional and negative powers, the same under logarithmic, exponential and Gaussian kernels, Hermite functions about the center, windows along the abscissa and along the value, the crossings of levels of the height, the moments of the derivative and of the running integral, and the moments of the square, the root, the logarithm, the reciprocal and the exponential of the value, which weigh the top of a line and its foot otherwise than the value does.

A term that the model lacks moves each family in its own way, and the pattern across the families is what the report reads. `python -m twinfit matrix` prints those readings for the model that the term files of a folder give, with no twin, in a few seconds:""")
    run("python -m twinfit matrix my_nmr")
    say(f"""Its first table is what each family reads: a chi square with its degrees of freedom, and its level at first order, which a family of many correlated statistics does not reach, so the report takes its levels from the twin. The chi square of every sample reads {every['chi_square']:.0f} on {every['degrees_of_freedom']} degrees of freedom, {every['chi_square'] / every['degrees_of_freedom']:.2f} for each, while {top} reads {readings[top]['chi_square']:.0f} on {readings[top]['degrees_of_freedom']}: a misfit of one shape stands in a few directions, which a family of few directions gathers and the sum over every sample spreads over all of them.""")
    nm = mx["noises"]["chance"]
    linear = [n for n in mx["noises"]["statistics"] if n in mx["readings"] and n != "every sample: the fit"]
    normal = list(nm["normal"].values())
    lit = [nm["correlated"][n] for n in linear]
    _claim(max(normal) <= 0.1 and nm["correlated"]["every sample: the fit"] < 0.15 and min(lit) > 0.5 and nm["correlated"]["a noise correlated from sample to sample"] > 0.9,
           "a normal noise holds every line at its level, and a correlated noise lights every family and not the plain fit")
    blind = max(nm[k][n] for k in ("heavy tails", "skewed") for n in linear)
    _claim(blind < 0.3 and nm["heavy tails"]["a noise with heavy tails"] > 0.9 and nm["skewed"]["a noise skewed to one side"] > 0.9, "no family reports heavy tails or a skewed noise, and the readings of the shape do")
    _claim(nm["grows with the signal"]["windows along the value: bands between levels of the height, each side"] > 0.9 and nm["grows with the signal"]["moments: powers 0 to 4"] < 0.1,
           "the windows along the value report a noise that grows with the signal, and the moments of powers do not")
    say(f"""The table of the noise asks the same of its law. Each law is drawn on these curves at the variance that the description states, so that only its law is wrong, and the line of each statistic is set where a normal noise of the same design leaves it on 95 draws of 100: the column of a normal noise drawn again reads from {min(normal) * 100:.0f} to {max(normal) * 100:.0f} cases of 100, the level.

A noise correlated from one sample to the next lights every family, in {min(lit) * 100:.0f} to {max(lit) * 100:.0f} cases of 100, and the plain fit in {nm['correlated']['every sample: the fit'] * 100:.0f}. So on curves whose noise is correlated, a family past its line is no sign of a missing term, and the reading of the neighbors, which reports that noise in {nm['correlated']['a noise correlated from sample to sample'] * 100:.0f}, says which it is.

Heavy tails and a skewed noise keep the covariance of the noise, so no family, linear in the residuals, reports them in more than {blind * 100:.0f} cases of 100, and the readings of the shape of the noise do. A noise that grows with the signal changes the covariance along the height, so the windows along the value report it, and the moments of the powers of the abscissa do not.""")
    say("""## 5. A term

A term is a file in `terms/`, `01.py`, then `02.py`, and so on, which you write. This prints the terms that a term file can name, each with the one parameter that it adds, and the form of a term file:""")
    run("python -m twinfit terms")
    say("The report points to the core and the flanks changing along the temperature, and the physics of a resonance names the width. `along` makes a shared parameter change along a setting, here the width along the temperature. A term file adds `why`, the reason in six words or more:")
    write("my_nmr/terms/01.py", 'from twinfit import mechanisms\nterm = mechanisms.along("width", "temperature")\nwhy = "the report: the core changes along temperature, and the width of a resonance changes with the temperature of the sample"\n')
    say("## 6. The second run")
    run("python -m twinfit loop my_nmr")
    best = _json(folder, "my_nmr/loop/best.json")
    _claim(best["best"] == "model 1", "the loop keeps the term and the best model is model 1")
    truth = fields.in_the_units_of_the_fit("nmr")
    rows = {name: (best["parameters"][name]["value"], best["parameters"][name]["bar_of_the_fit"], truth[name]) for name in truth}
    off = {name: abs(v - t) / b for name, (v, b, t) in rows.items()}
    _claim(all(o < 2.0 for o in off.values()), "the fit returns each true value within two bars")
    say("Two words of these lines are defined here, since section 9 uses them later: your list is `list_of_terms.py` of the project, the terms that you suspect, empty at the start, and matched is the question of the report that reads the residuals against each term of that list. With the list empty, no answer of matched is read. The package holds a list of its own, read only where you run the loop with `--dictionary`: it names terms, and the search for a term that you do not suspect runs without it.")          # read as the terms of terms/ and as a typo (the blind red team of 0.1.0rc16, mislead 3)
    say("The loop kept the term. Its output is `my_nmr/loop/best.md`, the best model with the five parts of the uncertainty of each parameter: the bar of the fit, the coverage and the bias on the twin, the bar of theory and the spread over the models. Their words open the page:")
    show("my_nmr/loop/best.md", kind="markdown")
    levels = fields.FIELDS["nmr"]["levels"]
    say(f"""The curves were drawn from a known truth, so the output can be checked. At the mean temperature, {sum(levels) / len(levels):g} K, the width was {truth['width']:.4g} Hz, and the fit gives {' ± '.join(_pm(*rows['width'][:2]))}, {off['width']:.1f} bars away. Its growth along the temperature, in the units of the term, was {truth['width along temperature']:.3g}, and the fit gives {' ± '.join(_pm(*rows['width along temperature'][:2]))}, {off['width along temperature']:.1f} bars away.""")
    _claim(any(t["own"] == "place" for t in _json(folder, "my_nmr/loop/best.json").get("own_trends") or []), "best.md of the worked spectra reads the places moving with the temperature")
    say("""The heights and the places of the worked spectra were drawn to rise with the temperature, each spectrum its own, as an instrument whose gain and lock drift with it: `best.md` names them as moving with the temperature, with their slopes, and they stay own here, since the question of this example is the width. Section 9 says how such a value is written into the model where it is the question.

A proposal of the last report that stands within the twin is one that the records do not ask for: its line of what it would cost if it were real names a systematic to say beside the value where the physics allows such a term, and adds nothing to the bar. The worked descriptions do not say in which order their curves were taken, so `best.md` says that it is not known.""")
    sd_t = float((sum((t - sum(levels) / len(levels)) ** 2 for t in levels) / len(levels)) ** 0.5)
    per_kelvin = rows["width along temperature"][0] * rows["width"][0] / sd_t
    true_per_kelvin = fields.TRUTH["nmr"]["width"] * fields.TRUTH["nmr"]["growth"]
    _claim(abs(per_kelvin - true_per_kelvin) < 0.1 * true_per_kelvin, "the growth per kelvin is within a tenth of the truth")
    cover30 = best["parameters"]["width"]["coverage_on_the_twin"]
    within30 = (0.95 * 0.05 / 30) ** 0.5
    say(f"""In hertz per kelvin the growth is the term's size times the width over the standard deviation of the temperatures, {sd_t:.4g} K, taken over the four records with their number: {per_kelvin:.3g} Hz per kelvin, where {true_per_kelvin:.3g} drew the curves.

Each run moves the `loop/` of the run before it to `loop.before/` while it works, puts it back where the run stops on a fault, and removes it where the run ends: so `loop/` is always the last whole run.

`loop/` holds, for each model K, its report `model_K.md`, its page for an agent `brief_K.md`, its numbers `model_K.json`, its residuals `residuals_K.csv` and what each of its terms does, `effects_K.npz`, then `equations.json`, `best.md` and `best.json`, `STATUS.json`, what the run is doing, `notebook/`, the record of every step of the run, and the visual log, `log.pdf`, its slides in `log/` and their list in `log.md`. A run that is refused writes `REFUSED.md` there, beside the files of the run before it.

`loop.before/` is removed where a run ends, so copy `loop/` where a run is to be kept. `loop.reuse/` keeps the twin, the report and the proposals of each model and the twin of each term: a run reads back what a run before it computed from the same files, records and fit, so a run after one more term file computes only what that file adds, and a run stopped at its budget goes on where it stopped.

Which number to quote: the value with its bar of the fit, when three things hold. The noise passes its check, as the line above the table says here. The coverage on the twin lies in the band of the rule of the twin at the number of sets read, the band that admits the fit: three standard deviations of a coverage about 0.95 with few sets, and 0.92 to 0.98 from 475 sets on, so that one rule says both whether the fit is admitted and whether its bar is quoted. And the bias on the twin is under one bar.

With 30 synthetic sets, the default, a coverage is known within {within30:.2f}, one standard deviation, so the coverage of the width, {cover30:.2f}, is {abs(cover30 - 0.95) / within30:.1f} of them from 0.95. The report of model 1 reads the same coverage on other sets, so its number differs: neither settles it on 30 sets, whose band is wide. Read it on more sets, and quote the run on 500 sets where the time allows, whose band is the one stated.

The section "To quote" of `best.md` writes each parameter that way, with the bar of theory and the half range over the models as separate terms where they are not zero, and says when a part forbids the quote:""")
    run("python -m twinfit loop my_nmr --sets 100")
    more = _json(folder, "my_nmr/loop/best.json")
    covers = {name: more["parameters"][name]["coverage_on_the_twin"] for name in truth}
    _claim(all(abs(c - 0.95) < 2.0 * (0.95 * 0.05 / 100) ** 0.5 + 1e-9 for c in covers.values()), "on 100 sets each coverage is within two of its standard deviations of 0.95")
    show("my_nmr/loop/best.md", _part((folder / "my_nmr/loop/best.md").read_text(), "| parameter |", "## The parameters of each curve"), "markdown")
    say(f"""On 100 sets the coverages read {' and '.join(f'{c:.2f}' for c in covers.values())}, each within two of their standard deviations of 0.95: the bar of the fit holds, and it is the uncertainty to quote.

When the noise fails its check, the coverage and the bias on the twin are marked not to be read, since the twin draws its sets at the noise that the description gives: quote no bar until the noise is given right, by a column of the curves or a value, as section 9 says, and the check passes.

`best.md` also gives the own parameters of each curve, its height, its place and its baseline, each with its bar. Each run of the loop writes `loop/` again at its own number of sets, so the sections after this one, run at the default of 30, read 30 sets.""")
    m1 = _json(folder, "my_nmr/loop/model_1.json")
    _claim(not any(a["beyond_the_twin"] for a in m1["where"]["answers"]), "the report of model 1 has no answer beyond the twin")
    say("And the report of the new model finds nothing more beyond its twin:")
    show("my_nmr/loop/model_1.md", _part((folder / "my_nmr/loop/model_1.md").read_text(), "## Where to search", "Every answer"), "markdown")
    from .log import log as draw_log
    slides = draw_log(folder / "my_nmr", dpi=60)          # the visual log of this run, his message 309: the loop drew it at 110 dots an inch, and the tutorial carries it smaller
    _claim(len(slides) == 2 * 9 + 1 and slides[-1].name == "99_best.png", "the visual log holds a slide for each phase of each of the two turns, two for the model, and the best model")
    import numpy as np
    eq = json.loads((folder / "my_nmr/loop/equations.json").read_text())
    temps = [r["settings"]["temperature"] for r in _json(folder, "my_nmr/loop/best.json")["own_parameters"]]
    with np.load(folder / "my_nmr/loop/effects_1.npz") as held:
        hot, cold = held[f"e0_{temps.index(max(temps))}"], held[f"e0_{temps.index(min(temps))}"]
    _claim(eq[1]["formula"] == r"$p \rightarrow p\,(1 + a\,z)$" and eq[1]["effects"][0]["value"] > 0 and hot.min() > -1e-12 and hot.max() > 0 and cold.max() < 1e-12 and cold.min() < 0,
           "the term of model 1 makes the hot records wider and the cold ones narrower, its formula one plus a times z")
    say(f"""The loop drew its visual log too, `my_nmr/loop/log.pdf`, and `python -m twinfit log my_nmr` draws it again from what the run left, fitting nothing. For each model it gives one slide for each of the eight phases of its turn: the records, the model with its equation and what each term does to the curves, its twin, the joint fit, what the fit leaves, where to search, the next term tried on the twin, and the decision. At the end comes the best model with the parts of its bars.

Here are its {len(slides)} slides, model 0 and then model 1. On the second slide of model 1 the equation reads the width of the first model times one plus a times z, z the temperature in standard deviations from its mean, and the third shows what that term does: the fitted curves less the same model with a at zero, wider on the hot records and narrower on the cold ones.""")
    out.append("\n\n".join((f"![model {int(p.name[:2])}, phase {p.name.split('_')[1]}](examples/figures/log_{p.name})" if p.name != "99_best.png" else f"![the best model](examples/figures/log_{p.name})") + (f"\n\n*{CAPTIONS[p.name]}*" if p.name in CAPTIONS else "") for p in slides))
    if figures is not None:
        import shutil
        for old in figures.glob("log_*.png"):
            old.unlink()
        for p in slides:
            shutil.copyfile(p, figures / f"log_{p.name}")
    say("""## 7. A term not asked for by the records

The loop is not a machine that keeps whatever it is given. Here is a term that no answer of the report points to, a ripple along the abscissa:""")
    write("my_nmr/terms/02.py", 'from twinfit import mechanisms\nterm = mechanisms.periodic(3.0, "cosine")\nwhy = "a test: no answer of the report points to a periodic function"\n')
    run("python -m twinfit loop my_nmr")
    best = _json(folder, "my_nmr/loop/best.json")
    _claim(best["best"] == "model 1", "the loop refuses the ripple and the best model stays model 1")
    say("It refused the ripple, and says why in `best.md`:")
    show("my_nmr/loop/best.md", _part((folder / "my_nmr/loop/best.md").read_text(), "Tried and refused", "The loop stops"), "markdown")
    say(f"""A term is kept when three things hold: its fit ends inside its bounds, the Bayesian criterion falls by more than {TIED:g}, and the twin admits the model with it. The loop asks for no next term where no proposal beyond the twin gains, at first order, what a term of one parameter costs to be kept: {TIED:g} and the logarithm of the number of samples.

A refusal by the twin, or for too few of its searches that ended, read on the 30 synthetic sets of the default, can turn on more: `best.md` says so beside each, and `--sets 100` reads it again before the term is set aside, since a right term was refused on 30 sets and kept on 100.

Kept means that the records are better told with it than without it, and not that it is the term of the physics: on a model that lacks another term, a term of another shape can be kept, since it takes part of that misfit. Its reason, in `why`, is what makes it a term of physics. Such a term can move a shared parameter away from its truth: while the last report still has an answer beyond its twin, `best.md` says that no value is to be quoted yet.

To take a kept term out, copy the folder, leave out its file in the copy and rename those after it so that the numbers run on, and run the loop on the copy. A term refused behind a wrong one is tried alone the same way, in a copy whose `terms/` holds only it, as `01.py`.

## 8. The proposals for the experiment

Each report ends with three proposals for the experiment, each of one of three kinds: a better setup, a better characterization of the setup, or a better data collection.

Each comes from the description, which may name the changes that you could make, each a lever `[addenda.lever.NAME]` of one kind: noise, samples, span or repeats with a factor, covariate with new values of a setting, or characterization with an ansatz that it measures better. Every lever is measured, and those that do not rank among the three are said below them.

Where the description names fewer than three, the package adds changes that every experiment has: a quieter setup, more samples, more repeats, a wider span. It draws two more from the description itself: for each setting along which the curves vary and that no lever names, curves half its range past the values measured, and for each item of the setup known only as assumed or not recorded, its measurement. Each is measured on the twin: the twin is drawn again as the change would give it, and fitted, and the error of the known parameters, the shared ones of the model, is set against their error under the records as they are:""")
    show("my_nmr/loop/model_1.md", _part((folder / "my_nmr/loop/model_1.md").read_text(), "## Three proposals for the experiment", None), "markdown")
    say("""## 9. The researcher's own curves

This writes the starting files of a project into a new folder:""")
    run("python -m twinfit init my_project")
    say("Its last line gives the folder of the package and the Python that runs it: here they are written as `/the/folder/that/holds/twinfit` and `/the/python/that/runs/it`, and on your machine they are the paths of yours. "
        "Fill the places marked `TO FILL` in the description: the loop refuses a description whose `name` or `what` still reads `TO FILL`, since every report quotes them. "
        "A description with no model under `[addenda]` starts from the generic model of the package, one peak or dip of one width on its own height, place and baseline, and says so. Here is the description to fill:")
    show("my_project/description.toml", kind="toml")
    say("Write your first model in `model.py`, from this start, and its `EQUATION`, the model in LaTeX between dollar signs: the visual log typesets it at each turn, and shows the Python line of `model` where it is left empty:")
    show("my_project/model.py", kind="python")
    say("""Put each curve in `traces/`, one file a curve: a first line that names the columns as `columns` names them, then one line a sample, the abscissa and the value separated by a comma, with no line of comments, as in section 3.

**The noise of each sample**, in four steps:

1. Where you know the noise of a sample, give it, `from = "value"` with its standard deviation, or `from = "column"` with a column of each file.
2. Where you do not know it, leave `[noise]` out, or write `from = "neighbours"`, the default: the package reads the variance of the noise as a law in the height of the curve above its foot, a + b h + c h^2, a floor, a counter and a noise in proportion to the signal, from the second differences of neighboring samples before any fit, which a curve smooth over three samples barely reaches.

   Where a curve is sharper than that, write `from = "residuals"`: the same law, read from the second differences of what the first model leaves, and the model is fitted again with it. Either way a missing term smooth over three samples is not read as noise, on a peak as on a decay.

3. A law of your own, as one from the instrument, goes in a column of each file: the residuals of each model are in `loop/residuals_K.csv`, whose columns are the curve, the abscissa, the value, the fitted value, the residual and the noise that each sample was given.
4. Make each change of the noise in a copy of the folder, made without its `loop/` and `agent/`, as `rsync -a --exclude loop --exclude loop.before --exclude loop.reuse --exclude agent --exclude .twinfit.lock --exclude .twinfit.command my_project/ my_copy/` makes it, run the loop on the copy, and quote the copy's `best.md`: the first folder keeps its run as the record of what the noise changed.

The report leads with a kind where that kind is the first answer of what beyond the twin and its excess passes that of every answer of the mean beyond the twin: a misfit of the mean far larger leaves part of itself in the reading of the noise, so there the term comes first, and the noise is read again once the mean holds.

Where it leads with a noise correlated from sample to sample, as a filter of the detector leaves it, its page gives the number of samples over which the correlation falls: keep one sample in every that many, in a copy of the folder, and run the loop on several such copies, each from another first sample, since a term that one keeps and the others refuse is not to be believed.

`python -m twinfit thin my_project N` writes up to three such copies beside the folder, fewer where N is small and two starts fall on one sample, each named `DIR_thin_N_S` with S the first sample that it keeps. Quote the values of one copy with the half range of the three beside them: its bars, about the root of N times those of the whole folder, are the ones to quote, since the whole folder counts each sample as apart.

Where the report leads with glitches, samples far from the curve in one record alone, the page of the model names them, and no term file holds them: in a copy of the folder, give them a larger noise in a column of the curves, `from = "column"`, or leave them out of the files, and run the loop on the copy. In the column, each such sample takes a noise as large as its residual, the whitened residual that the page gives times the noise of the description, so that it weighs as one sample of noise, and every other sample the noise of the description.

**Curves that are not peaks.** A decay that starts at a time, a step or an edge is read as a peak is. `centre` gives the place about which the report reads, as the start of a decay: the core is the part near it, the wings the parts far from it, and the asymmetry the difference of its two sides.

The width in the names of the proposals is the width of the fitted curve at half its height above its median value, the unit of their distances. The kinds are shapes of the residuals and not causes: a tail on the upper side of a decay can be a second, slower decay, which no term of the package writes, since `second` is a copy of the curve displaced along the abscissa. Write it in `model.py`, its share and its lifetime as shared parameters, as the component of section 10, and a law of a setting, as the lifetime 1/(1 + K Q) of a quencher Q, there too.

**The loop.** Then run the loop, `python -m twinfit loop my_project`, read the report, write a term, and run it again. A term file that was kept may be written again later: each run fits every file afresh. A value that you hold from theory or from another measurement enters as an ansatz, `[addenda.ansatz.NAME]` with a line `value` and a line `uncertainty`, `NAME` a setting of every curve or a key of a dictionary `HELD` that `model.py` reads: its uncertainty is then the bar of theory of each parameter that it moves.

**A model that grows.** Start from a model small enough that you understand every part of it. The report reads what that model lacks, and each term that you add is one answer. Four things are written in `model.py` or in the description, and not in a term file. Each is a change of the folder: make it in a copy, made without `loop/` and `agent/`, so that the run before it stays as its record, as the guide of an agent asks too:

1. A dependence that the physics already gives, as a width that grows with the power, can be written in `model.py` itself, a shared parameter that reads `design["power"]`: it needs no term file, and its value is quoted as it is written.
2. An own value that moves with a setting, as the place of a line that the temperature shifts, is read by no question of the report while it is own to each curve, and `best.md` names it with its slope: write it as shared parameters. Here the place of the line of section 3, written to move along the temperature:

```python
import numpy as np

SHARED = ("width", "place_at_300", "shift_per_kelvin")   # the place at 300 K, and its change for each kelvin
OWN = ("height", "baseline")                              # of each spectrum
START_SHARED = [2.0, 0.0, 0.0]
CENTRE = "height"                                         # a name and no more: `centre` below gives the place
NAME = "one Lorentzian line of one width, its place along the temperature"


def place(shared, design):
    return shared[1] + shared[2] * (design["temperature"] - 300.0)


def model(x, shared, own, design):
    return own[0] / (1.0 + (2.0 * (x - place(shared, design)) / shared[0]) ** 2) + own[1]


def start_own(record):
    return [float(record.y.max() - record.y.min()), float(np.median(record.y))]


def centre(shared, own, design):
    return place(shared, design)
```

This is the whole file: it takes the place of the `model.py` of section 3, whose `BOUNDS_SHARED` held one bound for its one shared parameter, and here holds none, since a bound is given for each shared parameter or for none, and a start for each, in `START_SHARED` or by a function `start_shared(records)` that reads them from the curves.

`CENTRE` names one of `OWN`, any one, since the package asks that it be one of them: here it is a name and no more, and `centre` returns the place from the shared parameters, so that the report reads the core and the wings about it, and a pedestal or a periodic term, placed about it, can still be fitted. On the spectra of section 3 this model gives the width that the model of own places gives, and the shift for each kelvin with its bar.

3. A setting taken in the order of the time cannot be told from a drift of the instrument: write the order under `[problem]` as `order = "..."`, and `best.md` says it beside the values.
4. Where the curves sit far from zero against their span, the span from the least abscissa of the curves to the largest, their middle more than three spans from zero, as lines near 587.6 nm a few hundredths of a nanometer wide, write `abscissa_offset = 587.6` under `[problem]`, and a refusal of the loop names the value: the abscissa is read less it everywhere, in `model.py`, in `start_own`, in the term files, in every place that the output gives and in `residuals_K.csv`, since the fit moves the place of a curve by steps that grow with the size of the place.

   So the starts and the bounds of a place, in `START_SHARED` and `BOUNDS_SHARED`, are written less it too, and `best.md` gives each place less it too, with a line that says so: add the offset back to quote a place.

The keys of the description beyond those of section 3, each with its form. The noise of each sample is given one of five ways. The output is for people, for agents or both, and the search of the literature is none, local, online or both, run by an agent of the model and the `effort` given.

An ansatz is a value held before the fit, with its uncertainty. A lever is a change of the experiment that the proposals weigh: `noise`, `samples`, `span` or `repeats` with a factor, `covariate` with the setting and its new values, `characterization` with the ansatz that it measures better and a factor on its uncertainty, or `setup` with `of`, the name of an item of `[setup]` that it changes, and a factor, a half where none is given. A lever of kind noise or setup is proposed as a better setup, one of kind characterization as a better characterization of the setup, and the other four as a better data collection:""")
    out.append(f"```toml\n{REFERENCE.rstrip()}\n```")
    _claim(_reads_as_a_description(folder), "the keys above are read by the package without a refusal")
    say("""An ansatz on a setting holds it at one value, which every curve gives: a setting that changes from curve to curve is fitted along, and not held. An ansatz on a key of `HELD` sets that key before each fit, and `model.py` reads it where it computes the curve:""")
    out.append("```python\nHELD = {\"resolution\": 0.5}          # the values that the model reads and does not fit, which an ansatz sets\n\n\ndef model(x, shared, own, design):\n    width = (shared[0] ** 2 + HELD[\"resolution\"] ** 2) ** 0.5\n    return own[0] / (1.0 + ((x - own[1]) / width) ** 2) + own[2]\n```")
    say("""Its uncertainty is then the bar of theory of each parameter that it moves, and a lever of kind characterization of it says what measuring it better would buy. A number written in a term file, as the distance of a second component, is held by nothing: where a measurement could fix it, move it into `model.py` as a key of `HELD`, with its ansatz.

## 10. A term of unknown distance or period

Some terms carry a number that is chosen and not fitted: the distance of a second component, the stretch of a tail, the period of a ripple. The worked project of chromatography has one:""")
    run("python -m twinfit example chromatography my_chrom")
    run("python -m twinfit loop my_chrom")
    c0 = _json(folder, "my_chrom/loop/model_0.json")["where"]["first_of_each_question"]
    _claim(c0["where"] == "the asymmetry" and c0["what"] is None, "the report of the chromatograms puts the misfit in their asymmetry, and no law of the noise")
    say("The report puts the misfit of the chromatograms in their asymmetry:")
    show("my_chrom/loop/model_0.md", _part((folder / "my_chrom/loop/model_0.md").read_text(), "## Where to search", "Every answer"), "markdown")
    chrom0 = _json(folder, "my_chrom/loop/model_0.json")["proposals"]["for_the_model"][0]
    _claim(chrom0["stands_for"] == "the asymmetry" and chrom0.get("as_a_term_file") is None, "the first proposal for the chromatograms is the asymmetry, and names no term")
    fit0 = _json(folder, "my_chrom/loop/model_0.json")
    chi, dof = fit0["chi_square"], fit0["degrees_of_freedom"]
    _claim(abs(chi / dof - 1.0) < 3.0 * (2.0 / dof) ** 0.5, "the chi square of the first model of the chromatograms agrees with the noise within three of its standard deviations")
    say(f"""The chi square of this first model, {chi:.1f} on {dof} degrees of freedom, agrees with the noise, and the report finds the second component all the same: a term that holds a small part of the curve moves the whole chi square by less than its own scatter, while a statistic aimed at one shape reads it. That is what the report is for.""")
    say("""In a chromatogram an asymmetry of the peaks is a compound that elutes close to the first, an impurity, or a tail of the column. The first proposal for the model is the asymmetry, measured with the shape that the residuals hold there:""")
    show("my_chrom/loop/model_0.md", _part((folder / "my_chrom/loop/model_0.md").read_text(), "## Three proposals for the model", "## Three proposals for the experiment"), "markdown")
    unit = _json(folder, "my_chrom/loop/model_0.json")["proposals"]["width_of_the_curves"]
    sigma = _json(folder, "my_chrom/loop/best.json")["parameters"]["width"]["value"]
    say(f"""Try the impurity first, a second component of the curve, `second(distance)`, with its distance in the unit of the abscissa. In the page of the proposals a width is the full width of the fitted curves at half their height above their median value, {unit:.2g} s here, while the parameter `width` of this model is the standard deviation of its Gaussian, {sigma:.2g} s. The distance is not read off the report, which names no term.

`scan` fits the term at each distance that you give, beside the model, and says which lowers the Bayesian criterion most. The distances are in the unit of the abscissa, here seconds, and `x` stands for the number that is scanned: `"second(x)"` is `second(distance)` at each value given. It chooses the stretch of `tail` or the period of `periodic` the same way, with the other argument in quotes, `"tail(x, 'upper')"` or `"periodic(x, 'sine')"`.

It refuses the term at a number that a kept term file of the folder already holds, and fits a second term of the same kind beside it at any other number: the files that the loop refused are not in the model that it scans beside. To scan the number of a term that a file holds, scan in a copy of the folder without that file.

Where the range that it gives is under the step of the values scanned, a parabola through the three falls about the largest gives it, and the top of that parabola is the value to write. That range is the bar of the chosen number where the noise is white: under a noise correlated over N samples it is about the root of N too narrow, so read it after the noise is thinned.

A distance has a side, and a term on the wrong side can still lower the criterion with a negative size, which takes its shape away from the curve: scan both signs, as -2 -1 1 2, keep the side of the larger fall and of a positive size, and read a kept term of negative size as the command line and `best.md` flag it. `displaced` points the other way: its long side is toward smaller abscissa where its reach is positive, as `python -m twinfit terms` says.""")
    run('python -m twinfit scan my_chrom "second(x)" 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.4')
    from .loop import scan
    scanned = scan(folder / "my_chrom", "second(x)", [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4])
    picked, (lo, hi) = scanned["best"]["x"], scanned["within_one"]
    near = [r["x"] for r in scanned["rows"] if r["how"] != "ended"]
    if near:
        _claim(near == sorted(r["x"] for r in scanned["rows"])[:len(near)], "the fits that did not end are those of the nearest distances")
        say(f"""At {' and '.join(f'{x:g}' for x in near)} s, the nearest distances, the second component lies so close to the peak that the fit cannot tell its size from the shape of the peak, and its search does not end. A fit that does not end is no answer, and the scan reads the others.""")
    say(f"""The distances from {lo:g} to {hi:g} s lower the criterion within 1 of the largest: the records do not tell them apart by the criterion, and the loop, with its twin, decides. Write the term at the distance that the scan gives, and run the loop:""")
    write("my_chrom/terms/01.py", f'from twinfit import mechanisms\nterm = mechanisms.second({picked:g})\nwhy = "the report: a second component in the asymmetry. An impurity that elutes later, at the distance that the scan gives"\n')
    run("python -m twinfit loop my_chrom")
    cb = _json(folder, "my_chrom/loop/best.json")
    _claim(cb["best"] == "model 1", "the loop keeps the second component")
    ct = fields.TRUTH["chromatography"]
    size = cb["parameters"][f"second component at {picked:g}"]
    after = ct["after"] * ct["width"]
    _claim(lo <= after <= hi, "the distance that drew the curves lies in the range that the scan gives")
    at = {r["x"]: r for r in scanned["rows"]}
    say(f"""The curves were drawn with an impurity {after:g} s after the peak, of {ct['impurity']:g} of its height. The scan gives the distance from {lo:g} to {hi:g} s, and {after:g} s is in that range. The size moves with the distance across it, from {at[lo]['size']:.2g} at {lo:g} s to {at[hi]['size']:.2g} at {hi:g} s: these records fix the distance and the size together and not each.

At {picked:g} s the loop kept the term with a size of {' ± '.join(_pm(size['value'], size['bar_of_the_fit']))}, and the width of the peak, {' ± '.join(_pm(cb['parameters']['width']['value'], cb['parameters']['width']['bar_of_the_fit']))} s with the term, was {sigma:.4g} s without it, where {ct['width']:g} s drew the curves. The width moves across the range of the distance too, as the last line of the scan says: half of its move is a part of its uncertainty that no column of `best.md` holds, since the distance is written and not fitted, and it is quoted beside them.

A measurement of the retention time of the impurity alone, a better characterization of the setup, would fix the distance, and the size and the width with it.

To carry the distance into the uncertainty of every other value, hold it: a key of `HELD` in `model.py`, an ansatz with half the range that the scan gives as its uncertainty, and a term that reads it at each call. That range comes from these records, so this ansatz narrows nothing: it carries the uncertainty of a chosen number into the other values, where an ansatz that narrows a value must come from outside the records.

A term of the package takes its number when its file is read, so `held.reading` builds it again at each call, from the value that `HELD` holds then. The new `01.py` takes the place of the one above, and it reads `HELD` by its name with no import: the loop gives every term file the names of `model.py`, its functions and its constants as well as `HELD`, so that a term reads them with no import:""")
    (folder / "my_chrom" / "model.py").write_text((folder / "my_chrom" / "model.py").read_text() + f"\n\nHELD = {{\"distance\": {picked:g}}}          # the distance of the impurity, held by the ansatz of the description\n")
    (folder / "my_chrom" / "description.toml").write_text((folder / "my_chrom" / "description.toml").read_text() + f"\n[addenda.ansatz.distance]\nvalue = {picked:g}\nuncertainty = {(hi - lo) / 2.0:g}\n")
    out.append(f"`my_chrom/model.py`, at its end:\n\n```python\nHELD = {{\"distance\": {picked:g}}}          # the distance of the impurity, held by the ansatz of the description\n```\n\n`my_chrom/description.toml`, at its end:\n\n```toml\n[addenda.ansatz.distance]\nvalue = {picked:g}\nuncertainty = {(hi - lo) / 2.0:g}\n```")
    write("my_chrom/terms/01.py", 'from twinfit import held, mechanisms\nterm = held.reading(lambda: mechanisms.second(HELD["distance"]))\nwhy = "the report: a second component in the asymmetry, at the distance that the scan gives and that the ansatz holds"\n')
    run("python -m twinfit loop my_chrom")
    hc = _json(folder, "my_chrom/loop/best.json")
    _claim(hc["parameters"]["width"]["bar_of_theory"]["together"] > 0.0 and "distance" in hc["parameters"]["width"]["bar_of_theory"]["of_each"], "the held distance puts a bar of theory on the width")
    show("my_chrom/loop/best.md", _part((folder / "my_chrom/loop/best.md").read_text(), "## To quote", "The bias on the twin here"), "markdown")
    say("""The bar of theory of each value is now the move that the distance makes across its uncertainty, and it is quoted beside the bar of the fit. The column says of theory whatever the source of the ansatz: here it is the range of a number that these records chose.

Where the records fix the distance well, as the coupling of a doublet whose two lines stand apart, a third way gives it a bar of its own: write the second line in `model.py`, with the distance as a shared parameter, and the fit gives the distance its bar and the twin its coverage, as every other shared value.

The term `second` places its copy from the place of the curve, so the place stays that of the first line: the middle of an equal doublet is the place plus half the distance, and the report reads the core and the wings about the first line. The size of a second component, a pedestal or a tail is its area over that of the curve, and never a distance.

Where the curves are doublets, start every curve on the same line: the `start_own` of the template takes the first sample from the left that reaches four fifths of the largest distance from the baseline, so that each curve starts on its lower line, and a page whose places differ by more than a width while the fit leaves much says that the starts may sit on different lines.

A component of its own width at its own distance, as the outflow under the emission line of a galaxy, is no term of the package, which adds one parameter a term: where the report puts a tail and a second component tied, write the whole line in `model.py`, its share, its distance and its width as shared parameters:

```python
SHARED = ("width", "share", "distance", "width_of_the_second")   # the share is the area of the second over that of the first
OWN = ("height", "place", "baseline")
START_SHARED = [50.0, 0.1, -50.0, 100.0]                          # in the unit of the abscissa, here km/s
BOUNDS_SHARED = ([1.0, 0.0, -500.0, 1.0], [500.0, 2.0, 500.0, 1000.0])


def model(x, shared, own, design):
    width, share, distance, second = shared
    first = np.exp(-0.5 * ((x - own[1]) / width) ** 2)
    other = share * (width / second) * np.exp(-0.5 * ((x - own[1] - distance) / second) ** 2)
    return own[0] * (first + other) + own[2]
```

with `start_own`, `CENTRE` and `centre` as in section 3. On spectra whose second component was two and a half times wider than the line, 120 km/s below it with a fifth of its area, this model gave each of the four values back within its bar of the fit.

Where the second component has the width of the first, as the impurity of the chromatograms above, its own width cannot be told from its share and its distance, and the term `second` is the one to keep. Any function of the settings is written the same way, as a height that falls as the exponential of a radius or a place that follows the tanh of it, and not only a line.""")
    say("""## 11. A term written by the researcher

When the term that the report points to is not among those of `python -m twinfit terms`, write it, in five steps.

1. **The term.** A term is a `Mechanism(name, kind, parameter, start, wrap, lower, upper, formula=...)`: a name, a word for its kind, any word, as `baseline`, which the package keeps with the term and reads in no result, the name of the one parameter that it adds, where the search of that parameter starts, `wrap`, the bounds of the parameter, none where `lower` and `upper` are not given, and `formula`, the term in LaTeX between dollar signs with `f` the model before it and `a` its parameter, as `r"$f(x) + a\\,x$"`, which the visual log typesets: without it the log says that its formula is not written. `python -m twinfit card` prints a term file of each kind in one page.
2. **Its wrap.** `wrap(model, problem, n)` takes the model as it stands, the problem, which holds every curve, and the number of shared parameters before the term, and returns the model with the term, a function `new(x, shared, own, design)`, where `shared[n]` is the new parameter and `shared[:n]` the parameters before it. The loop refuses a `wrap` that does not take these three arguments. `design` holds the settings of the curve, as the description names them, `own[k]` and `shared[k]` the parameters of the model in the order of `OWN` and `SHARED` of `model.py`, and `x`, `shared` and `own` are arrays of numpy.
3. **One of three forms.** `new` may add a term to the model, as a step under a peak whose size is a share of its height, `model(x, shared[:n], own, design) + shared[n] * own[0] * 0.5 * (1.0 + erf((x - own[1]) / (2.0 ** 0.5 * shared[0])))`, with `erf` of `scipy.special`. It may multiply the model, as a gain of the detector that changes along the abscissa, `model(x, shared[:n], own, design) * (1.0 + shared[n] * (x - middle))`, which multiplies the terms kept before it too, about a middle that you choose as the zero of the change.

   Or it may replace the line and return a curve of its own, with no call of `model`: then it returns the whole curve, its baseline too, since what it does not add again is gone. The loop refuses such a term where its curve no longer reads an own parameter that the model read, and says which.

4. **A curve that replaces the line**, written out whole, as `01.py` in a copy of the folder of section 10 that holds no other term file: a Gaussian with an exponential tail, for the chromatograms of section 10, whose model holds the width `s = shared[0]`, the height `own[0]`, the place `own[1]` and the baseline `own[2]`:

```python
import numpy as np
from scipy.special import erfc, erfcx
from twinfit.mechanisms import Mechanism


def wrap(model, problem, n):
    def new(x, shared, own, design):
        s, t = shared[0], shared[n]                    # the width and the time constant of the tail
        d = x - own[1]                                 # the distance from the place: own[1] - x turns the tail toward smaller abscissa
        z = (s / t - d / s) / np.sqrt(2.0)
        tail = own[0] * (s / t) * np.sqrt(np.pi / 2.0) * np.where(
            z > 0, np.exp(-0.5 * (d / s) ** 2) * erfcx(np.maximum(z, 0.0)), np.exp(0.5 * (s / t) ** 2 - d / t) * erfc(np.minimum(z, 0.0)))
        return tail + own[2]                           # the line is replaced, so its baseline is added again
    return new


term = Mechanism("an exponential tail of the line", "tail", "time constant", 1.0, wrap, 1e-3, 100.0)
why = "the report: the asymmetry of the chromatograms, and a detector that answers with a time constant"
```

   `own[0]` is the height that the curve keeps as `t` falls to zero, and its area is that of the Gaussian, so a tail that holds a share `a` of a curve is `(1 - a)` times the Gaussian and `a` times this curve, with `a` a shared parameter. With `erfcx` and `erfc` it gives a number at every distance: the product of an exponential and `erfcx` alone gives none past about 38 widths. Where the time constant follows a setting of the curves named load, write `t = shared[n] * design["load"]`.

5. **Its start and its place among the files.** Start the new parameter at a size that the residuals give, and a parameter that the curve divides by, as `t`, near the width of the curve and away from zero. A term reads the shared parameters as the terms before it leave them: a tail in `02.py` after `along("width", ...)` in `01.py` reads the width that `along` makes at each setting. A term that replaces the line drops the terms kept before it, so it goes first, as `01.py`.

   A term must be smooth in its parameters, since one whose curve jumps as a parameter moves keeps the searches from ending, and a term whose curve holds a value that is not a number is refused with its file named: compute its `wrap` on the curves alone, with the starts of the model, before the loop. The loop reads the term files in their order up to the first that is missing, so `03.py` is read only beside a `02.py`, and it reads the `why` of each file before its first fit.

Here is a slope of the baseline that grows with the temperature, written in the project of section 3:""")
    write("my_nmr/terms/03.py", '''import numpy as np
from twinfit.mechanisms import Mechanism


def wrap(model, problem, n):
    def new(x, shared, own, design):
        middle = 0.5 * (x[0] + x[-1])
        slope = shared[n] * (design["temperature"] - 305.0)
        return model(x, shared[:n], own, design) + slope * (x - middle)
    return new


term = Mechanism("a slope of the baseline that grows with the temperature", "baseline", "slope per kelvin", 0.0, wrap)
why = "a test of a term of my own: a drift of the baseline that a warmer probe could leave"
''')
    said = run("python -m twinfit loop my_nmr")
    _claim(_json(folder, "my_nmr/loop/best.json")["best"] == "model 1" and "a slope of the baseline that grows with the temperature" in said.split("Refused:", 1)[1].split(". The best", 1)[0],
           "the loop refuses the slope of the baseline, which the curves do not hold")
    say("""The curves of this project were drawn with a flat baseline, so the loop refuses the slope as it refused the ripple. A refused term stays where it is: the loop reads `01.py`, `02.py` and on, and stops at the first number that is missing, so a file taken out would hide every file after it. Each run tries a refused term again and refuses it again, with its reason, in `best.md`.

The terms of the package are made of copies of the curve. A shape that the physics names otherwise, as the exponential tail of a detector, or a ripple fixed to the abscissa, the same on every curve whatever its center, as the fringe of an etalon fixed to the abscissa, is a term of your own: the report gives its place and its kind, and you give its form.

`periodic` is placed about the center of each curve, and where the curves have centers apart, a fringe fixed to the abscissa read by it leaves a misfit that changes from curve to curve, which the report then reads along a setting. Here is that fringe, with its period as a chosen number. A term file of your own that defines `term_at(x)`, the term at the number x, is scanned by naming the file, and a term file reads the dictionary `HELD` of `model.py` by its name, so that a number held by an ansatz reaches the term:""")
    write("my_nmr/terms/04.py", '''import numpy as np
from twinfit.mechanisms import Mechanism


def term_at(period):
    def wrap(model, problem, n):
        def new(x, shared, own, design):
            return model(x, shared[:n], own, design) + shared[n] * np.cos(2.0 * np.pi * x / period)
        return new
    return Mechanism(f"a fringe fixed to the abscissa, of period {period:g}", "periodic", f"fringe of period {period:g}", 0.0, wrap)


term = term_at(3.0)
why = "a test of a term of my own with a chosen number: a fringe that the optics could leave"
''')
    fringe = run("python -m twinfit scan my_nmr terms/04.py 2 3 4")
    _claim("No value of x lowers the Bayesian criterion by more than" in fringe, "the scan of the fringe finds no period that the curves ask for")
    say("""No period lowers the criterion: these curves hold no fringe of this phase. A ripple of unknown phase is a cosine and a sine, two terms, and a scan of the cosine alone can find no period where the sine would: scan each. The scan reads the model of the other kept term files, and leaves the scanned file out of it.""")
    say("""## 12. With an agent

`init` also wrote `CLAUDE.md` and a skill for Claude Code. Open a fresh Claude Code chat in the project folder and ask it, for example: "Run the loop, read the report, and propose the next term." Where you will not be there to choose at step 5, add: "choose on your own at step 5, and say why." It reads `CLAUDE.md` first. At each iteration it follows ten steps, which `python -m twinfit guide` prints with the rest of its protocol, its section 0 and its part What you give back among them, which the steps name:""")
    agent = (HERE / "agent" / "AGENT.md").read_text()
    out.append("```\n" + _part(agent, "**Each iteration**", "**What you give back") + "\n```")
    say("""At step 5 you choose among its three proposals, or let it choose. Every term that it adds carries its reason in `why`: the physics that could leave the structure, offered as a candidate, which the loop and its twin then bear out or refuse. The package names no cause, and the agent names one only as a candidate.

## 13. More examples, of other fields

`python -m twinfit guide examples` prints the gallery, `twinfit/EXAMPLES.md`: projects of other fields, each written by `python -m twinfit example NAME DIR` with a truth that its folder does not hold: a pedestal under the lines, a satellite line, a width that saturates with the power, a noise correlated by a filter, a Fano profile, photon counts, an atomic clock whose center moves with the power of its probe, a band of electrophoresis with a tail, the transit of an exoplanet over a star darker at its limb, the fringes of an interferometer whose scan is bent, and the daily returns of a stock index with fat tails.

- For each it gives what the first report put first, what to do and what not to do, and what the loop then gave beside the truth.
- Find the one whose first report leads as yours does.
- Each is shown too at three times its noise, with none drawn, and with a second kind of noise, correlated over neighboring samples or that of a counter: `python -m twinfit example clock DIR --noise 3 --kind correlated` writes one so, to see how the first report and the outcome move with the noise of an example, and the page says what the noise changes in its tables.

## 14. Limits of the package

- It names a place, a kind and a setting of what is missing, never its cause.
- The twin draws its synthetic curves at the noise that the description gives, so it takes that noise as right. The check of the residuals against the height of the fitted curve, on the first page of each report, is what says whether it is.
- With 30 synthetic sets, the default, an answer can stand beyond all of them and no further: its chance is then the smallest that 30 sets give. More sets, `python -m twinfit loop DIR --sets 100`, read smaller levels and take longer.
- A term much smaller than the noise can be missing while no answer stands beyond the twin. When one answer alone stands just past its line, the loop says so and stops: a model that lacks nothing shows one such answer on some of its synthetic sets.
- Glitches and a noise that grows with the signal are kinds of the report that no term of the mean fills. The report says so when one leads, and how to give the noise of each sample instead.
- A proposal for the experiment is read on the synthetic sets of the twin, so its ratio moves from run to run: the page says within how much, about one over the root of the number of sets.
- A term in the wrong frame leaves a misfit of its own. A ripple read about each center where it is fixed to the abscissa changes from curve to curve, and the report then points along a setting: read a new answer along a setting against the frame of the terms kept before it.
- Two forms of one effect, as a Lorentzian added to a line and a Voigt line in its place, are not set side by side by one run of the loop, which fits each term file on top of the ones before it: compare them in two copies of the folder, each with one form, by their criterion and their twin.
- One misfit can light several answers at once: a slope of the baseline lights a slow trend over the span of each curve and a periodic structure along the abscissa together. The report says when kinds tie, and then their order means little: try each, and keep the one that the twin and the known parameters bear out.
- The bars of the own parameters of each curve count what the shared ones leave uncertain, but no twin reads them, as it reads the shared ones.
- `along` multiplies a parameter by one plus its size times the setting, so it cannot move a parameter that sits near zero, as a place read less the offset of its abscissa can. A place that moves with a setting is written in `model.py`, as the model of section 9 writes the shift for each kelvin, and as the gallery's clock writes the shift of its line along the power.
- `second` takes one distance and copies the curve there, at its own width: a second line whose distance or share follows a setting is written in `model.py`, as a function of the settings, as a place is.
- A lever along a setting reads the model at settings past those of the records, and holds as far as the model holds there. A lever of samples takes them as independent: under a noise correlated over N samples, twice the samples is less than twice the information.
- The loop can end in words that the runs above do not show. Where a noise correlated from sample to sample leads the last report, a term file cannot hold it, and the report of that model says what to do in a copy of the folder. Where the check of the noise of the best model fails, what was kept is not to be believed before that is read. Where the twin does not admit the fit of the best model, its bars are not to be read, and where a statistic misses its band by little, `--sets 500` reads it again.""")
    page = "\n\n".join(out) + "\n"
    return page.replace("{SECONDS}", f"ran in {time.time() - began:.0f} seconds")


REFERENCE = """[noise]                          # one of five ways
from = "value"                   # one standard deviation for every sample
value = 0.02
# from = "column"                # or a column of the curves that gives the standard deviation of each sample
# column = "sd"
# from = "neighbours"            # or the default, a law in the height of the curve, read from the differences of neighboring samples before any fit
# from = "residuals"             # or the same law, read from what the first model leaves, for a curve sharper than three samples
# from = "model"                 # or a law of yours: variance(mean, record, own) in model.py, the variance of each sample

[output]                         # who reads the output of the loop: both by default
for = ["people", "agents"]       # people: the visual log and model_K.md, agents: brief_K.md, one page a model

[literature]                     # a search of the literature by an agent, which `python -m twinfit literature DIR` briefs
search = "none"                  # none, local, online or both
folder = "papers"                # the papers held, for a local search
model = "sonnet"                 # the model of the agent that searches, and its effort
effort = "high"

[addenda.ansatz.resolution]      # a key of the dictionary HELD of model.py, or a setting that every curve gives at this value
value = 0.5
uncertainty = 0.05

[addenda.lever.quieter]          # noise, samples, span or repeats: with a factor
kind = "noise"
factor = 0.5
words = "a detector of half the noise"

[addenda.lever.hotter]           # covariate: new values of one setting
kind = "covariate"
covariate = "temperature"
values = [330.0, 340.0]
words = "two spectra at higher temperatures"

[addenda.lever.resolution]       # characterization: an ansatz measured better, its uncertainty times the factor
kind = "characterization"
of = "resolution"
factor = 0.5
words = "the resolution measured before the run"
"""


def _reads_as_a_description(folder: Path) -> bool:
    """Whether the keys of REFERENCE, added to the description of the worked project, are read by the package without a
    refusal: the reference is shown, and a key that the package refuses would send a researcher to a refusal (round 3 of
    the cold tests asked for the keys, which only refusals had taught)."""
    from .description import read
    text = (folder / "my_nmr" / "description.toml").read_text()
    trial = folder / "my_nmr" / "reference_trial.toml"
    trial.write_text(re.sub(r"(?ms)^\[addenda\.lever\..*", "", re.sub(r"(?ms)^\[noise\][^\[]*", "", text)) + "\n" + REFERENCE)
    try:
        d = read(trial)
        return {lever["name"] for lever in d.levers()} >= {"quieter", "hotter", "resolution"} and d.noise["from"] == "value"
    finally:
        trial.unlink()


def read_as(text: str) -> str:
    """The digest of a page of the tutorial as a tester reads it: all but the seconds that its run took. A round of the
    cold tests names the digest that its testers read, and a page that changed since is a page that no tester read."""
    import hashlib
    return hashlib.sha256(SECONDS.sub("", text).encode()).hexdigest()


WORDS = (          # a word that readers of the tutorial met before its definition: the pattern of its use, and the words that define it
    ("an ansatz", r"\bansatz", "**An ansatz**"), ("a lever", r"\blevers?\b", "**A lever**"), ("a bar", r"\bbars?\b", "**A bar**"),
    ("the five parts of the uncertainty", r"\bfive parts\b", "**The five parts of the uncertainty**"), ("the coverage", r"\bcoverage\b", "**The coverage and the bias on the twin**"),
    ("a scatter", r"\bscatters?\b", "**A scatter**"), ("the excess", r"\bexcess\b", "**The excess**"), ("a tie", r"\bti(?:e|es|ed)\b", "**Two kinds tie**"),
    ("the gain in chi square", r"\bgain in chi square\b", "**The gain in chi square**"), ("a profile", r"\bprofile\b", "**A profile**"),
    ("a Mechanism", r"\bMechanisms?\b", "**A Mechanism**"), ("a rung", r"\brungs?\b", "**A rung**"), ("a kernel", r"\bkernels?\b", "**A kernel**"),
    ("the Bayesian criterion", r"\bBayesian criterion\b", "**The Bayesian criterion**"))


def before_defined(text: str, pointer: str | None = None) -> list:
    """The words of WORDS that a page uses before the sentence that defines them, each said with the lines of its use
    and of its definition. A page that points to the definitions, by `pointer`, before its first use of any of them
    passes for them all, as the guide that points to section 1 of the tutorial. A word that the page never defines and
    never points to is said too. Testers met these words before their definitions round after round."""
    import re
    line = lambda at: text.count("\n", 0, at) + 1          # noqa: E731
    uses = {name: re.search(use, text, 0 if name == "a Mechanism" else re.I) for name, use, _ in WORDS}
    used = [m.start() for m in uses.values() if m]
    if pointer is not None and pointer in text and used and text.index(pointer) < min(used):
        return []
    found = []
    for name, _, mark in WORDS:
        first, at = uses[name], text.find(mark)
        if first is None:
            continue
        if at < 0:
            found.append(f"{name}: used on line {line(first.start())} and never defined")
        elif first.start() < at:
            found.append(f"{name}: used on line {line(first.start())} before its definition on line {line(at)}")
    return found


AGENTS = HERE / "agent" / "TUTORIAL.md"          # the tutorial for an agent, written from this one by `for_agents`

RULES = """The rules that decide each turn:

- The output of the loop for an agent is `DIR/loop/brief_K.md`, one page a model K, with `[output] for` holding "agents",
  as it does by default. `model_K.md` is the page for people and `model_K.json` holds every number.
- A term is asked only where an answer of the report stands beyond the twin: a right model passes the line of an answer
  on 5 sets of 100, so one answer just past it is weak evidence.
- The loop keeps a fitted term when its fit ends inside its bounds, the Bayesian criterion falls by more than 2 and the
  twin admits the larger model, and best.md gives the reason of a refusal.
- A term that the physics asks and the records cannot show is held at its theory, with hold = True in its file.
- Quote a value only from a model that its twin admits, with the parts of its bar that best.md gives.
- `python -m twinfit wait DIR` waits for a running loop within ten minutes and prints the last brief: never poll a run
  in another way, and never start two in one folder."""


def for_agents(people: str) -> str:
    """The tutorial for an agent, from the tutorial for people (his message 316 of the repository it is made in: "one
    only intended for humans and one only intended for AI agents, so optimised and efficient for them, without visual
    explanation unless necessary"). The facts are one source: each numbered section keeps its heading, its first
    sentence, every command and every file to write, as the tutorial for people gives them. What it leaves out is what a
    person needs and an agent does not: the pictures, the outputs that the commands print, which an agent reads when it
    runs them, and the reasons. The rules that decide each turn come first."""
    out = ["# twinfit for an agent: the tutorial in commands", "",
           "Written by `python -m twinfit.tutorial --write` from TUTORIAL.md, the tutorial for people, which holds the "
           "reasons, the outputs and the pictures: every command and every file here is taken from it.", "",
           "It builds on the card, `agent/AGENT.md`, and sets out every command and every file of the two worked projects. Not covered here: the reasons, the outputs and the pictures, which `TUTORIAL.md` gives.", "",
           "> Section 1 of `TUTORIAL.md` defines every word of the package that this page uses.", "", RULES, ""]
    lines = people.splitlines()
    i, first = 0, False
    while i < len(lines):
        line = lines[i]
        if re.match(r"## \d+\. ", line):
            out += [line, ""]
            first = True
            i += 1
            continue
        if line.startswith("```"):
            lang = line[3:].strip()
            j = i + 1
            while j < len(lines) and not lines[j].startswith("```"):
                j += 1
            if lang in ("bash", "python", "toml") and any(o.startswith("## ") for o in out[5:]):
                label = next((lines[k] for k in range(i - 1, max(i - 3, -1), -1) if lines[k].strip()), "")
                if lang != "bash" and re.fullmatch(r"`[^`]+`:", label.strip()):
                    out += [label, ""]
                out += lines[i:j + 1] + [""]
            i = j + 1
            continue
        if first and line.strip() and not line.startswith(("#", "!", "|", "`", "-", "*", ">")) and any(o.startswith("## ") for o in out[5:]):
            j = i
            while j < len(lines) and lines[j].strip() and not lines[j].startswith("```"):
                j += 1
            paragraph = " ".join(x.strip() for x in lines[i:j])
            if len(paragraph.split()) > 60:          # a long paragraph gives its first sentence, and the next where the first ends on a colon
                parts = re.split(r"(?<=[.:])\s", paragraph)
                paragraph = " ".join(parts[:2] if parts[0].endswith(":") else parts[:1])
            if paragraph.endswith(":"):          # a colon that introduced an output, which this page drops, would introduce the next command: "The report opens with the words that it uses:" stood over `python -m twinfit matrix` (the blind red team of 0.1.0rc13)
                block = next((x[3:].strip() for x in lines[j:] if x.startswith("```")), "")
                if block not in ("bash", "python", "toml"):
                    paragraph = paragraph[:-1] + "."
            paragraph = paragraph.replace("write it, in five steps.", "write it: `python -m twinfit card` gives its form, and section 11 of `TUTORIAL.md` its five steps. The file below is `terms/01.py` of a copy of the folder of section 10.")          # it promised five steps and showed none (the red team of his message 437)
            out += [paragraph, ""]
            first = False
            i = j
            continue
        i += 1
    return "\n".join(out).rstrip() + "\n"


def main(argv: list) -> int:
    if argv[1:2] == ["--check-words"]:          # the words of the tutorial and of the guide, each defined before its use: seconds, where --check takes minutes
        here = Path(__file__).resolve().parent
        found = before_defined((here / "TUTORIAL.md").read_text()) + ["the guide: " + f for f in before_defined((here / "agent" / "AGENT.md").read_text(), "Its section 1 defines the words that this page uses")]
        print("\n".join(found) if found else "the tutorial and the guide: every word of the list is defined before its use")
        return 1 if found else 0
    if argv[1:] not in (["--check"], ["--write"]):
        print("python -m twinfit.tutorial --check    runs every command of TUTORIAL.md again, in a temporary folder, and says whether any output differs\n"
              "python -m twinfit.tutorial --write    runs them and writes TUTORIAL.md again, for whoever changes the package\n"
              "python -m twinfit.tutorial --check-words    says whether the tutorial and the guide use a word of the list before its definition")
        return 1
    with tempfile.TemporaryDirectory() as tmp:
        page = build(Path(tmp), HERE / "examples" / "figures" if argv[1:] == ["--write"] else None)
    if argv[1:] == ["--check"]:
        if not PAGE.exists():
            print(f"REFUSED: {PAGE} is not written: run python -m twinfit.tutorial --write")
            return 1
        if SECONDS.sub("", page) != SECONDS.sub("", PAGE.read_text()):
            print(f"REFUSED: the outputs of the tutorial differ from {PAGE}: run python -m twinfit.tutorial --write and read what changed")
            return 1
        if not AGENTS.exists() or AGENTS.read_text() != for_agents(PAGE.read_text()):          # the tutorial for an agent is that of people in its form (his message 316)
            print(f"REFUSED: {AGENTS} is not what the tutorial for people gives: run python -m twinfit.tutorial --write")
            return 1
        words = main(["", "--check-words"])          # a word used before its definition fails the check too
        if words:
            return words
        print(f"the tutorial: every output of {PAGE.name} is what its command prints now")
        return 0
    PAGE.write_text(page)
    AGENTS.write_text(for_agents(page))
    print(f"written: {PAGE} and {AGENTS}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
