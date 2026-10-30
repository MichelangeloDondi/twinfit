# twinfit

Where a model lacks a term: the structure of the residuals of a fit, read against its twin.

In five minutes, the clock below. In an hour, sections 1 to 6 of `twinfit/TUTORIAL.md`. Whole, that tutorial, `twinfit/EXAMPLES.md` and, for an agent, the card of `twinfit/agent/AGENT.md`.

This page builds on a fit of curves that a researcher already makes, and sets out what the package reads in its residuals, how to install it and a first run. Not covered here: the method step by step, which the tutorial gives, and the projects of other fields, which `twinfit/EXAMPLES.md` gives.

> Section 1 of `twinfit/TUTORIAL.md` defines every word of the package that this page uses.

## Purpose

A model fitted to measured curves can pass the usual test of fit and still miss an effect, a broadening, a background or a drift, which biases the quoted parameters while the chi square stays near one. twinfit finds the structure that the effect leaves, says what it looks like and measures what it does to the parameters. It names no term of a field. What it was measured to do, under its limits below, was measured on synthetic curves, the twin of one model of the author's own and eleven projects of other fields. Its code is that of version 0.1.0, whose wheel two blind agents used with no block, as two others had used its pre-release 0.1.0rc22. The author read it on his own recorded traces, and no other person has yet used it on curves of their own, nor has anyone outside validated it.

## Continuing this work

The source and the bench on which every share below was counted are public, at https://github.com/MichelangeloDondi/twinfit. The bench holds the model of five terms whose synthetic curves gave the shares below, the programs of the counts with their tables and, for three of the five, the kept row of each draw, the reader of the criteria, and the pages of the validation, the roadmap and the decisions. The projects of the gallery are drawn by the package itself. `twinfit/CONTINUING.md`, beside this page in the folder that `python -m twinfit guide readme` names, is the page for a researcher who takes the work over: what is where and what to run first.

## Method

1. It fits one model to every curve of a design at once, with parameters shared by all the curves and parameters of each curve.
2. It reads what the fit leaves through families of statistics, generalized moments under many kernels, windows along the abscissa and along the value, and the same of transformations of the value, and sets each against a twin: synthetic curves drawn from the fitted model with the law of noise that the description gives the records, a normal law of each sample's deviation, fitted and read in the same way.
3. Its report matches the residuals across every family against each shape built from the fitted model, widths among them, and says where the model misses and along which setting.
4. A researcher, or an agent in the researcher's place, writes the physical effect as a term. twinfit keeps it only if the Bayesian criterion falls by more than 2 and the twin admits the new model, and ends each turn with three proposals for the model and three for the experiment.

![A fit that looks right, read against its twin](https://raw.githubusercontent.com/MichelangeloDondi/twinfit/main/twinfit/examples/figures/against_the_twin.png)

*A blind run of an agent on eight synthetic curves of the line of `RUBIDIUM.md`. (a, b) Turn 2, model 1, fits with a chi square of 1.014 a degree of freedom. (c, d) What it leaves, summed in pairs of windows at six scales along the density, a grid drawn for this figure where the report of the package reads six bands: noise about zero when nothing is missing, as on a set of the twin. (e, f) Model 1 stands beyond all 100 sets of its twin, and model 3 within them.*

## Output

The best model, and each parameter with five parts of its uncertainty, never folded into one number: the bar of the fit, its coverage on the twin with the bias of the estimate there, the bar of theory, the spread over models the records cannot tell apart, and none where they fix only a combination.

![Turn by turn](https://raw.githubusercontent.com/MichelangeloDondi/twinfit/main/twinfit/examples/figures/turn_by_turn.png)

*The same run, turn by turn. With model 1 the size of its term reads 31 bars too large at a chi square of 1.014. The width along the density takes it to 6.9 bars, and the width of every curve to 0.7, within its bars. The best of three blind runs on two draws of one synthetic model: it found three of its four terms, missed the smallest, and was stopped at two hours.*

## Its limits

An effect that a parameter of the model can take the place of leaves nothing in the residuals. The report says what a miss looks like, never what causes it. A model that its twin admits is not proven complete, and a term chosen after looking at the records keeps its stated level only once synthetic records with it have shown what the whole procedure returns.

What it was measured to do, each ± the standard error of its share:

- On synthetic sets drawn from one model of five terms, two of whose values were fitted to the author's own recorded curves, two held at their theory and one taken from the literature, each draw fitted with one of four of its terms left out, at the size that this model gives it, the report puts first, or tied for first, the place that it gives for that term on the same records without noise, on 0.780 ± 0.021 of 400 draws. A tie counts as a hit. At half that size, on 0.128 ± 0.017 of 400 draws.
- For the two of the four that change along a setting, it names that setting first on 0.600 ± 0.035 of 200 draws, at the floor of 0.60 that this version set itself.
- For one of the four, a width that grows steeply along a setting, as a density does, the share of the place is 0.440 ± 0.050 of 100 draws, below that floor.
- On a model with no term missing, read against the twin of each draw's own fit at the default of 30 synthetic sets, the report sees something beyond the twin, a false alarm, on 0.188 ± 0.025 of 250 draws, above the ceiling of 0.12 that this version set itself. At 100 sets, `python -m twinfit loop DIR --sets 100`, the same draws read 0.076 ± 0.017 of 250 draws.
- On the eleven synthetic projects of other fields of the gallery, each drawn again with a misfit that its folder does not hold, the report puts first where the misfit shows on 0.588 ± 0.055 of 80 draws, the project pedestal on none of its 20 draws. It names along which setting it changes, of which kind it is and its shape on all of their 120 draws, and the shape of a noise with heavy tails on 0 of 20. Pooled over every question, 0.759 ± 0.029 of 220 draws.

## By hand, or with an agent

`twinfit/TUTORIAL.md` takes a newcomer to the best model, and `twinfit/agent/TUTORIAL.md` is the same tutorial for an agent, its commands and its rules without the pictures. `twinfit/agent/AGENT.md` lets an agent, an AI assistant that runs commands, take the place of the researcher, from the card at its head. `twinfit/EXAMPLES.md` holds eleven projects of other fields.

## Requirements and installation

Python 3.12 or later, with numpy 2.5 or later and scipy 1.16 or later, and matplotlib for the visual log of a run. Nothing else. It has run on macOS, and it needs a system of the Unix kind: on Windows the loop cannot start, since it locks its files with the module `fcntl`, which Python holds on those systems alone.

numpy 2.5 is the floor of the author's choice: on a Mac with an arm64 processor, one thread, it convolves a curve of 9000 points with a kernel of 2001 in 0.80 ms, against 10.0 ms for numpy 1.26.4, while on Linux with an x86_64 processor the two took about the same time, 2.1 and 1.7 ms. numpy 2.5 asks Python 3.12 and scipy 1.16.

The package has run on Python 3.12 with numpy 2.5.0 and scipy 1.16.0, and on Python 3.14 with numpy 2.5.2 and scipy 1.18.0, each in a clean environment: on both the five minutes below and the whole gallery give the same numbers and the same words, all but their run times.

```
pip install twinfit
```

That installs the package from PyPI, with numpy and scipy. `pip install "twinfit[log]"` installs matplotlib with it, which the visual log, `python -m twinfit log DIR`, needs. A candidate of a release, before its upload, is not on PyPI: its wheel is installed with the same extra, `pip install "/path/to/twinfit-VERSION-py3-none-any.whl[log]"`, `VERSION` its number. From a copy of its source:

```
pip install numpy scipy
pip install /path/to/the/source
```

where `/path/to/the/source` is the folder that holds `pyproject.toml` and `twinfit/`. The second line installs the package. Without it, the package is used from that folder by its path:

```
export PYTHONPATH=/path/to/the/source
python -m twinfit guide
```

When pip installs the package it installs numpy 2.5 and scipy 1.16 or later in place of older ones, so a virtual environment of its own keeps other work apart, and on a Python under 3.12 it passes over this release and may take an older one of the name, where `pip install "twinfit>=0.1"` refuses. A source folder put on the path of Python passes no check of pip: there a Python under 3.11, or one without numpy or scipy, is refused by an ImportError whose last line says what is missing, and a numpy under 2.5 runs, on a Mac with an arm64 processor about ten times slower on a convolution and on Linux with an x86_64 processor about as fast.

On a computer with no Python and no administrator password, the standalone `uv` installs both in folders of your own, with no password. These are the steps:

```
mkdir -p uv && curl -fsSL https://github.com/astral-sh/uv/releases/latest/download/uv-aarch64-apple-darwin.tar.gz | tar -xz -C uv --strip-components=1
uv/uv python install 3.14
uv/uv venv --python 3.14 venv
uv/uv pip install --python venv/bin/python "twinfit[log]"
```

The last line installs the release from PyPI. A candidate, before its upload, is not there: install its wheel in its place, `uv/uv pip install --python venv/bin/python "/path/to/twinfit-VERSION-py3-none-any.whl[log]"`, or the line resolves an older release.

On another machine, the file of uv's releases that names it. uv keeps its Python and its cache in your home folder: set `UV_PYTHON_INSTALL_DIR` and `UV_CACHE_DIR` to a folder of your own to keep them there, and `MPLCONFIGDIR` where matplotlib cannot write its cache in your home folder. uv may warn that it cannot create `~/.local/bin`: the Python still installs where `UV_PYTHON_INSTALL_DIR` says. Then `source venv/bin/activate` gives the Python that holds the package, and the five minutes below begin at their second line.

## Five minutes: an atomic clock at five probe powers

```
pip install twinfit
python -m twinfit example clock clock
python -m twinfit loop clock
```

The first line installs the package, and the second writes a project into the folder `clock`: five resonances of a clock transition against the detuning of the probe, at 2 to 10 uW, the powers read by a photodiode calibrated within 2 per cent, and a first model of one Lorentzian line of one width and one center, the same at every power. `clock/description.toml` holds that calibration as an ansatz, a value held before the fit with its uncertainty, on `power scale`, which `clock/model.py` keeps in its dictionary `HELD`. The truth of the project is in the gallery, `python -m twinfit guide examples`, and the folder does not hold it.

The loop fits the five curves at once, reads their residuals against 30 synthetic sets drawn from the fit, and writes `clock/loop/model_0.md`. Its section "Where to search" gives the first answer to each question, and several stand beyond the twin. The question seen by says how the misfit is read and names no part of the curves, so read first the largest excess of the others, here the question changes along:

```
| changes along | which part of the residuals changes from curve to curve along a setting, read on the residuals themselves | the asymmetry along power | 18388.1 |
```

The report names no term: it says where the residuals are and how they change. The physics names the term. An asymmetry of a line that grows with the power of the probe is a center that moves with the power: a light shift, which the first model lacks. Write it into `clock/model.py`, whose whole text is then this, the center at a power P being the clock frequency plus the light shift times P, with P read through its calibration:

```python
"""A model of one's own: the centre of the line moves with the power of the probe, the light shift, from the frequency of
the clock at no power, which the model then measures."""
import numpy as np

SHARED = ("width", "clock frequency", "light shift")      # the full width, in Hz, the centre at no power less the nominal frequency, in Hz, and its shift, in Hz for each microwatt
OWN = ("height", "background")                           # of each record
START_SHARED = [1.0, 0.0, 0.0]
BOUNDS_SHARED = ([0.05, -3.0, -1.0], [10.0, 3.0, 1.0])
HELD = {"power scale": 1.0}                              # the powers read times this, calibrated within 2 per cent, which description.toml holds
NAME = "one Lorentzian line whose centre moves with the power"
EQUATION = r"$\frac{h}{1 + (2\,(x - \nu_0 - \delta\,P)/w)^2} + b$"          # the model in LaTeX, which the visual log typesets
LETTERS = {"h": "the height", r"\nu_0": "the clock frequency", r"\delta": "the light shift for each microwatt", "P": "the power", "w": "the width", "b": "the background"}          # each letter of EQUATION, which the visual log writes under it


def model(x, shared, own, design):
    centre = shared[1] + shared[2] * design["power"] * HELD["power scale"]
    return own[0] / (1.0 + (2.0 * (x - centre) / shared[0]) ** 2) + own[1]


def start_own(record):
    base = float(np.median(record.y))
    return [float(record.y.max() - base), base]
```

Then `python -m twinfit loop clock` again. Its report finds nothing beyond the twin, and it prints "kept: nothing. Refused: nothing", since the model of the file is the one that it fits and no term file was written. `clock/loop/best.md` quotes:

```
- clock frequency: 0.12000 ± 0.00038 of the fit
- light shift: 0.025051 ± 0.000058 of the fit ± 0.00050 of theory
```

The bar of theory is the 2 per cent of the calibration carried through the fit: the light shift is read against the power, so it carries the error of the scale of the power, and the clock frequency at no power carries none. Set both beside the truth in the gallery.

A term file `clock/terms/01.py` that holds `term = mechanisms.along("centre", "power")` would give the same light shift as the size of a term, but names no clock frequency as a parameter and reads no calibration: where the physics gives the law, write it in `model.py`. The model names no `CENTRE`, since its center is no parameter of one curve, so the proposals that place a term about the center of a curve, a tail or a period, are not fitted, and the report says so.

Each run takes from ten seconds to a minute on a laptop, the first the longer while matplotlib builds its font cache. The loop writes `clock/loop/` again at each run, so copy it to keep a report. Before a value is quoted, `best.md` may ask for `--sets 500`, which takes minutes. `python -m twinfit guide tutorial` prints the tutorial, which goes through every page that the loop writes, and `python -m twinfit guide examples` the gallery, ten more projects of other fields. `python -m twinfit guide rubidium` shows a real case: recorded traces of the two photon line of rubidium, as plots, with what the package reads on them.

## A project from measured curves

`python -m twinfit init DIR` writes the starting files of a project: its description, `description.toml`, and a first model, `model.py`. Each curve is one comma-separated file whose first line names its columns, with no line of comment and no byte order mark, which some spreadsheets write before the first column's name: a file of the lab written otherwise is converted first. The noise of the samples is one value, a column of each file, or read from neighboring samples. For a fraction of N atoms read by projection, the noise is a column of the root of p(1 - p)/N, p the fraction of the fit. Section 9 of the tutorial writes such a project whole.

## Its commands

```
python -m twinfit example NAME DIR [--noise L] [--kind K]
                                       writes a worked project into DIR: chromatography or nmr, an example of twinfit/EXAMPLES.md, or rubidium, the case of `guide rubidium`,
                                       whose noise is L times its own, 0 for none drawn, and of kind K, white, correlated or counts
python -m twinfit init DIR             writes the starting files of a project into DIR
python -m twinfit loop DIR [--n N] [--sets S] [--max-seconds T] [--workers W] [--dictionary]
                                       fits, reads and reports each model of the project in DIR, and fits the terms written in DIR/terms,
                                       up to N term files, 10 by default, each twin on S synthetic sets, 30 by default and 30 at least.
                                       With T, a run past T seconds of the clock stops at the end of its step and is refused, and the
                                       last whole run stays in DIR/loop. DIR/loop.reuse keeps each costly step, so the next run reads
                                       back what the files, records and fit have not changed, and a stopped run goes on where it stopped. With W, the twins of each model's proposals run side by side
                                       on W cores, every core by default, with the same numbers
                                       With --dictionary, the package's own list of generic terms is set against the residuals beside
                                       your list, list_of_terms.py: asked by you alone, since the package is for the term nobody suspects
python -m twinfit terms                lists the terms that a term file can name, the name that best.md gives each, and the form of a term file
python -m twinfit card                 how to write a term file, of the package or of your own, in one page
python -m twinfit log DIR              draws the visual log of the last run in DIR
python -m twinfit wait DIR [--timeout S]  waits up to S seconds for the loop that runs in DIR, 540 by default, and prints its last page for an agent
python -m twinfit literature DIR       the brief of the search of the literature that the description asks, for an agent to run
python -m twinfit scan DIR TERM X...   fits a term whose distance, stretch, period or power is chosen, at each value X
python -m twinfit matrix DIR           what every family reads on the model of DIR, and which statistic sees which law of the noise
python -m twinfit counts DIR [--floor F] writes a copy of DIR whose curves carry the noise of a counter, from the fitted values of the best model
python -m twinfit thin DIR N           writes up to three copies of DIR beside it, each keeping one sample in every N, for a noise correlated over N samples
python -m twinfit guide                prints the protocol for an agent that drives the package
python -m twinfit guide manual         prints the older manual path, for a model that one file model.py cannot hold
python -m twinfit guide examples       prints the gallery, twinfit/EXAMPLES.md: projects of other fields with their truth
python -m twinfit guide tutorial       prints the tutorial, twinfit/TUTORIAL.md, and `guide readme` this page, each with its path
python -m twinfit guide agent-tutorial prints the tutorial for an agent, twinfit/agent/TUTORIAL.md, with its path
python -m twinfit --version            prints the number of the package
python -m twinfit guide rubidium       prints a real case, recorded traces of a line of rubidium as plots, with the folder of its figures
python -m twinfit walk [HIDDEN] [DIR]  walks one rung of that older manual path, on synthetic records with a hidden term
python -m twinfit.tutorial --check     runs every command of the tutorial again, with matplotlib, and says whether an output changed
```

The tutorial, from nothing to the best model, is `twinfit/TUTORIAL.md`. Its first project is `python -m twinfit example nmr my_nmr`.
