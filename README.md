# twinfit

Where a model lacks a term: the structure of the residuals of a fit, read against its twin.

## What it is for

A model fitted to measured curves can pass the usual test of fit and still miss an effect, a broadening, a background or a drift, which biases the quoted parameters while the chi square stays near one. twinfit finds the structure that the effect leaves, says what it looks like and measures what it does to the parameters, for curves of any field.

## How it works

1. It fits one model to every curve of a design at once, with parameters shared by all the curves and parameters of each curve.
2. It reads what the fit leaves through families of statistics, among them sums over pairs of windows about the centre of the curve, and sets each against a twin: synthetic curves drawn from the fitted model with the noise of the records, fitted and read in the same way.
3. Its report says where the model misses, which shape among those that it builds from the fitted model the miss takes, and along which setting of the design it changes.
4. A researcher, or an agent in the researcher's place, writes the physical effect as a term. twinfit keeps it only if the Bayesian criterion falls by more than 2 and the twin admits the new model, and ends each turn with three proposals for the model and three for the experiment.

![A fit that looks right, read against its twin](https://raw.githubusercontent.com/MichelangeloDondi/twinfit/main/twinfit/examples/figures/against_the_twin.png)

*A blind run of an agent on eight synthetic curves of the line of `RUBIDIUM.md`. (a, b) Turn 2, model 1, fits with a chi square of 1.014 a degree of freedom. (c, d) What it leaves, summed in pairs of windows at six scales along the density, a grid drawn for this figure where the report of the package reads six bands: noise about zero when nothing is missing, as on a set of the twin. (e, f) Model 1 stands beyond all 100 sets of its twin, and model 3 within them.*

## What it gives

The best model, and each parameter with five parts of its uncertainty, never folded into one number: the bar of the fit, its coverage on the twin, the bar of theory, the spread over models the records cannot tell apart, and none where they fix only a combination.

![Turn by turn](https://raw.githubusercontent.com/MichelangeloDondi/twinfit/main/twinfit/examples/figures/turn_by_turn.png)

*The same run, turn by turn. With model 1 the size of its term reads 31 bars too large at a chi square of 1.014. The width along the density takes it to 6.9 bars, and the width of every curve to 0.7, within its bars. The best of three blind runs on two draws of one synthetic model: it found three of its four terms, missed the smallest, and was stopped at two hours.*

## Its limits

An effect that a parameter of the model can take the place of leaves nothing in the residuals. The report says what a miss looks like, never what causes it. A model that its twin admits is not proven complete, and a term chosen after looking at the records keeps its stated level only once synthetic records with it have shown what the whole procedure returns.

## By hand, or with an agent

`twinfit/TUTORIAL.md` takes a newcomer to the best model. `twinfit/agent/AGENT.md` lets an agent, an AI assistant that runs commands, take the place of the researcher. `twinfit/EXAMPLES.md` holds nine projects of other fields.

## What it needs

Python 3.12 or later, with numpy 2.5 or later and scipy 1.16 or later. Nothing else. numpy 2.5 is the floor of the author's choice: on a Mac with an arm64 processor, one thread, it convolves a curve of 9000 points with a kernel of 2001 in 0.80 ms, against 10.0 ms for numpy 1.26.4, while on Linux with an x86_64 processor the two took about the same time, 2.1 and 1.7 ms. numpy 2.5 asks Python 3.12 and scipy 1.16. The package has run on Python 3.12 with numpy 2.5.0 and scipy 1.16.0, and on Python 3.14 with numpy 2.5.2 and scipy 1.18.0, each in a clean environment: on both the five minutes below and the whole gallery give the same numbers and the same words, all but their run times.

```
pip install twinfit
```

That installs the package from PyPI, with numpy and scipy. From a copy of its source:

```
pip install numpy scipy
pip install /path/to/the/source
```

where `/path/to/the/source` is the folder that holds `pyproject.toml` and `twinfit/`. The second line installs the package. Without it, the package is used from that folder by its path:

```
export PYTHONPATH=/path/to/the/source
python -m twinfit guide
```

pip refuses a Python under 3.12, a numpy under 2.5 and a scipy under 1.16 when it installs the package. A source folder put on the path of Python passes no check of pip: there a Python under 3.11, or one without numpy or scipy, is refused by an ImportError whose last line says what is missing, and a numpy under 2.5 runs, about ten times slower on a convolution.

## Five minutes: an atomic clock read at five powers of the probe

```
pip install twinfit
python -m twinfit example clock clock
python -m twinfit loop clock
```

The first command writes a project into the folder `clock`: five resonances of a clock transition against the detuning of the probe, at 2 to 10 uW, the powers read by a photodiode calibrated within 2 per cent, and a first model of one Lorentzian line of one width and one centre, the same at every power. `clock/description.toml` holds that calibration as an ansatz, a value held before the fit with its uncertainty, on `power scale`, which `clock/model.py` keeps in its dictionary `HELD`. The truth of the project is in the gallery, `python -m twinfit guide examples`, and the folder does not hold it.

The loop fits the five curves at once, reads their residuals against 30 synthetic sets drawn from the fit, and writes `clock/loop/model_0.md`. Its section "Where to search" gives the first answer to each question, and several stand beyond the twin. Read first the one of the largest excess, here the question along:

```
| along | the setting along which a shared parameter of the model would change, read on the fit of each curve | power, read on centre | 52264.2 |
```

The centre moves with the power: a light shift, which the first model lacks. Write it into `clock/model.py`, whose whole text is then this, the centre at a power P being the clock frequency plus the light shift times P, with P read through its calibration:

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

The first proposal of the report, a term file `clock/terms/01.py` that holds `term = mechanisms.along("centre", "power")`, gives the same light shift as the size of a term, but names no clock frequency as a parameter and reads no calibration: where the physics gives the law, write it in `model.py`. The model names no `CENTRE`, since its centre is no parameter of one curve, so the proposals that place a term about the centre of a curve, a tail or a period, are not fitted, and the report says so.

Each run takes 10 to 20 seconds on a laptop. The loop writes `clock/loop/` again at each run, so copy it to keep a report. Before a value is quoted, `best.md` may ask for `--sets 500`, which takes minutes. `python -m twinfit guide tutorial` prints the tutorial, which goes through every page that the loop writes, and `python -m twinfit guide examples` the gallery, eight more projects of other fields. `python -m twinfit guide rubidium` shows a real case: recorded traces of the two photon line of rubidium, as plots, with what the package reads on them.

## Your own curves

`python -m twinfit init DIR` writes the starting files of a project: its description, `description.toml`, and a first model, `model.py`. Each curve is one comma-separated file whose first line names its columns, with no line of comment: a file of the lab written otherwise is converted first. The noise of the samples is one value, a column of each file, or read from neighbouring samples. For a fraction of N atoms read by projection, the noise is a column of the root of p(1 - p)/N, p the fraction of the fit. Section 9 of the tutorial writes such a project whole.

## Its commands

```
python -m twinfit example NAME DIR [--noise L] [--kind K]
                                       writes a worked project into DIR: chromatography or nmr, or an example of twinfit/EXAMPLES.md,
                                       whose noise is L times its own, 0 for none drawn, and of kind K, white, correlated or counts
python -m twinfit init DIR             writes the starting files of a project into DIR
python -m twinfit loop DIR [--n N] [--sets S] [--max-seconds T]
                                       fits, reads and reports each model of the project in DIR, and fits the terms written in DIR/terms,
                                       up to N term files, 10 by default, each twin on S synthetic sets, 30 by default and 30 at least.
                                       With T, a run past T seconds of the clock stops at the end of its step and is refused, and the
                                       last whole run stays in DIR/loop
python -m twinfit terms                lists the terms of the package, the name by which a report prints each, and the form of a term file
python -m twinfit scan DIR TERM X...   fits a term whose distance, stretch, period or power is chosen, at each value X
python -m twinfit counts DIR [--floor F] writes a copy of DIR whose curves carry the noise of a counter, from the fitted values of the best model
python -m twinfit thin DIR N           writes up to three copies of DIR beside it, each keeping one sample in every N, for a noise correlated over N samples
python -m twinfit guide                prints the protocol for an agent that drives the package
python -m twinfit guide manual         prints the older manual path, for a model that one file model.py cannot hold
python -m twinfit guide examples       prints the gallery, twinfit/EXAMPLES.md: projects of other fields with their truth
python -m twinfit guide tutorial       prints the tutorial, twinfit/TUTORIAL.md, and `guide readme` this page, each with its path
python -m twinfit guide rubidium       prints a real case, recorded traces of a line of rubidium as plots, with the folder of its figures
python -m twinfit walk [HIDDEN] [DIR]  walks one rung of that older manual path, on synthetic records with a hidden term
python -m twinfit.tutorial --check     runs every command of the tutorial again and says whether an output changed
```

The tutorial, from nothing to the best model, is `twinfit/TUTORIAL.md`. Its first project is `python -m twinfit example nmr my_nmr`.
