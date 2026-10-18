# twinfit

Where a model lacks a term: the structure of the residuals of a fit, read against its twin.

## What it is for

A model fitted to measured curves can pass the usual test of fit and still miss an effect: a broadening, a background, a drift, a distortion of the instrument. The missing effect biases the parameters that are quoted, and a chi square near one does not show it, since it sums over every sample what is often a small structure in one part of the curves. twinfit finds that structure, says what it looks like, and measures what it does to the parameters, for curves of any field: spectra, chromatograms, light curves, resonances, bands of a gel.

## How it works

1. It fits one model to every curve of a design at once, curves recorded at several temperatures, powers, flows or filters, with parameters shared by all the curves and parameters of each curve.
2. It reads what the fit leaves through families of statistics: moments of the residuals with several kernels and their ratios, moments of their derivative and of their running integral, windows along the abscissa, bands along the height of the curve, crossings of levels, and the chi of every sample.
3. It sets every reading against a twin: synthetic curves drawn from the fitted model with the noise of the records, fitted and read in the same way. A reading that the twin does not reach is evidence that the model misses something, and the twin says how large a reading of a right model may be.
4. It compares what the fit leaves with generic shapes built from the fitted model itself: a shifted copy of the curve, a wider pedestal under it, a tail on one side, a slow polynomial of the baseline, a ripple, glitches, a noise that grows with the signal or is correlated from sample to sample. Its report says where the model misses, in the core, the flanks, the wings or the asymmetry of the curves, which of those shapes the miss resembles, and along which setting of the design it changes.
5. A researcher, or an agent in the researcher's place, writes the physical effect that would leave that structure as a term of the model. twinfit fits it, keeps or refuses it on its twin, and reads the new model again. That is the loop.
6. Each turn of the loop ends with three proposals for the model and three for the experiment, each measured on the twin. A proposal for the model comes with the bias of every known parameter on the twin without the proposed term and with it.

![A fit that looks right, read against its twin](https://raw.githubusercontent.com/MichelangeloDondi/twinfit/main/twinfit/examples/figures/against_the_twin.png)

*A blind run of an agent with the package, on eight synthetic curves of the line of the real case of `RUBIDIUM.md`, at four temperatures and five powers, drawn from a model whose terms the agent was never told. (a, b) Its second model, the line widened by the time that each atom spends in the beam, fits with a chi square of 1.014 a degree of freedom, and what it leaves looks like noise sample by sample. (c) The distances from the centre of the line cut into 2, 4, 8, 16, 32 and 64 windows, each a pair at equal distances on the two sides, and in each the sum of what the fit leaves, the curves weighted by the density of the vapour less its mean, in its standard deviation after the fit: a number that is noise about zero when nothing is missing, as (d) shows on a synthetic set of the twin. (e) The moments of orders 0, 2 and 4 within windows of growing half width, along the density, beside each of 100 sets of the twin. (f) Each scale read together: where the windows are wide, what the fit leaves adds up and the noise does not, and the second model stands beyond all 100 sets of its twin, while the fourth, with the width that grows with the density and the width of every curve that the agent wrote next, stands within it.*

## What it gives

The best model after the turns, and each of its parameters with five parts of its uncertainty, never folded into one number: the bar of the fit, the coverage of its interval on the twin with the bias beside it, the bar of the values held from theory, the spread over the models that the records cannot tell apart from the best, and no bar at all where the records fix a combination of two parameters and not each.

![Turn by turn](https://raw.githubusercontent.com/MichelangeloDondi/twinfit/main/twinfit/examples/figures/turn_by_turn.png)

*The same run, turn by turn. Top, what each model leaves, every sample. Middle, its windows along the density, read together against 100 sets of its twin. Bottom, each parameter less the value that made the curves, which the agent never saw, in its own bars, and under each turn what made the agent write the next term. With the time in the beam alone the size of that term reads 31 bars too large while the chi square is 1.014: the terms that the model lacks hide in it. The width along the density, which the windows asked for, takes it to 6.9 bars. The width of every curve, asked because that size still stood beyond the theory of the beam that the description gives, takes it to 0.7, and every parameter of the fourth model lies within its bars of the truth.*

## Its limits

If a parameter already in the model can take the place of the missing effect completely, the effect leaves nothing in the residuals, and no statistic can see it. A structure unlike every generic shape is described by the nearest of them, so the report says what the miss looks like and never what causes it: the physics is the researcher's. And a term chosen because the report ranked it first was chosen after looking at the records, so its interval keeps its stated level only once synthetic records with the term have shown what the whole procedure returns.

## By hand, or with an agent

`twinfit/TUTORIAL.md` takes a newcomer from nothing to the best model. `twinfit/agent/AGENT.md` is a protocol that lets an agent, an AI assistant that runs commands, take the place of the researcher: read the report, choose a term from it and from the physics, have it checked by a second agent, fit it and report. `twinfit/EXAMPLES.md` holds nine projects of other fields, each drawn with a truth that its folder does not hold, and what the first report of each led with.

## In more words

You give it curves recorded across a design, such as spectra at several temperatures, chromatograms at several flows or light curves of several stars, with the model that should describe them. It fits the model to every curve at once, reads the structure of what the fit leaves, curve by curve, with families of statistics of the curve, of its values and of their transformations, and reads how that structure changes from one setting of the design to the next. It sets every reading against a twin: synthetic curves drawn from the fitted model with the noise of the records, fitted and read the same way. It then says where a term is missing, in the core of the curves, in their wings or in their asymmetry, of which kind, and along which setting. It fits the terms that you, or an agent, write in answer, keeps or refuses each, and ends each iteration with three proposals for the model and three for the experiment. Its output is the best model, with five parts of the uncertainty of each parameter.

## What it needs

Python 3.12 or later, with numpy 2.5 or later and scipy 1.16 or later. Nothing else. numpy 2.5 is the floor because it convolves about ten times faster than the versions before it, which a term with a kernel uses on every curve, and it asks Python 3.12 and scipy 1.16. The package has run on Python 3.12 with numpy 2.5.0 and scipy 1.16.0, and on Python 3.14 with numpy 2.5.2 and scipy 1.18.0, each in a clean environment: on both the five minutes below and the whole gallery give the same numbers and the same words, all but their run times.

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
