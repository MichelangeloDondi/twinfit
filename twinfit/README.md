# twinfit

Find what your model is missing. twinfit reads the residuals of a fit and tells you where the model misses, what the miss looks like, and what it would do to the parameters you share across the curves.

In five minutes, the clock below. In an hour, sections 1 to 6 of `twinfit/TUTORIAL.md`. The projects of other fields are in `twinfit/EXAMPLES.md`, and the protocol for an agent in `twinfit/agent/AGENT.md`.

> Section 1 of `twinfit/TUTORIAL.md` defines the terms used on this page.

## Purpose

A model fitted to measured curves can pass the usual test of fit and still lack an effect: a broadening, a background, a drift. The missing effect biases the fitted parameters while the reduced chi square stays near one.

twinfit starts from the data alone. It computes statistics of the residuals and compares each with the same statistic on synthetic curves drawn from the fitted model, its twin. Where the residuals stand out, it reports what its statistics can read of them: in which region of the curve the excess lies (in the core or in the wings, for example, or as an asymmetry between the two sides), what structure they see (a slow trend or a periodic structure, for example, or a departure of the noise from its law), the size of the excess, the setting of the experiment along which it varies, and what each shared parameter would move by if the excess were a real term of that size. A structure that none of its statistics reads still shows as an excess, without a name. The physics comes next, and it is yours: you write the term that explains the structure, and twinfit keeps it if the records support it. It suggests no term of its own unless you ask for its generic list with `--dictionary`.

Version 0.1.0 was measured on synthetic curves of the author's own experiment, the real case of `python -m twinfit guide rubidium`, and on eleven synthetic projects of other fields. The results are summarized below and given in full in `docs/VALIDATION.md` of the repository.

## Continuing this work

The source and the bench on which the package was measured are public at https://github.com/MichelangeloDondi/twinfit: the model that drew the synthetic curves, the programs of the counts with their tables and, for three of the counts, their kept rows, and the pages of the validation, the roadmap and the decisions. `twinfit/CONTINUING.md`, in the folder that the last line of `python -m twinfit guide readme` names, is written for a researcher who takes the work over.

## Method

1. One model is fitted to all the curves at once, with parameters shared by the curves and parameters of each curve.
2. The residuals are read through families of statistics: generalized moments under several kernels, sums over windows along the abscissa and along the value, and the same for transformations of the value. Each statistic is compared with its value on the twin, synthetic curves drawn from the fitted model with the law of noise of the records, fitted and read in the same way.
3. The report reads, across all the families, what the residuals hold once what the model's own parameters can produce is removed, a slow trend or a periodic structure among them, and states where the model misses and along which setting.
4. You, or an agent in your place, write the physical effect as a term. twinfit keeps it if the Bayesian information criterion falls by more than 2 and the twin admits the new model. Each turn ends with three proposals for the model and three for the experiment.

![A fit that looks right, read against its twin](https://raw.githubusercontent.com/MichelangeloDondi/twinfit/main/twinfit/examples/figures/against_the_twin.png)

*A blind run of an agent on eight synthetic curves of the line of `RUBIDIUM.md`. (a, b) Turn 2, model 1, which fits with a chi square of 1.014 a degree of freedom. (c, d) The residuals, summed in pairs of windows at six scales along the density, on a grid drawn for this figure where the report reads six bands. When nothing is missing, as on a set of the twin, the sums are noise about zero. (e, f) Model 1 lies beyond all 100 sets of its twin, and model 3 within them.*

## Output

The best model and, for each parameter, its uncertainty in five parts kept separate: the standard error of the fit, its coverage on the twin with the bias of the estimate there, the uncertainty carried from theory, the spread over the models that the records cannot tell apart, and no value where the records fix only a combination of parameters.

![Turn by turn](https://raw.githubusercontent.com/MichelangeloDondi/twinfit/main/twinfit/examples/figures/turn_by_turn.png)

*The same run, turn by turn. With model 1 the size of its term reads 31 bars too large at a chi square of 1.014. A width that varies along the density takes it to 6.9 bars, and a width of each curve to 0.7, within its bars. This was the best of three blind runs on two draws of one synthetic model. It found three of the four terms, missed the smallest, and was stopped after two hours.*

## What was measured

Each share is given with its standard error. Two stand outside their thresholds and are marked. The full table is in `docs/VALIDATION.md`.

- The place of a hidden term, put first or tied for first: 0.780 ± 0.021 of 400 draws, on a model of five terms fitted with one of four terms left out at its size. First alone: 0.432 ± 0.025. At half that size: 0.128 ± 0.017 of 400.
- The setting along which the term varies, named first: 0.600 ± 0.035 of 200 draws.
- The place of a width that grows steeply along a setting: 0.440 ± 0.050 of 100 draws, below its threshold of 0.60.
- False alarms on a model with nothing missing, at the default of 30 sets: 0.188 ± 0.025 of 250 draws, above its ceiling of 0.12, and 0.076 ± 0.017 at 100 sets.
- Eleven projects of other fields, pooled over every question: 0.759 ± 0.029 of 220 draws, the place of the misfit alone 0.588 ± 0.055 of 80, with the project pedestal and the shape of a noise with heavy tails on 0 of 20 each.

An effect that a parameter of the model can absorb leaves nothing in the residuals, and a model that its twin admits is not thereby complete.

## By hand or with an agent

`twinfit/TUTORIAL.md` takes you from a first project to the best model. `twinfit/agent/TUTORIAL.md` is the same for an agent, and `twinfit/agent/AGENT.md` lets an agent take your place. `twinfit/EXAMPLES.md` holds eleven projects of other fields.

## Requirements and installation

Python 3.12 or later, numpy 2.5 or later, scipy 1.16 or later, and matplotlib for the visual log. On Windows the loop cannot start, because it locks its files with `fcntl`.

```
pip install twinfit
```

`pip install "twinfit[log]"` adds matplotlib. A candidate of a release is not on PyPI: install its wheel in its place, `pip install "/path/to/twinfit-VERSION-py3-none-any.whl[log]"`. From a copy of the source, `pip install /path/to/the/source`, or set `PYTHONPATH` to that folder.

numpy 2.5 is required because it convolves about ten times faster than numpy 1.26 on an arm64 Mac. pip replaces older versions of numpy and scipy, so a virtual environment of its own is advisable.

On a computer with no Python and no administrator password, `uv` installs everything in your own folders:

```
mkdir -p uv && curl -fsSL https://github.com/astral-sh/uv/releases/latest/download/uv-aarch64-apple-darwin.tar.gz | tar -xz -C uv --strip-components=1
uv/uv python install 3.14
uv/uv venv --python 3.14 venv
uv/uv pip install --python venv/bin/python "twinfit[log]"
```

The last line installs the release from PyPI. For a candidate, install its wheel with `uv/uv pip install --python venv/bin/python "/path/to/twinfit-VERSION-py3-none-any.whl[log]"`. On another machine, use the uv archive for that machine. `UV_PYTHON_INSTALL_DIR`, `UV_CACHE_DIR` and `MPLCONFIGDIR` move the caches of uv and matplotlib out of your home folder when it cannot be written. uv also tries to link its Python under `~/.local/bin` and warns when it cannot. The install works regardless. Then `source venv/bin/activate`, and start at the second line of the five minutes.

## Five minutes: an atomic clock at five probe powers

```
pip install twinfit
python -m twinfit example clock clock
python -m twinfit loop clock
```

The second line writes a project into `clock`: five resonances of a clock transition against the detuning of the probe, at 2 to 10 uW, with a first model of one Lorentzian line of one width and one center. The powers were read by a photodiode calibrated within 2 per cent, which `clock/description.toml` holds as an ansatz under `power scale`. The truth of the project is in the gallery, `python -m twinfit guide examples`.

The loop fits the five curves at once, reads their residuals against 30 synthetic sets drawn from the fit, and writes `clock/loop/model_0.md`. Its section "Where to search" gives the first answer to each question. The largest excess that names a part of the curves is here the question `changes along`:

```
| changes along | which part of the residuals changes from curve to curve along a setting, read on the residuals themselves | the asymmetry along power | 18388.1 |
```

The report names no term. The physics does: an asymmetry that grows with the power of the probe is a center that moves with the power, a light shift, which the first model lacks. Write it into `clock/model.py`. The center at a power P becomes the clock frequency plus the light shift times P:

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

Run `python -m twinfit loop clock` again. The report finds nothing beyond the twin, and `clock/loop/best.md` reports:

```
- clock frequency: 0.12000 ± 0.00038 of the fit
- light shift: 0.025051 ± 0.000058 of the fit ± 0.00050 of theory
```

The uncertainty from theory is the 2 per cent of the calibration, carried through the fit to the light shift and not to the clock frequency. Compare both with the truth in the gallery.

A term file `clock/terms/01.py` with `term = mechanisms.along("centre", "power")` would give the same light shift as the size of a term, but it names no clock frequency and reads no calibration. Where the physics gives the law, write it in `model.py`.

Each run takes between ten seconds and a minute on a laptop. The loop writes `clock/loop/` again at each run, so copy it to keep a report. `best.md` may ask for `--sets 500` before a value is quoted. `python -m twinfit guide tutorial` prints the tutorial, `guide examples` the gallery, and `guide rubidium` a real case: recorded traces of the two photon line of rubidium, as plots, with what the package reads on them.

## A project from measured curves

`python -m twinfit init DIR` writes the starting files of a project: `description.toml` and a first `model.py`. Each curve is one comma-separated file whose first line names its columns, without a comment line or a byte order mark. The noise of the samples is one value, a column of each file, or read from neighboring samples. Section 9 of the tutorial writes such a project in full.

## Commands

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
python -m twinfit terms                lists the terms of its table, the name that best.md gives each, and the form of a term file. The table lacks exponential_width until 0.1.1 (docs/ROADMAP.md, section 4)
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

The tutorial, from a first project to the best model, is `twinfit/TUTORIAL.md`. Its first project is written by `python -m twinfit example nmr my_nmr`.
