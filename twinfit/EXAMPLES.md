# Examples of the package

11 projects of other fields, each written by `python -m twinfit example NAME DIR`, each drawn with a truth that its folder does not hold. The projects pedestal, satellite, saturation, filter, fano and counts were each first worked by a fresh agent that was not told the truth, before they joined the gallery, and clock, gel, transit, fringes and returns were written for the gallery and worked by no agent yet.

This page builds on sections 1 to 6 of `TUTORIAL.md` and sets out what the loop finds on each project against its truth, at its own noise and at four others. Not covered here: the run of a project of one's own, which section 9 of the tutorial gives.

> Section 1 of `TUTORIAL.md` defines every word of the package that this page uses.

For each the loop ran three times: on the first model alone, after what to do, and after what not to do. Then at four other noises, on the first model and after what to do: three times its own, none drawn, and a second kind of noise as large as its own and three times as large. `python -m twinfit example NAME DIR --noise L --kind K` writes each: `L` is a multiple of the example's own noise, 0 for none drawn, and `K` one of `white`, `correlated` and `counts`, the example's own kind where `K` is not given.

Every number below is read from those runs, which `python -m twinfit.examples.gallery --write` makes again for those who build the package.

Find the example whose first report leads as yours does, and read what to do there. Each answer of a report is given with its excess, 0 at the median of its statistic over the synthetic sets of the twin and 1 at their 95th percentile, and each report with its line, the excess that the largest answer of a synthetic set passes on 5 sets of 100: an answer past it stands beyond the twin. The excess says how clearly a term shows against the noise, and not how large the term is.

| example | the field | what the first report puts first beyond its twin |
|---|---|---|
| pedestal | gas cell absorption lines at five pressures | where: the wings, excess 52.2 · where: the flanks, excess 45.4 · its line 2.73 |
| satellite | gas cell absorption lines with weak satellites | where: the asymmetry, excess 9.12e+03 · where: the wings, excess 4.77e+03 · where: the flanks, excess 4.49e+03 · its line 2.73 |
| saturation | laser-driven resonance fluorescence at five powers | changes along: the flanks along power, excess 1.17e+04 · changes along: the core along power, excess 8.97e+03 · where: the core, excess 5.08e+03 · its line 3.08 |
| filter | absorption dips behind a detector's low-pass filter | what: a noise correlated from sample to sample, excess 73.3 · changes along: the core along temperature, excess 26.8 · where: the flanks, excess 12.8 · its line 3.28 |
| fano | photoionization across an autoionizing resonance at four pressures | where: the asymmetry, excess 5.84e+05 · changes along: the asymmetry along pressure, excess 3.29e+04 · where: the core, excess 2.29e+04 · its line 3.18 |
| counts | emission-line photon counts at four pump powers | what: a noise that grows with the signal, excess 347 · what: glitches, excess 66.3 · changes along: the core along power, excess 26.3 · its line 2.37 |
| clock | atomic clock resonance at five probe powers | changes along: the asymmetry along power, excess 1.84e+04 · changes along: the core along power, excess 1.49e+04 · changes along: the flanks along power, excess 1.61e+03 · its line 2.47 |
| gel | gel electrophoresis band at four voltages | where: the flanks, excess 3.01e+03 · where: the asymmetry, excess 3e+03 · where: the wings, excess 544 · its line 2.48 |
| transit | transit of an exoplanet in four filters | where: the core, excess 1.28e+03 · changes along: the core along wavelength, excess 164 · what: a noise correlated from sample to sample, excess 22.6 · its line 2.74 |
| fringes | piezo-scanned Michelson interferometer at five wavelengths | what: a noise correlated from sample to sample, excess 44.8 · where: the asymmetry, excess 14.6 · what: glitches, excess 6.94 · its line 2.42 |
| returns | daily stock index returns in five volatility periods | where: the wings, excess 47.5 · changes along: the flanks along volatility, excess 40.2 · where: the core, excess 37.9 · its line 2.82 |

## The tables at other noises

Each row is the project written with `--noise L --kind K`, its first model run, and then what to do, whose term the third column gives. For a term written at the top of a scan, it gives the range of the scan too: the values whose fall is within 1 of its top.

At `L` times its noise the description states `L` times the noise, so the bars grow about `L` times. With none drawn it states the example's own noise: the bars stay those of the own row, the chi square falls near 0, and the check of the noise says that the noise is stated too large, which is why those values are not yet to be quoted.

The rows of white and of correlated noise draw one stream of random numbers, scaled: a row at three times its noise holds the draw of the row at its own noise, three times larger. Their chi squares are then equal, unless the project also draws something of a fixed size, as the glitches of the saturation, and the two rows are not independent experiments.

A counter's noise at three times is drawn as nine times fewer counts, scaled up nine times, so the root of the counts is a third of the noise. What to do for a counter, `python -m twinfit counts DIR`, states that root: there it states the noise three times too small, and the chi square comes near 9.

A parameter of a term written at the top of a scan is named by that term, as "as pedestal 8": its size is read at that stretch, distance or period, and a size at another one is another quantity. So the sizes of two rows are compared at one stretch only. Where the range of the scan is wide, the records do not fix the stretch, and the bar of the size, read at one stretch, is too small.

## pedestal: gas cell absorption lines at five pressures

`python -m twinfit example pedestal DIR`. The truth, which the folder does not hold: a pedestal 6 times wider, of relative area 0.04, under a Gaussian of width 4 MHz.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 1 | model 0 |
| its terms | one Gaussian line of one width | one Gaussian line of one width, a pedestal 6 times wider than the curve | one Gaussian line of one width |
| chi square for each degree of freedom | 1.124 | 0.994 | 1.124 |
| the last report, first beyond its twin | where: the wings, excess 52.2 · where: the flanks, excess 45.4 · its line 2.73 | nothing | where: the wings, excess 52.2 · where: the flanks, excess 45.4 · its line 2.73 |
| width, true 4 | 4.0128 ± 0.0011, not yet to be quoted | 3.9986 ± 0.0012 | 4.0128 ± 0.0011, not yet to be quoted |
| pedestal 6, true 0.04 | not in the model | 0.0418 ± 0.0018 | not in the model |

The first report proposes for the model: a slow trend over the span of each curve · the shape that the residuals hold in the wings · a periodic structure of 2.5 cycles over the span. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.76 · a wider span, a better data collection: the error of the known parameters times 0.76.

**What to do**: the report puts the misfit in the wings and the flanks: in a gas cell, wings that the line lacks are a wider component under it, a pedestal, so scan its stretch and write the term at the top of the scan.

Here `python -m twinfit scan DIR "pedestal(x)" 3 4 6 8 12` puts its top at 6: the term file holds `term = mechanisms.pedestal(6)`.

**What not to do**: make the width change along the pressure, a setting that the report does not name. The loop refused it on 30 synthetic sets: width changes along pressure.

**At other noises**: the first model, then what to do. Its own noise is white, and its second kind correlated over neighbouring samples. The rows are written by `--noise 3`, `--noise 0`, `--kind correlated` and `--noise 3 --kind correlated`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | width, true 4 | pedestal 6, true 0.04 |
|---|---|---|---|---|---|
| white, 3 times as large as its own noise | where: the flanks, excess 6.11 · where: the wings, excess 4.84 · its line 2.86 | model 1, 0.994, `term = mechanisms.pedestal(8)`, its scan within 1 of its top from 8 to 26 | nothing | 3.9971 ± 0.0036 | 0.083 ± 0.010, as pedestal 8 |
| no noise drawn | where: the wings, excess 45.6 · where: the flanks, excess 38.6 · where: the core, excess 2.86 · its line 2.69 | model 1, 0.000, `term = mechanisms.pedestal(6.73)`, its scan within 1 of its top from 5.36 to 6.73 | nothing | 4.0005 ± 0.0012, not yet to be quoted | 0.0499 ± 0.0023, as pedestal 6.73, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | where: the flanks, excess 66.7 · where: the wings, excess 42.7 · what: a noise correlated from sample to sample, excess 17.9 · its line 2.91 | model 1, 0.948, `term = mechanisms.pedestal(10.1)`, its scan within 1 of its top from 8.26 to 12 | what: a noise correlated from sample to sample, excess 17.9 · changes along: the core along pressure, excess 17.2 · changes along: the flanks along pressure, excess 7.92 · its line 2.65 | 3.9985 ± 0.0012, not yet to be quoted | 0.1420 ± 0.0058, as pedestal 10.1, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | what: a noise correlated from sample to sample, excess 18.4 · changes along: the core along pressure, excess 15.3 · where: the flanks, excess 13.5 · its line 2.68 | model 1, 0.948, `term = mechanisms.pedestal(56)`, its scan within 1 of its top from 36 to 60 | what: a noise correlated from sample to sample, excess 18.4 · changes along: the core along pressure, excess 16 · changes along: the flanks along pressure, excess 9.22 · its line 2.63 | 3.9931 ± 0.0035, not yet to be quoted | 23.2 ± 2.3, as pedestal 56, not yet to be quoted |

## satellite: gas cell absorption lines with weak satellites

`python -m twinfit example satellite DIR`. The truth, which the folder does not hold: a second component 12 MHz above, of relative area 0.08, beside a Gaussian of width 4 MHz.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 1 | model 0 |
| its terms | one Gaussian line of one width | one Gaussian line of one width, a second component 12.015 away | one Gaussian line of one width |
| chi square for each degree of freedom | 38.051 | 0.995 | 38.051 |
| the last report, first beyond its twin | where: the asymmetry, excess 9.12e+03 · where: the wings, excess 4.77e+03 · where: the flanks, excess 4.49e+03 · its line 2.73 | nothing | where: the asymmetry, excess 9.12e+03 · where: the wings, excess 4.77e+03 · where: the flanks, excess 4.49e+03 · its line 2.73 |
| width, true 4 | 4.1163 ± 0.0011, not yet to be quoted | 3.9989 ± 0.0011 | 4.1163 ± 0.0011, not yet to be quoted |
| second component at 12, true 0.08 | not in the model | 0.08004 ± 0.00021, as second component at 12.015 | not in the model |

The first report proposes for the model: the shape that the residuals hold in the asymmetry · a slow trend over the span of each curve · a periodic structure of 2 cycles over the span. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · pressure past its range, a better data collection: the error of the known parameters times 0.90 · more samples, a better data collection: the error of the known parameters times 0.95.

**What to do**: the report puts the misfit in the asymmetry: in a gas cell, a weak line beside each line is a second component, so scan its distance, `python -m twinfit scan DIR "second(x)" ...`, and write the term at the top of the scan.

Here `python -m twinfit scan DIR "second(x)" 6 8 10 12 14 16 18` puts its top at 11.95, and a finer scan about it, at 11.8522 11.8845 11.9167 11.949 11.9813 12.0135 12.0458, at 12.02: the term file holds `term = mechanisms.second(12.015)`.

**What not to do**: widen the line with a kernel, as if the satellite were wings of the line. The loop refused it on 30 synthetic sets: a kernel with exponential wings, by a hair, which the loop itself asks to run again at 500 sets before the term is set aside: at more sets it may be kept, and that a wrong term is kept says that the noise check and the chi square, not the twin alone, judge it.

**At other noises**: the first model, then what to do. Its own noise is white, and its second kind correlated over neighbouring samples. The rows are written by `--noise 3`, `--noise 0`, `--kind correlated` and `--noise 3 --kind correlated`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | width, true 4 | second component at 12, true 0.08 |
|---|---|---|---|---|---|
| white, 3 times as large as its own noise | where: the asymmetry, excess 1.15e+03 · where: the flanks, excess 529 · where: the wings, excess 527 · its line 2.93 | model 1, 0.995, `term = mechanisms.second(12.046)`, its scan within 1 of its top from 12 to 12.1 | nothing | 3.9967 ± 0.0032 | 0.08011 ± 0.00062, as second component at 12.046 |
| no noise drawn | where: the asymmetry, excess 9.02e+03 · where: the wings, excess 4.8e+03 · where: the flanks, excess 4.54e+03 · its line 2.75 | model 1, 0.000, `term = mechanisms.second(12)`, its scan within 1 of its top from 12 to 12 | nothing | 4.0000 ± 0.0011, not yet to be quoted | 0.08000 ± 0.00021, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | where: the asymmetry, excess 1.02e+04 · where: the flanks, excess 4.86e+03 · where: the wings, excess 4.73e+03 · its line 2.91 | model 1, 0.939, `term = mechanisms.second(12.045)`, its scan within 1 of its top from 12 to 12.1 | what: a noise correlated from sample to sample, excess 18.6 · changes along: the asymmetry along pressure, excess 13 · changes along: the core along pressure, excess 9.35 · its line 2.9 | 3.9965 ± 0.0011, not yet to be quoted | 0.08002 ± 0.00021, as second component at 12.045, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | where: the asymmetry, excess 1.33e+03 · where: the flanks, excess 619 · where: the wings, excess 429 · its line 2.63 | model 1, 0.939, `term = mechanisms.second(12.136)`, its scan within 1 of its top from 12.1 to 12.2 | what: a noise correlated from sample to sample, excess 19.3 · changes along: the asymmetry along pressure, excess 13.8 · changes along: the core along pressure, excess 9.93 · its line 3.01 | 3.9895 ± 0.0032, not yet to be quoted | 0.08007 ± 0.00061, as second component at 12.136, not yet to be quoted |

## saturation: laser-driven resonance fluorescence at five powers

`python -m twinfit example saturation DIR`. The truth, which the folder does not hold: full width 1.0 MHz times the root of 1 + P / 50 mW, the broadening by saturation, and height 1.5 P / (P + 50 mW), place 0.3 MHz in every record, background 0.05, three glitches of +0.08 in the record at 40 mW, at -9.3, 2.1 and 7.7 MHz, noise 0.004.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 0 |
| its terms | one Lorentzian line of one width | one Lorentzian line whose width saturates with the power | one Lorentzian line of one width |
| chi square for each degree of freedom | 12.573 | 1.150 | 12.573 |
| the last report, first beyond its twin | changes along: the flanks along power, excess 1.17e+04 · changes along: the core along power, excess 8.97e+03 · where: the core, excess 5.08e+03 · its line 3.08 | what: glitches, excess 111 · its line 2.69 | changes along: the flanks along power, excess 1.17e+04 · changes along: the core along power, excess 8.97e+03 · where: the core, excess 5.08e+03 · its line 3.08 |
| width at no power, true 1 | not in the model | 0.9969 ± 0.0028, not yet to be quoted | not in the model |
| saturation power, true 50 | not in the model | 49.36 ± 0.42, not yet to be quoted | not in the model |

The first report proposes for the model: the shape that the residuals hold in the core · the shape that the residuals hold in the flanks, changing along power · the shape that the residuals hold in the core, changing along power. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.62 · more repeats, a better data collection: the error of the known parameters times 0.67.

**What to do**: write the law that the physics gives, the width saturating with the power, as shared parameters that read the power, in model.py.

The `model.py` that it writes in the place of the first:

```python
"""A model of one's own: the width of the line grows with the power by the law of saturation, its two numbers shared."""
import numpy as np

SHARED = ("width at no power", "saturation power")      # the full width at half maximum at no power, in MHz, and the power of saturation, in mW
OWN = ("height", "place", "background")                 # of each record
START_SHARED = [1.5, 30.0]
BOUNDS_SHARED = ([0.05, 1.0], [20.0, 1000.0])
CENTRE = "place"
NAME = "one Lorentzian line whose width saturates with the power"
EQUATION = r"$\frac{h}{1 + (2\,(x - x_0)/w)^2} + b,\quad w = w_0\,\sqrt{1 + P/P_s}$"          # the model in LaTeX, which the visual log typesets
LETTERS = {"h": "the height", "x_0": "the place", "w_0": "the width at no power", "P": "the power", "P_s": "the power of saturation", "b": "the background"}          # each letter of EQUATION, which the visual log writes under it


def model(x, shared, own, design):
    width = shared[0] * np.sqrt(1.0 + design["power"] / shared[1])
    return own[0] / (1.0 + (2.0 * (x - own[1]) / width) ** 2) + own[2]


def start_own(record):
    base = float(np.median(record.y))
    return [float(record.y.max() - base), float(record.x[np.argmax(record.y)]), base]


def centre(shared, own, design):
    return own[1]
```

**What not to do**: a width that grows in proportion to the power, the generic term where the physics gives a law. The loop refused it on 30 synthetic sets: width changes along power, by a hair, which the loop itself asks to run again at 500 sets before the term is set aside: at more sets it may be kept, and that a wrong term is kept says that the noise check and the chi square, not the twin alone, judge it.

**At other noises**: the first model, then what to do. Its own noise is white, and its second kind correlated over neighbouring samples. The rows are written by `--noise 3`, `--noise 0`, `--kind correlated` and `--noise 3 --kind correlated`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | width at no power, true 1 | saturation power, true 50 |
|---|---|---|---|---|---|
| white, 3 times as large as its own noise | changes along: the flanks along power, excess 1.24e+03 · changes along: the core along power, excess 1.06e+03 · where: the core, excess 606 · its line 3 | model 0, 1.008, `the model.py of the gallery` | what: glitches, excess 8.16 · its line 2.65 | 0.9894 ± 0.0086, not yet to be quoted | 47.9 ± 1.2, not yet to be quoted |
| no noise drawn | changes along: the flanks along power, excess 1.15e+04 · changes along: the core along power, excess 9.16e+03 · where: the core, excess 5.19e+03 · its line 2.88 | model 0, 0.160, `the model.py of the gallery` | what: glitches, excess 118 · its line 2.64 | 1.0007 ± 0.0028, not yet to be quoted | 50.09 ± 0.43, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | changes along: the core along power, excess 9.62e+03 · changes along: the flanks along power, excess 9.45e+03 · where: the core, excess 5.15e+03 · its line 2.74 | model 0, 1.156, `the model.py of the gallery` | what: glitches, excess 115 · what: a noise correlated from sample to sample, excess 25.7 · where: the flanks, excess 8.55 · its line 2.64 | 0.9892 ± 0.0029, not yet to be quoted | 47.92 ± 0.41, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | changes along: the flanks along power, excess 1.12e+03 · changes along: the core along power, excess 1.08e+03 · where: the core, excess 602 · its line 2.48 | model 0, 1.006, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 27.1 · what: glitches, excess 7.82 · changes along: the asymmetry along power, excess 6.55 · its line 2.59 | 0.9663 ± 0.0087, not yet to be quoted | 43.9 ± 1.1, not yet to be quoted |

## filter: absorption dips behind a detector's low-pass filter

`python -m twinfit example filter DIR`. The truth, which the folder does not hold: the model is whole, width 1.2 pm, and the noise is correlated from sample to sample, 0.8 from one to the next, of spread 0.002 as stated.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 0 |
| its terms | one Gaussian dip of one width | one Gaussian dip of one width | one Gaussian dip of one width |
| chi square for each degree of freedom | 1.028 | 1.051 | 1.028 |
| the last report, first beyond its twin | what: a noise correlated from sample to sample, excess 73.3 · changes along: the core along temperature, excess 26.8 · where: the flanks, excess 12.8 · its line 3.28 | nothing | what: a noise correlated from sample to sample, excess 73.3 · changes along: the core along temperature, excess 26.8 · where: the flanks, excess 12.8 · its line 3.28 |
| width, true 1.2 | 1.20024 ± 0.00033, not yet to be quoted | 1.2007 ± 0.0013 | 1.20024 ± 0.00033, not yet to be quoted |

The first report proposes for the model: the shape that the residuals hold in the flanks · a periodic structure of 21 cycles over the span · the shape that the residuals hold in the core, changing along temperature. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.70 · more repeats, a better data collection: the error of the known parameters times 0.71.

**What to do**: read the correlated noise as noise: keep one sample in every few with `python -m twinfit thin DIR N`, at the span that the report gives, and run the loop on the copy.

Here `python -m twinfit thin DIR 15` writes 3 copies, and the loop runs on the first, do_thin_15_0.

**What not to do**: hold the correlation with a term of the curve, a pedestal. The loop refused it on 30 synthetic sets: a pedestal 6 times wider than the curve.

**At other noises**: the first model, then what to do. Its own noise is correlated over neighbouring samples, and its second kind white. The rows are written by `--noise 3`, `--noise 0`, `--kind white` and `--noise 3 --kind white`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | width, true 1.2 |
|---|---|---|---|---|
| correlated over neighbouring samples, 3 times as large as its own noise | what: a noise correlated from sample to sample, excess 71.8 · changes along: the core along temperature, excess 25.8 · where: the flanks, excess 13.7 · its line 2.86 | model 0, 1.051, `python -m twinfit thin DIR 15` | nothing | 1.2021 ± 0.0038 |
| no noise drawn | nothing | model 0, 0.000, `python -m twinfit thin DIR 15` | nothing | 1.2000 ± 0.0013, not yet to be quoted |
| white, as large as its own noise | changes along: the core along temperature, excess 3.05 · its line 3 | model 0, 1.051, `python -m twinfit thin DIR 15` | nothing | 1.2000 ± 0.0013 |
| white, 3 times as large as its own noise | nothing | model 0, 1.051, `python -m twinfit thin DIR 15` | nothing | 1.2001 ± 0.0038 |

## fano: photoionization across an autoionizing resonance at four pressures

`python -m twinfit example fano DIR`. The truth, which the folder does not hold: a Fano profile, q = 3 and half width 2 meV at every pressure, its top near one, and the first model is a Lorentzian, which lacks the interference.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 0 |
| its terms | one Lorentzian line of one half width | one Fano profile | one Lorentzian line of one half width |
| chi square for each degree of freedom | 475.375 | 1.001 | 475.375 |
| the last report, first beyond its twin | where: the asymmetry, excess 5.84e+05 · changes along: the asymmetry along pressure, excess 3.29e+04 · where: the core, excess 2.29e+04 · its line 3.18 | nothing | where: the asymmetry, excess 5.84e+05 · changes along: the asymmetry along pressure, excess 3.29e+04 · where: the core, excess 2.29e+04 · its line 3.18 |
| half_width, true 2 | 2.2254 ± 0.0011, not yet to be quoted | 1.99966 ± 0.00088 | 2.2254 ± 0.0011, not yet to be quoted |
| q, true 3 | not in the model | 2.9971 ± 0.0017 | not in the model |

The first report proposes for the model: the shape that the residuals hold in the asymmetry · a slow trend over the span of each curve · the shape that the residuals hold in the asymmetry, changing along pressure. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more repeats, a better data collection: the error of the known parameters times 0.64 · more samples, a better data collection: the error of the known parameters times 0.65.

**What to do**: the asymmetry is the line's own shape, which no term of the package holds: write the profile that the physics gives, a Fano profile, as the model.

The `model.py` that it writes in the place of the first:

```python
"""A model of one's own: a Fano profile, the line of a discrete state that interferes with a continuum, its half width
and its asymmetry q shared."""
import numpy as np

SHARED = ("half_width", "q")                 # the half width of the resonance, in meV, and the asymmetry of the profile
OWN = ("height", "place", "background")      # of each record
START_SHARED = [1.5, 2.0]
BOUNDS_SHARED = ([0.05, 0.1], [20.0, 50.0])
CENTRE = "place"
NAME = "one Fano profile"
EQUATION = r"$h\,\frac{(q + \epsilon)^2}{(1 + \epsilon^2)(1 + q^2)} + b,\quad \epsilon = (x - x_0)/\Gamma$"          # the model in LaTeX, which the visual log typesets
LETTERS = {"h": "the height", "q": "the asymmetry", "x_0": "the place", r"\Gamma": "the half width", "b": "the background"}          # each letter of EQUATION, which the visual log writes under it


def model(x, shared, own, design):
    e = (x - own[1]) / shared[0]
    q = shared[1]
    return own[0] * ((q + e) ** 2 / (1.0 + e ** 2)) / (q * q + 1.0) + own[2]


def start_own(record):
    base = float(np.median(record.y))
    top = int(np.argmax(np.abs(record.y - base)))
    return [float(record.y[top] - base), float(record.x[top]), base]


def centre(shared, own, design):
    return own[1]
```

**What not to do**: a tail on one side, the generic term of an asymmetry. The loop refused it on 30 synthetic sets: a tail on the upper side, 2 times wider.

**At other noises**: the first model, then what to do. Its own noise is white, and its second kind correlated over neighbouring samples. The rows are written by `--noise 3`, `--noise 0`, `--kind correlated` and `--noise 3 --kind correlated`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | half_width, true 2 | q, true 3 |
|---|---|---|---|---|---|
| white, 3 times as large as its own noise | where: the asymmetry, excess 6.53e+04 · changes along: the asymmetry along pressure, excess 3.3e+03 · where: the core, excess 2.42e+03 · its line 3.18 | model 0, 1.001, `the model.py of the gallery` | nothing | 1.9990 ± 0.0026 | 2.9912 ± 0.0051 |
| no noise drawn | where: the asymmetry, excess 5.68e+05 · changes along: the asymmetry along pressure, excess 3.26e+04 · where: the core, excess 2.23e+04 · its line 3.14 | model 0, 0.000, `the model.py of the gallery` | nothing | 2.00000 ± 0.00088, not yet to be quoted | 3.0000 ± 0.0017, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | where: the asymmetry, excess 5.8e+05 · changes along: the asymmetry along pressure, excess 3.39e+04 · where: the core, excess 2.33e+04 · its line 3.06 | model 0, 1.000, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 69.1 · where: the asymmetry, excess 13 · changes along: the asymmetry along pressure, excess 5.1 · its line 2.64 | 1.99870 ± 0.00088, not yet to be quoted | 2.9911 ± 0.0017, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | where: the asymmetry, excess 6.27e+04 · changes along: the asymmetry along pressure, excess 3.7e+03 · where: the core, excess 2.42e+03 · its line 3.12 | model 0, 1.000, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 71 · where: the asymmetry, excess 10.9 · changes along: the asymmetry along pressure, excess 4.86 · its line 2.75 | 1.9961 ± 0.0026, not yet to be quoted | 2.9735 ± 0.0050, not yet to be quoted |

## counts: emission-line photon counts at four pump powers

`python -m twinfit example counts DIR`. The truth, which the folder does not hold: the model is whole, half width 1.5 GHz, and the noise is that of a counter, the root of the counts, 6.3 on the background and 40 at the top of the strongest line.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 0 |
| its terms | one Lorentzian line of one half width | one Lorentzian line of one half width | one Lorentzian line of one half width |
| chi square for each degree of freedom | 4.509 | 0.984 | 4.509 |
| the last report, first beyond its twin | what: a noise that grows with the signal, excess 347 · what: glitches, excess 66.3 · changes along: the core along power, excess 26.3 · its line 2.37 | nothing | what: a noise that grows with the signal, excess 347 · what: glitches, excess 66.3 · changes along: the core along power, excess 26.3 · its line 2.37 |
| half_width, true 1.5 | 1.50179 ± 0.00094, not yet to be quoted | 1.5009 ± 0.0032 | 1.50179 ± 0.00094, not yet to be quoted |

The first report proposes for the model: a periodic structure of 272.5 cycles over the span · the shape that the residuals hold in the core, changing along power · the shape that the residuals hold in the flanks, changing along power. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.71 · more repeats, a better data collection: the error of the known parameters times 0.75.

**What to do**: give each sample the noise of a counter, the root of the fitted counts, with `python -m twinfit counts DIR`, and run the loop on the copy.

Here `python -m twinfit counts DIR` writes the copy, and the loop runs on it.

**What not to do**: hold the misfit at the top of the lines with a term of the curve, a pedestal. The loop refused it on 30 synthetic sets: a pedestal 6 times wider than the curve.

**At other noises**: the first model, then what to do. Its own noise is that of a counter, and its second kind white. The rows are written by `--noise 3`, `--noise 0`, `--kind white` and `--noise 3 --kind white`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | half_width, true 1.5 |
|---|---|---|---|---|
| that of a counter, 3 times as large as its own noise | what: a noise that grows with the signal, excess 320 · what: glitches, excess 106 · changes along: the core along power, excess 20.7 · its line 2.45 | model 0, 8.853, `python -m twinfit counts DIR` | what: glitches, excess 30.8 · changes along: the core along power, excess 15.1 · where: the asymmetry, excess 13.8 · its line 2.69 | 1.5030 ± 0.0032, not yet to be quoted |
| no noise drawn | nothing | model 0, 0.000, `python -m twinfit counts DIR` | nothing | 1.5000 ± 0.0032, not yet to be quoted |
| white, as large as its own noise | nothing | model 0, 0.458, `python -m twinfit counts DIR` | what: a noise that grows with the signal, excess 8.27 · its line 2.49 | 1.5008 ± 0.0032, not yet to be quoted |
| white, 3 times as large as its own noise | nothing | model 0, 4.127, `python -m twinfit counts DIR` | what: a noise that grows with the signal, excess 77.4 · what: glitches, excess 20.7 · where: the asymmetry, excess 4.76 · its line 2.56 | 1.5023 ± 0.0032, not yet to be quoted |

## clock: atomic clock resonance at five probe powers

`python -m twinfit example clock DIR`. The truth, which the folder does not hold: the center moves with the power of the probe, 0.12 Hz from the nominal frequency at no power and 0.025 Hz for each microwatt, the light shift, and the width is 1.5 Hz at every power.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 0 |
| its terms | one Lorentzian line of one width and one centre | one Lorentzian line whose centre moves with the power | one Lorentzian line of one width, its centre its own in each curve |
| chi square for each degree of freedom | 31.812 | 0.972 | 0.972 |
| the last report, first beyond its twin | changes along: the asymmetry along power, excess 1.84e+04 · changes along: the core along power, excess 1.49e+04 · changes along: the flanks along power, excess 1.61e+03 · its line 2.47 | nothing | nothing |
| width, true 1.5 | 1.51614 ± 0.00058, not yet to be quoted | 1.49970 ± 0.00058 | 1.49970 ± 0.00058 |
| clock frequency, true 0.12 | not in the model | 0.12000 ± 0.00038 | not in the model |
| light shift, true 0.025 | not in the model | 0.025051 ± 0.000058 ± 0.00050 of theory | not in the model |

The first report proposes for the model: the shape that the residuals hold in the flanks, changing along power · the shape that the residuals hold in the core, changing along power · the shape that the residuals hold in the asymmetry, changing along power. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.66 · more repeats, a better data collection: the error of the known parameters times 0.68.

**What to do**: write the light shift that the physics gives, the center moving in proportion to the power read through its calibration, which the description holds within 2 per cent, so that the model measures the frequency at no power and the light shift carries the bar of that calibration.

The `model.py` that it writes in the place of the first:

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

**What not to do**: free the center of each curve: the fit is right, and the frequency of the clock is then in no parameter, the shift hidden in five centers.

**At other noises**: the first model, then what to do. Its own noise is white, and its second kind correlated over neighbouring samples. The rows are written by `--noise 3`, `--noise 0`, `--kind correlated` and `--noise 3 --kind correlated`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | width, true 1.5 | clock frequency, true 0.12 | light shift, true 0.025 |
|---|---|---|---|---|---|---|
| white, 3 times as large as its own noise | changes along: the asymmetry along power, excess 2.04e+03 · changes along: the core along power, excess 1.7e+03 · changes along: the flanks along power, excess 178 · its line 2.28 | model 0, 0.972, `the model.py of the gallery` | nothing | 1.4991 ± 0.0017 | 0.1200 ± 0.0012 | 0.02515 ± 0.00017 ± 0.00050 of theory |
| no noise drawn | changes along: the asymmetry along power, excess 1.84e+04 · changes along: the core along power, excess 1.46e+04 · changes along: the flanks along power, excess 1.62e+03 · its line 2.47 | model 0, 0.000, `the model.py of the gallery` | nothing | 1.50000 ± 0.00058, not yet to be quoted | 0.12000 ± 0.00038, not yet to be quoted | 0.025000 ± 0.000058 ± 0.00050 of theory, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | changes along: the asymmetry along power, excess 1.83e+04 · changes along: the core along power, excess 1.53e+04 · changes along: the flanks along power, excess 1.6e+03 · its line 2.48 | model 0, 0.930, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 28.7 · changes along: the core along power, excess 14.8 · where: the asymmetry, excess 8.66 · its line 3.79 | 1.49900 ± 0.00058, not yet to be quoted | 0.12006 ± 0.00038, not yet to be quoted | 0.025138 ± 0.000058 ± 0.00050 of theory, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | changes along: the asymmetry along power, excess 2.04e+03 · changes along: the core along power, excess 1.87e+03 · changes along: the flanks along power, excess 168 · its line 2.21 | model 0, 0.930, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 28.7 · changes along: the core along power, excess 14.7 · where: the asymmetry, excess 8.66 · its line 3.78 | 1.4970 ± 0.0017, not yet to be quoted | 0.1202 ± 0.0012, not yet to be quoted | 0.02542 ± 0.00017 ± 0.00051 of theory, not yet to be quoted |

## gel: gel electrophoresis band at four voltages

`python -m twinfit example gel DIR`. The truth, which the folder does not hold: a tail on the side of slower migration, a copy of the band three times wider below its place, of size 0.15 as the tail of the package counts it, beside a Gaussian band of standard deviation 1 mm.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 1 | model 1 |
| its terms | one Gaussian band of one width | one Gaussian band of one width, a tail on the lower side, 3.03 times wider | one Gaussian band of one width, a pedestal 3 times wider than the curve |
| chi square for each degree of freedom | 7.304 | 0.999 | 5.371 |
| the last report, first beyond its twin | where: the flanks, excess 3.01e+03 · where: the asymmetry, excess 3e+03 · where: the wings, excess 544 · its line 2.48 | nothing | where: the asymmetry, excess 2.53e+03 · changes along: the asymmetry along voltage, excess 58.3 · what: a noise that grows with the signal, excess 56.7 · its line 2.85 |
| width, true 1 | 1.02129 ± 0.00024, not yet to be quoted | 1.00023 ± 0.00027 | 0.99926 ± 0.00033, not yet to be quoted |
| tail lower 3, true 0.15 | not in the model | 0.14918 ± 0.00086, as tail lower 3.03 | not in the model |

The first report proposes for the model: the shape that the residuals hold in the flanks · the shape that the residuals hold in the asymmetry · a slow trend over the span of each curve. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.62 · voltage past its range, a better data collection: the error of the known parameters times 0.81.

**What to do**: the report puts the misfit in the flanks and the asymmetry: in a gel a band trails on the side of slower migration, a tail, so scan its stretch and write the term at the top of the scan.

Here `python -m twinfit scan DIR "tail(x, 'lower')" 2 3 4 5 6` puts its top at 3.29, and a finer scan about it, at 3.20106 3.23071 3.26035 3.29 3.31965 3.34929 3.37894, at 3.201, and a finer scan about it, at 3.11212 3.14177 3.17142 3.20106 3.23071 3.26035 3.29, at 3.112, and a finer scan about it, at 3.02318 3.05283 3.08248 3.11212 3.14177 3.17141 3.20106, at 3.023, and a finer scan about it, at 2.93424 2.96389 2.99354 3.02318 3.05283 3.08248 3.11212, at 3.03: the term file holds `term = mechanisms.tail(3.03, 'lower')`.

**What not to do**: a pedestal under the band, wider on both sides alike.

**At other noises**: the first model, then what to do. Its own noise is white, and its second kind correlated over neighbouring samples. The rows are written by `--noise 3`, `--noise 0`, `--kind correlated` and `--noise 3 --kind correlated`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | width, true 1 | tail lower 3, true 0.15 |
|---|---|---|---|---|---|
| white, 3 times as large as its own noise | where: the flanks, excess 425 · where: the asymmetry, excess 310 · where: the wings, excess 65.3 · its line 2.41 | model 1, 0.999, `term = mechanisms.tail(3.06464, 'lower')`, its scan within 1 of its top from 3.06 to 3.15 | nothing | 1.00054 ± 0.00080 | 0.1478 ± 0.0026, as tail lower 3.06464 |
| no noise drawn | where: the flanks, excess 3.2e+03 · where: the asymmetry, excess 3e+03 · where: the wings, excess 476 · its line 2.27 | model 1, 0.000, `term = mechanisms.tail(3, 'lower')`, its scan within 1 of its top from 2.98 to 3.02 | nothing | 1.00000 ± 0.00027, not yet to be quoted | 0.15000 ± 0.00086, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | where: the flanks, excess 3.38e+03 · where: the asymmetry, excess 2.52e+03 · where: the wings, excess 577 · its line 2.4 | model 1, 0.986, `term = mechanisms.tail(3.101, 'lower')`, its scan within 1 of its top from 3.08 to 3.12 | what: a noise correlated from sample to sample, excess 27.9 · changes along: the core along voltage, excess 13.6 · changes along: the asymmetry along voltage, excess 5.18 · its line 2.67 | 1.00091 ± 0.00027, not yet to be quoted | 0.14717 ± 0.00086, as tail lower 3.101, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | where: the flanks, excess 411 · where: the asymmetry, excess 273 · where: the wings, excess 76 · its line 2.35 | model 1, 0.986, `term = mechanisms.tail(3.37339, 'lower')`, its scan within 1 of its top from 3.28 to 3.37 | what: a noise correlated from sample to sample, excess 28.9 · changes along: the core along voltage, excess 13.4 · changes along: the asymmetry along voltage, excess 5.27 · its line 2.34 | 1.00298 ± 0.00079, not yet to be quoted | 0.1423 ± 0.0025, as tail lower 3.37339, not yet to be quoted |

## transit: transit of an exoplanet in four filters

`python -m twinfit example transit DIR`. The truth, which the folder does not hold: a star darker at its limb, by a coefficient u = 0.6 - 0.3 (wavelength - 450 nm) / 350 nm, 0.6 in the blue and 0.3 at 800 nm, under a planet of depth 0.010 and half duration 1.2 h.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 1 |
| its terms | the dip of a uniform star | the dip of a star darker at its limb | the dip of a uniform star, a power 2 of the abscissa under every curve |
| chi square for each degree of freedom | 2.397 | 0.984 | 2.342 |
| the last report, first beyond its twin | where: the core, excess 1.28e+03 · changes along: the core along wavelength, excess 164 · what: a noise correlated from sample to sample, excess 22.6 · its line 2.74 | nothing | where: the core, excess 1.03e+03 · changes along: the core along wavelength, excess 149 · where: the flanks, excess 29.6 · its line 2.62 |
| half_duration, true 1.2 | 1.14558 ± 0.00052, not yet to be quoted | 1.1997 ± 0.0011 | 1.14304 ± 0.00055, not yet to be quoted |
| limb darkening at 450 nm, true 0.6 | not in the model | 0.5967 ± 0.0072 | not in the model |
| its fall to 800 nm, true 0.3 | not in the model | 0.3091 ± 0.0098 | not in the model |

The first report proposes for the model: the shape that the residuals hold in the core · the shape that the residuals hold in the core, changing along wavelength · a periodic structure of 3.5 cycles over the span. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.75 · more repeats, a better data collection: the error of the known parameters times 0.77.

**What to do**: write the darkening of the limb of the star that the physics gives, its coefficient falling with the wavelength, as the model.

The `model.py` that it writes in the place of the first:

```python
"""A model of one's own: the dip of a star darker at its limb, by a coefficient that falls in proportion to the wavelength."""
import numpy as np

SHARED = ("half_duration", "limb darkening at 450 nm", "its fall to 800 nm")      # hours, and the coefficient u of the limb and its fall from 450 to 800 nm
OWN = ("depth", "centre", "baseline")                                            # of each record
START_SHARED = [1.0, 0.4, 0.1]
BOUNDS_SHARED = ([0.2, -1.0, -1.0], [2.5, 1.5, 1.5])
CENTRE = "centre"
NAME = "the dip of a star darker at its limb"
EQUATION = r"$b - d\,\frac{1 - u\,(1 - \mu)}{1 - u/3},\quad u = u_{450} - \Delta u\,\frac{\lambda - 450}{350}$"          # the model in LaTeX, which the visual log typesets
LETTERS = {"b": "the baseline", "d": "the depth", "u": "the darkening of the limb", r"\mu": "the cosine at the limb", r"\lambda": "the wavelength, nm"}          # each letter of EQUATION, which the visual log writes under it
EDGE = 0.05                                 # the time over which the planet enters the disc of the star, in hours


def model(x, shared, own, design):
    u = shared[1] - shared[2] * (design["wavelength"] - 450.0) / 350.0
    t = x - own[1]
    inside = 1.0 / (1.0 + np.exp(-(shared[0] - np.abs(t)) / EDGE))
    rim = 1.0 - (t / shared[0]) ** 2
    mu = np.sqrt(0.5 * (rim + np.sqrt(rim ** 2 + 1e-4)))          # the cosine of the angle at the limb, smooth at the edge of the star, where a root clipped at zero has a corner that a search cannot cross
    return own[2] - own[0] * inside * (1.0 - u * (1.0 - mu)) / (1.0 - u / 3.0)


def start_own(record):
    base = float(np.median(record.y))
    return [float(base - record.y.min()), float(record.x[np.argmin(record.y)]), base]


def centre(shared, own, design):
    return own[1]
```

**What not to do**: a curvature of the baseline under every curve, the generic term of a round bottom.

**At other noises**: the first model, then what to do. Its own noise is white, and its second kind correlated over neighbouring samples. The rows are written by `--noise 3`, `--noise 0`, `--kind correlated` and `--noise 3 --kind correlated`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | half_duration, true 1.2 | limb darkening at 450 nm, true 0.6 | its fall to 800 nm, true 0.3 |
|---|---|---|---|---|---|---|
| white, 3 times as large as its own noise | where: the core, excess 131 · changes along: the core along wavelength, excess 23.2 · what: a noise correlated from sample to sample, excess 3.19 · its line 2.62 | model 0, 0.984, `the model.py of the gallery` | nothing | 1.2008 ± 0.0032 | 0.591 ± 0.022 | 0.324 ± 0.029 |
| no noise drawn | where: the core, excess 1.23e+03 · changes along: the core along wavelength, excess 150 · what: a noise correlated from sample to sample, excess 18.5 · its line 2.61 | model 0, 0.001, `the model.py of the gallery` | nothing | 1.1992 ± 0.0011, not yet to be quoted | 0.5993 ± 0.0073, not yet to be quoted | 0.3018 ± 0.0098, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | where: the core, excess 1.06e+03 · changes along: the core along wavelength, excess 220 · what: a noise correlated from sample to sample, excess 63.6 · its line 2.61 | model 0, 0.997, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 45.5 · where: the core, excess 3.41 · its line 2.82 | 1.2000 ± 0.0011, not yet to be quoted | 0.5869 ± 0.0072, not yet to be quoted | 0.3199 ± 0.0098, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | where: the core, excess 138 · what: a noise correlated from sample to sample, excess 45.9 · changes along: the core along wavelength, excess 35.6 · its line 2.49 | model 0, 0.996, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 48.2 · its line 2.53 | 1.2014 ± 0.0030, not yet to be quoted | 0.561 ± 0.021, not yet to be quoted | 0.358 ± 0.030, not yet to be quoted |

## fringes: piezo-scanned Michelson interferometer at five wavelengths

`python -m twinfit example fringes DIR`. The truth, which the folder does not hold: the piezo moves the mirror 1.012 nm for each nm of the scan, and its displacement is bent by a square of the position, 1.012 (x + 0.003 x^2 / 2000 nm), 6 nm off a linear scan at both ends of the scan.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 0 |
| its terms | one cosine fringe of period half the wavelength | one cosine fringe of period half the wavelength, the scan bent | one cosine fringe of period half the wavelength |
| chi square for each degree of freedom | 2.425 | 1.008 | 2.425 |
| the last report, first beyond its twin | what: a noise correlated from sample to sample, excess 44.8 · where: the asymmetry, excess 14.6 · what: glitches, excess 6.94 · its line 2.42 | nothing | what: a noise correlated from sample to sample, excess 44.8 · where: the asymmetry, excess 14.6 · what: glitches, excess 6.94 · its line 2.42 |
| piezo scale, true 1.012 | 1.012018 ± 0.000030, not yet to be quoted | 1.012020 ± 0.000030 | 1.012018 ± 0.000030, not yet to be quoted |
| scan bend, true 0.003 | not in the model | 0.003001 ± 0.000057 | not in the model |

The first report proposes for the model: a periodic structure of 16.5 cycles over the span · the shape that the residuals hold in the asymmetry, changing along wavelength · the shape that the residuals hold in the asymmetry. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · wavelength past its range, a better data collection: the error of the known parameters times 0.56 · a wider span, a better data collection: the error of the known parameters times 0.81.

**What to do**: read the first answers, a noise correlated from sample to sample with the asymmetry, as the smooth ripple that a fringe leaves where its phase drifts, and write the physics of the instrument that makes it drift, a piezo that does not move the mirror in proportion to the scan: its displacement bent by a square of the position, as the model.

The `model.py` that it writes in the place of the first:

```python
"""A model of one's own: a cosine fringe whose mirror displacement is the scan times a calibration of the piezo, bent by a square of the position."""
import numpy as np

SHARED = ("piezo scale", "scan bend")      # the displacement of the mirror for each nm of the scan, and the bend of the response: the displacement is scale (x + bend x^2 / HALF)
OWN = ("mean", "visibility", "phase")       # of each record
START_SHARED = [1.0, 0.0]
BOUNDS_SHARED = ([0.9, -0.05], [1.1, 0.05])
NAME = "one cosine fringe of period half the wavelength, the scan bent"
EQUATION = r"$m\,(1 + v\cos(4\pi d/\lambda + \phi)),\quad d = k\,(x + c\,x^2/H)$"          # the model in LaTeX, which the visual log typesets
LETTERS = {"m": "the mean", "v": "the visibility", "d": "the displacement", r"\lambda": "the wavelength", r"\phi": "the phase", "k": "the scale of the piezo", "c": "the bend", "H": "half the scan"}          # each letter of EQUATION, which the visual log writes under it
HALF = 2000.0                               # half the span of the scan, in nm


def model(x, shared, own, design):
    displacement = shared[0] * (x + shared[1] * x ** 2 / HALF)
    return own[0] * (1.0 + own[1] * np.cos(4.0 * np.pi * displacement / design["wavelength"] + own[2]))


def start_own(record):
    k = 4.0 * np.pi / record.design["wavelength"]
    basis = np.column_stack([np.ones_like(record.x), np.cos(k * record.x), np.sin(k * record.x)])
    a, b, c = np.linalg.lstsq(basis, record.y, rcond=None)[0]
    return [float(a), float(np.hypot(b, c) / a), float(np.arctan2(-c, b))]
```

**What not to do**: let the calibration of the piezo change with the wavelength, a setting that the report does not name. The loop refused it on 30 synthetic sets: piezo scale changes along wavelength.

**At other noises**: the first model, then what to do. Its own noise is white, and its second kind correlated over neighbouring samples. The rows are written by `--noise 3`, `--noise 0`, `--kind correlated` and `--noise 3 --kind correlated`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | piezo scale, true 1.012 | scan bend, true 0.003 |
|---|---|---|---|---|---|
| white, 3 times as large as its own noise | where: the asymmetry, excess 4.44 · what: a noise correlated from sample to sample, excess 3.95 · its line 2.42 | model 0, 1.008, `the model.py of the gallery` | nothing | 1.012059 ± 0.000089 | 0.00300 ± 0.00017 |
| no noise drawn | what: a noise correlated from sample to sample, excess 46.5 · where: the asymmetry, excess 5.41 · where: the flanks, excess 2.49 · its line 2.42 | model 0, 0.000, `the model.py of the gallery` | nothing | 1.012000 ± 0.000030, not yet to be quoted | 0.003000 ± 0.000057, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | what: a noise correlated from sample to sample, excess 65.2 · where: the asymmetry, excess 9.83 · changes along: the asymmetry along wavelength, excess 8.75 · its line 2.42 | model 0, 0.944, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 22.4 · where: the core, excess 3.85 · changes along: the asymmetry along wavelength, excess 3.42 · its line 2.34 | 1.012003 ± 0.000030, not yet to be quoted | 0.002995 ± 0.000056, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | what: a noise correlated from sample to sample, excess 25.8 · changes along: the asymmetry along wavelength, excess 4.77 · where: the core, excess 4.14 · its line 2.42 | model 0, 0.944, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 22.4 · where: the core, excess 3.85 · changes along: the asymmetry along wavelength, excess 3.42 · its line 2.39 | 1.012012 ± 0.000089, not yet to be quoted | 0.00298 ± 0.00017, not yet to be quoted |

## returns: daily stock index returns in five volatility periods

`python -m twinfit example returns DIR`. The truth, which the folder does not hold: the returns follow a Student law whose tail index, its degrees of freedom, is 4: fat tails, with a standard deviation 0.9 of the quoted volatility of the period.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 0 |
| its terms | one normal law of the returns | a Student law of the returns | one normal law of the returns |
| chi square for each degree of freedom | 3.723 | 0.837 | 3.723 |
| the last report, first beyond its twin | where: the wings, excess 47.5 · changes along: the flanks along volatility, excess 40.2 · where: the core, excess 37.9 · its line 2.82 | nothing | where: the wings, excess 47.5 · changes along: the flanks along volatility, excess 40.2 · where: the core, excess 37.9 · its line 2.82 |
| volatility ratio, true 0.9 | 0.7407 ± 0.0038, not yet to be quoted | 0.898 ± 0.012 | 0.7407 ± 0.0038, not yet to be quoted |
| tail index, true 4 | not in the model | 4.03 ± 0.14 | not in the model |

The first report proposes for the model: the shape that the residuals hold in the wings · the shape that the residuals hold in the core · the shape that the residuals hold in the flanks, changing along volatility. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.52 · more repeats, a better data collection: the error of the known parameters times 0.61.

**What to do**: write the law that the market gives for returns, a Student law whose tail index, its degrees of freedom, is shared by the periods and whose standard deviation stays a fraction of the quoted volatility, as the model.

The `model.py` that it writes in the place of the first:

```python
"""A model of one's own: a Student law of the daily returns, the fat tails of a market, its standard deviation a fraction of the quoted volatility and its tail index shared by all the periods."""
import numpy as np
from scipy.special import gammaln

SHARED = ("volatility ratio", "tail index")             # the standard deviation of the returns over the volatility that the options quoted, and the tail index of the law, its degrees of freedom: the fewer, the fatter the tails
OWN = ("returns", "mean")                               # of each record
START_SHARED = [0.9, 8.0]
BOUNDS_SHARED = ([0.2, 2.5], [3.0, 200.0])
CENTRE = "mean"
NAME = "a Student law of the returns"
EQUATION = r"$N\,\frac{\Gamma((\nu + 1)/2)}{\Gamma(\nu/2)\,\sqrt{\nu\pi}\,s}\,(1 + t^2/\nu)^{-(\nu + 1)/2},\quad t = (x - m)/s$"          # the model in LaTeX, which the visual log typesets
LETTERS = {"N": "the count", r"\nu": "the tail index", "s": "the scale", "m": "the mean"}          # each letter of EQUATION, which the visual log writes under it
BIN = 0.2                                               # the width of a bin of the histogram, in per cent


def model(x, shared, own, design):
    nu = shared[1]
    scale = shared[0] * design["volatility"] / np.sqrt(nu / (nu - 2.0))          # the scale that gives the law the standard deviation of the model
    t = (x - own[1]) / scale
    log = gammaln(0.5 * (nu + 1.0)) - gammaln(0.5 * nu) - 0.5 * np.log(nu * np.pi) - 0.5 * (nu + 1.0) * np.log1p(t ** 2 / nu)
    return own[0] * BIN * np.exp(log) / scale


def start_own(record):
    return [float(record.y.sum()), float(np.average(record.x, weights=np.maximum(record.y, 0.0)))]


def centre(shared, own, design):
    return own[1]
```

**What not to do**: let the ratio of the standard deviation to the quoted volatility change along the volatility, the generic term where the law of the market gives the shape of every period. The loop refused it on 30 synthetic sets: volatility ratio changes along volatility.

**At other noises**: the first model, then what to do. Its own noise is that of a counter, and its second kind white. The rows are written by `--noise 3`, `--noise 0`, `--kind white` and `--noise 3 --kind white`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | volatility ratio, true 0.9 | tail index, true 4 |
|---|---|---|---|---|---|
| that of a counter, 3 times as large as its own noise | where: the wings, excess 11.2 · changes along: the wings along volatility, excess 4.69 · changes along: the flanks along volatility, excess 4.4 · its line 3.95 | model 0, 0.979, `the model.py of the gallery` | nothing | 0.943 ± 0.042 | 3.68 ± 0.37 |
| no noise drawn | where: the wings, excess 47.9 · where: the core, excess 36.3 · changes along: the wings along volatility, excess 28.3 · its line 2.89 | model 0, 0.000, `the model.py of the gallery` | nothing | 0.900 ± 0.012, not yet to be quoted | 4.00 ± 0.14, not yet to be quoted |
| white, as large as its own noise | where: the wings, excess 49.6 · changes along: the flanks along volatility, excess 37.4 · where: the core, excess 31.4 · its line 2.69 | model 0, 0.921, `the model.py of the gallery` | nothing | 0.884 ± 0.011 | 4.19 ± 0.15 |
| white, 3 times as large as its own noise | changes along: the flanks along volatility, excess 10.1 · where: the wings, excess 4.07 · its line 3.47 | model 0, 0.920, `the model.py of the gallery` | nothing | 0.855 ± 0.028 | 4.62 ± 0.55 |
