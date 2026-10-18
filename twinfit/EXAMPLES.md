# Examples of the package

9 projects of other fields, each written by `python -m twinfit example NAME DIR`, each drawn with a truth that its folder does not hold. The projects pedestal, satellite, saturation, filter, fano and counts were each first worked by a fresh agent that was not told the truth, before they joined the gallery, and clock, gel and transit were written for the gallery and worked by no agent yet. For each the loop ran three times: on the first model alone, after what to do, and after what not to do. Then at four other noises, on the first model and after what to do: three times its own, none drawn, and a second kind of noise as large as its own and three times as large. `python -m twinfit example NAME DIR --noise L --kind K` writes each: L is a multiple of the example's own noise, 0 for none drawn, and K one of `white`, `correlated` and `counts`, the example's own kind where K is not given. Every number below is read from those runs, which `python -m twinfit.examples.gallery --write` makes again for those who build the package.

Find the example whose first report leads as yours does, and read what to do there. Each answer of a report is given with its excess, 0 at the median of its statistic over the synthetic sets of the twin and 1 at their 95th percentile, and each report with its line, the excess that the largest answer of a synthetic set passes on 5 sets of 100: an answer past it stands beyond the twin. The excess says how clearly a term shows against the noise, and not how large the term is.

| example | the field | what the first report puts first beyond its twin |
|---|---|---|
| pedestal | absorption lines of a gas cell at five pressures | what: a wider pedestal under the curve, excess 229 · what: a second component of the curve, excess 131 · what: a polynomial of the level, excess 122 · its line 2.83 |
| satellite | absorption lines of a gas cell with a weak line beside each | what: a second component of the curve, excess 5.45e+04 · what: a tail on one side, excess 2.41e+04 · what: a wider pedestal under the curve, excess 1.84e+04 · its line 2.73 |
| saturation | fluorescence of a resonance driven at five powers of a laser | along: power, excess 3.87e+04 · changes along: the flanks along power, excess 1.17e+04 · changes along: the core along power, excess 8.97e+03 · its line 4.24 |
| filter | absorption dips recorded behind the low-pass filter of a detector | what: a noise correlated from sample to sample, excess 64.7 · changes along: the core along temperature, excess 26.8 · what: a periodic function, excess 21 · its line 3.33 |
| fano | photoionization across an autoionizing resonance at four pressures | what: a tail on one side, excess 1.29e+06 · where: the asymmetry, excess 5.84e+05 · what: a periodic function, excess 4.85e+05 · its line 3.94 |
| counts | photon counts of an emission line at four powers of the pump | what: a noise that grows with the signal, excess 323 · what: glitches, excess 53.2 · changes along: the core along power, excess 26.3 · its line 2.76 |
| clock | the resonance of an atomic clock read at five powers of the probe | along: power, excess 5.23e+04 · changes along: the asymmetry along power, excess 1.84e+04 · changes along: the core along power, excess 1.49e+04 · its line 2.48 |
| gel | a band of electrophoresis in a gel at four voltages | what: a second component of the curve, excess 7.72e+03 · what: a tail on one side, excess 7.55e+03 · where: the flanks, excess 3.01e+03 · its line 2.76 |
| transit | the transit of an exoplanet seen through four filters | where: the core, excess 1.28e+03 · what: a second component of the curve, excess 970 · what: a periodic function, excess 735 · its line 2.92 |

## How to read the tables at other noises

Each row is the project written with `--noise L --kind K`, its first model run, and then what to do, whose term the third column gives. For a term written at the top of a scan, it gives the range of the scan too: the values whose fall is within 1 of its top.

At L times its noise the description states L times the noise, so the bars grow about L times. With none drawn it states the example's own noise: the bars stay those of the own row, the chi square falls near 0, and the check of the noise says that the noise is stated too large, which is why those values are not yet to be quoted.

The rows of white and of correlated noise draw one stream of random numbers, scaled: a row at three times its noise holds the draw of the row at its own noise, three times larger. Their chi squares are then equal, unless the project also draws something of a fixed size, as the glitches of the saturation, and the two rows are not independent experiments.

A counter's noise at three times is drawn as nine times fewer counts, scaled up nine times, so the root of the counts is a third of the noise. What to do for a counter, `python -m twinfit counts DIR`, states that root: there it states the noise three times too small, and the chi square comes near 9.

A parameter of a term written at the top of a scan is named by that term, as "as pedestal 8": its size is read at that stretch, distance or period, and a size at another one is another quantity. So the sizes of two rows are compared at one stretch only. Where the range of the scan is wide, the records do not fix the stretch, and the bar of the size, read at one stretch, is too small.

## pedestal: absorption lines of a gas cell at five pressures

`python -m twinfit example pedestal DIR`. The truth, which the folder does not hold: a pedestal 6 times wider, of relative area 0.04, under a Gaussian of width 4 MHz.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 1 | model 0 |
| its terms | one Gaussian line of one width | one Gaussian line of one width, a pedestal 6 times wider than the curve | one Gaussian line of one width |
| chi square for each degree of freedom | 1.124 | 0.994 | 1.124 |
| the last report, first beyond its twin | what: a wider pedestal under the curve, excess 229 · what: a second component of the curve, excess 131 · what: a polynomial of the level, excess 122 · its line 2.83 | nothing | what: a wider pedestal under the curve, excess 229 · what: a second component of the curve, excess 131 · what: a polynomial of the level, excess 122 · its line 2.83 |
| width, true 4 | 4.0128 ± 0.0011, not yet to be quoted | 3.9986 ± 0.0012 | 4.0128 ± 0.0011, not yet to be quoted |
| pedestal 6, true 0.04 | not in the model | 0.0418 ± 0.0018 | not in the model |

The first report proposes for the model: a pedestal 8 times wider than the curve, `term = mechanisms.pedestal(8)` · a curvature of the level, its own in each record, `term = mechanisms.level(2)` · a second component 1.41 widths away on both sides. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.73 · a wider span, a better data collection: the error of the known parameters times 0.76.

**What to do**: scan the stretch of the shape that the report names first, a wider pedestal under the line, and write the term at the top of the scan. Here `python -m twinfit scan DIR "pedestal(x)" 3 4 6 8 12` puts its top at 6: the term file holds `term = mechanisms.pedestal(6)`.

**What not to do**: make the width change along the pressure, a setting that the report does not name. The loop refused it on 30 synthetic sets: width changes along pressure.

**At other noises**: the first model, then what to do. Its own noise is white, and its second kind correlated over neighbouring samples. The rows are written by `--noise 3`, `--noise 0`, `--kind correlated` and `--noise 3 --kind correlated`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | width, true 4 | pedestal 6, true 0.04 |
|---|---|---|---|---|---|
| white, 3 times as large as its own noise | what: a wider pedestal under the curve, excess 30.1 · what: a second component of the curve, excess 16.5 · what: a polynomial of the level, excess 15.2 · its line 3.08 | model 1, 0.994, `term = mechanisms.pedestal(8)`, its scan within 1 of its top from -10 to 26 | nothing | 3.9971 ± 0.0036 | 0.0827 ± 0.0100, as pedestal 8 |
| no noise drawn | what: a wider pedestal under the curve, excess 208 · what: a second component of the curve, excess 117 · what: a polynomial of the level, excess 110 · its line 2.83 | model 1, 0.000, `term = mechanisms.pedestal(6.73)`, its scan within 1 of its top from 5.36 to 6.73 | nothing | 4.0005 ± 0.0012, not yet to be quoted | 0.0499 ± 0.0023, as pedestal 6.73, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | what: a wider pedestal under the curve, excess 260 · what: a second component of the curve, excess 158 · what: a polynomial of the level, excess 141 · its line 3.13 | model 1, 0.948, `term = mechanisms.pedestal(10.1)`, its scan within 1 of its top from 8.26 to 12 | what: a noise correlated from sample to sample, excess 31.2 · changes along: the core along pressure, excess 17.2 · what: a polynomial of the level, excess 8.74 · its line 2.72 | 3.9985 ± 0.0012, not yet to be quoted | 0.1420 ± 0.0058, as pedestal 10.1, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | what: a wider pedestal under the curve, excess 42.1 · what: a noise correlated from sample to sample, excess 32.1 · what: a second component of the curve, excess 28.9 · its line 2.74 | model 1, 0.948, `term = mechanisms.pedestal(56)`, its scan within 1 of its top from 36 to 60 | what: a noise correlated from sample to sample, excess 31.2 · changes along: the core along pressure, excess 16 · changes along: the flanks along pressure, excess 9.22 · its line 2.58 | 3.9931 ± 0.0035, not yet to be quoted | 23.2 ± 2.3, as pedestal 56, not yet to be quoted |

## satellite: absorption lines of a gas cell with a weak line beside each

`python -m twinfit example satellite DIR`. The truth, which the folder does not hold: a second component 12 MHz above, of relative area 0.08, beside a Gaussian of width 4 MHz.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 1 | model 0 |
| its terms | one Gaussian line of one width | one Gaussian line of one width, a second component 12.015 away | one Gaussian line of one width |
| chi square for each degree of freedom | 38.051 | 0.995 | 38.051 |
| the last report, first beyond its twin | what: a second component of the curve, excess 5.45e+04 · what: a tail on one side, excess 2.41e+04 · what: a wider pedestal under the curve, excess 1.84e+04 · its line 2.73 | nothing | what: a second component of the curve, excess 5.45e+04 · what: a tail on one side, excess 2.41e+04 · what: a wider pedestal under the curve, excess 1.84e+04 · its line 2.73 |
| width, true 4 | 4.1163 ± 0.0011, not yet to be quoted | 3.9989 ± 0.0011 | 4.1163 ± 0.0011, not yet to be quoted |
| second component at 12, true 0.08 | not in the model | 0.08004 ± 0.00021, as second component at 12.015 | not in the model |

The first report proposes for the model: a second component 1.41 widths above the centre, `term = mechanisms.second(13.5)` · a tail on the upper side, 4 times wider, `term = mechanisms.tail(4, "upper")` · a pedestal 2 times wider than the curve, `term = mechanisms.pedestal(2)`. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.92 · a wider span, a better data collection: the error of the known parameters times 0.96.

**What to do**: scan the distance of a second component, `python -m twinfit scan DIR "second(x)" ...`, and write the term at the top of the scan. Here `python -m twinfit scan DIR "second(x)" 6 8 10 12 14 16 18` puts its top at 11.95, and a finer scan about it, at 11.8522 11.8845 11.9167 11.949 11.9813 12.0135 12.0458, at 12.02: the term file holds `term = mechanisms.second(12.015)`.

**What not to do**: widen the line with a kernel, as if the satellite were wings of the line. The loop refused it on 30 synthetic sets: a kernel with exponential wings, by a hair, which the loop itself asks to run again at 500 sets before the term is set aside: at more sets it may be kept, and that a wrong term is kept says that the noise check and the chi square, not the twin alone, judge it.

**At other noises**: the first model, then what to do. Its own noise is white, and its second kind correlated over neighbouring samples. The rows are written by `--noise 3`, `--noise 0`, `--kind correlated` and `--noise 3 --kind correlated`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | width, true 4 | second component at 12, true 0.08 |
|---|---|---|---|---|---|
| white, 3 times as large as its own noise | what: a second component of the curve, excess 6.08e+03 · what: a tail on one side, excess 3.08e+03 · what: a wider pedestal under the curve, excess 2.02e+03 · its line 2.93 | model 1, 0.995, `term = mechanisms.second(12.046)`, its scan within 1 of its top from 12 to 12.1 | nothing | 3.9967 ± 0.0032 | 0.08011 ± 0.00062, as second component at 12.046 |
| no noise drawn | what: a second component of the curve, excess 5.44e+04 · what: a tail on one side, excess 2.4e+04 · what: a wider pedestal under the curve, excess 1.85e+04 · its line 2.75 | model 1, 0.000, `term = mechanisms.second(12)`, its scan within 1 of its top from 12 to 12 | nothing | 4.0000 ± 0.0011, not yet to be quoted | 0.08000 ± 0.00021, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | what: a second component of the curve, excess 5.46e+04 · what: a tail on one side, excess 2.81e+04 · what: a wider pedestal under the curve, excess 1.81e+04 · its line 2.91 | model 1, 0.939, `term = mechanisms.second(12.045)`, its scan within 1 of its top from 12 to 12.1 | what: a noise correlated from sample to sample, excess 30.8 · changes along: the asymmetry along pressure, excess 13 · changes along: the core along pressure, excess 9.35 · its line 2.9 | 3.9965 ± 0.0011, not yet to be quoted | 0.08002 ± 0.00021, as second component at 12.045, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | what: a second component of the curve, excess 6.29e+03 · what: a tail on one side, excess 3.19e+03 · what: a wider pedestal under the curve, excess 1.92e+03 · its line 2.63 | model 1, 0.939, `term = mechanisms.second(12.136)`, its scan within 1 of its top from 12.1 to 12.2 | what: a noise correlated from sample to sample, excess 30.7 · changes along: the asymmetry along pressure, excess 13.8 · changes along: the core along pressure, excess 9.93 · its line 3.01 | 3.9895 ± 0.0032, not yet to be quoted | 0.08007 ± 0.00061, as second component at 12.136, not yet to be quoted |

## saturation: fluorescence of a resonance driven at five powers of a laser

`python -m twinfit example saturation DIR`. The truth, which the folder does not hold: full width 1.0 MHz times the root of 1 + P / 50 mW, the broadening by saturation, and height 1.5 P / (P + 50 mW), place 0.3 MHz in every record, background 0.05, three glitches of +0.08 in the record at 40 mW, at -9.3, 2.1 and 7.7 MHz, noise 0.004.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 0 |
| its terms | one Lorentzian line of one width | one Lorentzian line whose width saturates with the power | one Lorentzian line of one width |
| chi square for each degree of freedom | 12.573 | 1.150 | 12.573 |
| the last report, first beyond its twin | along: power, excess 3.87e+04 · changes along: the flanks along power, excess 1.17e+04 · changes along: the core along power, excess 8.97e+03 · its line 4.24 | what: glitches, excess 96 · its line 2.85 | along: power, excess 3.87e+04 · changes along: the flanks along power, excess 1.17e+04 · changes along: the core along power, excess 8.97e+03 · its line 4.24 |
| width at no power, true 1 | not in the model | 0.9969 ± 0.0028, not yet to be quoted | not in the model |
| saturation power, true 50 | not in the model | 49.36 ± 0.42, not yet to be quoted | not in the model |

The first report proposes for the model: width changes along power, `term = mechanisms.along("width", "power")` · a curvature of the level, its own in each record, `term = mechanisms.level(2)` · a periodic term of period 3.528, `term = mechanisms.periodic(3.528, "cosine")`. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more repeats, a better data collection: the error of the known parameters times 0.68 · more samples, a better data collection: the error of the known parameters times 0.76.

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
| white, 3 times as large as its own noise | along: power, excess 4.41e+03 · changes along: the flanks along power, excess 1.24e+03 · changes along: the core along power, excess 1.06e+03 · its line 4.22 | model 0, 1.008, `the model.py of the gallery` | what: glitches, excess 8.34 · its line 2.87 | 0.9894 ± 0.0086, not yet to be quoted | 47.9 ± 1.2, not yet to be quoted |
| no noise drawn | along: power, excess 3.8e+04 · changes along: the flanks along power, excess 1.15e+04 · changes along: the core along power, excess 9.16e+03 · its line 4.06 | model 0, 0.160, `the model.py of the gallery` | what: glitches, excess 93.4 · its line 2.77 | 1.0007 ± 0.0028, not yet to be quoted | 50.09 ± 0.43, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | along: power, excess 4e+04 · changes along: the core along power, excess 9.62e+03 · changes along: the flanks along power, excess 9.45e+03 · its line 4.25 | model 0, 1.156, `the model.py of the gallery` | what: glitches, excess 97 · what: a noise correlated from sample to sample, excess 54.7 · what: a periodic function, excess 22.8 · its line 2.87 | 0.9892 ± 0.0029, not yet to be quoted | 47.92 ± 0.41, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | along: power, excess 4.86e+03 · changes along: the flanks along power, excess 1.12e+03 · changes along: the core along power, excess 1.08e+03 · its line 3.49 | model 0, 1.006, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 55.5 · what: a periodic function, excess 23.7 · what: a second component of the curve, excess 19.6 · its line 3.08 | 0.9663 ± 0.0087, not yet to be quoted | 43.9 ± 1.1, not yet to be quoted |

## filter: absorption dips recorded behind the low-pass filter of a detector

`python -m twinfit example filter DIR`. The truth, which the folder does not hold: the model is whole, width 1.2 pm, and the noise is correlated from sample to sample, 0.8 from one to the next, of spread 0.002 as stated.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 0 |
| its terms | one Gaussian dip of one width | one Gaussian dip of one width | one Gaussian dip of one width |
| chi square for each degree of freedom | 1.028 | 1.051 | 1.028 |
| the last report, first beyond its twin | what: a noise correlated from sample to sample, excess 64.7 · changes along: the core along temperature, excess 26.8 · what: a periodic function, excess 21 · its line 3.33 | nothing | what: a noise correlated from sample to sample, excess 64.7 · changes along: the core along temperature, excess 26.8 · what: a periodic function, excess 21 · its line 3.33 |
| width, true 1.2 | 1.20024 ± 0.00033, not yet to be quoted | 1.2007 ± 0.0013 | 1.20024 ± 0.00033, not yet to be quoted |

The first report proposes for the model: a periodic term of period 0.7728, `term = mechanisms.periodic(0.7728, "cosine")` · a second component 1.41 widths below the centre, `term = mechanisms.second(-3.98)` · a curvature of the level, its own in each record, `term = mechanisms.level(2)`. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.70 · more repeats, a better data collection: the error of the known parameters times 0.71.

**What to do**: read the correlated noise as noise: keep one sample in every few with `python -m twinfit thin DIR N`, at the span that the report gives, and run the loop on the copy. Here `python -m twinfit thin DIR 15` writes 3 copies, and the loop runs on the first, do_thin_15_0.

**What not to do**: hold the correlation with a term of the curve, a pedestal. The loop refused it on 30 synthetic sets: a pedestal 6 times wider than the curve.

**At other noises**: the first model, then what to do. Its own noise is correlated over neighbouring samples, and its second kind white. The rows are written by `--noise 3`, `--noise 0`, `--kind white` and `--noise 3 --kind white`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | width, true 1.2 |
|---|---|---|---|---|
| correlated over neighbouring samples, 3 times as large as its own noise | what: a noise correlated from sample to sample, excess 66.8 · changes along: the core along temperature, excess 25.8 · where: the flanks, excess 13.7 · its line 2.89 | model 0, 1.051, `python -m twinfit thin DIR 15` | nothing | 1.2021 ± 0.0038 |
| no noise drawn | nothing | model 0, 0.000, `python -m twinfit thin DIR 15` | nothing | 1.2000 ± 0.0013, not yet to be quoted |
| white, as large as its own noise | what: a noise that grows with the signal, excess 3.12 · changes along: the core along temperature, excess 3.05 · its line 3.04 | model 0, 1.051, `python -m twinfit thin DIR 15` | nothing | 1.2000 ± 0.0013 |
| white, 3 times as large as its own noise | nothing | model 0, 1.051, `python -m twinfit thin DIR 15` | nothing | 1.2001 ± 0.0038 |

## fano: photoionization across an autoionizing resonance at four pressures

`python -m twinfit example fano DIR`. The truth, which the folder does not hold: a Fano profile, q = 3 and half width 2 meV at every pressure, its top near one, and the first model is a Lorentzian, which lacks the interference.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 0 |
| its terms | one Lorentzian line of one half width | one Fano profile | one Lorentzian line of one half width |
| chi square for each degree of freedom | 475.375 | 1.001 | 475.375 |
| the last report, first beyond its twin | what: a tail on one side, excess 1.29e+06 · where: the asymmetry, excess 5.84e+05 · what: a periodic function, excess 4.85e+05 · its line 3.94 | nothing | what: a tail on one side, excess 1.29e+06 · where: the asymmetry, excess 5.84e+05 · what: a periodic function, excess 4.85e+05 · its line 3.94 |
| half_width, true 2 | 2.2254 ± 0.0011, not yet to be quoted | 1.99966 ± 0.00088 | 2.2254 ± 0.0011, not yet to be quoted |
| q, true 3 | not in the model | 2.9971 ± 0.0017 | not in the model |

The first report proposes for the model: a second component 1 widths below the centre, `term = mechanisms.second(-4.22)` · a periodic term of period 6.493, `term = mechanisms.periodic(6.493, "cosine")` · a tail on the lower side, 4 times wider, `term = mechanisms.tail(4, "lower")`. And for the experiment: more repeats, a better data collection: the error of the known parameters times 0.48 · a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.58.

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
| white, 3 times as large as its own noise | what: a tail on one side, excess 1.4e+05 · what: a periodic function, excess 7.32e+04 · where: the asymmetry, excess 6.53e+04 · its line 3.93 | model 0, 1.001, `the model.py of the gallery` | nothing | 1.9990 ± 0.0026 | 2.9912 ± 0.0051 |
| no noise drawn | what: a tail on one side, excess 1.29e+06 · what: a periodic function, excess 5.9e+05 · where: the asymmetry, excess 5.68e+05 · its line 3.9 | model 0, 0.000, `the model.py of the gallery` | nothing | 2.00000 ± 0.00088, not yet to be quoted | 3.0000 ± 0.0017, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | what: a tail on one side, excess 1.29e+06 · where: the asymmetry, excess 5.8e+05 · what: a periodic function, excess 4.69e+05 · its line 3.94 | model 0, 1.000, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 65 · what: a polynomial of the level, excess 16.6 · what: a periodic function, excess 13.1 · its line 2.73 | 1.99870 ± 0.00088, not yet to be quoted | 2.9911 ± 0.0017, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | what: a tail on one side, excess 1.36e+05 · where: the asymmetry, excess 6.27e+04 · what: a periodic function, excess 5.84e+04 · its line 3.77 | model 0, 1.000, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 63.2 · what: a polynomial of the level, excess 16.5 · what: a periodic function, excess 14.7 · its line 2.92 | 1.9961 ± 0.0026, not yet to be quoted | 2.9735 ± 0.0050, not yet to be quoted |

## counts: photon counts of an emission line at four powers of the pump

`python -m twinfit example counts DIR`. The truth, which the folder does not hold: the model is whole, half width 1.5 GHz, and the noise is that of a counter, the root of the counts, 6.3 on the background and 40 at the top of the strongest line.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 0 |
| its terms | one Lorentzian line of one half width | one Lorentzian line of one half width | one Lorentzian line of one half width |
| chi square for each degree of freedom | 4.509 | 0.984 | 4.509 |
| the last report, first beyond its twin | what: a noise that grows with the signal, excess 323 · what: glitches, excess 53.2 · changes along: the core along power, excess 26.3 · its line 2.76 | nothing | what: a noise that grows with the signal, excess 323 · what: glitches, excess 53.2 · changes along: the core along power, excess 26.3 · its line 2.76 |
| half_width, true 1.5 | 1.50179 ± 0.00094, not yet to be quoted | 1.5009 ± 0.0032 | 1.50179 ± 0.00094, not yet to be quoted |

The first report proposes for the model: a periodic term of period 0.8074, `term = mechanisms.periodic(0.8074, "cosine")` · a second component 0.5 widths away on both sides · a tail on the lower side, 4 times wider, `term = mechanisms.tail(4, "lower")`. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.71 · more repeats, a better data collection: the error of the known parameters times 0.75.

**What to do**: give each sample the noise of a counter, the root of the fitted counts, with `python -m twinfit counts DIR`, and run the loop on the copy. Here `python -m twinfit counts DIR` writes the copy, and the loop runs on it.

**What not to do**: hold the misfit at the top of the lines with a term of the curve, a pedestal. The loop refused it on 30 synthetic sets: a pedestal 6 times wider than the curve.

**At other noises**: the first model, then what to do. Its own noise is that of a counter, and its second kind white. The rows are written by `--noise 3`, `--noise 0`, `--kind white` and `--noise 3 --kind white`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | half_width, true 1.5 |
|---|---|---|---|---|
| that of a counter, 3 times as large as its own noise | what: a noise that grows with the signal, excess 352 · what: glitches, excess 91.8 · what: a periodic function, excess 30.2 · its line 3.13 | model 0, 8.853, `python -m twinfit counts DIR` | what: glitches, excess 33.8 · what: a periodic function, excess 20.3 · what: a polynomial of the level, excess 15.8 · its line 2.52 | 1.5030 ± 0.0032, not yet to be quoted |
| no noise drawn | nothing | model 0, 0.000, `python -m twinfit counts DIR` | nothing | 1.5000 ± 0.0032, not yet to be quoted |
| white, as large as its own noise | nothing | model 0, 0.458, `python -m twinfit counts DIR` | what: a noise that grows with the signal, excess 8.44 · its line 2.48 | 1.5008 ± 0.0032, not yet to be quoted |
| white, 3 times as large as its own noise | nothing | model 0, 4.127, `python -m twinfit counts DIR` | what: a noise that grows with the signal, excess 85 · what: glitches, excess 19.3 · what: a periodic function, excess 17.8 · its line 2.5 | 1.5023 ± 0.0032, not yet to be quoted |

## clock: the resonance of an atomic clock read at five powers of the probe

`python -m twinfit example clock DIR`. The truth, which the folder does not hold: the centre moves with the power of the probe, 0.12 Hz from the nominal frequency at no power and 0.025 Hz for each microwatt, the light shift, and the width is 1.5 Hz at every power.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 0 |
| its terms | one Lorentzian line of one width and one centre | one Lorentzian line whose centre moves with the power | one Lorentzian line of one width, its centre its own in each curve |
| chi square for each degree of freedom | 31.812 | 0.972 | 0.972 |
| the last report, first beyond its twin | along: power, excess 5.23e+04 · changes along: the asymmetry along power, excess 1.84e+04 · changes along: the core along power, excess 1.49e+04 · its line 2.48 | nothing | nothing |
| width, true 1.5 | 1.51614 ± 0.00058, not yet to be quoted | 1.49970 ± 0.00058 | 1.49970 ± 0.00058 |
| clock frequency, true 0.12 | not in the model | 0.12000 ± 0.00038 | not in the model |
| light shift, true 0.025 | not in the model | 0.025051 ± 0.000058 ± 0.00050 of theory | not in the model |

The first report proposes for the model: centre changes along power, `term = mechanisms.along("centre", "power")` · a second component 0.35 widths away on both sides · a slope of the level, its own in each record, `term = mechanisms.level(1)`. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.66 · more repeats, a better data collection: the error of the known parameters times 0.67.

**What to do**: write the light shift that the physics gives, the centre moving in proportion to the power read through its calibration, which the description holds within 2 per cent, so that the model measures the frequency at no power and the light shift carries the bar of that calibration.

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


def model(x, shared, own, design):
    centre = shared[1] + shared[2] * design["power"] * HELD["power scale"]
    return own[0] / (1.0 + (2.0 * (x - centre) / shared[0]) ** 2) + own[1]


def start_own(record):
    base = float(np.median(record.y))
    return [float(record.y.max() - base), base]
```

**What not to do**: free the centre of each curve: the fit is right, and the frequency of the clock is then in no parameter, the shift hidden in five centres.

**At other noises**: the first model, then what to do. Its own noise is white, and its second kind correlated over neighbouring samples. The rows are written by `--noise 3`, `--noise 0`, `--kind correlated` and `--noise 3 --kind correlated`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | width, true 1.5 | clock frequency, true 0.12 | light shift, true 0.025 |
|---|---|---|---|---|---|---|
| white, 3 times as large as its own noise | along: power, excess 5.86e+03 · changes along: the asymmetry along power, excess 2.04e+03 · changes along: the core along power, excess 1.7e+03 · its line 2.48 | model 0, 0.972, `the model.py of the gallery` | nothing | 1.4991 ± 0.0017 | 0.1200 ± 0.0012 | 0.02515 ± 0.00017 ± 0.00050 of theory |
| no noise drawn | along: power, excess 5.2e+04 · changes along: the asymmetry along power, excess 1.84e+04 · changes along: the core along power, excess 1.46e+04 · its line 2.48 | model 0, 0.000, `the model.py of the gallery` | nothing | 1.50000 ± 0.00058, not yet to be quoted | 0.12000 ± 0.00038, not yet to be quoted | 0.025000 ± 0.000058 ± 0.00050 of theory, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | along: power, excess 5.27e+04 · changes along: the asymmetry along power, excess 1.83e+04 · changes along: the core along power, excess 1.53e+04 · its line 2.48 | model 0, 0.930, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 51.8 · what: a tail on one side, excess 21.4 · changes along: the core along power, excess 15.2 · its line 2.85 | 1.49900 ± 0.00058, not yet to be quoted | 0.12006 ± 0.00038, not yet to be quoted | 0.025138 ± 0.000058 ± 0.00050 of theory, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | along: power, excess 6.02e+03 · changes along: the asymmetry along power, excess 2.02e+03 · changes along: the core along power, excess 1.87e+03 · its line 2.49 | model 0, 0.930, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 51.8 · what: a tail on one side, excess 21.3 · changes along: the core along power, excess 15.2 · its line 2.85 | 1.4970 ± 0.0017, not yet to be quoted | 0.1202 ± 0.0012, not yet to be quoted | 0.02542 ± 0.00017 ± 0.00051 of theory, not yet to be quoted |

## gel: a band of electrophoresis in a gel at four voltages

`python -m twinfit example gel DIR`. The truth, which the folder does not hold: a tail on the side of slower migration, a copy of the band three times wider below its place, of size 0.15 as the tail of the package counts it, beside a Gaussian band of standard deviation 1 mm.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 1 | model 1 |
| its terms | one Gaussian band of one width | one Gaussian band of one width, a tail on the lower side, 3.03 times wider | one Gaussian band of one width, a pedestal 3 times wider than the curve |
| chi square for each degree of freedom | 7.304 | 0.999 | 5.371 |
| the last report, first beyond its twin | what: a second component of the curve, excess 7.72e+03 · what: a tail on one side, excess 7.55e+03 · where: the flanks, excess 3.01e+03 · its line 2.76 | nothing | what: a tail on one side, excess 6.56e+03 · what: a second component of the curve, excess 4.33e+03 · where: the asymmetry, excess 2.53e+03 · its line 2.97 |
| width, true 1 | 1.02129 ± 0.00024, not yet to be quoted | 1.00023 ± 0.00027 | 0.99926 ± 0.00033, not yet to be quoted |
| tail lower 3, true 0.15 | not in the model | 0.14918 ± 0.00086, as tail lower 3.03 | not in the model |

The first report proposes for the model: a second component 1.41 widths below the centre, `term = mechanisms.second(-3.37)` · a tail on the lower side, 4 times wider, `term = mechanisms.tail(4, "lower")` · a pedestal 4 times wider than the curve, `term = mechanisms.pedestal(4)`. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.67 · more repeats, a better data collection: the error of the known parameters times 0.86.

**What to do**: scan the stretch of the tail on the side of slower migration that the report points to, and write the term at the top of the scan. Here `python -m twinfit scan DIR "tail(x, 'lower')" 2 3 4 5 6` puts its top at 3.29, and a finer scan about it, at 3.20106 3.23071 3.26035 3.29 3.31965 3.34929 3.37894, at 3.201, and a finer scan about it, at 3.11212 3.14177 3.17142 3.20106 3.23071 3.26035 3.29, at 3.112, and a finer scan about it, at 3.02318 3.05283 3.08248 3.11212 3.14177 3.17141 3.20106, at 3.023, and a finer scan about it, at 2.93424 2.96389 2.99354 3.02318 3.05283 3.08248 3.11212, at 3.03: the term file holds `term = mechanisms.tail(3.03, 'lower')`.

**What not to do**: a pedestal under the band, wider on both sides alike.

**At other noises**: the first model, then what to do. Its own noise is white, and its second kind correlated over neighbouring samples. The rows are written by `--noise 3`, `--noise 0`, `--kind correlated` and `--noise 3 --kind correlated`.

| noise | the first report, first beyond its twin | after what to do: the best model, its chi square, the term written | its last report, first beyond its twin | width, true 1 | tail lower 3, true 0.15 |
|---|---|---|---|---|---|
| white, 3 times as large as its own noise | what: a second component of the curve, excess 899 · what: a tail on one side, excess 816 · where: the flanks, excess 425 · its line 2.72 | model 1, 0.999, `term = mechanisms.tail(3.06464, 'lower')`, its scan within 1 of its top from 3.06 to 3.15 | nothing | 1.00054 ± 0.00080 | 0.1478 ± 0.0026, as tail lower 3.06464 |
| no noise drawn | what: a tail on one side, excess 7.57e+03 · what: a second component of the curve, excess 7.57e+03 · where: the flanks, excess 3.2e+03 · its line 2.67 | model 1, 0.000, `term = mechanisms.tail(3, 'lower')`, its scan within 1 of its top from 2.98 to 3.02 | nothing | 1.00000 ± 0.00027, not yet to be quoted | 0.15000 ± 0.00086, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | what: a second component of the curve, excess 8.19e+03 · what: a tail on one side, excess 7.59e+03 · where: the flanks, excess 3.38e+03 · its line 2.79 | model 1, 0.986, `term = mechanisms.tail(3.101, 'lower')`, its scan within 1 of its top from 3.08 to 3.12 | what: a noise correlated from sample to sample, excess 63.9 · changes along: the core along voltage, excess 13.6 · along: voltage, excess 13.6 · its line 2.82 | 1.00091 ± 0.00027, not yet to be quoted | 0.14717 ± 0.00086, as tail lower 3.101, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | what: a second component of the curve, excess 928 · what: a tail on one side, excess 848 · where: the flanks, excess 411 · its line 2.7 | model 1, 0.986, `term = mechanisms.tail(3.37339, 'lower')`, its scan within 1 of its top from 3.28 to 3.37 | what: a noise correlated from sample to sample, excess 62.5 · along: voltage, excess 14 · changes along: the core along voltage, excess 13.4 · its line 2.8 | 1.00298 ± 0.00079, not yet to be quoted | 0.1423 ± 0.0025, as tail lower 3.37339, not yet to be quoted |

## transit: the transit of an exoplanet seen through four filters

`python -m twinfit example transit DIR`. The truth, which the folder does not hold: a star darker at its limb, by a coefficient u = 0.6 - 0.3 (wavelength - 450 nm) / 350 nm, 0.6 in the blue and 0.3 at 800 nm, under a planet of depth 0.010 and half duration 1.2 h.

| | the first model | after what to do | after what not to do |
|---|---|---|---|
| the best model | model 0 | model 0 | model 1 |
| its terms | the dip of a uniform star | the dip of a star darker at its limb | the dip of a uniform star, a power 2 of the abscissa under every curve |
| chi square for each degree of freedom | 2.397 | 0.984 | 2.342 |
| the last report, first beyond its twin | where: the core, excess 1.28e+03 · what: a second component of the curve, excess 970 · what: a periodic function, excess 735 · its line 2.92 | nothing | what: a second component of the curve, excess 1.6e+03 · what: a periodic function, excess 1.11e+03 · where: the core, excess 1.03e+03 · its line 3.03 |
| half_duration, true 1.2 | 1.14558 ± 0.00052, not yet to be quoted | 1.1997 ± 0.0011 | 1.14304 ± 0.00055, not yet to be quoted |
| limb darkening at 450 nm, true 0.6 | not in the model | 0.5967 ± 0.0072 | not in the model |
| its fall to 800 nm, true 0.3 | not in the model | 0.3091 ± 0.0098 | not in the model |

The first report proposes for the model: a second component 0.71 widths away on both sides · a periodic term of period 1.844, `term = mechanisms.periodic(1.844, "cosine")` · half_duration changes along wavelength, `term = mechanisms.along("half_duration", "wavelength")`. And for the experiment: a quieter setup, a better setup: the error of the known parameters times 0.50 · more samples, a better data collection: the error of the known parameters times 0.74 · more repeats, a better data collection: the error of the known parameters times 0.76.

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
| white, 3 times as large as its own noise | what: a periodic function, excess 139 · where: the core, excess 131 · what: a second component of the curve, excess 105 · its line 2.91 | model 0, 0.984, `the model.py of the gallery` | nothing | 1.2008 ± 0.0032 | 0.591 ± 0.022 | 0.324 ± 0.029 |
| no noise drawn | where: the core, excess 1.23e+03 · what: a second component of the curve, excess 985 · what: a periodic function, excess 589 · its line 2.91 | model 0, 0.001, `the model.py of the gallery` | nothing | 1.1992 ± 0.0011, not yet to be quoted | 0.5993 ± 0.0073, not yet to be quoted | 0.3018 ± 0.0098, not yet to be quoted |
| correlated over neighbouring samples, as large as its own noise | what: a periodic function, excess 1.22e+03 · where: the core, excess 1.06e+03 · what: a second component of the curve, excess 954 · its line 2.94 | model 0, 0.997, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 71.9 · what: a periodic function, excess 17.7 · what: a second component of the curve, excess 6.1 · its line 2.73 | 1.2000 ± 0.0011, not yet to be quoted | 0.5869 ± 0.0072, not yet to be quoted | 0.3199 ± 0.0098, not yet to be quoted |
| correlated over neighbouring samples, 3 times as large as its own noise | where: the core, excess 138 · what: a periodic function, excess 95.4 · what: a second component of the curve, excess 86 · its line 2.76 | model 0, 0.996, `the model.py of the gallery` | what: a noise correlated from sample to sample, excess 72 · what: a periodic function, excess 10.3 · what: a second component of the curve, excess 6.17 · its line 2.78 | 1.2014 ± 0.0030, not yet to be quoted | 0.561 ± 0.021, not yet to be quoted | 0.358 ± 0.030, not yet to be quoted |
