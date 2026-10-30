# twinfit for an agent: the tutorial in commands

Written by `python -m twinfit.tutorial --write` from TUTORIAL.md, the tutorial for people, which holds the reasons, the outputs and the pictures: every command and every file here is taken from it.

It builds on the card, `agent/AGENT.md`, and sets out every command and every file of the two worked projects. Not covered here: the reasons, the outputs and the pictures, which `TUTORIAL.md` gives.

> Section 1 of `TUTORIAL.md` defines every word of the package that this page uses.

The rules that decide each turn:

- The output of the loop for an agent is `DIR/loop/brief_K.md`, one page a model K, with `[output] for` holding "agents",
  as it does by default. `model_K.md` is the page for people and `model_K.json` holds every number.
- A term is asked only where an answer of the report stands beyond the twin: a right model passes the line of an answer
  on 5 sets of 100, so one answer just past it is weak evidence.
- The loop keeps a fitted term when its fit ends inside its bounds, the Bayesian criterion falls by more than 2 and the
  twin admits the larger model, and best.md gives the reason of a refusal.
- A term that the physics asks and the records cannot show is held at its theory, with hold = True in its file.
- Quote a value only from a model that its twin admits, with the parts of its bar that best.md gives.
- `python -m twinfit wait DIR` waits for a running loop within ten minutes and prints the last brief: never poll a run
  in another way, and never start two in one folder.

## 1. Purpose and vocabulary of the package

You give it curves recorded across a design, such as spectra at several temperatures, and the model that should describe them.

## 2. Requirements and installation

Python 3.12 or later, with numpy 2.5 or later and scipy 1.16 or later, and matplotlib for the visual log of a run: numpy 2.5 convolves about ten times faster than numpy 1.26 on a Mac with an arm64 processor, and about as fast on Linux with an x86_64 processor, and the oldest versions that the package has run on are Python 3.12, numpy 2.5.0 and scipy 1.16.0.

```bash
python -c "import twinfit; print(twinfit.__version__)"
```

## 3. A worked project

This writes a project of nuclear magnetic resonance into a new folder, `my_nmr`:

```bash
python -m twinfit example nmr my_nmr
```

`my_nmr/description.toml`:

```toml
[problem]
name = "nmr"
what = "one resonance of a sample measured at four temperatures. Each spectrum has its own height, frequency and baseline. The first model is one Lorentzian line of one width for all four"
abscissa = "frequency"
abscissa_unit = "Hz"
value = "intensity"
value_unit = "a.u."

[covariates.temperature]
unit = "K"

[[trace]]
file = "traces/temperature_290.csv"
columns = ["frequency", "intensity"]
temperature = 290

[[trace]]
file = "traces/temperature_300.csv"
columns = ["frequency", "intensity"]
temperature = 300

[[trace]]
file = "traces/temperature_310.csv"
columns = ["frequency", "intensity"]
temperature = 310

[[trace]]
file = "traces/temperature_320.csv"
columns = ["frequency", "intensity"]
temperature = 320

[noise]
from = "value"
value = 0.01

[addenda]
model = "model.py"

[addenda.lever.quieter]
kind = "noise"
factor = 0.5
words = "a quieter detector, which halves the noise of each sample"

[addenda.lever.repeats]
kind = "repeats"
factor = 3
words = "each run made three times"

[addenda.lever.slower]
kind = "samples"
factor = 2
words = "twice the samples over the same span"
```

`my_nmr/model.py`:

```python
"""The first model of the spectra: one Lorentzian line of one width, on its own height, place and baseline."""
import numpy as np

SHARED = ("width",)                     # the full width of the line at half its height, in Hz, the same for every spectrum
OWN = ("height", "place", "baseline")   # of each spectrum
START_SHARED = [2.0]
BOUNDS_SHARED = ([0.05], [20.0])
CENTRE = "place"
NAME = "one Lorentzian line of one width"
EQUATION = r"$\frac{h}{1 + (2\,(x - x_0)/w)^2} + b$"          # the model in mathtext, which the visual log writes at each turn
LETTERS = {"h": "the height", "x_0": "the place", "w": "the full width at half height", "b": "the baseline"}          # each letter of EQUATION, which the visual log writes under it


def model(x, shared, own, design):
    return own[0] / (1.0 + (2.0 * (x - own[1]) / shared[0]) ** 2) + own[2]


def start_own(record):
    return [float(record.y.max() - record.y.min()), float(record.x[np.argmax(record.y)]), float(np.median(record.y))]


def centre(shared, own, design):
    return own[1]
```

## 4. The first run

The loop runs here in the foreground, since each run of this section takes a minute or two.

```bash
python -m twinfit loop my_nmr
```

```bash
python -m twinfit matrix my_nmr
```

## 5. A term

A term is a file in `terms/`, `01.py`, then `02.py`, and so on, which you write. This prints the terms that a term file can name, each with the one parameter that it adds, and the form of a term file:

```bash
python -m twinfit terms
```

`my_nmr/terms/01.py`:

```python
from twinfit import mechanisms
term = mechanisms.along("width", "temperature")
why = "the report: the core changes along temperature, and the width of a resonance changes with the temperature of the sample"
```

## 6. The second run

```bash
python -m twinfit loop my_nmr
```

Two words of these lines are defined here, since section 9 uses them later: your list is `list_of_terms.py` of the project, the terms that you suspect, empty at the start, and matched is the question of the report that reads the residuals against each term of that list.

```bash
python -m twinfit loop my_nmr --sets 100
```

## 7. A term not asked for by the records

The loop is not a machine that keeps whatever it is given. Here is a term that no answer of the report points to, a ripple along the abscissa:

`my_nmr/terms/02.py`:

```python
from twinfit import mechanisms
term = mechanisms.periodic(3.0, "cosine")
why = "a test: no answer of the report points to a periodic function"
```

```bash
python -m twinfit loop my_nmr
```

## 8. The proposals for the experiment

Each report ends with three proposals for the experiment, each of one of three kinds: a better setup, a better characterization of the setup, or a better data collection.

## 9. The researcher's own curves

This writes the starting files of a project into a new folder:

```bash
python -m twinfit init my_project
```

`my_project/description.toml`:

```toml
# The description of the data collection, which `python -m twinfit loop` reads. Fill the places marked `TO FILL`.

[problem]
name = "my problem"                      # `TO FILL`: a short name
what = "TO FILL: the physical problem in words: what was measured, how, and what the first model holds"
abscissa = "x"                           # `TO FILL`: the name of the abscissa of the curves, and its unit below
abscissa_unit = "s"                      # `TO FILL`: the unit of the abscissa
# abscissa_offset = 587.6                # optional: taken from the abscissa of every curve, where the curves sit far from zero against their span
value = "signal"                         # `TO FILL`: the name of what was measured, and its unit below
value_unit = "V"                         # `TO FILL`: the unit of the value, or "1" where it has none

[covariates.setting]                     # `TO FILL`: each setting under which a curve was taken, with its unit, one table each
unit = "none"                            # `TO FILL`

[[trace]]                                # `TO FILL`: one table a curve: its file in traces/, its two columns, its settings
file = "traces/one.csv"
columns = ["x", "signal"]
setting = 1

[noise]                                  # the standard deviation of the noise of each sample, given in one of five ways:
from = "neighbours"                      # "neighbours", the default where [noise] is left out: a law in the height of the curve above its foot, a + b h + c h^2, read from the second differences of neighboring samples before any fit
                                         # "residuals": the same law, read from what the first model leaves, for a curve sharper than three samples, whose second differences hold the curve
                                         # "value", with a line value = 0.01: one number for every sample, in the unit of the value, where you know it
                                         # "column", with a line column = "sd": a third column of each file
                                         # "model": a law of yours, a function variance(mean, record, own) in model.py that gives the variance of each sample

# [setup.detector]                       # optional, one table each: an instrument or a setting of the data collection, in words, with its value,
# words = "the gain of the detector"     # its uncertainty and its unit where they are known, and how it is known: "measured", "set", "from its
# value = 1e6                            # maker", "assumed" or "not recorded". The package reads no physics in it: whoever chooses the terms
# unit = "V/A"                           # starts from it, and a lever of kind "setup" names one of them
# known = "set"

[addenda]                                # all optional: with no model the loop starts from the package's generic model
model = "model.py"                       # the first model
list_of_terms = "list_of_terms.py"       # optional: the terms that you suspect, which the report sets against the residuals. Empty at the start

[addenda.lever.quieter]                  # optional: a change of the experiment, which the three proposals for the experiment weigh against the others
kind = "noise"
factor = 0.5
words = "a quieter detector, which halves the noise of each sample"
```

`my_project/model.py`:

```python
"""The first model of this project, which `python -m twinfit loop` reads. Fill the places marked TO FILL.

Start with the small model that is understood in every part. The loop fits it to every curve at once, reads what it
leaves against its twin, and says where a term is missing: each term that you then write in terms/ is fitted beside it.
"""
import numpy as np

SHARED = ("width",)                     # TO FILL: the parameters that every curve shares, which the analysis is after
OWN = ("height", "place", "baseline")   # TO FILL: the parameters of each curve alone
START_SHARED = [1.0]                    # TO FILL: a start for each shared parameter, in the order of SHARED
BOUNDS_SHARED = ([1e-9], [1e6])         # optional: the lower bounds and the upper bounds of the shared parameters, a width above zero since a curve divides by it
CENTRE = "place"                        # an own parameter, which `centre` below goes with: the report reads the core and the wings of a curve about the place that `centre` returns. Where the place is shared, as one that a setting moves, name any own parameter here, as the height, and let `centre` return the place from `shared`
NAME = "one peak of one width"          # TO FILL: the model in a few words, as the report names it
EQUATION = ""                           # the model in LaTeX between dollar signs, as r"$h\,e^{-(x - x_0)^2/(2 w^2)} + b$", which the visual log typesets at each turn, `python -m twinfit log`. Left empty, the log writes the Python line that `model` returns
LETTERS = {}                            # each letter of EQUATION with its meaning, as {"h": "the height", "x_0": "the place", "w": "the width", "b": "the baseline"}, which the log writes under it


def model(x, shared, own, design):
    """TO FILL: the mean of one curve. `x` is its abscissa, `shared` and `own` hold the parameters in the order of
    SHARED and OWN, and `design` the settings of the curve, as the description gives them."""
    return own[0] * np.exp(-0.5 * ((x - own[1]) / shared[0]) ** 2) + own[2]


def start_own(record):
    """TO FILL: a start for the own parameters of one curve, read from its values, in the order of OWN. This one reads a
    peak or a dip on a baseline that may slope: the baseline as a straight line through the means of the first and the
    last twentieth of the curve, the place as the first sample from the left whose distance from that line reaches four
    fifths of the largest, that distance as the height, negative for a dip, and the line there as the baseline. So every
    curve of a doublet starts on its lower line, and a dip on a slope is not taken for an end of the slope, where the
    median as a baseline took the end of a slope as steep as the dip. A curve started away
    from its peak or its dip is moved by no shared parameter, and the first fit is then refused."""
    order = np.argsort(record.x)
    x, y = np.asarray(record.x, dtype=float)[order], np.asarray(record.y, dtype=float)[order]
    n = max(len(y) // 20, 1)
    left, right = float(np.mean(y[:n])), float(np.mean(y[-n:]))
    line = left + (right - left) * (x - x[0]) / ((x[-1] - x[0]) or 1.0)
    up = np.abs(y - line)
    top = int(np.flatnonzero(up >= 0.8 * up.max())[0])
    return [float(y[top] - line[top]), float(x[top]), float(line[top])]


def start_shared(records):
    """optional, and where it is defined it takes the place of START_SHARED: a start for each shared parameter, read
    from the curves. This one reads the width of the peak or the dip of each curve where its distance from the straight
    line through its ends is half the largest, and gives their mean as a standard deviation. Change it or remove it when
    SHARED changes."""
    widths = []
    for record in records:
        order = np.argsort(record.x)
        x, y = np.asarray(record.x, dtype=float)[order], np.asarray(record.y, dtype=float)[order]
        n = max(len(y) // 20, 1)
        left, right = float(np.mean(y[:n])), float(np.mean(y[-n:]))
        line = left + (right - left) * (x - x[0]) / ((x[-1] - x[0]) or 1.0)
        up = np.abs(y - line)
        above = x[up >= 0.5 * up.max()] if up.max() > 0 else x
        widths.append((above.max() - above.min()) / 2.3548 if above.size > 1 and above.max() > above.min() else float(np.ptp(x)) / 10.0)
    return [float(np.mean(widths))]


def centre(shared, own, design):
    return own[1]
```

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

```toml
[noise]                          # one of five ways
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
```

```python
HELD = {"resolution": 0.5}          # the values that the model reads and does not fit, which an ansatz sets


def model(x, shared, own, design):
    width = (shared[0] ** 2 + HELD["resolution"] ** 2) ** 0.5
    return own[0] / (1.0 + ((x - own[1]) / width) ** 2) + own[2]
```

## 10. A term of unknown distance or period

Some terms carry a number that is chosen and not fitted: the distance of a second component, the stretch of a tail, the period of a ripple. The worked project of chromatography has one:

```bash
python -m twinfit example chromatography my_chrom
```

```bash
python -m twinfit loop my_chrom
```

```bash
python -m twinfit scan my_chrom "second(x)" 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.4
```

`my_chrom/terms/01.py`:

```python
from twinfit import mechanisms
term = mechanisms.second(1.4)
why = "the report: a second component in the asymmetry. An impurity that elutes later, at the distance that the scan gives"
```

```bash
python -m twinfit loop my_chrom
```

```python
HELD = {"distance": 1.4}          # the distance of the impurity, held by the ansatz of the description
```

```toml
[addenda.ansatz.distance]
value = 1.4
uncertainty = 0.2
```

`my_chrom/terms/01.py`:

```python
from twinfit import held, mechanisms
term = held.reading(lambda: mechanisms.second(HELD["distance"]))
why = "the report: a second component in the asymmetry, at the distance that the scan gives and that the ansatz holds"
```

```bash
python -m twinfit loop my_chrom
```

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

## 11. A term written by the researcher

When the term that the report points to is not among those of `python -m twinfit terms`, write it: `python -m twinfit card` gives its form, and section 11 of `TUTORIAL.md` its five steps. The file below is `terms/01.py` of a copy of the folder of section 10.

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

`my_nmr/terms/03.py`:

```python
import numpy as np
from twinfit.mechanisms import Mechanism


def wrap(model, problem, n):
    def new(x, shared, own, design):
        middle = 0.5 * (x[0] + x[-1])
        slope = shared[n] * (design["temperature"] - 305.0)
        return model(x, shared[:n], own, design) + slope * (x - middle)
    return new


term = Mechanism("a slope of the baseline that grows with the temperature", "baseline", "slope per kelvin", 0.0, wrap)
why = "a test of a term of my own: a drift of the baseline that a warmer probe could leave"
```

```bash
python -m twinfit loop my_nmr
```

`my_nmr/terms/04.py`:

```python
import numpy as np
from twinfit.mechanisms import Mechanism


def term_at(period):
    def wrap(model, problem, n):
        def new(x, shared, own, design):
            return model(x, shared[:n], own, design) + shared[n] * np.cos(2.0 * np.pi * x / period)
        return new
    return Mechanism(f"a fringe fixed to the abscissa, of period {period:g}", "periodic", f"fringe of period {period:g}", 0.0, wrap)


term = term_at(3.0)
why = "a test of a term of my own with a chosen number: a fringe that the optics could leave"
```

```bash
python -m twinfit scan my_nmr terms/04.py 2 3 4
```

## 12. With an agent

At step 5 you choose among its three proposals, or let it choose. Every term that it adds carries its reason in `why`: the physics that could leave the structure, offered as a candidate, which the loop and its twin then bear out or refuse. The package names no cause, and the agent names one only as a candidate.

## 13. More examples, of other fields

## 14. Limits of the package
