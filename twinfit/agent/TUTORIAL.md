# twinfit for an agent: the tutorial in commands

Written by `python -m twinfit.tutorial --write` from TUTORIAL.md, the tutorial for people, which holds the reasons, the outputs and the pictures: every command and every file here is taken from it.

It builds on the card, `agent/AGENT.md`, and sets out every command in a block and every file that the page writes, those that `example` writes standing in the folders that it makes, of the two worked projects. Not covered here: the reasons, the outputs and the pictures, which `TUTORIAL.md` gives.

> Section 1 of `TUTORIAL.md` defines every word of the package that this page uses.

The rules that decide each turn:

- The output of the loop for an agent is `DIR/loop/brief_K.md`, one page a model K, with `[output] for` holding "agents",
  as it does by default. `model_K.md` is the page for people and `model_K.json` holds every number.
- A term is asked only where an answer of the report stands beyond the twin. The line is set so that a right model passes
  it on 5 sets of 100, and at the default of 30 sets a run of a model that lacks nothing reads some answer past it more
  often, as section 14 of `TUTORIAL.md` measures: one answer just past it is weak evidence.
- A gain in chi square ranks the directions of the report and is no forecast of what a term brings. The loop asks for
  no term in four cases, which the line "Next: stop" of the last brief names: a kind of the noise leads, one answer
  alone stands just past its line, the last term file was refused and no proposal beyond the twin gains what a term
  costs, or no answer stands beyond the twin.
- After each run, read in this order: its exit code, 1 where a term file was refused before any fit, with
  no new `best.md`, and `DIR/loop/REFUSED.md` beside the files of the run before it where one stood. Then its line `kept: ...`, the lines "Answers beyond the twin" and "Next"
  of the last `brief_K.md`, whether its `model_K.md` says "The twin admits the fit" and the check of the noise passes,
  and each "Not yet to be quoted" of `best.md`, which writes no line that a value may be quoted: the brief's "where
  best.md says that it may be quoted" means a value with no such mark.
- Of the five parts, the bar of the fit, the bar of theory and the half range over the models are the uncertainty,
  added in squares where one number is asked. The coverage and the bias on the twin are checks of the bar of the fit.
- An own value that moves with a setting is read first as what a missing shared term leaves. Write the shared term
  that the report asks first: an own value that still moves after it is a term of its own.
- `loop` takes `--n N`, `--sets S`, `--max-seconds T`, `--workers W` and `--dictionary`. A run that names no `--sets`
  reads 30 sets: the run whose numbers you quote is the last, at the sets that you quote it at. A change of the noise,
  of `model.py`, or a kept term taken out, is made in a copy made without the runs, `rsync -a --exclude loop --exclude
  loop.before --exclude loop.reuse --exclude loop.costs.tsv --exclude agent --exclude .twinfit.lock --exclude
  .twinfit.command DIR/ COPY/`. `example` writes no `agent/`: make it with `mkdir -p DIR/agent`.
- The loop keeps a fitted term when its fit ends inside its bounds, the Bayesian criterion falls by more than 2 and the
  twin admits the larger model, and best.md gives the reason of a refusal.
- A term that the physics asks and the records cannot show is held at its theory, with hold = True in its file.
- Quote a value only from a model that its twin admits, with the parts of its bar that best.md gives.
- `python -m twinfit wait DIR` waits for a running loop within ten minutes and prints the last brief: never poll a run
  in another way, and never start two in one folder.

## 1. Purpose and vocabulary of the package

You give it curves recorded across a design, such as spectra at several temperatures, and the model that should describe them.

| the word of this page | its name in the literature |
|---|---|
| a shared parameter, an own parameter | a global parameter of a joint fit, a parameter of one curve, a nuisance parameter where the analysis is not after it |
| a bar, first order | the standard error from the inverse of the curvature matrix of the chi square, under the linearization of the model about the best fit |
| the twin | a parametric bootstrap from the fitted model (Efron and Tibshirani 1993) |
| a scatter, the bias on the twin | the bootstrap standard deviation of an estimator, and its bootstrap estimate of the bias |
| the coverage on the twin | the coverage probability of the interval of 1.96 standard errors, read on the parametric bootstrap |
| a generalized moment, its kernel, a window | a moment of the residuals under a weight function, read as a diagnostic and not as the estimator of the generalized method of moments, the weight function, a sum truncated to a range of the abscissa or of the fitted value |
| the excess, the line | a statistic standardized on its bootstrap distribution by its median and its 95th percentile, and the 95th percentile, over the bootstrap sets, of the largest standardized statistic of each set over every answer of the report, a familywise critical value of the maximum statistic with the whole report as its family (Westfall and Young 1993) |
| the Bayesian criterion | the Bayesian information criterion (Schwarz 1978), a fall of more than 2 being positive evidence on the scale of Kass and Raftery (1995) |
| the gain in chi square | the likelihood ratio statistic of the model extended by one parameter along the shape that the residuals hold, for a normal noise of known variance |
| a profile | the profile of the criterion over a parameter held at each of its values, the values within 1 of its minimum giving one standard deviation for a normal noise (Wilks 1938) |
| an ansatz, its bar of theory | a parameter held at a value from outside the records, and its uncertainty carried as a systematic one, by the move of each value across that uncertainty |
| a chi function | the loss function of an M estimator (Huber 1964) |
| glitches, a noise correlated from sample to sample | outliers, a serially correlated noise |
| the shapes of the residuals: a periodic structure about each centre, fixed to the abscissa or on one side | the periodogram of the residuals by least squares on cosines and sines (Lomb 1976, Scargle 1982), in the frame of each curve, in the frame of the abscissa and on half of each curve |
| a part in quadrature with the line | the Hilbert transform of the fitted curve, the imaginary part of its analytic signal (Bracewell 2000) |
| a structure in one place and at one scale | the continuous wavelet transform of the residuals, with a wavelet of the kind of Morlet's (Torrence and Compo 1998) |
| each curve predicted by the fit of the others | the leave one out cross validation of the joint fit, read at first order (Stone 1974) |
| the chi functions in windows of position, of the value | a scan statistic over windows of position (Kulldorff 1997), and a test of the law of the noise by bands of its values, as the chi square test of fit (Pearson 1900), each calibrated by the parametric bootstrap |
| a correlated twin | a parametric bootstrap with a first order autoregressive noise, its coefficient read on the differences of the residuals |
| a rung | one model of a nested sequence, each term reducing to the model before it at a size of zero: a term that replaces the line, as the tail of section 11, is not nested |
| the rule of the twin | bands of acceptance for the bootstrap coverage and bias of each parameter |
| a chance | a p value read by Monte Carlo simulation |
| the five parts of the uncertainty | the standard error, its two bootstrap checks, a systematic by displacement and a systematic of the choice of model |
| the report as a whole | a test of the specification by weighted moments of the residuals, of the family of the conditional moment tests (Bierens 1990), calibrated by the parametric bootstrap |

Every statistical tool of the package, what it reads and where its answer stands:

| the tool | what it reads | where it stands |
|---|---|---|
| the twin | synthetic curves drawn from the fitted model with the noise law of the records, fitted and read as the records are | beside every answer of the report and every bar of `best.md` |
| every sample | the chi square of the fit, the misfit as a whole | the answer "seen by" |
| the generalized moments, of the curve, of its derivative and of its running integral, and their ratios | the residuals under powers, fractional powers, inverse and logarithmic kernels, cosines and sines and Hermite functions | the answer "seen by" |
| the windows along the abscissa and along the value, the two sides added and subtracted, the crossings | the residuals in pairs of intervals about the centre, and in bands of the fitted value | the answers "where" and "seen by" |
| the transformations of the value | the same families on the logarithm and other transformations of the value | the answer "seen by" |
| the records together | every record added in bins of the abscissa | the answer "seen by" |
| the chi functions over all the residuals, and the plane of position and residual value | the cost of each residual under eight laws, and the table of position against value | the answer "seen by" |
| the chi functions in windows of position and in bands of the value | thirteen costs, the square or the count beside the twelve of the chi functions in windows, under smooth windows about each centre or fixed to the abscissa, and in bands of the residual's value | the answer "seen by" |
| the chi functions in cells of position and of the level or of the value | every cost of the residual in each cell of a smooth window of position times a window of the curve's level or of the residual's value, on a ladder of scales | the answer "seen by" |
| the fit under every chi function | the model fitted again under seven costs, and each shared parameter's move from the fit by the square in the spread that a right model gives it | the answer "seen by" |
| the running sum | the largest excursion of each curve's running sum of its residuals | the answer "seen by" |
| the Fourier transform and the differences of the residuals | every cost of the residuals' Fourier coefficients in windows of frequency, and of their steps from sample to sample, `where.fourier_chi` and `where.difference_chi`, beside the report | a researcher's call |
| each curve predicted by the fit of the others | the rise of each curve's chi square when the shared parameters are fitted without it | the answer "seen by" |
| the shapes of the residuals | a slow trend, periodic structures in two frames and on each side, a part in quadrature with the line and a wavelet | the answer "shape" and the proposals |
| the shape of the noise | skewness, heavy tails, glitches, a noise that grows with the signal and a noise correlated from sample to sample | the answers "what" and "noise" |
| what changes along a setting | each place of the residuals read against each setting of the design | the answer "changes along" |
| the noise matrix | the chance that each statistic reports each of five laws of the noise on your records | `python -m twinfit matrix DIR` |
| the proposals | three directions for the model, each with its gain in chi square and its cost on the twin, and three for the experiment | the report of each model |
| the known parameters before the term is named | the estimates of each shared parameter that a term of the class that leads does not spoil | the report of each model |
| the Bayesian criterion and the profile | whether a term is kept, and a number that is chosen and not fitted | `best.md`, the brief of each model and the log, and the profile by `python -m twinfit scan DIR FILE VALUES` |
| the five parts of the uncertainty | the bar of the fit, its coverage and its bias on the twin, the bar of theory and the spread over the models | `best.md` |
| the combined estimate | each shared parameter estimated by every family that gives one, beside the fit, a check of the model | `best.md` |
| what the records tell apart | the combinations of parameters correlated beyond 0.9 that the records hold, each with the bar to quote | `best.md` |
| the values in windows | each value read in windows about the centre and off the peak | `best.md` |
| the own values along a setting | an own parameter that moves from curve to curve along a setting of the design | `best.md` |
| a correlated twin | the twin drawing each curve's noise with the correlation that its residuals give | `correlated = true` under `[noise]` |
| the traces for a target bar | how many times each record is taken so that a parameter's bar reaches a target, checked on the twin | `python -m twinfit traces DIR name=bar` |

## 2. Requirements and installation

Python 3.12 or later, with numpy 2.5 or later and scipy 1.16 or later, and matplotlib for the visual log of a run.

```bash
pip install "twinfit[log]"
```

```bash
python -c "import twinfit; print(twinfit.__version__)"
```

Its first line, on this project: `0.1.3`

## 3. A worked project

This writes a project of nuclear magnetic resonance into a new folder, `my_nmr`:

```bash
python -m twinfit example nmr my_nmr
```

Its first line, on this project: `written into my_nmr: description.toml, model.py, an empty terms/, and 4 curves in traces/. Next: python -m twinfit loop my_nmr`

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

The loop runs here in the foreground, since each run of this section takes about a minute.

```bash
python -m twinfit loop my_nmr
```

Its line of what it kept, on this project: `kept: nothing. Refused: nothing. The best model: model 0.`

```bash
python -m twinfit matrix my_nmr
```

## 5. A term

A term is a file in `terms/`, `01.py`, then `02.py`, and so on, which you write. This prints the terms that a term file can name, each but `exponential_width()` described with the one parameter that it adds, and the form of a term file:

```bash
python -m twinfit terms
```

Its first line, on this project: `The terms of the package, each with the one parameter a that it adds:`

Save this as `my_nmr/terms/01.py`:

```python
from twinfit import mechanisms
term = mechanisms.along("width", "temperature")
why = "the report: the core and the flanks change along temperature, and the width of a resonance changes with the temperature of the sample"
```

## 6. The second run

```bash
python -m twinfit loop my_nmr
```

Its line of what it kept, on this project: `kept: 01.py, width changes along temperature. Refused: nothing. The best model: model 1.`

The run's output uses two words: your list, `list_of_terms.py`, the terms that you suspect, which `init` writes empty and section 9 shows, and which a project of `example` does not hold, so that its list is empty, and matched, the question of the report that reads the residuals against each term of that list.

```bash
python -m twinfit loop my_nmr --sets 100
```

Its line of what it kept, on this project: `kept: 01.py, width changes along temperature. Refused: nothing. The best model: model 1.`

```bash
python -m twinfit log my_nmr
```

Its first line, on this project: `The visual log: my_nmr/loop/log.pdf, 29 slides, each phase of each of 2 turns first as one plot, then the equations, the effects of the terms, the curves in a row, the report and the criteria, the pictures in my_nmr/loop/log.`

1. **The curves**: the curves and their samples, the rows of `loop/residuals_K.csv`, one plot a sweep of the description, three sweeps or six levels of a grid at most.
2. **The model**: the bar of each of its first ten shared parameters in per cent of its value, `per_cent_of_its_value` of the phase `identify`, and the largest of its `correlations`.
3. **The twin**: each of its first eight shared parameters' `bias_over_spread`, `spread_over_bound` and `coverage` against the `bands` of the rule, under `admission` of the phase `twin`.
4. **The fit**: `chi_square` over `degrees_of_freedom` of `loop/model_K.json`, a right model 1 within 1.96 times the root of 2 over the degrees, times the root of `samples_of_the_correlation` of the phase `read` where its `correlation_of_neighbours` passes 0.3.
5. **What the fit leaves**: the mean square of each curve's `residual` over its `noise` in `loop/residuals_K.csv`, a right model 1.
6. **Where to search**: `where.first_of_each_question` of `loop/model_K.json` names the first answer of each question, its excess is `excess` of the same answer in `where.answers`, against `where.largest_excess_on_the_twin`. Where `where.where_read_along` names a place along a setting, as on model 0 here, no place stands beyond the twin alone, and the slide reads that answer of the question changes along and its excess.
7. **The next term**: `gain` of the phase `candidate on the twin`, in the notebook `NNN_model_K_candidate_on_the_twin.json`, `on_the_records` beside its `mean` and `scatter` on the twin's sets.
8. **The decision**: the chi square over its degrees before and after the next term, from `loop/model_K.json` and the next model's, and on its second slide each model's `criterion`, read from its own `loop/model_K.json`, `loop/best.json` comparing them from a kept term on.

## 7. A term not asked for by the records

The loop is not a machine that keeps whatever it is given. Here is a term that no answer of the report points to, a ripple along the abscissa:

Save this as `my_nmr/terms/02.py`:

```python
from twinfit import mechanisms
term = mechanisms.periodic(3.0, "cosine")
why = "a test: no answer of the report points to a periodic function"
```

```bash
python -m twinfit loop my_nmr
```

Its line of what it kept, on this project: `kept: 01.py, width changes along temperature. Refused: 02.py, a cosine of period 3. The best model: model 1. best.md gives the reason of each refusal.`

## 8. The proposals for the experiment

Each report ends with three proposals for the experiment, each of one of three kinds: a better setup, a better characterization of the setup, or a better data collection.

```bash
python -m twinfit traces my_nmr width=0.0025
```

Its first line, on this project: `The model that the term files give: 01.py, on 4 records.`

## 9. The researcher's own curves

This writes the starting files of a project into a new folder:

```bash
python -m twinfit init my_project
```

Its first line, on this project: `written: my_project/description.toml`

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

`my_project/list_of_terms.py`:

```python
"""Your list of terms: the terms that you suspect, which the report sets against the residuals, each against its twin,
in the question matched. It is empty at the start, and the package names none: with it empty the report reads the
residuals alone. A term is written here as a term file writes it, a term of the package that `python -m twinfit terms`
lists or a Mechanism of your own, and the names of model.py are read here as a term file reads them."""
from twinfit import mechanisms  # noqa: F401

TERMS = []
```

A `model.py` whose place moves along the temperature, written again whole for `my_nmr_place/`, a copy of `my_nmr`:

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

The keys of `description.toml`, each table with its form:

```toml
[noise]                          # one of five ways
from = "value"                   # one standard deviation for every sample
value = 0.02
# from = "column"                # or a column of the curves that gives the standard deviation of each sample
# column = "sd"
# from = "neighbours"            # or the default, a law in the height of the curve, read from the differences of neighboring samples before any fit
# from = "residuals"             # or the same law, read from what the first model leaves, for a curve sharper than three samples
# from = "model"                 # or a law of yours: variance(mean, record, own) in model.py, the variance of each sample
# correlated = true              # with any of them, where the instrument joins neighbouring samples: the twin draws each curve's noise correlated as its residuals' differences say. A ripple of a few samples reads as such a correlation, so seek it first

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

A `model.py` that reads a held value, its `HELD` and its model:

```python
HELD = {"resolution": 0.5}          # the values that the model reads and does not fit, which an ansatz sets


def model(x, shared, own, design):
    width = shared[0] + HELD["resolution"]          # a Lorentzian seen through a Lorentzian resolution: their widths add, where a Gaussian resolution makes a Voigt line
    return own[0] / (1.0 + (2.0 * (x - own[1]) / width) ** 2) + own[2]
```

## 10. A term of unknown distance or period

Some terms carry a number that is chosen and not fitted: the distance of a second component, the stretch of a tail, the period of a ripple. The worked project of chromatography has one:

```bash
python -m twinfit example chromatography my_chrom
```

Its first line, on this project: `written into my_chrom: description.toml, model.py, an empty terms/, and 4 curves in traces/. Next: python -m twinfit loop my_chrom`

```bash
python -m twinfit loop my_chrom
```

Its line of what it kept, on this project: `kept: nothing. Refused: nothing. The best model: model 0.`

```bash
python -m twinfit scan my_chrom "second(x)" 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.4
```

```bash
python -m twinfit scan my_chrom "second(x)" 1 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9
```

Save this as `my_chrom/terms/01.py`:

```python
from twinfit import mechanisms
term = mechanisms.second(1.45)
why = "the report: a second component in the asymmetry. An impurity that elutes later, at the distance that the scan gives"
```

```bash
python -m twinfit loop my_chrom
```

Its line of what it kept, on this project: `kept: 01.py, a second component 1.45 away. Refused: nothing. The best model: model 1.`

Add at the end of `my_chrom/model.py`:

```python
HELD = {"distance": 1.45}          # the distance of the impurity, held by the ansatz of the description
```

Add at the end of `my_chrom/description.toml`:

```toml
[addenda.ansatz.distance]
value = 1.45
uncertainty = 0.25
```

Save this as `my_chrom/terms/01.py`:

```python
from twinfit import held, mechanisms
term = held.reading(lambda: mechanisms.second(HELD["distance"]))
why = "the report: a second component in the asymmetry, at the distance that the scan gives and that the ansatz holds"
```

```bash
python -m twinfit loop my_chrom
```

Its line of what it kept, on this project: `kept: 01.py, a second component 1.45 away. Refused: nothing. The best model: model 1.`

A part of a `model.py` of your own for spectra of galaxies in km/s, which no example holds:

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

When the term that the report points to is not among those of `python -m twinfit terms`, write it: `python -m twinfit card` gives its form, and section 11 of `TUTORIAL.md` its five steps.

1. **The term.** A term is a `Mechanism(name, family, parameter, start, wrap, lower, upper, formula=...)`: a name, a word for its family, any word, as `baseline`, which the package keeps with the term and reads in no result, the name of the one parameter that it adds, where the search of that parameter starts, `wrap`, the bounds of the parameter, none where `lower` and `upper` are not given, and `formula`, the term in LaTeX between dollar signs with `f` the model before it and `a` its parameter, as `r"$f(x) + a\,x$"`, which the visual log typesets: without it the log says that its formula is not written. `python -m twinfit card` prints a term file of each kind in one page. The name of the parameter is printed as it is written, so it carries its unit, as "shift of the centre per microwatt, in Hz per microwatt".
2. **Its wrap.** `wrap(model, problem, n)` takes the model as it stands, the problem, which holds every curve, and the number of shared parameters before the term. It returns the model with the term, a function `new(x, shared, own, design)`, where `shared[n]` is the new parameter and `shared[:n]` the parameters before it. The loop refuses a `wrap` that does not take these three arguments. `design` holds the settings of the curve, as the description names them, `own[k]` and `shared[k]` the parameters of the model in the order of `OWN` and `SHARED` of `model.py`, and `x`, `shared` and `own` are arrays of numpy. A term that moves the abscissa calls the model at a moved `x`, as `model(x - shared[n] * design["power"], shared[:n], own, design)`, which moves each curve by `shared[n]` times its power. A held number is read inside `new` as `HELD["name"]`, since the loop gives every term file the `HELD` of `model.py`. A term that moves the place leaves `centre` where `model.py` puts it, so a place that moves is written in `model.py`, as the `model.py` of section 9 whose place moves along the temperature does.
3. **One of three forms.** `new` may add a term to the model, as a step under a peak whose size is a share of its height, `model(x, shared[:n], own, design) + shared[n] * own[0] * 0.5 * (1.0 + erf((x - own[1]) / (2.0 ** 0.5 * shared[0])))`, with `erf` of `scipy.special`. It may multiply the model, as a gain of the detector that changes along the abscissa, `model(x, shared[:n], own, design) * (1.0 + shared[n] * (x - middle))`. That multiplies the terms kept before it too, about a middle that you choose as the zero of the change. Or it may replace the line and return a curve of its own, with no call of `model`: then it returns the whole curve, its baseline too, since what it does not add again is gone. The loop refuses such a term where its curve no longer reads an own parameter that the model read, and says which.

4. **A curve that replaces the line**, written out whole, as `01.py` of a fresh project of chromatography, `my_chrom_tail`, which holds no other term file. It is a Gaussian with an exponential tail, for the chromatograms of section 10, whose model holds the width `s = shared[0]`, the height `own[0]`, the place `own[1]` and the baseline `own[2]`. It is the other form of the asymmetry that section 10 read, a tail of the column where section 10 wrote an impurity, and it is run here beside it:

```bash
python -m twinfit example chromatography my_chrom_tail
```

Its first line, on this project: `written into my_chrom_tail: description.toml, model.py, an empty terms/, and 4 curves in traces/. Next: python -m twinfit loop my_chrom_tail`

Save this as `my_chrom_tail/terms/01.py`:

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

```bash
python -m twinfit loop my_chrom_tail
```

Its line of what it kept, on this project: `kept: 01.py, an exponential tail of the line. Refused: nothing. The best model: model 1.`

5. **Its start and its place among the files.** Start the new parameter at a size that the residuals give, and a parameter that the curve divides by, as `t`, near the width of the curve and away from zero. A term reads the shared parameters as the terms before it leave them: a tail in `02.py` after `along("width", ...)` in `01.py` reads the width that `along` makes at each setting. A term that replaces the line drops the terms kept before it, so it goes first, as `01.py`.

Save this as `my_nmr/terms/03.py`:

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

Its line of what it kept, on this project: `kept: 01.py, width changes along temperature. Refused: 02.py, a cosine of period 3; 03.py, a slope of the baseline that grows with the temperature. The best model: model 1. best.md gives the reason of each refusal.`

Save this as `my_nmr/terms/04.py`:

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

The protocol of an agent's turn is `agent/AGENT.md`, which `python -m twinfit guide` prints.

## 13. More examples, of other fields

`python -m twinfit guide examples` prints the gallery, `twinfit/EXAMPLES.md`: projects of other fields, each written by `python -m twinfit example NAME DIR` with a truth that its folder does not hold.

```bash
python -m twinfit example satellite my_satellite
```

Its first line, on this project: `written into my_satellite, its noise white: description.toml, model.py, an empty terms/ and the curves in traces/. Its truth and what the loop gives are in twinfit/EXAMPLES.md. Next: python -m twinfit loop my_satellite`

```bash
python -m twinfit loop my_satellite --n 0
```

Its line of what it kept, on this project: `kept: nothing. Refused: nothing. The best model: model 0.`

## 14. Limits of the package

- It names a place, a kind and a setting of what is missing, never its cause.
- The twin draws its synthetic curves at the noise that the description gives, so it takes that noise as right.
- With 30 synthetic sets, the default, an answer can stand beyond all of them and no further: its chance is then the smallest that 30 sets give.
- On models that lack nothing, at the default of 30 sets with its twin's seed drawn from the records, a run of 0.1.3 read some answer past its line on 0.080 ± 0.017 of 250 draws, as `README.md` of the package says, where the line is set for 0.05 by its construction.
- A term much smaller than the noise can be missing while no answer stands beyond the twin.
- Glitches and a noise that grows with the signal are kinds of the report that no term of the mean fills.
- A proposal for the experiment is read on the synthetic sets of the twin, so its ratio moves from run to run: the page says within how much, about one over the root of the number of sets.
- A term in the wrong frame leaves a misfit of its own.
- Two forms of one effect, as a Lorentzian added to a line and a Voigt line in its place, are not set side by side by one run of the loop, which fits each term file on top of the ones before it.
- One misfit can light several answers at once: a slope of the baseline lights a slow trend over the span of each curve and a periodic structure along the abscissa together.
- The bars of the own parameters of each curve count what the shared ones leave uncertain, but no twin reads them, as it reads the shared ones.
- `along` multiplies a parameter by one plus its size times the setting, so it cannot move a parameter that sits near zero, as a place read less the offset of its abscissa can.
- `second` takes one distance and copies the curve there, at its own width: a second line whose distance or share follows a setting is written in `model.py`, as a function of the settings, as a place is.
- A lever along a setting reads the model at settings past those of the records, and holds as far as the model holds there.
- The loop can end in words that the runs above do not show.
