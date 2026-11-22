# The tutorial of twinfit

This tutorial takes you from nothing to the best model of a set of curves, with the uncertainty of each of its parameters and six proposals, three for the model and three for the experiment. It runs on two worked projects that the package writes, one of nuclear magnetic resonance and one of chromatography, so nothing of yours is needed to follow it. Their curves were drawn from a known truth, so you can see whether the package finds it.

If you already fit your curves with a model that you trust for a precision measurement, section 9 runs the package on your own curves and your own `model.py`: the report reads that model against its twin and says where its residuals hold a structure that the twin does not, its region, kind and setting, how far each shared parameter moves under other chi functions, and which parameters the curves cannot tell apart, the checks that its chi square does not make. The worked projects below teach the words and the reports first.

It builds on the five minutes of the README and sets out the whole path of the package, its words, its reports and its choices. Not covered here: the protocol of an agent, which `agent/AGENT.md` gives, and the projects of other fields, which `EXAMPLES.md` gives.

> Section 1 below says what the package is for and defines its words: in a first list the words that the first run needs, then the others in a glossary, and a table sets those that the literature of statistics names otherwise beside their names there.

How to read this page: read section 1 to the end of its first list, then sections 2 to 6, the shortest path through the package. They run the worked project of nuclear magnetic resonance, write its first term in section 5 and read its best model in section 6. Look a word up in the glossary of section 1 when you meet it. The other sections answer one question each.

| to | read |
|---|---|
| look up a word of the package | the glossary of section 1 |
| install the package and check its version | section 2 |
| see the whole path once, on curves drawn from a known truth | sections 3 to 6 |
| see the loop refuse a term that the records do not ask for | section 7 |
| read the proposals for the experiment | section 8 |
| run the loop on your own curves | section 9 |
| fit a term whose number is chosen and not fitted, a distance or a period | section 10 |
| write a term of your own in Python | section 11 |
| hand the loop to an agent | section 12, and `agent/AGENT.md` |
| find a project of another field whose first report reads as yours | section 13 |
| know what the package cannot do | section 14 |

Every command in a block of its own below, but the installation of section 2, was run as it is written, in an empty folder, and every output is what it printed. A command named inside a sentence, as `nohup` or `rsync`, is not run by the check. `python -m twinfit.tutorial --check`, with matplotlib installed by `pip install "twinfit[log]"`, runs them all again, in a temporary folder, and says whether any output changed. Here the whole of it ran in 266 seconds. A machine busy with other work takes longer.

## 1. Purpose and vocabulary of the package

You give it curves recorded across a design, such as spectra at several temperatures, and the model that should describe them. It fits the model to every curve at once. It reads the structure of what the fit leaves, the residuals, with families of statistics: generalized moments of the residuals and of their transformations, and windows along the abscissa and along the value. It reads how that structure changes from one setting of the design to the next.

It sets every reading against a twin: synthetic curves drawn from the fitted model with the noise of the records, fitted and read the same way. The twin draws the noise of each sample from a normal law of the deviation that the description gives it, apart from the others, unless the description writes `correlated = true` under `[noise]`. Where a counter gives few counts, its noise is skewed, and the twin does not hold that.

The loop, `python -m twinfit loop`, does all of it on a folder. It fits the first model, reads it against its twin and writes a report of where to search. Then it fits each term that you write in turn, keeps or refuses it, and writes the best model with the uncertainty of each parameter. Where the records stand beyond the twin, the model lacks something, and the report says where in the curves, of which kind, and along which setting. You write a term in answer. The package fits it, keeps or refuses it, and says what it does to the other parameters.

The package names a place, a kind and a setting. It never names a cause: which term of physical origin fills them is your choice, or the choice of an agent in your place.

The words that the first run needs, the others of sections 4 to 11 in the glossary below:

- **A curve** is one recording, a spectrum or a chromatogram, as a column of abscissas and one of values. The description calls its file a trace, and the package calls a curve with its noise a record: the three words name one thing.
- **A setting** is a quantity under which a curve was taken and that the description gives for it, as the temperature of the sample. The package calls a setting that changes from curve to curve a covariate. The settings of one curve, by name, are its `design`, the dictionary that the model receives.
- **A bar** is the standard deviation of a parameter that a fit gives from the curvature of its chi square. **First order** is the model taken as straight in its parameters about the fit, so that bar is a first order bar.
- **A shared parameter** is one value for every curve, as a width. **An own parameter** is one value for each curve, as its height. The value of an own parameter in one curve is an own value, and the pages say an own value moves where it changes from curve to curve beyond its bars. The analysis is after the shared parameters, and the reports call them the known parameters: known in the sense that the model names them, and not that their values are known.
- **The excess** of an answer is its statistic on the records less the median of the same statistic on the synthetic sets, over the span from that median to their 95th percentile. An answer is beyond the twin where its excess passes a line, the largest excess that the synthetic sets of a right model reach on 95 of 100 of them. Each report prints its line, read on the twin of its model. The line is itself read on the synthetic sets of the run: measured at 0.1.2 on 250 draws of models that lack nothing, before the chi functions were read, a run at the default of 30 sets with a seed for each draw, as 0.1.3 draws one from the records, read some answer past it on 0.068 ± 0.016 of them, as `README.md` of the package says, the count of 0.1.3 to come, where the line is set for 0.05 by its construction. So one answer just past its line is weak evidence.
- **The Bayesian criterion** of a model is its chi square plus the number of its fitted parameters times the logarithm of the number of samples. A model with one term more is kept when its fit ends inside its bounds, its criterion falls by more than 2, and its twin admits it, as section 7 says.

The glossary, the other words of the package, to look up when you meet them:

- **An ansatz** is a value held before the fit, with its uncertainty, from outside the records: section 9 writes one. **A lever** is a change of the experiment that the proposals weigh on the twin, as a quieter detector: section 8 reads them.
- **The coverage and the bias on the twin** of a parameter are read on synthetic sets drawn from the fitted model. The coverage is the share of the sets on which its interval of 1.96 bars of the fit holds the value that made them, near 0.95 where the bar is right. The bias is the mean of its fitted values less that value, in bars. Section 6 reads them.
- **A scatter** is the standard deviation of a parameter over the synthetic sets of the twin, the spread of the estimate. A report gives a bias on the twin in scatters, and `best.md` in bars of the fit.
- **The rule of the twin** admits the fit of a model where, on its synthetic sets, the coverage of each shared parameter lies in a band about 0.95 and its bias is small beside its scatter, as the report of each model states it. The band is three standard deviations of a coverage of 0.95 at the number of sets read, and 0.92 to 0.98 from 475 sets on. The bars of a model that its twin does not admit are not to be read.
- **A chance** is the share of the synthetic sets of a right model on which a reading goes at least as far as on the records: a p value read by simulation.
- **The level** of a curve is its baseline. In a table of readings, as those that `python -m twinfit matrix` prints, the level of a reading is its chance at first order, and a test at the level 0.05 is one set to pass a right model on 95 cases of 100.
- **A generalized moment** is a sum over the samples of a residual times a weight of the abscissa. The weight is a power of the distance from the center of the curve, a fractional power, an inverse or a logarithmic weight, a cosine or a sine. It is also taken of the derivative of the curve, of its running integral, of its logarithm, and as a ratio of two moments.
- **A window** keeps the samples of a pair of intervals at equal distances from the center, a window along the abscissa, or the samples whose fitted value lies between two levels of its height, a window along the value.
- **A kernel** is the weight of a generalized moment. The report names each family of statistics by its kernels and its windows, as "moments: powers 0 to 8", "inverse and logarithmic kernels", "the two sides subtracted: pairs of intervals" or "the crossings of nine levels". The last reads the abscissa where the curve crosses each of nine levels of its height.
- **Two kinds tie** where their excesses lie within a factor of 1.5 of each other, and **two places tie** where their excesses differ by less than 1.96 standard deviations of their difference, read on 30 sets of a twin that carries the misfit, the fit plus what it leaves averaged over 9 samples. The records then do not tell them apart, the report names both, and a term of either is as asked.
- **The five parts of the uncertainty** of a value are its bar of the fit, its coverage on the twin and its bias on the twin, with two more. Its bar of theory comes where values are held, and the spread over the models comes from the models that the records cannot tell apart. The bar of the fit, the bar of theory and the half range over the models are the uncertainty, quoted side by side as `best.md` writes them and added in squares where one number is asked. The coverage and the bias are checks of the bar of the fit. Section 6 reads them.
- **The gain in chi square** of a proposal is how far the chi square falls when the shape that the residuals hold there is added to the model with one free size and every parameter fitted again. It ranks the directions beyond the twin, and it is no forecast of the fall that a term of your own brings, which can be far larger: section 7 sets the two side by side.
- **A profile** of a number that is chosen and not fitted, as the distance of a second component, is the best criterion at each value of it: the values within 1 of the best run from about one standard deviation below it to one above, where the records fix it well. `python -m twinfit scan` prints it as a table of the criterion at each value, section 10.
- **A Mechanism** is a term written in Python, section 11. **A rung** is one model of the loop: model 0, the first, and one rung more for each term kept on it.
- **Glitches** are a few samples far from the curve in one record. **A noise correlated from sample to sample** is a noise that moves neighboring samples together. The report names both among the kinds of a misfit, and no term of the curve fills them: they are read in the noise.
- **A chi function** is the cost given to each residual over its noise: its square in a plain fit, another function of it in a fit that weighs large residuals less, or more. Where the residuals are those of the noise, the fits agree on average, and the fit of the square scatters least. The report averages eight chi functions over the residuals of the fit and maps each so that, for residuals near the normal law, each grows about as their square or slower, each a statistic read against the twin.
- **The shapes of the residuals** that the report reads are a slow trend over the span of each curve, a periodic structure about each centre, one fixed to the abscissa and in phase from curve to curve, one on one side of the curves, a part in quadrature with the line, the Hilbert transform of the fitted curve, and a structure in one place and at one scale, a cosine and a sine of one cycle under a Gaussian. Each is scored on what the fit cannot absorb and set against the twin, and a proposal of the model is measured with the one that leads.
- **Each curve predicted by the fit of the others** is the rise of each curve's chi square when the shared parameters are fitted without it, at first order, summed over the curves: a term that changes from curve to curve is read there, and one that every curve holds alike is not.
- **The chi functions in windows** are the square, or the count in bands of the value, beside each chi function and four more, the powers 0.5 and 6 and the redescending costs of Tukey and of Welsch, thirteen costs, summed under a smooth window of position, about the centres or fixed to the abscissa, against the cost's mean over all the residuals, or over the residuals of one band of their value against what a normal law gives in their robust scale: a cost gathered in one place, or a law of the noise wrong in one band of its values, stands out, which a sum over every residual spreads thin. Each family is one answer, its largest over every window, read against the twin's largest.
- **A correlated twin**, `correlated = true` under `[noise]`, draws each curve's noise correlated from sample to sample, with the coefficient that the differences of its residuals give. A ripple of a few samples reads as such a correlation, so seek one first.
- **A budget** is the seconds of the clock that `--max-seconds` gives a run of the loop: a run past it stops at the end of the step in course, and the next run goes on where it stopped.

The words of this page beside their names in the literature of statistics, for a reader who knows those:

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

Bierens, H. J. (1990), A consistent conditional moment test of functional form, Econometrica 58, 1443. Bracewell, R. N. (2000), The Fourier Transform and Its Applications, third edition, McGraw-Hill. Efron, B. and Tibshirani, R. J. (1993), An Introduction to the Bootstrap, Chapman and Hall. Huber, P. J. (1964), Robust estimation of a location parameter, Annals of Mathematical Statistics 35, 73. Kass, R. E. and Raftery, A. E. (1995), Bayes factors, Journal of the American Statistical Association 90, 773. Kulldorff, M. (1997), A spatial scan statistic, Communications in Statistics: Theory and Methods 26, 1481. Lomb, N. R. (1976), Least-squares frequency analysis of unequally spaced data, Astrophysics and Space Science 39, 447. Pearson, K. (1900), On the criterion that a given system of deviations from the probable in the case of a correlated system of variables is such that it can be reasonably supposed to have arisen from random sampling, Philosophical Magazine 50, 157. Scargle, J. D. (1982), Studies in astronomical time series analysis. II. Statistical aspects of spectral analysis of unevenly spaced data, Astrophysical Journal 263, 835. Schwarz, G. (1978), Estimating the dimension of a model, Annals of Statistics 6, 461. Stone, M. (1974), Cross-validatory choice and assessment of statistical predictions, Journal of the Royal Statistical Society B 36, 111. Torrence, C. and Compo, G. P. (1998), A practical guide to wavelet analysis, Bulletin of the American Meteorological Society 79, 61. Westfall, P. H. and Young, S. S. (1993), Resampling Based Multiple Testing, Wiley. Wilks, S. S. (1938), The large sample distribution of the likelihood ratio for testing composite hypotheses, Annals of Mathematical Statistics 9, 60.

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

Python 3.12 or later, with numpy 2.5 or later and scipy 1.16 or later, and matplotlib for the visual log of a run. numpy 2.5 convolves 12.5 times faster than numpy 1.26 on a Mac with an arm64 processor, 0.80 ms against 10.0 ms, and 1.2 times slower on Linux with an x86_64 processor, 2.11 ms against 1.74 ms, the times that the comment of `pyproject.toml` holds. The oldest versions that the package has run on are Python 3.12, numpy 2.5.0 and scipy 1.16.0. Those floors are set by the speed of numpy on the Mac: the package itself starts on Python 3.11 with numpy and scipy, and refuses an older Python.

The loop runs on macOS and on Linux. It locks its files with `fcntl`, which Windows does not have, so on Windows run it under the Windows Subsystem for Linux.

Install it first, with matplotlib for the visual log, by a command that this page does not run:

```bash
pip install "twinfit[log]"
```

Or `pip install twinfit` without matplotlib, or from a folder that holds the package, `pip install /the/folder/that/holds/twinfit`, or put that folder on the path of Python, `export PYTHONPATH=/the/folder/that/holds/twinfit`. Where the Python of your machine is called `python3`, or lives at a path of its own, type that in the place of `python` in every command of this page. Then this prints the version of the package:

```bash
python -c "import twinfit; print(twinfit.__version__)"
```

```
0.1.3rc4
```

## 3. A worked project

This writes a project of nuclear magnetic resonance into a new folder, `my_nmr`:

```bash
python -m twinfit example nmr my_nmr
```

```
written into my_nmr: description.toml, model.py, an empty terms/, and 4 curves in traces/. Next: python -m twinfit loop my_nmr
```

It writes what a researcher writes before any fit. These files are on disk now, and the blocks of this section show them: you type nothing. A block under the words `Save this as` is a file that you write, from section 5 on:

```
my_nmr/
my_nmr/description.toml
my_nmr/model.py
my_nmr/terms/
my_nmr/traces/
my_nmr/traces/temperature_290.csv
my_nmr/traces/temperature_300.csv
my_nmr/traces/temperature_310.csv
my_nmr/traces/temperature_320.csv
```

**The curves**, one file a curve in `traces/`, with a line of names and then two columns, the abscissa and the value. The first lines of one:

`my_nmr/traces/temperature_290.csv`:

```
frequency,intensity
-20,0.052952
-19.9,0.066234
-19.8,0.064911
```

**The description of the data collection**, `description.toml`. It says what was measured, the abscissa and the value with their units, and each setting under which a curve was taken, here the temperature of the sample. Then one table a curve, with its file, its columns and its setting, and the noise of each sample, here one standard deviation for every sample. Then the addenda, all optional. One names the first model, which the loop takes from the package's generic model where a description names none. Here three more are changes of the experiment, which the proposals for the experiment weigh against each other.

A key that the package does not read is refused, so that a misspelled one is not dropped in silence. The order in which the curves were taken goes under `[problem]` as `order`, since a setting taken in the order of the time cannot be told from a drift:

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

**The first model**, `model.py`: the curve as a function of the abscissa `x`. Its parameters are `SHARED` by every curve or `OWN` to each. Here it is one Lorentzian line of one width for all four spectra, each with its own height, place and baseline. `CENTRE` names the own parameter that places a curve, and `centre` returns it, so that the report can tell the core of a curve from its wings. `START_SHARED` and `start_own` give the fit its starts, and `BOUNDS_SHARED` the bounds of the shared parameters, a list of the lower bounds and one of the upper.

`model`, `SHARED`, `OWN` and `start_own` are required, with a start of the shared parameters, `START_SHARED` or a function `start_shared(records)` in its place, `CENTRE` and `centre` come together or not at all, and the rest is optional. `NAME` says the model in a few words, as the reports name it, and `HELD` holds the values that it reads and does not fit, section 9. `record` holds a curve, `record.x` its abscissa and `record.y` its values, and `design` its settings by name, as the description gives them:

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

The loop runs here in the foreground, since each run of this section takes about a minute. A run of your own records may take an hour. Start it in the background with `nohup`, `nohup python -m twinfit loop DIR > DIR/run.log 2>&1 &`, and read its end with `python -m twinfit wait DIR`. Then a closed window or the limit of a tool does not stop it, and the card of the guide for agents runs it in the background too, with its log in `DIR/agent/run.log`.

```bash
python -m twinfit loop my_nmr
```

```
model 0: the first model is fitted, and its twin and its report are read on 30 synthetic sets
best.md: the five parts of each value of the best model, on 30 synthetic sets, and the combined estimate of the families
kept: nothing. Refused: nothing. The best model: model 0.
The report of each model and the output are in my_nmr/loop. Write my_nmr/terms/01.py with the next term, then run the loop again. An own value moves from curve to curve, which no question of the report reads: best.md says which, and how to write it.
The visual log: my_nmr/loop/log.pdf, 15 slides, each phase of each of 1 turn first as one plot, then the equations, the effects of the terms, the curves in a row, the report and the criteria, the pictures in my_nmr/loop/log.
```

Its last lines name the folder of the reports and the visual log, `my_nmr/loop/log.pdf`, a slide for each phase of each model, which section 6 shows.

The loop fitted the first model, `model 0`, and wrote its report to `my_nmr/loop/model_0.md`. The report opens with the words that it uses:

`my_nmr/loop/model_0.md`:

```markdown
## The words of this page

- **The twin**: synthetic records drawn from this fitted model with the noise that the description gives, each fitted again as the records were. What the records show beyond what the twin shows is what the model lacks.
- **The bar of the fit**: the standard deviation of a parameter that the fit gives from the curvature of its chi square, at first order. The first order bound is that bar.
- **A scatter**: the standard deviation of a parameter over the synthetic sets of the twin, which the rule of the twin calls the spread of the estimate. A bias of 3 scatters is a value three of its own spreads from the truth of the twin.
- **The excess** of an answer: its statistic on the records less the median of the same statistic on the synthetic sets, over the span from that median to their 95th percentile. An excess of 1 is the 95th percentile of the twin, and an excess under 0 an answer that the records show less than the median of the twin does.
- **The chance** of an answer: the share of the synthetic sets on which its statistic is as large as on the records. With 30 sets the smallest chance is 0.032, and the answers that share it are ordered by their excess.
- **The level** of a curve, in a polynomial of the level, is its baseline.
- **A search**: the fit of one synthetic set, which ends when its steps no longer lower its chi square. A **Mechanism** is a term written in Python, as section 11 of the tutorial shows.
- **Seen by**: the family of statistics that sees the misfit most. It says how the misfit is read, and not which term to add.
- **The known parameters**: the shared parameters of the model, those of every curve. **The rule of the twin**: the bands within which the fit of this model must give back, on its synthetic sets, the values that made them, its searches ending on nearly all sets, its bias under a share of its spread and its spread near its bar, for its bars to be read.
- **Beyond the twin**: an excess larger than 4.00, the largest excess that a synthetic set of a right model reaches on 95 of 100 of them. Some answer is beyond it on 5 of 100 sets where the model lacks nothing. The line is read on the twin of this model, so it differs a little from one model to the next.
- **The size** of a term: the one parameter that it adds, which `python -m twinfit terms` defines for each. For a term that adds a part of the curve it is the area of that part over that of the curve without it, so that its share of the whole is a / (1 + a), and never a distance.
```

Then the fit. Its chi square is 1923.0 on 1591 degrees of freedom, 1.21 for each, where a right model with a right noise gives 1 with a standard deviation of 0.035. So something is missing. The page says it in its first line:

`my_nmr/loop/model_0.md`:

```markdown
## The fit

- **Chi square**: 1923.0 on 1591 degrees of freedom, 1.21 a degree of freedom, above 1 by more than three of its standard deviations, 0.035 each: the model lacks something, or the noise is stated too small.
- **The twin admits the fit**: the fit of this model gives back, on synthetic sets drawn from it, the values that made them, within the bands of its rule. That says the fit can be read, and not that the model is right.
- **What the twin measured**, each bias in spreads of the estimate on its own sets: searches ended on 1.00 of the sets. width: a bias of +0.13 of its spread, a spread of 0.87 of its bar of the fit, an interval that holds the truth on 0.97 of the sets.
- **The rule of the twin**: searches ended on 0.98 of the sets, a bias under 0.548 of the spread, a spread of 0.606 to 1.39 of the first order bound, an interval that holds the truth on 0.831 to 1 of the sets. With 30 sets each of these is known only within its standard error, so each band is widened to three of them. The bands as stated are a bias under 0.25, a spread of 0.85 to 1.15 and a coverage of 0.92 to 0.98, and from 475 sets on they are held as stated. The verdict above used the wider bands, since the twin read 30 sets.
- **The noise**: The twin draws its synthetic sets at the noise that the description gives, so it takes that noise as right: its coverage and its bias cannot see a noise stated wrong. The residuals against the height of the fitted curve say whether it is. The mean squared whitened residual, the residual over the noise of its sample squared, in 6 bands of the height of the fitted curve above its lowest value, from low to high: 1.08, 0.92, 1.10, 0.81, 1.07, 2.20, where 1 is right. Not every band is 1 within three standard errors: there the noise that the description gives is not that of the records, or the model misses. Where the report finds an answer of the curve beyond the twin, a place, a kind or a setting, read it first as the model that misses, and read this check again once the term is in. Where the answer that leads is one of the noise, a noise that grows with the curve or glitches, the noise comes first, as the section of the proposals says. Where nothing stands beyond the twin and a band still is not 1, the noise is wrong there: a band above 1 says that the stated noise is too small at that height, and the bars of the fit too small by about the square root of the band, and a band under 1 the reverse. `[noise] from = "residuals"` in the description reads it from the second differences of what the first model leaves, which a term that the model lacks, smooth over three samples, does not reach, or give a law of your own. The bands outside: 2.20 +- 0.09, 13.9 of its standard errors from 1.

| parameter | value | bar of the fit |
|---|---|---|
| width | 2.1367 | 0.0049 |
```

The twin admits the fit: it gives back the width that made it. But the check of the noise fails in its highest band, where the line is tallest. That check cannot tell a noise stated wrong from a model that misses. While the report finds an answer beyond the twin, as here, read it as the model that misses, and read the check again once the term is in, as section 6 does. Until then `best.md` marks the coverage and the bias on the twin as not to be read.

Then the first answer of each question. An answer stands beyond the twin past a line. An excess of 1 is the 95th percentile of one answer on the synthetic sets. The report reads many answers at once, so its line is higher: the largest excess of all the answers together that a set of a right model reaches on 95 of 100 sets. Some answer then passes the line on 5 of 100 sets of a right model, and not each answer on 5 of 100:

`my_nmr/loop/model_0.md`:

```markdown
## Where to search

Read on 30 synthetic sets. The first answer beyond the twin of each question, beyond meaning an excess above 4.00 here:

| question | what it asks | first answer beyond the twin | its excess |
|---|---|---|---|
| where | the core of the curves, their flanks, their wings or their asymmetry | the flanks, read along a setting, the flanks along temperature: no place stands beyond the twin alone | 34.2 |
| shape | how the residuals are shaped: a slow trend over the span of each curve, a periodic structure, a part in quadrature with the line, or a structure in one place and at one scale | none beyond the twin. The largest: a periodic structure fixed to the abscissa, in phase from curve to curve | 3.6 |
| what | glitches, a noise that grows with the signal, or a noise correlated from sample to sample | none beyond the twin. The largest: a noise that grows with the signal | 2.2 |
| changes along | which part of the residuals changes from curve to curve along a setting, read on the residuals themselves | the flanks along temperature | 34.2 |
| seen by | the family of statistics that sees the misfit most | each curve predicted by the fit of the others | 62.9 |
| noise | the shape of the noise: skewed to one side, or with heavy tails | none beyond the twin. The largest: a noise skewed to one side | 0.6 |
```

Read it row by row. No place stands beyond the twin on its own, no shape and no law of the noise. But what the fit leaves changes along the temperature, in the core of the lines and in their flanks, so the row of where names the flanks, the first of the two, read along that setting. The package says no more than that: it names no term. The term is yours, from the physics of the curves.

Where two places stand beyond the twin and the records scatter their excesses alike, the row names both, tied. The report reads the scatter of each excess on sets drawn from the fit plus what the fit leaves, and the records do not tell which of the two leads.

In a resonance, what changes the core and the flanks of a line together with the temperature of the sample is its width, which the motion of the molecules sets. The proposals for the model are the directions of the report, each measured with the shape that the residuals hold there, fitted beside the model. On their twin, the shape that the residuals hold in a place, freed of the fit's own directions, costs the width nothing at first order, since the residuals hold none of what the fit took into its own parameters, while a slow trend or a periodic structure, not freed, reads a cost: what a term of that place and that setting would cost lies in that part, which only a term fitted beside the model shows.

`my_nmr/loop/model_0.md`:

```markdown
## Three proposals for the model

Each is a direction read on the residuals, for an answer of where, of shape or of what changes along a setting: those beyond the twin first, then by the gain in chi square that the shape that the residuals hold there brings when it is fitted beside the model. The package names no term: which term of physical origin has that place, shape and setting is the researcher's to say, in a term file. On its twin, synthetic curves drawn with the term at its fitted size are fitted without it and with it: the bias without it is what the term would cost the known parameters if it were real, and says nothing of whether it is. The shape that the residuals hold in a place is freed of the fit's own directions, so it costs them nothing at first order, and its line says so: a shape built apart from the residuals, a slow trend, a periodic structure, a part in quadrature with the line or a wavelet, is not freed, and its cost is read as any other.

In the names of the proposals a width is the full width of the fitted curves at half their height above their median value, here 2.1 in the unit of the abscissa: where the wings of a curve fill its window, it is narrower than the width at half its peak, and it is not a parameter of the model. Above the centre is toward larger abscissa.

- the shape that the residuals hold in the flanks, changing along temperature, for the residuals in the flanks change along temperature, beyond the twin at an excess of 34.2: a gain in chi square of 4.0
  - read on the residuals, this shape costs width nothing at first order: on its twin width is biased by -0.4 scatters, -0.00133 in its own unit, without it and by -0.4 scatters, -0.00133 in its own unit, with it, the twin's own bias, since the residuals hold none of what the fit took into its own parameters. What a term in this place would cost lies in that part, which no residual shows: fit a term of the physics here and read how far width moves
- the shape that the residuals hold in the core, changing along temperature, for the residuals in the core change along temperature, beyond the twin at an excess of 25.2: a gain in chi square of 2.4
  - read on the residuals, this shape costs width nothing at first order: on its twin width is biased by -0.4 scatters, -0.00133 in its own unit, without it and by -0.4 scatters, -0.00133 in its own unit, with it, the twin's own bias, since the residuals hold none of what the fit took into its own parameters. What a term in this place would cost lies in that part, which no residual shows: fit a term of the physics here and read how far width moves
- the shape that the residuals hold in the wings, changing along temperature, for the residuals in the wings change along temperature, beyond the twin at an excess of 6.8: a gain in chi square of 1.9
  - read on the residuals, this shape costs width nothing at first order: on its twin width is biased by -0.4 scatters, -0.00133 in its own unit, without it and by -0.4 scatters, -0.00133 in its own unit, with it, the twin's own bias, since the residuals hold none of what the fit took into its own parameters. What a term in this place would cost lies in that part, which no residual shows: fit a term of the physics here and read how far width moves

## The known parameters before the term is named

Where an answer stands beyond the twin, each shared parameter is read five ways while the missing term has no name: the fit as it is, the combination of the families that is free of the leading class of the slow trend, the one basis fixed before the residuals are read, the model with the generic term of the largest gain, width changes along temperature, fitted and at first order, and the combination free of that term. The generic terms are those that the model and the design alone give, each shared parameter along each covariate, and no list of terms is read. Where the five agree, the missing term does not move the parameter in the directions tried. Where they part, the spread says how far it may. None of them sees a part of the missing term that lies in the fit's own directions and in no basis tried.

| parameter | the fit | free of the class | with the generic term | that term at first order | free of that term |
|---|---|---|---|---|---|
| width | 2.1367 ± 0.0049 | none | 2.1182 ± 0.0050 | 2.1158 ± 0.0051 | 3.455 ± 0.030 |
```

A width in the page of the proposals is a unit of their own: the full width of the fitted curves at half their height above their median value, as the page says above its list. It is not the parameter `width` of `model.py`, which here is a full width too, and in section 10 a standard deviation.

The first direction is the one whose shape lowers the chi square most, and a term that takes part of another's misfit can lead. On records that lack a slope of the baseline, a second component written first can be kept, and fall to nothing once the slope is written. A term kept first is not looked at again when a later one arrives. So `best.md` names a kept term that has fallen within two bars of none, and the model without it is tried in a copy, as section 7 says.

The gain of a proposal ranks the directions and is no forecast of the fall that a term of yours brings: the shape fitted at first order holds a little of the misfit, and a term that holds it whole lowers the chi square much further. Section 7 sets the two side by side.

The row seen by names the family of statistics that sees the misfit most. A family turns the residuals into a few numbers. Its moments are of the abscissa, and of the distance from the center to whole, fractional and negative powers, and the same under logarithmic, exponential and Gaussian kernels. Others are Hermite functions about the center, windows along the abscissa and along the value, and the crossings of levels of the height. Others again are the moments of the derivative and of the running integral, and the moments of the square, the root, the logarithm, the reciprocal and the exponential of the value, which weigh the top of a line and its foot otherwise than the value does. Beside the families stand the eight chi functions, each averaged over the residuals, and the plane of position and residual value, which weigh large residuals, and where they fall, otherwise than the square does.

### In depth, optional on a first reading: what each family reads

Section 5 goes on without this part, which shows what each family of statistics reads, and how each law of the noise lights them.

A term that the model lacks moves each family in its own way, and the pattern across the families is what the report reads. `python -m twinfit matrix` prints those readings for the model that the term files of a folder give, with no twin, in a few seconds:

```bash
python -m twinfit matrix my_nmr
```

```
The families of statistics on the model of my_nmr, its first model, each term file kept by the Bayesian criterion alone or by its file's keep = True. No twin is run: a level here is that of a chi square under a right model at first order, which the families of many statistics do not reach, and the report of the loop gives the levels of the twin.

## What each family reads

Each family turns the residuals into a few statistics. Its reading is their chi square in the covariance that they have under a right model, with as many degrees of freedom as the family has directions that the model does not hold. The smallest levels first.

| family | chi square | degrees of freedom | level |
|---|---|---|---|
| moments: negative orders of the distance | 336.0 | 16 | 8.6e-62 |
| the exponential of the value: power kernels 0 to 8 | 379.1 | 36 | 7.7e-59 |
| the square of the value: power kernels 0 to 8 | 378.2 | 36 | 1.2e-58 |
| moments: inverse and logarithmic kernels | 352.4 | 27 | 2.2e-58 |
| ratios of moments of the distance: orders 0.5 to 4 over order 0 | 357.8 | 32 | 1.1e-56 |
| moments: fractional powers of the distance | 330.6 | 24 | 1.1e-55 |
| moments: exponential and Gaussian kernels | 355.3 | 36 | 3.9e-54 |
| the two sides added: pairs of intervals, the centre held | 333.2 | 28 | 5.8e-54 |
| Hermite functions about the centre: orders 3 to 8 | 314.6 | 28 | 3e-50 |
| the reciprocal of the value: power kernels 0 to 8 | 318.5 | 36 | 5.8e-47 |
| the square root of the value: power kernels 0 to 8 | 314.9 | 36 | 2.9e-46 |
| windows along the value: bands between levels of the height, each side | 415.9 | 84 | 1.9e-45 |
| the logarithm of the value: power kernels 0 to 8 | 308.8 | 36 | 4.4e-45 |
| windows along the abscissa: pairs of intervals, the two sides added and subtracted | 345.4 | 52 | 6.3e-45 |
| moments: powers 0 to 8 | 244.2 | 32 | 1.6e-34 |
| moments of the running integral: powers 0 to 8 | 244.9 | 36 | 6.7e-33 |
| moments of the derivative: powers 0 to 8 | 232.2 | 36 | 1.6e-30 |
| the abscissa against the value: the crossings of nine levels, each side | 296.2 | 72 | 5.7e-29 |
| moments: powers 0 to 4 | 131.6 | 16 | 3.2e-20 |
| moments: cosines and sines | 99.5 | 32 | 7.7e-09 |
| every sample: the fit | 1923.0 | 1591 | 1.6e-08 |
| the two sides subtracted: pairs of intervals | 20.0 | 28 | 0.86 |

## The shape of the noise

The whitened residuals of every curve as they stand, their skewness over its standard error under a normal law and their excess kurtosis over its own, each a standard normal deviate at first order where the noise is normal. A misfit of the mean leaves tails of its own, and the report reads them once the shape that leads is taken out: a noise skewed to one side 1.16, a noise with heavy tails 2.12.

## Which statistic sees which law of the noise

Each law of the noise is drawn 400 times on these records, whitened and of unit variance on every sample, so that the noise that the description states stays right on average and only its law is wrong, and taken through the fit as the residuals take it. The line of each statistic is where it stands on 95 of 100 draws of a normal noise of the same design, so each cell is a chance by simulation at the level 0.05, and owes nothing to first order. After the families come the five readings of the noise that the report gives, the records together, and two opponents that practitioners read, the runs test and Durbin and Watson, which the report does not read. The last row reads every other statistic at one line, as the report reads its answers: the chance that the largest excess of them all passes where a normal noise stands on 95 of 100 draws. The laws:

- normal: a normal noise as the description states it, drawn again: each of its cells is the level that the line holds
- correlated: a noise correlated from one sample to the next by 0.5 in the order of the abscissa
- heavy tails: Student's t of 3 degrees
- skewed: a chi square of 4 degrees less its mean
- grows with the signal: a variance in proportion to the height above the foot of the fitted curve plus a tenth of its span
- glitches: one sample in 500 at six deviations

| statistic | normal | correlated | heavy tails | skewed | grows with the signal | glitches |
|---|---|---|---|---|---|---|
| moments: powers 0 to 4 | 0.06 | 0.93 | 0.05 | 0.04 | 0.00 | 0.06 |
| moments: powers 0 to 8 | 0.04 | 0.99 | 0.07 | 0.06 | 0.00 | 0.08 |
| moments: fractional powers of the distance | 0.04 | 0.95 | 0.06 | 0.04 | 0.49 | 0.04 |
| moments: inverse and logarithmic kernels | 0.04 | 0.97 | 0.07 | 0.03 | 0.45 | 0.04 |
| moments: negative orders of the distance | 0.04 | 0.88 | 0.07 | 0.06 | 0.81 | 0.05 |
| moments: exponential and Gaussian kernels | 0.06 | 0.99 | 0.10 | 0.05 | 0.74 | 0.06 |
| Hermite functions about the centre: orders 3 to 8 | 0.06 | 0.99 | 0.07 | 0.05 | 0.76 | 0.05 |
| moments: cosines and sines | 0.05 | 0.98 | 0.08 | 0.06 | 0.01 | 0.09 |
| windows along the abscissa: pairs of intervals, the two sides added and subtracted | 0.06 | 0.99 | 0.09 | 0.08 | 0.96 | 0.07 |
| the two sides added: pairs of intervals, the centre held | 0.07 | 0.95 | 0.09 | 0.07 | 0.75 | 0.06 |
| the two sides subtracted: pairs of intervals | 0.05 | 0.95 | 0.06 | 0.10 | 0.83 | 0.08 |
| windows along the value: bands between levels of the height, each side | 0.04 | 0.97 | 0.10 | 0.10 | 1.00 | 0.07 |
| the abscissa against the value: the crossings of nine levels, each side | 0.04 | 0.80 | 0.12 | 0.12 | 1.00 | 0.09 |
| moments of the derivative: powers 0 to 8 | 0.05 | 1.00 | 0.08 | 0.07 | 0.00 | 0.08 |
| moments of the running integral: powers 0 to 8 | 0.03 | 1.00 | 0.07 | 0.04 | 0.00 | 0.08 |
| ratios of moments of the distance: orders 0.5 to 4 over order 0 | 0.04 | 0.98 | 0.07 | 0.06 | 0.74 | 0.04 |
| the logarithm of the value: power kernels 0 to 8 | 0.03 | 0.94 | 0.09 | 0.06 | 0.78 | 0.05 |
| the reciprocal of the value: power kernels 0 to 8 | 0.03 | 0.90 | 0.09 | 0.06 | 0.68 | 0.04 |
| the square of the value: power kernels 0 to 8 | 0.05 | 0.97 | 0.07 | 0.06 | 0.78 | 0.04 |
| the square root of the value: power kernels 0 to 8 | 0.03 | 0.94 | 0.07 | 0.05 | 0.76 | 0.04 |
| the exponential of the value: power kernels 0 to 8 | 0.05 | 0.97 | 0.08 | 0.08 | 0.77 | 0.06 |
| every sample: the fit | 0.04 | 0.04 | 0.20 | 0.11 | 0.03 | 0.12 |
| glitches | 0.05 | 0.00 | 0.99 | 0.96 | 0.99 | 0.97 |
| a noise skewed to one side | 0.05 | 0.07 | 0.85 | 1.00 | 0.39 | 0.46 |
| a noise with heavy tails | 0.04 | 0.06 | 1.00 | 1.00 | 1.00 | 0.95 |
| a noise that grows with the signal | 0.06 | 0.00 | 0.42 | 0.18 | 1.00 | 0.16 |
| a noise correlated from sample to sample | 0.04 | 1.00 | 0.03 | 0.04 | 0.15 | 0.05 |
| the records together: every record added in bins of the abscissa | 0.05 | 0.94 | 0.07 | 0.06 | 0.01 | 0.07 |
| the opponent: the runs test of the signs, record by record | 0.07 | 1.00 | 0.04 | 0.04 | 0.07 | 0.04 |
| the opponent: Durbin and Watson, the lag one correlation, record by record | 0.05 | 1.00 | 0.03 | 0.04 | 0.05 | 0.04 |
| every statistic at one line, as the report reads them | 0.08 | 1.00 | 1.00 | 1.00 | 1.00 | 0.96 |
```

Its first table is what each family reads: a chi square with its degrees of freedom, and its level at first order. A family of many correlated statistics does not hold that first order level, so the report takes its levels from the twin. The chi square of every sample reads 1923 on 1591 degrees of freedom, 1.21 for each, while the exponential of the value: power kernels 0 to 8 reads 379 on 36: a misfit of one shape stands in a few directions, which a family of few directions gathers and the sum over every sample spreads over all of them.

The table of the noise asks the same of its law. Each law is drawn on these curves at the variance that the description states, so that only its law is wrong, and the line of each statistic is set where a normal noise of the same design leaves it on 95 draws of 100: the column of a normal noise drawn again reads from 0.03 to 0.08, as the table prints it, the level.

A noise correlated from one sample to the next lights every family, in 80 to 100 cases of 100, and the plain fit in 4. So on curves whose noise is correlated, a family past its line is no sign of a missing term, and the reading of the neighbors, which reports that noise in 100, says which it is.

Heavy tails and a skewed noise keep the covariance of the noise, so no family, linear in the residuals, reports them in more than 12 cases of 100, and the readings of the shape of the noise do. A noise that grows with the signal changes the covariance along the height, so the windows along the value report it, and the moments of the powers of the abscissa do not.

## 5. A term

A term is a file in `terms/`, `01.py`, then `02.py`, and so on, which you write. This prints the terms that a term file can name, each but `exponential_width()` described with the one parameter that it adds, and the form of a term file:

```bash
python -m twinfit terms
```

```
The terms of the package, each with the one parameter a that it adds:

    widened(size)            the curve convolved with a two sided exponential kernel of scale size(design) * sqrt(a), size a
                             number or a function of the settings of a curve, as lambda design: 0.5
    displaced(reach)         a mixture of copies of the curve displaced by reach(design) * sqrt(a) * v, v from 0 to 1 of
                             density 2 v, reach a number or a function of the settings, as for widened. About its mean,
                             its long side points toward smaller abscissa where reach is positive, toward larger where
                             it is negative: the opposite of second, whose positive distance is after the curve
    along(parameter, key, power=1)
                             a shared parameter, named as SHARED names it, that changes along a setting of the description,
                             named as its [covariates] names it: p (1 + a z), z the setting less its mean over the
                             records, over its standard deviation over them, taken with their number, not one less.
                             power=2 takes z squared less its mean
    gaussian_width(along=None)
                             the curve convolved with a Gaussian of standard deviation a, in the unit of the abscissa, on
                             the steps of the records, which need equal steps: a width that the line lacks, as the noise
                             of a source or the response of an instrument. along, the name of a setting or a function of the settings of a curve, scales a in it
    lorentzian_width(along=None)
                             the same with a Lorentzian of full width a at half its height: a width of a lifetime, and
                             along="concentration" makes it in proportion to that setting over its mean over the
                             records, a the width at that mean, and names it
    pedestal(stretch)        a copy of the curve `stretch` times wider under it, of relative area a, 0 or more unless
                             `lower` says otherwise, stretch a number
    bent()                   the abscissa is not linear: x + a (x - middle)^2 / half its span
    second(distance)         a second copy of the curve, `distance` away along the abscissa, in its unit, after the curve
                             where it is positive and before it where negative, of relative area a, 0 or more unless
                             `lower` says otherwise
    tail(stretch, side)      a copy of the curve `stretch` times wider on one side of its centre alone, of relative area a,
                             side "upper" or "lower". A wider copy, not the exponential tailing of a column or of a
                             detector, which a Mechanism of your own writes
    level(power)             a power of the abscissa under every curve, of one size for all: a ((x - middle) / half)^power,
                             the middle and the half span those of the abscissas of all the records together, power 1, 2
                             or 3 for a slope, a curvature or a cubic term of the baseline
    periodic(period, part)   a cosine or a sine of the abscissa of one period, in its unit, part "cosine" by default or
                             "sine": a cos or a sin(2 pi (x - centre) / period), about the centre of each curve. A ripple
                             fixed to the abscissa, the same on every curve whatever its centre, as the fringe of an
                             etalon, is a Mechanism of your own with x in the place of x - centre:
                             section 11 of the tutorial writes one

Glitches and a change of the noise are no term of the mean: they are read in the weights of the fit, the noise of each
sample that the description gives. Each term here is made of copies of the curve, wider, displaced or on one side: a shape
that the physics names otherwise, as the exponential tail of a detector or of a column, is a Mechanism of your own. The report
gives the place, the shape and the setting of a misfit, and never the name of a term: which term of physical origin fills it
is the researcher's to say. A Mechanism is Mechanism(name, family, parameter, start, wrap,
lower=-inf, upper=inf): family is a word that the package keeps with the term and reads in no result, lower and upper bound
the new parameter, and wrap(model, problem, n) returns any curve, the line replaced as well as a term added to it.

The name that a report and best.md give a term, and the call of a term file that writes it:

    a kernel with exponential wings               widened(0.5)
    a displacement of linear density              displaced(0.5)
    a pedestal 8 times wider than the curve       pedestal(8)
    a quadratic map of the abscissa               bent()
    a second component 1.5 away                   second(1.5)
    a tail on the upper side, 2 times wider       tail(2, 'upper')
    a power 1 of the abscissa under every curve   level(1)
    a cosine of period 3                          periodic(3.0)
    a sine of period 3                            periodic(3.0, 'sine')
    width changes along concentration             along('width', 'concentration')
    a Gaussian width beside the line              gaussian_width()
    a Lorentzian width beside the line            lorentzian_width()
    an exponential width beside the line          exponential_width()

A term file is DIR/terms/01.py, then 02.py, and so on, one term each, read by the loop in that order:

    from twinfit import mechanisms
    term = mechanisms.along("width", "concentration")
    why = "the report: the core changes along concentration, and a width grows with the concentration of the sample"

`term` is one of the terms above, or a mechanism of your own, and `why` says why you add it, in six words or more. A
file that is kept may be written again later: the loop fits every file again at each run. The loop fits the
model with the term and keeps it when three things hold: its fit ends inside its bounds, its Bayesian criterion falls
by more than 2, and the twin admits the larger model. The Bayesian criterion is the chi square plus the number of
fitted parameters times the logarithm of the number of samples, so a fall of more than 2 is a gain that one more
parameter does not buy alone. Otherwise it refuses the term, with the reason, but where the file says `keep = True`,
for a term of solid physical ground, which it keeps whatever the criterion.

A term that the physics asks for and the records cannot show is held at its value of theory and not fitted: give the
term its theory, as mechanisms.pedestal(4.0, theory=(0.05, 0.01)), its size and the bar of that size, and write a line
hold = True in its file. The loop keeps it at that size whatever the criterion says, and carries its bar into the bar of
theory of every other value.
```

The report points to the core and the flanks changing along the temperature, and the physics of a resonance names the width. `along` makes a shared parameter change along a setting, here the width along the temperature. A term file adds `why`, the reason in six words or more:

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

```
model 0: the first model is fitted, and its twin and its report are read on 30 synthetic sets
01.py: width changes along temperature is fitted, and read on 30 synthetic sets: each term costs more than the one before it, since the model grows
best.md: the five parts of each value of the best model, on 30 synthetic sets, and the combined estimate of the families
kept: 01.py, width changes along temperature. Refused: nothing. The best model: model 1.
The report of each model and the output are in my_nmr/loop. No answer of where, shape, what or changes along in the last report, nor a term of your list, stands beyond its twin: the records ask for no more, and a further term rests on theory alone. An own value moves from curve to curve, which no question of the report reads: best.md says which, and how to write it.
The visual log: my_nmr/loop/log.pdf, 29 slides, each phase of each of 2 turns first as one plot, then the equations, the effects of the terms, the curves in a row, the report and the criteria, the pictures in my_nmr/loop/log.
```

The run's output uses two words: your list, `list_of_terms.py`, the terms that you suspect, which `init` writes empty and section 9 shows, and which a project of `example` does not hold, so that its list is empty, and matched, the question of the report that reads the residuals against each term of that list. With the list empty, no answer of matched is read. The package holds a list of its own, read only where you run the loop with `--dictionary`: it names terms, and the search for a term that you do not suspect runs without it.

The loop kept the term. Its output is `my_nmr/loop/best.md`, the best model with the five parts of the uncertainty of each parameter: the bar of the fit, the coverage and the bias on the twin, the bar of theory and the spread over the models. Their words stand under its heading `The words of the table`. Where shared parameters are correlated beyond 0.9, its section `What the records tell apart` gives the combinations of them that the records hold, each with the bar to quote, read on the twin, and beyond 0.999 it quotes none, and its section `The values in windows` gives each value read in windows about the centre and off the peak:

`my_nmr/loop/best.md`:

```markdown
# The best model: model 1

Its terms: one Lorentzian line of one width, width changes along temperature.

- width changes along temperature, in 01.py, width along temperature: the report: the core and the flanks change along temperature, and the width of a resonance changes with the temperature of the sample

The description does not say in which order the curves were taken, `order` under [problem]: a setting taken in the order of the time cannot be told from a drift.

The best model is the last model of the loop, the first model where it kept no term: the one of the smallest Bayesian criterion, the chi square plus the number of fitted parameters times the logarithm of the number of samples, among the models that it compares, from the first model on, or from the last term held at its theory or kept by its file's keep = True, since a model that lacks such a term is not its rival. The models whose criterion is within 2 of it are those that the records cannot tell apart from it.

| model | its term | chi square | fitted parameters | Bayesian criterion | less that of the best |
|---|---|---|---|---|---|
| model 0 | the first model | 1922.97 | 13 | 2018.91 | 328.22 |
| model 1 | width changes along temperature | 1587.38 | 14 | 1690.70 | 0.00 |

## The words of the table

- **The value**: the fitted value of the parameter in the best model. A parameter that a term adds is defined by its term: `python -m twinfit terms` lists them.
- **The bar of the fit**: the standard deviation of the parameter that the fit gives at first order, from the curvature of its chi square.
- **The coverage on the twin**: the share of synthetic sets of the best model on which the interval of 1.96 bars of the fit holds the value that made them. 0.95 is right, and with the 30 sets read here a coverage is known within 0.04, one standard deviation. Well under 0.95, the bar of the fit is too small, and well over it, too large. The sets are drawn from seeds read from the records, so a run gives the same numbers again on the same records, and records of another project draw sets of their own: a coverage is one draw of its uncertainty, which more sets narrow.
- **The bias on the twin, in bars**: the mean of the values fitted on those sets less the value that made them, in bars of the fit. It is given beside the value and not taken out of it: `best.json` gives the value less it too.
- **The bar of theory**: where the model holds values of theory, half the move of the parameter when each held value is set one bar of theory below and above its centre. None held: the description holds no value of theory, which an ansatz in its addenda would, `[addenda.ansatz.NAME]` with its value and its uncertainty.
- **The spread over the models**: half the range of the parameter over the best model and the models that the loop fitted and that the records cannot tell apart from it. A term that was never written into a term file is in no model, so two terms that fit alike and were not both tried leave it at 0.
- **Not told apart**: a parameter correlated with another beyond 0.999 gets no bar of the fit, since the records fix a combination of the two and not each.

## The parameters

The twin draws its synthetic sets at the noise that the description gives, so it takes that noise as right: its coverage and its bias cannot see a noise stated wrong. The residuals against the height of the fitted curve say whether it is. The mean squared whitened residual, the residual over the noise of its sample squared, in 6 bands of the height of the fitted curve above its lowest value, from low to high: 1.04, 0.91, 1.02, 0.82, 1.03, 1.11, where 1 is right. Each is 1 within three standard errors.

| parameter | value | bar of the fit | coverage on the twin | bias on the twin, in bars | bar of theory | spread over the models |
|---|---|---|---|---|---|---|
| width | 2.1182 | 0.0050 | 0.97 | -0.01 | none held | 0 |
| width along temperature | 0.0433 | 0.0024 | 0.97 | 0.32 | none held | 0 |

## The combined estimate of the families

Each parameter estimated by the 5 families of statistics of the report that give an estimate, the fit among them, weighted by how they vary together on synthetic sets of the best model, with its standard deviation read on 60 further sets. It is a check of the model and not a better value: under a right model and normal errors no combination scatters less than the fit, and weights read on synthetic sets cost a little. Its interval of 1.96 standard deviations should hold the truth on 0.95 of the further sets, within 0.06. Its gap from the fit is read in the scatter of their difference on the same sets: beyond 3, the families and the fit read the records apart, which says that the model lacks a term, as the test of two estimators of Hausman (1978) reads it.

| parameter | combined estimate | its standard deviation over the fit's | its interval holds the truth on the twin | gap from the fit, in the scatter of their difference |
|---|---|---|---|---|
| width | 2.1152 ± 0.0058 | 1.23 | 0.90 of 60 | -1.0 |
| width along temperature | 0.0420 ± 0.0023 | 1.04 | 0.95 of 60 | -1.2 |

## The values in windows

The fit read again in windows of growing half width about the centre of each curve, in widths of the curve, and off the peak, every sample beyond a half width, each estimate set against the fit's by its contrast, read at first order. On a right model every window agrees with the fit within the scatter of their difference, and a term that the model lacks moves the windows where it lives first: the core, the wings, or both. A level under 0.05 is a window that disagrees.

| window | samples | width | width along temperature | level of its contrast with the fit |
|---|---|---|---|---|
| within 0.5 | 80 | 2.138 ± 0.058 | 0.039 ± 0.027 | 0.933 |
| within 1 | 160 | 2.112 ± 0.014 | 0.0375 ± 0.0068 | 0.584 |
| within 2 | 320 | 2.1171 ± 0.0076 | 0.0391 ± 0.0036 | 0.286 |
| within 4 | 640 | 2.1151 ± 0.0057 | 0.0416 ± 0.0027 | 0.203 |
| within 8 | 1280 | 2.1167 ± 0.0051 | 0.0432 ± 0.0024 | 0.352 |
| beyond 0.25 | 1564 | 2.1170 ± 0.0089 | 0.0464 ± 0.0042 | 0.671 |
| beyond 0.5 | 1524 | 2.137 ± 0.022 | 0.056 ± 0.010 | 0.171 |
| beyond 1 | 1444 | 2.18 ± 0.12 | 0.075 ± 0.055 | 0.585 |

Read on 30 synthetic sets, the coverage and the bias are provisional: `--sets 100` reads each coverage within 0.02.

## To quote

Each parameter as a value with the bar of the fit, then, as separate terms, the bar of theory and the half range over the models where they are not zero. Where one bar is asked, add them in squares, since each comes from another source. The bias on the twin is not taken out while it is under one bar. Where the noise, the coverage or the bias fails, the line says so.

- width: 2.1182 ± 0.0050 of the fit, at the mean of temperature over the curves, 305
- width along temperature: 0.0433 ± 0.0024 of the fit, a share of width for each deviation of temperature over the curves, 11.18: the line below gives it per unit
- height, place, baseline, of each curve: Each is in the section of the parameters of each curve with its bar of the fit, which counts what the shared parameters leave uncertain. The twin reads no coverage and no bias of them: quote each with its bar of the fit alone, and with the settings of its curve
- width grows by 0.00820 ± 0.00044 per unit of temperature: the size of its term times width over 11.1803, the deviation of temperature over the curves, with the bar that the covariance of the two gives. At each value of temperature: 290, 1.9952 ± 0.0089, 300, 2.0772 ± 0.0057, 310, 2.1592 ± 0.0052, 320, 2.2413 ± 0.0077

## Notes

- The bias on the twin here is in bars of the fit, on sets of its own: the report of each model reads it in spreads of the estimate, on other sets, so the two differ.
- The coverage on the twin here is read on sets of its own too, and a coverage counted on 30 sets scatters by about 0.04 from one draw of the sets to another: the coverage of a parameter in the report of its model can differ from this one by two or three times that, with no fault in either.
- The spread over the models reads the models of this folder alone. A form of a term tried in a copy, as a law of the physics in `model.py` in the place of a term of the package, is not in it: compare such forms by their Bayesian criterion and quote the range of their values beside the five parts.

## The parameters of each curve

Each with its bar of the fit, which counts what the shared parameters leave uncertain: the twin reads the shared ones alone. An own value that moves with a setting is a term of its own: read them down each column.

- height moves with temperature: a line through its values gives 0.00994 ± 0.00014 per unit of temperature, 69 bars from none, a chi square of 2.4 on 2 about the line. While height is own to each curve, no question of the report reads it: write it in `model.py` as shared parameters that read `design["temperature"]`, as section 9 of the tutorial shows, and the report reads what is left. Across 4 curves a line can also hold by chance, so ask whether the physics gives it.
- place moves with temperature: a line through its values gives 0.02983 ± 0.00014 per unit of temperature, 207 bars from none, a chi square of 1.1 on 2 about the line. While place is own to each curve, no question of the report reads it: write it in `model.py` as shared parameters that read `design["temperature"]`, as section 9 of the tutorial shows, and the report reads what is left. Across 4 curves a line can also hold by chance, so ask whether the physics gives it.

| curve | temperature | height | place | baseline |
|---|---|---|---|---|
| traces/temperature_290.csv | 290 | 0.9985 ± 0.0033 | -0.4462 ± 0.0036 | 0.05075 ± 0.00057 |
| traces/temperature_300.csv | 300 | 1.1033 ± 0.0030 | -0.1443 ± 0.0033 | 0.04963 ± 0.00056 |
| traces/temperature_310.csv | 310 | 1.1969 ± 0.0029 | 0.1496 ± 0.0031 | 0.05041 ± 0.00056 |
| traces/temperature_320.csv | 320 | 1.2991 ± 0.0032 | 0.4504 ± 0.0029 | 0.04973 ± 0.00059 |

Tried and refused: none.

The loop stops here: no answer of where, what or along stands beyond the twin, so the records ask for no more term. A term that the physics asks for is held at its value of theory, in terms/02.py, with hold = True.
```

The curves were drawn from a known truth, so the output can be checked. At the mean temperature, 305 K, the width was 2.12 Hz, and the fit gives 2.1182 ± 0.0050, 0.4 bars away. Its growth along the temperature, in the units of the term, was 0.0422, and the fit gives 0.0433 ± 0.0024, 0.5 bars away.

The heights and the places of the worked spectra were drawn to rise with the temperature, each spectrum its own, as an instrument whose gain and lock drift with it. `best.md` names them as moving with the temperature, with their slopes. They stay own here, since the question of this example is the width. Section 9 says how such a value is written into the model where it is the question.

A proposal of the last report that stands within the twin is one that the records do not ask for. Where its line says what it would cost if it were real, a direction defined apart from the residuals, it names a systematic to say beside the value where the physics allows such a term, and adds nothing to the bar. The shape that the residuals hold in a place, freed of the fit's own directions, costs nothing at first order, and its line says so: a slow trend or a periodic structure is not freed, and its line reads a cost. The worked descriptions do not say in which order their curves were taken, so `best.md` says that it is not known.

In hertz per kelvin the growth is the term's size times the width over the standard deviation of the temperatures, 11.18 K, taken over the four records with their number: 0.0082 Hz per kelvin, where 0.008 drew the curves.

Each run moves the `loop/` of the run before it to `loop.before/` while it works, puts it back where the run stops on a fault, and removes it where the run ends: so `loop/` is always the last whole run.

`loop/` holds, for each model K, its report `model_K.md`, its page for an agent `brief_K.md` and its numbers `model_K.json`. It holds its residuals `residuals_K.csv` and what each of its terms does, `effects_K.npz`. Then come `equations.json`, `best.md`, `best.json` and `STATUS.json`, which says what the run is doing. `notebook/` records every step of the run, and the visual log is `log.pdf`, its slides in `log/` and their list in `log.md`. A run that is refused writes `REFUSED.md` there, beside the files of the run before it.

`loop.before/` is removed where a run ends, so copy `loop/` where a run is to be kept. `loop.reuse/` keeps the twin, the report and the proposals of each model and the twin of each term. A run reads back what a run before it computed from the same files, records and fit. So a run after one more term file computes only what that file adds, and a run stopped at its budget goes on where it stopped.

Which number to quote. Of the five parts, the bar of the fit, the bar of theory and the half range over the models are the uncertainty of a value, and the coverage and the bias on the twin are checks of the bar of the fit. Where the twin admits the fit and the last report asks for nothing more, as section 7 says, the value with its bar may be quoted when three things hold. The noise passes its check, as the line above the table says here. The coverage on the twin lies in the band of the rule of the twin at the number of sets read, the band that admits the fit. That band is three standard deviations of a coverage about 0.95 with few sets, and 0.92 to 0.98 from 475 sets on. So one rule says both whether the fit is admitted and whether its bar is quoted. And the bias on the twin is under one bar. `best.md` writes no line that says a value may be quoted. It writes "Not yet to be quoted" beside a value that may not be, with the reason, so the brief's words "where best.md says that it may be quoted" name a value with no such mark. The report of a model that its twin admits says "The twin admits the fit".

With 30 synthetic sets, the default, a coverage is known within 0.04, one standard deviation, so the coverage of the width, 0.97, is 0.4 of them from 0.95. The report of model 1 reads the same coverage on other sets, so its number differs: neither settles it on 30 sets, whose band is wide. Read it on more sets, and quote the run on 500 sets where the time allows, whose band is the one stated.

The section "To quote" of `best.md` writes each parameter that way, with the bar of theory and the half range over the models as separate terms where they are not zero, and says when a part forbids the quote:

```bash
python -m twinfit loop my_nmr --sets 100
```

```
model 0: the first model is fitted, and its twin and its report are read on 100 synthetic sets
01.py: width changes along temperature is fitted, and read on 100 synthetic sets: each term costs more than the one before it, since the model grows
best.md: the five parts of each value of the best model, on 100 synthetic sets, and the combined estimate of the families
kept: 01.py, width changes along temperature. Refused: nothing. The best model: model 1.
The report of each model and the output are in my_nmr/loop. No answer of where, shape, what or changes along in the last report, nor a term of your list, stands beyond its twin: the records ask for no more, and a further term rests on theory alone. An own value moves from curve to curve, which no question of the report reads: best.md says which, and how to write it.
The visual log: my_nmr/loop/log.pdf, 29 slides, each phase of each of 2 turns first as one plot, then the equations, the effects of the terms, the curves in a row, the report and the criteria, the pictures in my_nmr/loop/log.
```

`my_nmr/loop/best.md`:

```markdown
| parameter | value | bar of the fit | coverage on the twin | bias on the twin, in bars | bar of theory | spread over the models |
|---|---|---|---|---|---|---|
| width | 2.1182 | 0.0050 | 0.93 | 0.02 | none held | 0 |
| width along temperature | 0.0433 | 0.0024 | 0.94 | 0.15 | none held | 0 |

## The combined estimate of the families

Each parameter estimated by the 5 families of statistics of the report that give an estimate, the fit among them, weighted by how they vary together on synthetic sets of the best model, with its standard deviation read on 60 further sets. It is a check of the model and not a better value: under a right model and normal errors no combination scatters less than the fit, and weights read on synthetic sets cost a little. Its interval of 1.96 standard deviations should hold the truth on 0.95 of the further sets, within 0.06. Its gap from the fit is read in the scatter of their difference on the same sets: beyond 3, the families and the fit read the records apart, which says that the model lacks a term, as the test of two estimators of Hausman (1978) reads it.

| parameter | combined estimate | its standard deviation over the fit's | its interval holds the truth on the twin | gap from the fit, in the scatter of their difference |
|---|---|---|---|---|
| width | 2.1152 ± 0.0058 | 1.23 | 0.90 of 60 | -1.0 |
| width along temperature | 0.0420 ± 0.0023 | 1.04 | 0.95 of 60 | -1.2 |

## The values in windows

The fit read again in windows of growing half width about the centre of each curve, in widths of the curve, and off the peak, every sample beyond a half width, each estimate set against the fit's by its contrast, read at first order. On a right model every window agrees with the fit within the scatter of their difference, and a term that the model lacks moves the windows where it lives first: the core, the wings, or both. A level under 0.05 is a window that disagrees.

| window | samples | width | width along temperature | level of its contrast with the fit |
|---|---|---|---|---|
| within 0.5 | 80 | 2.138 ± 0.058 | 0.039 ± 0.027 | 0.933 |
| within 1 | 160 | 2.112 ± 0.014 | 0.0375 ± 0.0068 | 0.584 |
| within 2 | 320 | 2.1171 ± 0.0076 | 0.0391 ± 0.0036 | 0.286 |
| within 4 | 640 | 2.1151 ± 0.0057 | 0.0416 ± 0.0027 | 0.203 |
| within 8 | 1280 | 2.1167 ± 0.0051 | 0.0432 ± 0.0024 | 0.352 |
| beyond 0.25 | 1564 | 2.1170 ± 0.0089 | 0.0464 ± 0.0042 | 0.671 |
| beyond 0.5 | 1524 | 2.137 ± 0.022 | 0.056 ± 0.010 | 0.171 |
| beyond 1 | 1444 | 2.18 ± 0.12 | 0.075 ± 0.055 | 0.585 |

## To quote

Each parameter as a value with the bar of the fit, then, as separate terms, the bar of theory and the half range over the models where they are not zero. Where one bar is asked, add them in squares, since each comes from another source. The bias on the twin is not taken out while it is under one bar. Where the noise, the coverage or the bias fails, the line says so.

- width: 2.1182 ± 0.0050 of the fit, at the mean of temperature over the curves, 305
- width along temperature: 0.0433 ± 0.0024 of the fit, a share of width for each deviation of temperature over the curves, 11.18: the line below gives it per unit
- height, place, baseline, of each curve: Each is in the section of the parameters of each curve with its bar of the fit, which counts what the shared parameters leave uncertain. The twin reads no coverage and no bias of them: quote each with its bar of the fit alone, and with the settings of its curve
- width grows by 0.00820 ± 0.00044 per unit of temperature: the size of its term times width over 11.1803, the deviation of temperature over the curves, with the bar that the covariance of the two gives. At each value of temperature: 290, 1.9952 ± 0.0089, 300, 2.0772 ± 0.0057, 310, 2.1592 ± 0.0052, 320, 2.2413 ± 0.0077

## Notes

- The bias on the twin here is in bars of the fit, on sets of its own: the report of each model reads it in spreads of the estimate, on other sets, so the two differ.
- The coverage on the twin here is read on sets of its own too, and a coverage counted on 100 sets scatters by about 0.02 from one draw of the sets to another: the coverage of a parameter in the report of its model can differ from this one by two or three times that, with no fault in either.
- The spread over the models reads the models of this folder alone. A form of a term tried in a copy, as a law of the physics in `model.py` in the place of a term of the package, is not in it: compare such forms by their Bayesian criterion and quote the range of their values beside the five parts.
```

On 100 sets the coverages read 0.93 and 0.94, each within two of their standard deviations of 0.95: the bar of the fit holds, and it is the uncertainty to quote.

When the noise fails its check, the coverage and the bias on the twin are marked not to be read, since the twin draws its sets at the noise that the description gives. Quote no bar until the noise is given right, by a column of the curves or a value, as section 9 says, and the check passes.

`best.md` also gives the own parameters of each curve, its height, its place and its baseline, each with its bar. Each run of the loop writes `loop/` again at its own number of sets, so the sections after this one, run at the default of 30, read 30 sets.

And the report of the new model finds nothing more beyond its twin:

`my_nmr/loop/model_1.md`:

```markdown
## Where to search

Read on 100 synthetic sets. The first answer beyond the twin of each question, beyond meaning an excess above 3.01 here:

| question | what it asks | first answer beyond the twin | its excess |
|---|---|---|---|
| where | the core of the curves, their flanks, their wings or their asymmetry | none beyond the twin. The largest: the wings | 0.1 |
| shape | how the residuals are shaped: a slow trend over the span of each curve, a periodic structure, a part in quadrature with the line, or a structure in one place and at one scale | none beyond the twin. The largest: a periodic structure on one side of the curves | 0.7 |
| what | glitches, a noise that grows with the signal, or a noise correlated from sample to sample | none beyond the twin. The largest: a noise that grows with the signal | 2.2 |
| changes along | which part of the residuals changes from curve to curve along a setting, read on the residuals themselves | none beyond the twin. The largest: the asymmetry along temperature | 0.2 |
| seen by | the family of statistics that sees the misfit most | none beyond the twin. The largest: the plane of position and residual value: intervals of position by bands of value | 1.3 |
| noise | the shape of the noise: skewed to one side, or with heavy tails | none beyond the twin. The largest: a noise skewed to one side | 0.6 |

No answer of where, shape, what or changes along, nor of a term of your list, stands beyond the twin: on these records the model leaves nothing that its own synthetic sets do not leave, though no answer of matched was read. A term much smaller than the noise can still be missing, and more sets, or records of more samples, read smaller ones.
```

The visual log of a run is drawn again from what the run left, fitting nothing, by one command:

```bash
python -m twinfit log my_nmr
```

```
The visual log: my_nmr/loop/log.pdf, 29 slides, each phase of each of 2 turns first as one plot, then the equations, the effects of the terms, the curves in a row, the report and the criteria, the pictures in my_nmr/loop/log.
```

The loop drew its visual log too, `my_nmr/loop/log.pdf`, and `python -m twinfit log my_nmr` draws it again from what the run left, fitting nothing. For each model it gives each of the eight phases of its turn first as one plot, with one line of words and, where the phase made one, one number to check, and beside four of them the slides that say more: thirteen kinds of slide a model, a kind on more than one page where its curves, its terms or its proposals ask it. They are the curves as recorded, the model with the bar of its first ten shared parameters, its equation and what each term does to the curves, its twin, the joint fit with what it leaves above it and then the curves in a row, what the fit leaves along the setting, where the residuals leave the twin's band and the answers of the report with the six proposals, the next term tried on the twin, and the decision, what the term changes and each model's criterion. A name longer than its place is cut at a word and ends in ..., and `loop/equations.json` and `loop/best.md` hold it whole. At the end comes the best model with the parts of its bars.

Here are its 29 slides, model 0 and then model 1. On the third slide of model 1 the equation reads the width of the first model times one plus a times z, z the temperature in standard deviations from its mean. The fourth shows what that term does: the fitted curves less the same model with a at zero, wider on the hot records and narrower on the cold ones.

![model 0, phase 1](examples/figures/log_00_1_records.png)

*The four spectra of the project as recorded, one plot along the temperature.*

![model 0, phase 2a](examples/figures/log_00_2a_model.png)

*The first model and the bar of each of its shared parameters, in per cent of its value.*

![model 0, phase 2b](examples/figures/log_00_2b_model.png)

*The equation of the first model, one Lorentzian line of one width for all four spectra.*

![model 0, phase 2c](examples/figures/log_00_2c_model.png)

*What each term of the first model does to the curves.*

![model 0, phase 3](examples/figures/log_00_3_twin.png)

*On synthetic sets drawn from its fit, the first model returns what drew them: each parameter stands inside the white bands.*

![model 0, phase 4a](examples/figures/log_00_4a_fit.png)

*The joint fit of the first model through the four spectra, and above them what it leaves, in noises.*

![model 0, phase 4b](examples/figures/log_00_4b_fit.png)

*The four spectra in a row, in the order of their sweep, each over its height with what the first model leaves under it.*

![model 0, phase 5](examples/figures/log_00_5_residuals.png)

*What the first model leaves on each spectrum, its mean square over the noise, along the temperature.*

![model 0, phase 6a](examples/figures/log_00_6a_where.png)

*The place where the residuals of the first model change along the temperature beyond the twin, shaded on them, the spectra in the order of the temperature.*

![model 0, phase 6b](examples/figures/log_00_6b_where.png)

*The answers of the first report beyond the twin, ordered by their excess, and its proposals, those for the experiment that it has no room for standing on the page after it.*

![model 0, phase 6b2](examples/figures/log_00_6b2_where.png)

*The proposals of the first report for the experiment that the slide before had no room for, whole.*

![model 0, phase 7](examples/figures/log_00_7_next_term.png)

*The width along the temperature, planted on the twin and found again: its gain on the records beside the twin's.*

![model 0, phase 8a](examples/figures/log_00_8a_decision.png)

*What the width along the temperature changes: what the fit leaves before it in grey, after it in colour.*

![model 0, phase 8b](examples/figures/log_00_8b_decision.png)

*The width along the temperature lowers the criterion and is kept.*

![model 1, phase 1](examples/figures/log_01_1_records.png)

*The same four spectra, read by model 1.*

![model 1, phase 2a](examples/figures/log_01_2a_model.png)

*Model 1 and the bar of each of its shared parameters.*

![model 1, phase 2b](examples/figures/log_01_2b_model.png)

*The equation of model 1, the width of the first model times one plus a times z.*

![model 1, phase 2c](examples/figures/log_01_2c_model.png)

*What the new term does: the curves wider on the hot records and narrower on the cold ones.*

![model 1, phase 3](examples/figures/log_01_3_twin.png)

*The twin admits model 1, as a term is kept only where it does.*

![model 1, phase 4a](examples/figures/log_01_4a_fit.png)

*The joint fit of model 1, and what it leaves above it.*

![model 1, phase 4b](examples/figures/log_01_4b_fit.png)

*The four spectra in a row under model 1, each with what it leaves under it.*

![model 1, phase 5](examples/figures/log_01_5_residuals.png)

*What model 1 leaves on each spectrum, along the temperature.*

![model 1, phase 6a](examples/figures/log_01_6a_where.png)

*After model 1 no place of the residuals stands beyond the twin.*

![model 1, phase 6b](examples/figures/log_01_6b_where.png)

*No answer stands beyond the twin after model 1, so no term is asked.*

![model 1, phase 6b2](examples/figures/log_01_6b2_where.png)

*The proposals of the report of model 1 for the experiment that the slide before had no room for, whole.*

![model 1, phase 7](examples/figures/log_01_7_next_term.png)

*No term is tried after model 1, and the loop ends.*

![model 1, phase 8a](examples/figures/log_01_8a_decision.png)

*No model is made after model 1.*

![model 1, phase 8b](examples/figures/log_01_8b_decision.png)

*Model 1 holds the lowest criterion of the run.*

![the best model](examples/figures/log_99_best.png)

*The best model, model 1, each value with its bar of the fit and its bias on the twin.*

Each phase's first slide holds one number to check where its phase made one, and here is where the run keeps it, `K` the model and the notebook `my_nmr/loop/notebook/NNN_model_K_PHASE.json`. Where a model leaves nothing beyond its twin, as model 1 here, the slide of phase 6 says so in its title with no excess, and on the last model of a run, model 1 here too, the slides of phases 7 and 8 say in words that no next term was tried, or that the term tried was refused and why, and that no model came after:

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

```
model 0: the first model is fitted, and its twin and its report are read on 30 synthetic sets
01.py: width changes along temperature is fitted, and read on 30 synthetic sets: each term costs more than the one before it, since the model grows
02.py: a cosine of period 3 is fitted, and read on 30 synthetic sets: each term costs more than the one before it, since the model grows
best.md: the five parts of each value of the best model, on 30 synthetic sets, and the combined estimate of the families
kept: 01.py, width changes along temperature. Refused: 02.py, a cosine of period 3. The best model: model 1. best.md gives the reason of each refusal.
The report of each model and the output are in my_nmr/loop. No answer of where, shape, what or changes along in the last report, nor a term of your list, stands beyond its twin: the records ask for no more, and a further term rests on theory alone. An own value moves from curve to curve, which no question of the report reads: best.md says which, and how to write it.
The visual log: my_nmr/loop/log.pdf, 29 slides, each phase of each of 2 turns first as one plot, then the equations, the effects of the terms, the curves in a row, the report and the criteria, the pictures in my_nmr/loop/log.
```

It refused the ripple, and says why in `best.md`:

`my_nmr/loop/best.md`:

```markdown
Tried and refused:
- a cosine of period 3, in 02.py: the Bayesian criterion rises by 7.20 with it, where it must fall by more than 2: these records do not tell the model with it from the model without it. A term of solid physical ground is better kept whatever the criterion says: write keep = True in its file. Its size, fitted: 0.00015 ± 0.00036. Fitted alone in the place of 01.py, its criterion is 334.78 higher than that of the model with 01.py

A term refused after others were kept may be refused because a kept term stands in its place. To try it alone, copy this folder, keep in the copy's terms/ only that file, renamed 01.py, and run the loop on the copy. To take a kept term out, copy the folder, leave out its file and rename those after it so that the numbers run on.
```

A term is kept when three things hold: its fit ends inside its bounds, the Bayesian criterion falls by more than 2, and the twin admits the model with it. The loop asks for no next term in four cases, in the order of `best.md`, and the line "Next: stop" of the last brief says which. A kind of the noise leads the last report, which a term file cannot hold. One answer alone stands just past its line, which more sets read again below 500 sets, and which at 500 says that the model may lack a term. The last term file tried was refused, and no proposal beyond the twin gains, at first order, what a term of one parameter costs to be kept, 2 and the logarithm of the number of samples. Or no answer stands beyond its twin. A gain is no forecast of the fall that a term brings. The first report of this project gave its proposals beyond the twin gains of 4.0, 2.4, 1.9, each under the cost of 9.4, and the width along the temperature then lowered the chi square by 335.6. So the loop asked for a term there, and it stops on the gains only once a term has been refused. A gain can be far off in either direction, and a proposal within its twin asks for nothing, whatever its gain.

This run names no `--sets`, so it reads 30 sets, and `best.md` gives its coverages and biases on 30 again, as does each run below that names none: the run whose numbers you quote is the last, at the number of sets that you quote it at, here `python -m twinfit loop my_nmr --sets 100` again.

The refused `02.py` stays in `terms/`. The loop reads the term files up to the first number that is missing, and ten at most unless `--n N` says more, so a file taken out would hide every file after it. Each run fits it again and refuses it again, which costs the time of its fit.

A refusal by the twin, or for too few of its searches that ended, read on the 30 synthetic sets of the default, can turn on more. `best.md` says so beside each, and `--sets 100` reads it again before the term is set aside, since a right term was refused on 30 sets and kept on 100.

Kept means that the records are better told with it than without it, and not that it is the term of the physics. On a model that lacks another term, a term of another shape can be kept, since it takes part of that misfit. Its reason, in `why`, is what makes it a term of physics. Such a term can move a shared parameter away from its truth: while the last report still has an answer beyond its twin, `best.md` says that no value is to be quoted yet.

To take a kept term out, copy the folder, leave out its file in the copy and rename those after it so that the numbers run on, and run the loop on the copy. A term refused behind a wrong one is tried alone the same way, in a copy whose `terms/` holds only it, as `01.py`.

## 8. The proposals for the experiment

Each report ends with three proposals for the experiment, each of one of three kinds: a better setup, a better characterization of the setup, or a better data collection.

Each comes from the description, which may name the changes that you could make, each a lever `[addenda.lever.NAME]` of one kind. The kinds are noise, samples, span or repeats with a factor, covariate with new values of a setting, characterization with an ansatz that it measures better, or setup with the item of `[setup]` that it changes, as section 9 gives them. Every lever is measured on the twin but a setup item that no value of the model reads, which the report prints as not measured with its reason, and those that do not rank among the three are said below them.

The package adds the change of each kind that the description does not name, changes that every experiment has: a quieter setup, more samples, more repeats, a wider span, so a description that names three gets the fourth. It draws levers of three more kinds from the description itself. For each setting along which the curves vary and that no lever names, it proposes curves half its range past the values measured. For each item of the setup known only as assumed or not recorded, it proposes its measurement. For each value that an ansatz holds with its bar and that no lever names, it proposes that value measured with half its bar, `NAME measured better`. Each is measured on the twin, but the setup item that no value of the model reads. The twin is drawn again as the change would give it, and fitted, and the error of the known parameters, the shared ones of the model, is set against their error under the records as they are:

`my_nmr/loop/model_1.md`:

```markdown
## Three proposals for the experiment

Each is of one of three kinds: a better setup, a better characterization of the setup, or a better data collection. Each is measured on the twin, with the model as it is, since no proposal for the model stands beyond the twin. The error of a value is the root of its squared bias and its squared scatter there, over the same under the records as they are: under 1, the change helps.

- A better setup, quieter: a quieter detector, which halves the noise of each sample. In the twin, the noise of every sample is multiplied by the factor. The error of the known parameters is 0.50 of that under the records
  - of each parameter, the number above being their geometric mean, the n-th root of the product of the n of them: width 0.50, width along temperature 0.50
- A better data collection, repeats: each run made three times. In the twin, each record is taken the factor times, each with its own noise. The error of the known parameters is 0.52 of that under the records, known within about 0.10 on 30 sets, since the change draws noise of its own
  - of each parameter, the number above being their geometric mean, the n-th root of the product of the n of them: width 0.61 ± 0.11, width along temperature 0.44 ± 0.08
- A better data collection, temperature past its range: records at temperature 335 and 275, half its range past the values measured, 290 to 320. In the twin, records are added at the new values of one covariate, each with the own parameters of the record nearest in that covariate. The error of the known parameters is 0.55 of that under the records, known within about 0.10 on 30 sets, since the change draws noise of its own
  - of each parameter, the number above being their geometric mean, the n-th root of the product of the n of them: width 0.75 ± 0.14, width along temperature 0.40 ± 0.07
- Also measured, not among the three, a better data collection, slower: twice the samples over the same span. The error is 0.75 of that under the records, known within about 0.14
- Also measured, not among the three, a better data collection, a wider span: each record spans one and a half times its span about its middle, with the same step. The error is 1.09 of that under the records, known within about 0.20
```

A proposal reads its change on the twin of that design. `traces` asks the other way round: how many times each record is taken for a bar that you name, by the law of independent records, read again on synthetic sets of that design, here for half the bar of the width, 0.0025 in the unit of the abscissa:

```bash
python -m twinfit traces my_nmr width=0.0025
```

```
The model that the term files give: 01.py, on 4 records.
width: its bar is 0.00497 on these records. To reach 0.0025, each record is taken 4 times, 16 traces, where the law gives 0.00249. On 50 sets of the twin of that design: bias +5e-05, scatter 0.0023, error 0.00231 ± 0.00023, which holds the target.
This holds under the model alone: a term that the records miss biases a parameter by as much whatever the number of traces. Read the report of the loop first.
```

The law takes the records as independent, and the twin of that design checks it: where the error that it reads does not hold the target, the fit is not linear enough there for the law, and the line says to take more. Neither sees a term that the records miss, which biases the width by as much whatever the number of traces, so `traces` comes after the report of the loop, as its last line says.

## 9. The researcher's own curves

This writes the starting files of a project into a new folder:

```bash
python -m twinfit init my_project
```

```
written: my_project/description.toml
written: my_project/model.py
written: my_project/list_of_terms.py
written: my_project/CLAUDE.md
written: my_project/.claude/skills/twinfit/SKILL.md
made: traces/, terms/, agent/, for the curves, the term files and the record of an agent
the package is imported from /the/folder/that/holds/twinfit. It needs numpy and scipy, which /the/python/that/runs/it has. Fill description.toml and model.py, put the curves in traces/, then: python -m twinfit loop my_project
```

Its last line gives the folder of the package and the Python that runs it: here they are written as `/the/folder/that/holds/twinfit` and `/the/python/that/runs/it`, and on your machine they are the paths of yours. Fill the places marked `TO FILL` in the description: the loop refuses a description whose `name` or `what` still reads `TO FILL`, since every report quotes them. A description with no model under `[addenda]` starts from the generic model of the package, one peak or dip of one width on its own height, place and baseline, and says so. Here is the description to fill:

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

Write your first model in `model.py`, from this start, and its `EQUATION`, the model in LaTeX between dollar signs: the visual log typesets it at each turn, and shows the Python line of `model` where it is left empty:

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

`init` also wrote your list of terms, `list_of_terms.py`, which the description names under `[addenda]`: the terms that you suspect, each written as a term file writes one, which the question matched of the report sets against the residuals. It is empty, and with it empty the report reads the residuals alone:

`my_project/list_of_terms.py`:

```python
"""Your list of terms: the terms that you suspect, which the report sets against the residuals, each against its twin,
in the question matched. It is empty at the start, and the package names none: with it empty the report reads the
residuals alone. A term is written here as a term file writes it, a term of the package that `python -m twinfit terms`
lists or a Mechanism of your own, and the names of model.py are read here as a term file reads them."""
from twinfit import mechanisms  # noqa: F401

TERMS = []
```

Put each curve in `traces/`, one file a curve. Its first line names the columns as `columns` names them, then come one line a sample, the abscissa and the value separated by a comma, with no line of comments, as in section 3.

### The noise of each sample

**The noise of each sample**, in four steps:

1. Where you know the noise of a sample, give it, `from = "value"` with its standard deviation, or `from = "column"` with a column of each file.
2. Where you do not know it, leave `[noise]` out, or write `from = "neighbours"`, the default. The package then reads the variance of the noise as a law in the height of the curve above its foot, a + b h + c h^2, a floor, a counter and a noise in proportion to the signal. It reads that law from the second differences of neighboring samples before any fit, which a curve smooth over three samples barely reaches.

   Where a curve is sharper than that, write `from = "residuals"`: the same law, read from the second differences of what the first model leaves, and the model is fitted again with it. Either way a missing term smooth over three samples is not read as noise, on a peak as on a decay.

3. A law of your own, as one from the instrument, goes in a column of each file, or in `model.py` as a function `variance(mean, record, own)` with `from = "model"`, the fifth way of the keys below. The residuals of each model are in `loop/residuals_K.csv`, whose columns are the curve, the abscissa, the value, the fitted value, the residual and the noise that each sample was given.
4. Make each change of the noise in a copy of the folder, made without its `loop/` and `agent/`, as `rsync -a --exclude loop --exclude loop.before --exclude loop.reuse --exclude loop.costs.tsv --exclude agent --exclude .twinfit.lock --exclude .twinfit.command my_project/ my_copy/` makes it. Run the loop on the copy, and quote the copy's `best.md`: the first folder keeps its run as the record of what the noise changed.

The report leads with a kind where that kind is the first answer of what beyond the twin and its excess passes that of every answer of the mean beyond the twin. A misfit of the mean far larger leaves part of itself in the reading of the noise, so there the term comes first, and the noise is read again once the mean holds.

Where it leads with a noise correlated from sample to sample, as a filter of the detector leaves it, its page gives the number of samples over which the correlation falls. Keep one sample in every that many, in a copy of the folder, and run the loop on several such copies, each from another first sample. A term that one keeps and the others refuse is not to be believed.

`python -m twinfit thin my_project N` writes up to three such copies beside the folder, fewer where N is small and two starts fall on one sample, each named `DIR_thin_N_S` with S the first sample that it keeps. Quote the values of one copy with the half range of the three beside them. Its bars, about the root of N times those of the whole folder, are the ones to quote, since the whole folder counts each sample as apart.

Where the report leads with glitches, samples far from the curve in one record alone, the page of the model names them, and no term file holds them. In a copy of the folder, give them a larger noise in a column of the curves, `from = "column"`, or leave them out of the files, and run the loop on the copy. In the column, each such sample takes a noise as large as its residual, the whitened residual that the page gives times the noise of the description, so that it weighs as one sample of noise. Every other sample takes the noise of the description.

### Curves that are not peaks

**Curves that are not peaks.** A decay that starts at a time, a step or an edge is read as a peak is. `centre` gives the place about which the report reads, as the start of a decay: the core is the part near it, the wings the parts far from it, and the asymmetry the difference of its two sides.

The width in the names of the proposals is the width of the fitted curve at half its height above its median value, the unit of their distances. The kinds are shapes of the residuals and not causes. A tail on the upper side of a decay can be a second, slower decay, which no term of the package writes, since `second` is a copy of the curve displaced along the abscissa. Write it in `model.py`, its share and its lifetime as shared parameters, as the component of section 10, and a law of a setting, as the lifetime 1/(1 + K Q) of a quencher Q, there too.

### The loop, and a model that grows

**The loop.** Then run the loop, `python -m twinfit loop my_project`, read the report, write a term, and run it again. A term file that was kept may be written again later: each run fits every file afresh. A value that you hold from theory or from another measurement enters as an ansatz, `[addenda.ansatz.NAME]` with a line `value` and a line `uncertainty`, `NAME` a setting of every curve or a key of a dictionary `HELD` that `model.py` reads. Its uncertainty is then the bar of theory of each parameter that it moves.

`loop` takes five options after the folder: `--n N`, the term files that it reads, at most N, 10 by default, `--sets S`, the synthetic sets of each twin, 30 by default, `--max-seconds T`, its budget on the clock, `--workers W`, the cores on which the twins of the proposals run, every core by default, and `--dictionary`, the package's own list of terms, which you alone ask for. Each run writes a line of what it cost, its seconds, its CPU and its workers, into `loop.costs.tsv` beside `loop/`.

**A model that grows.** Start from a model small enough that you understand every part of it. The report reads what that model lacks, and each term that you add is one answer. Four things are written in `model.py` or in the description, and not in a term file. Each is a change of the folder: make it in a copy, made without `loop/` and `agent/`, so that the run before it stays as its record, as the guide of an agent asks too:

1. A dependence that the physics already gives, as a width that grows with the power, can be written in `model.py` itself, a shared parameter that reads `design["power"]`. It needs no term file, and its value is quoted as it is written.
2. An own value that moves with a setting, as the place of a line that the temperature shifts, is read by no question of the report while it is own to each curve. `best.md` names it with its slope: write it as shared parameters, as the `model.py` below does. Read a moving own value first as what a missing shared term leaves. Where the place of a line moves along a setting and the model holds one place for all the curves, the heights and the baselines can move with it, and they stand still once the place is written. So write first the shared term that the report asks for: an own value that still moves after it is a term of its own.
3. A setting taken in the order of the time cannot be told from a drift of the instrument: write the order under `[problem]` in words, as `order = "from the coldest to the hottest, one an hour"` or `order = "not known"`, and `best.md` says it beside the values.
4. Where the curves sit far from zero against their span, the span from the least abscissa of the curves to the largest, their middle more than three spans from zero, as lines near 587.6 nm a few hundredths of a nanometer wide, write `abscissa_offset = 587.6` under `[problem]`, and a refusal of the loop names the value. The abscissa is then read less it everywhere, in `model.py`, in `start_own`, in the term files, in every place that the output gives and in `residuals_K.csv`. The fit moves the place of a curve by steps that grow with the size of the place, which the offset keeps small.

   So the starts and the bounds of a place, in `START_SHARED` and `BOUNDS_SHARED`, are written less it too, and `best.md` gives each place less it too, with a line that says so: add the offset back to quote a place.

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

This is the whole file, and in the copy it takes the place of the `model.py` of section 3. Make the copy as step 4 of the noise above makes one, `rsync -a --exclude loop --exclude loop.before --exclude loop.reuse --exclude loop.costs.tsv --exclude agent --exclude .twinfit.lock --exclude .twinfit.command my_nmr/ my_nmr_place/`, write this file in it, and run `python -m twinfit loop my_nmr_place`. The copy keeps the term files of `my_nmr`, `01.py` kept and `02.py` refused again, and its `best.md` gives `place_at_300` and `shift_per_kelvin` with their bars where the run of `my_nmr` named the places as moving. That file's `BOUNDS_SHARED` held one bound for its one shared parameter, and this one holds none, since a bound is given for each shared parameter or for none. A start is given for each, in `START_SHARED` or by a function `start_shared(records)` that reads them from the curves.

`CENTRE` names one of `OWN`, any one, since the package asks that it be one of them. Here it is a name and no more, and `centre` returns the place from the shared parameters. So the report reads the core and the wings about that place, and a pedestal or a periodic term, placed about it, can still be fitted. A term file that moves the place leaves `centre` where `model.py` puts it, so the report would read the core and the wings about the old place: a place that moves is written here, with `centre` returning it. On the spectra of section 3 this model gives the width that the model of own places gives, and the shift for each kelvin with its bar.

### The keys of the description, and a held value

The keys of the description beyond those of section 3, each with its form. The noise of each sample is given one of five ways. The output is for people, for agents or both, and the search of the literature is none, local, online or both, run by an agent of the model and the `effort` given.

An ansatz is a value held before the fit, with its uncertainty. A lever is a change of the experiment that the proposals weigh. Its kinds are `noise`, `samples`, `span` or `repeats` with a factor, and `covariate` with the setting and its new values. They are also `characterization` with the ansatz that it measures better and a factor on its uncertainty, and `setup` with `of`, the name of an item of `[setup]` that it changes, and a factor, a half where none is given. A lever of kind noise or setup is proposed as a better setup, one of kind characterization as a better characterization of the setup, and the other four as a better data collection.

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

An ansatz on a setting holds it at one value, which every curve gives: a setting that changes from curve to curve is fitted along, and not held. An ansatz on a key of `HELD` sets that key before each fit, and `model.py` reads it where it computes the curve.

A `model.py` that reads a held value, its `HELD` and its model:

```python
HELD = {"resolution": 0.5}          # the values that the model reads and does not fit, which an ansatz sets


def model(x, shared, own, design):
    width = shared[0] + HELD["resolution"]          # a Lorentzian seen through a Lorentzian resolution: their widths add, where a Gaussian resolution makes a Voigt line
    return own[0] / (1.0 + (2.0 * (x - own[1]) / width) ** 2) + own[2]
```

The fitted `width` is then the line's own, without the resolution, so on curves drawn with no resolution it reads 0.5 less than without `HELD`. Its uncertainty is then the bar of theory of each parameter that it moves, and a lever of kind characterization of it says what measuring it better would buy. A number written in a term file, as the distance of a second component, is held by nothing: where a measurement could fix it, move it into `model.py` as a key of `HELD`, with its ansatz.

## 10. A term of unknown distance or period

Some terms carry a number that is chosen and not fitted: the distance of a second component, the stretch of a tail, the period of a ripple. The worked project of chromatography has one:

```bash
python -m twinfit example chromatography my_chrom
```

```
written into my_chrom: description.toml, model.py, an empty terms/, and 4 curves in traces/. Next: python -m twinfit loop my_chrom
```

```bash
python -m twinfit loop my_chrom
```

```
model 0: the first model is fitted, and its twin and its report are read on 30 synthetic sets
best.md: the five parts of each value of the best model, on 30 synthetic sets, and the combined estimate of the families
kept: nothing. Refused: nothing. The best model: model 0.
The report of each model and the output are in my_chrom/loop. Write my_chrom/terms/01.py with the next term, then run the loop again. An own value moves from curve to curve, which no question of the report reads: best.md says which, and how to write it.
The visual log: my_chrom/loop/log.pdf, 15 slides, each phase of each of 1 turn first as one plot, then the equations, the effects of the terms, the curves in a row, the report and the criteria, the pictures in my_chrom/loop/log.
```

The report puts the misfit of the chromatograms in their asymmetry:

`my_chrom/loop/model_0.md`:

```markdown
## Where to search

Read on 30 synthetic sets. The first answer beyond the twin of each question, beyond meaning an excess above 3.73 here:

| question | what it asks | first answer beyond the twin | its excess |
|---|---|---|---|
| where | the core of the curves, their flanks, their wings or their asymmetry | the asymmetry | 4.0 |
| shape | how the residuals are shaped: a slow trend over the span of each curve, a periodic structure, a part in quadrature with the line, or a structure in one place and at one scale | a part in quadrature with the line, the Hilbert transform of the fitted curve | 8.3 |
| what | glitches, a noise that grows with the signal, or a noise correlated from sample to sample | none beyond the twin. The largest: a noise that grows with the signal | 0.6 |
| changes along | which part of the residuals changes from curve to curve along a setting, read on the residuals themselves | none beyond the twin. The largest: the asymmetry along temperature | 1.1 |
| seen by | the family of statistics that sees the misfit most | none beyond the twin. The largest: moments of the running integral: powers 0 to 8 | 3.6 |
| noise | the shape of the noise: skewed to one side, or with heavy tails | none beyond the twin. The largest: a noise skewed to one side | 0.0 |
```

The chi square of this first model, 1253.0 on 1191 degrees of freedom, agrees with the noise, while the report puts the misfit in the asymmetry at an excess of 4.0 against its line of 3.73, and the loop says above whether that asks for more: a term that holds a small part of the curve moves the whole chi square by about its own scatter, while a statistic aimed at one shape reads it. Whether that asymmetry is a second component is the physics of the chromatogram, which the package does not name.

In a chromatogram an asymmetry of the peaks is a compound that elutes close to the first, an impurity, or a tail of the column. The first proposal for the model is the asymmetry, measured with the shape that the residuals hold there:

`my_chrom/loop/model_0.md`:

```markdown
## Three proposals for the model

Each is a direction read on the residuals, for an answer of where, of shape or of what changes along a setting: those beyond the twin first, then by the gain in chi square that the shape that the residuals hold there brings when it is fitted beside the model. The package names no term: which term of physical origin has that place, shape and setting is the researcher's to say, in a term file. On its twin, synthetic curves drawn with the term at its fitted size are fitted without it and with it: the bias without it is what the term would cost the known parameters if it were real, and says nothing of whether it is. The shape that the residuals hold in a place is freed of the fit's own directions, so it costs them nothing at first order, and its line says so: a shape built apart from the residuals, a slow trend, a periodic structure, a part in quadrature with the line or a wavelet, is not freed, and its cost is read as any other.

In the names of the proposals a width is the full width of the fitted curves at half their height above their median value, here 2.4 in the unit of the abscissa: where the wings of a curve fill its window, it is narrower than the width at half its peak, and it is not a parameter of the model. Above the centre is toward larger abscissa.

- the shape that the residuals hold in the asymmetry, for the asymmetry, beyond the twin at an excess of 4.0: a gain in chi square of 283.8
  - read on the residuals, this shape costs width nothing at first order: on its twin width is biased by 0.1 scatters, 0.00027 in its own unit, without it and by 0.1 scatters, 0.000286 in its own unit, with it, the twin's own bias, since the residuals hold none of what the fit took into its own parameters. What a term in this place would cost lies in that part, which no residual shows: fit a term of the physics here and read how far width moves
- a structure of 2 widths at +0.79 widths from the centre, for a structure in one place and at one scale, a wavelet of the residuals, beyond the twin at an excess of 3.9: a gain in chi square of 34.0
  - if the term were real at that size, width would be biased on its twin by 0.9 scatters, 0.0025 in its own unit, without it and by 0.1 scatters, 0.000207 in its own unit, with it: what the term would cost, and not a sign that it is real. A scatter of a model that lacks a term can be many bars of its fit
- a part in quadrature with the line, the Hilbert transform of the fitted curve, for a part in quadrature with the line, the Hilbert transform of the fitted curve, beyond the twin at an excess of 8.3: a gain in chi square of 20.0
  - if the term were real at that size, width would be biased on its twin by 0.2 scatters, 0.00048 in its own unit, without it and by 0.1 scatters, 0.000266 in its own unit, with it: what the term would cost, and not a sign that it is real. A scatter of a model that lacks a term can be many bars of its fit

The asymmetry leads. Where no proposal above holds it once fitted, a shape that no term of the package writes may be missing, asymmetric about the place of the line: an exponential fall of the line on one side, its height times x over one plus x squared, x the distance from its place in its widths, or the derivative of the line. A Mechanism of your own writes it, as section 11 of the tutorial shows, and the physics of the records says which.

## The known parameters before the term is named

Where an answer stands beyond the twin, each shared parameter is read five ways while the missing term has no name: the fit as it is, the combination of the families that is free of the leading class of the slow trend, the one basis fixed before the residuals are read, the model with the generic term of the largest gain, width changes along temperature, fitted and at first order, and the combination free of that term. The generic terms are those that the model and the design alone give, each shared parameter along each covariate, and no list of terms is read. Where the five agree, the missing term does not move the parameter in the directions tried. Where they part, the spread says how far it may. None of them sees a part of the missing term that lies in the fit's own directions and in no basis tried.

| parameter | the fit | free of the class | with the generic term | that term at first order | free of that term |
|---|---|---|---|---|---|
| width | 1.0306 ± 0.0030 | none | 1.0311 ± 0.0031 | 1.0311 ± 0.0030 | 0.855 ± 0.012 |
```

Try the impurity first, a second component of the curve, `second(distance)`, with its distance in the unit of the abscissa. In the page of the proposals a width is the full width of the fitted curves at half their height above their median value, 2.4 s here, while the parameter `width` of this model is the standard deviation of its Gaussian, 1 s. The distance is not read off the report, which names no term.

`scan` fits the term at each distance that you give, beside the model, and says which lowers the Bayesian criterion most. The distances are in the unit of the abscissa, here seconds, and `x` stands for the number that is scanned: `"second(x)"` is `second(distance)` at each value given. It chooses the stretch of `tail` or the period of `periodic` the same way, with the other argument in quotes, `"tail(x, 'upper')"` or `"periodic(x, 'sine')"`.

It refuses the term at a number that a kept term file of the folder already holds, and fits a second term of the same kind beside it at any other number. The model that it scans beside holds each term file that says `hold`, at its theory, and each other whose fit ends inside its bounds and lowers the Bayesian criterion by more than 2, or whose file says `keep = True`, with no twin run, so a file that the loop refused on its twin alone is in it. To scan the number of a term that a file holds, scan in a copy of the folder without that file.

Where the range that it gives is under the step of the values scanned, a parabola through the three falls about the largest gives it, and the top of that parabola is the value to write. That range is the bar of the chosen number where the noise is white: under a noise correlated over N samples it is about the root of N too narrow, so read it after the noise is thinned.

A distance has a side, and a term on the wrong side can still lower the criterion with a negative size, which takes its shape away from the curve. Scan both signs, as -2 -1 1 2, and keep the side of the larger fall and of a positive size. Read a kept term of negative size as the command line and `best.md` flag it. `displaced` points the other way: its long side is toward smaller abscissa where its reach is positive, as `python -m twinfit terms` says.

```bash
python -m twinfit scan my_chrom "second(x)" 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.4
```

```
second(x) beside the model of my_chrom, its first model, at each x, in the unit that the term takes: a distance or a period in the unit of the abscissa, a stretch or a power as a number. The fall of the Bayesian criterion is from the model without the term, the size is that of the term and the bar its bar of the fit. No twin is run here: the loop decides.
| x | its fit | fall of the Bayesian criterion, negative where it rises | size | bar |
|---|---|---|---|---|
| 0.6 | did not end | -7.1 | 0.000 | 0.074 |
| 0.8 | did not end | 27.1 | 0.25 | 0.18 |
| 1 | ended | 38.4 | 0.158 | 0.038 |
| 1.2 | ended | 41.5 | 0.089 | 0.015 |
| 1.4 | ended | 42.3 | 0.0577 | 0.0087 |
| 1.6 | ended | 42.0 | 0.0414 | 0.0060 |
| 1.8 | ended | 41.0 | 0.0318 | 0.0046 |
| 2 | ended | 39.2 | 0.0256 | 0.0037 |
| 2.4 | ended | 34.2 | 0.0183 | 0.0028 |
The largest fall: at x = 1.4. The values of x whose fall is within 1 of it, which is one standard deviation of x as the profile gives it, the criterion of the best fit at each x, if the model is otherwise right: 1.20 to 1.60. Where the range is wide, these records fix the size and x together and not each. To try it, write `term = mechanisms.second(1.4)` in the next term file, then run the loop.
Across that range the other values of the model move: width from 0.9884 to 1.0029. Half of each move is the part of its uncertainty that the chosen x carries: quote it beside the five parts of best.md.
```

At 0.6 and 0.8 s, the nearest distances, the second component lies so close to the peak that the fit cannot tell its size from the shape of the peak, and its search does not end. A fit that does not end is no answer, and the scan reads the others.

The distances from 1.2 to 1.6 s lower the criterion within 1 of the largest. That range runs between values scanned, so it is as coarse as their step: scan again about it at a step finer than the range, before it is read or held:

```bash
python -m twinfit scan my_chrom "second(x)" 1 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9
```

```
second(x) beside the model of my_chrom, its first model, at each x, in the unit that the term takes: a distance or a period in the unit of the abscissa, a stretch or a power as a number. The fall of the Bayesian criterion is from the model without the term, the size is that of the term and the bar its bar of the fit. No twin is run here: the loop decides.
| x | its fit | fall of the Bayesian criterion, negative where it rises | size | bar |
|---|---|---|---|---|
| 1 | ended | 38.4 | 0.158 | 0.038 |
| 1.05 | ended | 39.6 | 0.135 | 0.028 |
| 1.1 | ended | 40.4 | 0.116 | 0.022 |
| 1.15 | ended | 41.0 | 0.101 | 0.018 |
| 1.2 | ended | 41.5 | 0.089 | 0.015 |
| 1.25 | ended | 41.8 | 0.079 | 0.013 |
| 1.3 | ended | 42.1 | 0.070 | 0.011 |
| 1.35 | ended | 42.2 | 0.0635 | 0.0097 |
| 1.4 | ended | 42.3 | 0.0577 | 0.0087 |
| 1.45 | ended | 42.3 | 0.0527 | 0.0078 |
| 1.5 | ended | 42.3 | 0.0484 | 0.0071 |
| 1.55 | ended | 42.2 | 0.0447 | 0.0065 |
| 1.6 | ended | 42.0 | 0.0414 | 0.0060 |
| 1.65 | ended | 41.8 | 0.0385 | 0.0056 |
| 1.7 | ended | 41.6 | 0.0360 | 0.0052 |
| 1.75 | ended | 41.3 | 0.0338 | 0.0049 |
| 1.8 | ended | 41.0 | 0.0318 | 0.0046 |
| 1.85 | ended | 40.6 | 0.0300 | 0.0043 |
| 1.9 | ended | 40.2 | 0.0283 | 0.0041 |
The largest fall: at x = 1.45. The values of x whose fall is within 1 of it, which is one standard deviation of x as the profile gives it, the criterion of the best fit at each x, if the model is otherwise right: 1.20 to 1.70. Where the range is wide, these records fix the size and x together and not each. To try it, write `term = mechanisms.second(1.45)` in the next term file, then run the loop.
Across that range the other values of the model move: width from 0.9884 to 1.0055. Half of each move is the part of its uncertainty that the chosen x carries: quote it beside the five parts of best.md.
```

At this step the distances from 1.2 to 1.7 s lower the criterion within 1 of the largest, against 1.2 to 1.6 s at the first. The records do not tell them apart by the criterion, and the loop, with its twin, decides. Write the term at the distance that the finer scan gives, and run the loop:

Save this as `my_chrom/terms/01.py`:

```python
from twinfit import mechanisms
term = mechanisms.second(1.45)
why = "the report: a second component in the asymmetry. An impurity that elutes later, at the distance that the scan gives"
```

```bash
python -m twinfit loop my_chrom
```

```
model 0: the first model is fitted, and its twin and its report are read on 30 synthetic sets
01.py: a second component 1.45 away is fitted, and read on 30 synthetic sets: each term costs more than the one before it, since the model grows
best.md: the five parts of each value of the best model, on 30 synthetic sets, and the combined estimate of the families
kept: 01.py, a second component 1.45 away. Refused: nothing. The best model: model 1.
The report of each model and the output are in my_chrom/loop. No answer of where, shape, what or changes along in the last report, nor a term of your list, stands beyond its twin: the records ask for no more, and a further term rests on theory alone. An own value moves from curve to curve, which no question of the report reads: best.md says which, and how to write it.
The visual log: my_chrom/loop/log.pdf, 29 slides, each phase of each of 2 turns first as one plot, then the equations, the effects of the terms, the curves in a row, the report and the criteria, the pictures in my_chrom/loop/log.
```

The curves were drawn with an impurity 1.4 s after the peak, of 0.06 of its height. The scan gives the distance from 1.2 to 1.7 s, and 1.4 s is in that range. The size moves with the distance across it, from 0.089 at 1.2 s to 0.036 at 1.7 s: these records fix the distance and the size together and not each.

At 1.45 s the loop kept the term with a size of 0.0527 ± 0.0078, and the width of the peak, 0.9984 ± 0.0055 s with the term, was 1.031 s without it, where 1 s drew the curves. The width moves across the range of the distance too, as the last line of the scan says. Half of its move is a part of its uncertainty that no column of `best.md` holds, since the distance is written and not fitted, and it is quoted beside them.

A measurement of the retention time of the impurity alone, a better characterization of the setup, would fix the distance, and the size and the width with it.

To carry the distance into the uncertainty of every other value, hold it. It takes a key of `HELD` in `model.py`, an ansatz with half the range that the scan gives as its uncertainty, and a term that reads it at each call. That range comes from these records, so this ansatz narrows nothing: it carries the uncertainty of a chosen number into the other values, where an ansatz that narrows a value must come from outside the records.

A term of the package takes its number when its file is read, so `held.reading` builds it again at each call, from the value that `HELD` holds then. The new `01.py` takes the place of the one above, and it reads `HELD` by its name with no import. The loop gives every term file the names of `model.py`, its functions and its constants as well as `HELD`, so that a term reads them with no import:

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

```
model 0: the first model is fitted, and its twin and its report are read on 30 synthetic sets
01.py: a second component 1.45 away is fitted, and read on 30 synthetic sets: each term costs more than the one before it, since the model grows
best.md: the five parts of each value of the best model, on 30 synthetic sets, and the combined estimate of the families
kept: 01.py, a second component 1.45 away. Refused: nothing. The best model: model 1.
The report of each model and the output are in my_chrom/loop. No answer of where, shape, what or changes along in the last report, nor a term of your list, stands beyond its twin: the records ask for no more, and a further term rests on theory alone. An own value moves from curve to curve, which no question of the report reads: best.md says which, and how to write it.
The visual log: my_chrom/loop/log.pdf, 29 slides, each phase of each of 2 turns first as one plot, then the equations, the effects of the terms, the curves in a row, the report and the criteria, the pictures in my_chrom/loop/log.
```

`my_chrom/loop/best.md`:

```markdown
## To quote

Each parameter as a value with the bar of the fit, then, as separate terms, the bar of theory and the half range over the models where they are not zero. Where one bar is asked, add them in squares, since each comes from another source. The bias on the twin is not taken out while it is under one bar. Where the noise, the coverage or the bias fails, the line says so.

- width: 0.9984 ± 0.0055 of the fit ± 0.0086 of theory
- second component at 1.45: 0.0527 ± 0.0078 of the fit ± 0.026 of theory
- height, place, baseline, of each curve: Each is in the section of the parameters of each curve with its bar of the fit, which counts what the shared parameters leave uncertain. The twin reads no coverage and no bias of them: quote each with its bar of the fit alone, and with the settings of its curve

## Notes

- The bias on the twin here is in bars of the fit, on sets of its own: the report of each model reads it in spreads of the estimate, on other sets, so the two differ.
- The coverage on the twin here is read on sets of its own too, and a coverage counted on 30 sets scatters by about 0.04 from one draw of the sets to another: the coverage of a parameter in the report of its model can differ from this one by two or three times that, with no fault in either.
- The spread over the models reads the models of this folder alone. A form of a term tried in a copy, as a law of the physics in `model.py` in the place of a term of the package, is not in it: compare such forms by their Bayesian criterion and quote the range of their values beside the five parts.
- second component at 1.45: the number in its name is chosen and not fitted, so it has no bar of the fit, where the row above gives the size of the term and its bar: quote that number with the range that `python -m twinfit scan` gives about it, and hold it with an ansatz of half that range where its uncertainty is to reach the other values, or write it in `model.py` as a shared parameter, which the fit then gives a bar and the twin a coverage. The size of a second component is its area over that of the curve, and the place stays that of the first line, so the middle of a doublet is the place plus half its distance.
- second component at 1.45: its size holds at the number written in its term file, or in `HELD` of the model where the term reads it there, a distance, a stretch or a period. `python -m twinfit scan` says how far the records fix that number, and the size moves with it.

## The parameters of each curve

Each with its bar of the fit, which counts what the shared parameters leave uncertain: the twin reads the shared ones alone. An own value that moves with a setting is a term of its own: read them down each column.

- height moves with temperature: a line through its values gives 0.01016 ± 0.00026 per unit of temperature, 39 bars from none, a chi square of 0.9 on 2 about the line. While height is own to each curve, no question of the report reads it: write it in `model.py` as shared parameters that read `design["temperature"]`, as section 9 of the tutorial shows, and the report reads what is left. Across 4 curves a line can also hold by chance, so ask whether the physics gives it.
- place moves with temperature: a line through its values gives 0.03016 ± 0.00040 per unit of temperature, 76 bars from none, a chi square of 0.4 on 2 about the line. While place is own to each curve, no question of the report reads it: write it in `model.py` as shared parameters that read `design["temperature"]`, as section 9 of the tutorial shows, and the report reads what is left. Across 4 curves a line can also hold by chance, so ask whether the physics gives it.

| curve | temperature | height | place | baseline |
|---|---|---|---|---|
| traces/temperature_30.csv | 30 | 1.0086 ± 0.0056 | 11.5549 ± 0.0094 | 0.0504 ± 0.0012 |
| traces/temperature_40.csv | 40 | 1.1038 ± 0.0057 | 11.8484 ± 0.0091 | 0.0507 ± 0.0012 |
| traces/temperature_50.csv | 50 | 1.2123 ± 0.0058 | 12.1537 ± 0.0087 | 0.0509 ± 0.0012 |
| traces/temperature_60.csv | 60 | 1.3114 ± 0.0060 | 12.4578 ± 0.0084 | 0.0491 ± 0.0012 |

Tried and refused: none.

The loop stops here: no answer of where, what or along stands beyond the twin, so the records ask for no more term. A term that the physics asks for is held at its value of theory, in terms/02.py, with hold = True.
```

The bar of theory of each value is now the move that the distance makes across its uncertainty, and it is quoted beside the bar of the fit. The column says of theory whatever the source of the ansatz: here it is the range of a number that these records chose.

Where the records fix the distance well, as the coupling of a doublet whose two lines stand apart, a third way gives it a bar of its own. Write the second line in `model.py`, with the distance as a shared parameter, and the fit gives the distance its bar and the twin its coverage, as every other shared value.

The term `second` places its copy from the place of the curve, so the place stays that of the first line. The middle of an equal doublet is the place plus half the distance, and the report reads the core and the wings about the first line. The size of a second component, a pedestal or a tail is its area over that of the curve, and never a distance.

Where the curves are doublets, start every curve on the same line. The `start_own` of the template takes the first sample from the left that reaches four fifths of the largest distance from the baseline, so that each curve starts on its lower line. A page whose places differ by more than a width while the fit leaves much says that the starts may sit on different lines.

A component of its own width at its own distance, as the outflow under the emission line of a galaxy, is no term of the package, which adds one parameter a term. Where the report puts a tail and a second component tied, write the whole line in `model.py`, its share, its distance and its width as shared parameters.

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

with `start_own`, `CENTRE` and `centre` as in section 3. It is a part to adapt, and no run of this page fits it.

Where the second component has the width of the first, as the impurity of the chromatograms above, its own width cannot be told from its share and its distance, and the term `second` is the one to keep. Any function of the settings is written the same way, as a height that falls as the exponential of a radius or a place that follows the tanh of it, and not only a line.

## 11. A term written by the researcher

When the term that the report points to is not among those of `python -m twinfit terms`, write it, in five steps.

1. **The term.** A term is a `Mechanism(name, family, parameter, start, wrap, lower, upper, formula=...)`: a name, a word for its family, any word, as `baseline`, which the package keeps with the term and reads in no result, the name of the one parameter that it adds, where the search of that parameter starts, `wrap`, the bounds of the parameter, none where `lower` and `upper` are not given, and `formula`, the term in LaTeX between dollar signs with `f` the model before it and `a` its parameter, as `r"$f(x) + a\,x$"`, which the visual log typesets: without it the log says that its formula is not written. `python -m twinfit card` prints a term file of each kind in one page. The name of the parameter is printed as it is written, so it carries its unit, as "shift of the centre per microwatt, in Hz per microwatt".
2. **Its wrap.** `wrap(model, problem, n)` takes the model as it stands, the problem, which holds every curve, and the number of shared parameters before the term. It returns the model with the term, a function `new(x, shared, own, design)`, where `shared[n]` is the new parameter and `shared[:n]` the parameters before it. The loop refuses a `wrap` that does not take these three arguments. `design` holds the settings of the curve, as the description names them, `own[k]` and `shared[k]` the parameters of the model in the order of `OWN` and `SHARED` of `model.py`, and `x`, `shared` and `own` are arrays of numpy. A term that moves the abscissa calls the model at a moved `x`, as `model(x - shared[n] * design["power"], shared[:n], own, design)`, which moves each curve by `shared[n]` times its power. A held number is read inside `new` as `HELD["name"]`, since the loop gives every term file the `HELD` of `model.py`. A term that moves the place leaves `centre` where `model.py` puts it, so a place that moves is written in `model.py`, as the `model.py` of section 9 whose place moves along the temperature does.
3. **One of three forms.** `new` may add a term to the model, as a step under a peak whose size is a share of its height, `model(x, shared[:n], own, design) + shared[n] * own[0] * 0.5 * (1.0 + erf((x - own[1]) / (2.0 ** 0.5 * shared[0])))`, with `erf` of `scipy.special`. It may multiply the model, as a gain of the detector that changes along the abscissa, `model(x, shared[:n], own, design) * (1.0 + shared[n] * (x - middle))`. That multiplies the terms kept before it too, about a middle that you choose as the zero of the change. Or it may replace the line and return a curve of its own, with no call of `model`: then it returns the whole curve, its baseline too, since what it does not add again is gone. The loop refuses such a term where its curve no longer reads an own parameter that the model read, and says which.

4. **A curve that replaces the line**, written out whole, as `01.py` of a fresh project of chromatography, `my_chrom_tail`, which holds no other term file. It is a Gaussian with an exponential tail, for the chromatograms of section 10, whose model holds the width `s = shared[0]`, the height `own[0]`, the place `own[1]` and the baseline `own[2]`. It is the other form of the asymmetry that section 10 read, a tail of the column where section 10 wrote an impurity, and it is run here beside it:

```bash
python -m twinfit example chromatography my_chrom_tail
```

```
written into my_chrom_tail: description.toml, model.py, an empty terms/, and 4 curves in traces/. Next: python -m twinfit loop my_chrom_tail
```

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

```
model 0: the first model is fitted, and its twin and its report are read on 30 synthetic sets
01.py: an exponential tail of the line is fitted, and read on 30 synthetic sets: each term costs more than the one before it, since the model grows
best.md: the five parts of each value of the best model, on 30 synthetic sets, and the combined estimate of the families
kept: 01.py, an exponential tail of the line. Refused: nothing. The best model: model 1.
The report of each model and the output are in my_chrom_tail/loop. No answer of where, shape, what or changes along in the last report, nor a term of your list, stands beyond its twin: the records ask for no more, and a further term rests on theory alone. An own value moves from curve to curve, which no question of the report reads: best.md says which, and how to write it.
The visual log: my_chrom_tail/loop/log.pdf, 29 slides, each phase of each of 2 turns first as one plot, then the equations, the effects of the terms, the curves in a row, the report and the criteria, the pictures in my_chrom_tail/loop/log.
```

The loop keeps the tail, and the twin admits it. Its model reads a Bayesian criterion of 1303.31, against 1302.92 for the impurity of section 10 at 1.45 s: within 2 of each other, so the records cannot tell the two forms apart. The width reads 0.9591 ± 0.0075 s under the tail and 0.9984 ± 0.0055 s under the impurity, where 1 s drew the curves: 5.5 and 0.3 bars of the fit from it. Neither `best.md` holds the gap, since each form stands in a folder of its own.

Where the physics does not choose the form, the value to quote is conditional on the term chosen, and half the gap of the two widths, 0.020 s, is a part of its uncertainty, quoted beside the bars as the spread over the models. Here the curves were drawn with an impurity, which the twin of neither form says: a measurement of the impurity alone, or of the answer of the detector, chooses.

   `own[0]` is the height that the curve keeps as `t` falls to zero, and its area is that of the Gaussian. So a tail that holds a share `a` of a curve is `(1 - a)` times the Gaussian and `a` times this curve, with `a` a shared parameter. With `erfcx` and `erfc` it gives a number at every distance: the product of an exponential and `erfcx` alone gives none past about 38 widths. Where the time constant follows a setting of the curves named load, write `t = shared[n] * design["load"]`.

5. **Its start and its place among the files.** Start the new parameter at a size that the residuals give, and a parameter that the curve divides by, as `t`, near the width of the curve and away from zero. A term reads the shared parameters as the terms before it leave them: a tail in `02.py` after `along("width", ...)` in `01.py` reads the width that `along` makes at each setting. A term that replaces the line drops the terms kept before it, so it goes first, as `01.py`.

   A term must be smooth in its parameters, since one whose curve jumps as a parameter moves keeps the searches from ending. A term whose curve holds a value that is not a number is refused with its file named: compute its `wrap` on the curves alone, with the starts of the model, before the loop. The loop reads the term files in their order up to the first that is missing, so `03.py` is read only beside a `02.py`, and it reads the `why` of each file before its first fit.

Here is a slope of the baseline that grows with the temperature, written in the project of section 3:

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

```
model 0: the first model is fitted, and its twin and its report are read on 30 synthetic sets
01.py: width changes along temperature is fitted, and read on 30 synthetic sets: each term costs more than the one before it, since the model grows
02.py: a cosine of period 3 is fitted, and read on 30 synthetic sets: each term costs more than the one before it, since the model grows
03.py: a slope of the baseline that grows with the temperature is fitted, and read on 30 synthetic sets: each term costs more than the one before it, since the model grows
best.md: the five parts of each value of the best model, on 30 synthetic sets, and the combined estimate of the families
kept: 01.py, width changes along temperature. Refused: 02.py, a cosine of period 3; 03.py, a slope of the baseline that grows with the temperature. The best model: model 1. best.md gives the reason of each refusal.
The report of each model and the output are in my_nmr/loop. No answer of where, shape, what or changes along in the last report, nor a term of your list, stands beyond its twin: the records ask for no more, and a further term rests on theory alone. An own value moves from curve to curve, which no question of the report reads: best.md says which, and how to write it.
The visual log: my_nmr/loop/log.pdf, 29 slides, each phase of each of 2 turns first as one plot, then the equations, the effects of the terms, the curves in a row, the report and the criteria, the pictures in my_nmr/loop/log.
```

The curves of this project were drawn with a flat baseline, so the loop refuses the slope as it refused the ripple. A refused term stays where it is, as section 7 says, and each run tries it again and refuses it again, with its reason, in `best.md`.

The terms of the package are made of copies of the curve. A shape that the physics names otherwise, as the exponential tail of a detector, or a ripple fixed to the abscissa, the same on every curve whatever its center, as the fringe of an etalon fixed to the abscissa, is a term of your own. The report gives its place and its kind, and you give its form.

`periodic` is placed about the center of each curve. Where the curves have centers apart, a fringe fixed to the abscissa read by it leaves a misfit that changes from curve to curve, which the report then reads along a setting. Here is that fringe, with its period as a chosen number. A term file of your own that defines `term_at(x)`, the term at the number x, is scanned by naming the file. A term file reads the dictionary `HELD` of `model.py` by its name, so that a number held by an ansatz reaches the term:

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

```
terms/04.py beside the model of my_nmr with 01.py, at each x, in the unit that the term takes: a distance or a period in the unit of the abscissa, a stretch or a power as a number. The fall of the Bayesian criterion is from the model without the term, the size is that of the term and the bar its bar of the fit. No twin is run here: the loop decides.
| x | its fit | fall of the Bayesian criterion, negative where it rises | size | bar |
|---|---|---|---|---|
| 2 | ended | -7.3 | -0.00010 | 0.00035 |
| 3 | ended | -6.8 | 0.00027 | 0.00036 |
| 4 | ended | -6.9 | 0.00025 | 0.00037 |
No value of x lowers the Bayesian criterion by more than 2, the least that the loop keeps: these records do not ask for this term at any x scanned.
```

No period lowers the criterion: these curves hold no fringe of this phase. A ripple of unknown phase is a cosine and a sine, two terms, and a scan of the cosine alone can find no period where the sine would: scan each. The scan reads the model of the other term files, a file that says `hold` held at its theory and each other kept when its fit ends inside its bounds and its criterion falls by more than 2, or when its file says `keep = True`, with no twin, and leaves the scanned file out of it.

## 12. With an agent

The protocol of an agent's turn is `agent/AGENT.md`, which `python -m twinfit guide` prints. It needs Claude Code, or an agent that reads files and runs commands: nothing in this section runs without one. `init` also wrote `CLAUDE.md` and a skill for Claude Code. Open a fresh Claude Code chat in the project folder and ask it, for example: "Run the loop, read the report, and propose the next term." Who chooses at step 5 is your setting, `DIR/agent.toml`: `chooses = "agent"` lets it choose at each turn and say why, `chooses = "researcher"` has it ask you at each turn, and with no such file it asks. It reads `CLAUDE.md` first. At each iteration it follows ten steps, which `python -m twinfit guide` prints whole with the rest of its protocol, its section 0 and its part What you give back among them. By their names:

1. Read the structure of the residuals.
2. Search the literature.
3. Estimate the size of each effect.
4. Choose the three most promising effects.
5. The scientist chooses.
6. Refine the model.
7. Ask a subagent to red team the refinement.
8. Run the joint fits with the new model.
9. Report the improvement to the scientist.
10. Present three ways to improve the design of the experiment.

At step 4 it proposes the three main improvements of the model, each with the non-ideality that the statistical tools read, the mathematical need that answers it with its gain on the twin, and its physical ground, a term, an effect or a refined approximation, and at step 5 one of them is chosen as `DIR/agent.toml` sets. Each turn leaves its page in `DIR/agent/iteration_K.md`, for a reader who did not watch it: its figure, what the statistical tools read, the three improvements, which was chosen, by whom and why over the others, on physical ground and on the statistics, and what its fit showed, and `DIR/agent/report.md` joins the turns. Every term that it adds carries its reason in `why`: the physics that could leave the structure, offered as a candidate, which the loop and its twin then bear out or refuse. The package names no cause, and the agent names one only as a candidate.

## 13. More examples, of other fields

`python -m twinfit guide examples` prints the gallery, `twinfit/EXAMPLES.md`: projects of other fields, each written by `python -m twinfit example NAME DIR` with a truth that its folder does not hold. They hold a pedestal under the lines, a satellite line, a width that saturates with the power, a noise correlated by a filter, a Fano profile, photon counts, an atomic clock whose center moves with the power of its probe, a band of electrophoresis with a tail, the transit of an exoplanet over a star darker at its limb, the fringes of an interferometer whose scan is bent, and the daily returns of a stock index with fat tails.

- For each it gives what the first report put first, what to do and what not to do, and what the loop then gave beside the truth.
- Find the one whose first report leads as yours does. The two worked projects of this page, nuclear magnetic resonance and chromatography, are not in it.
- Each is shown too at three times its noise, with none drawn, and with a second kind of noise, correlated over neighboring samples or that of a counter: `python -m twinfit example clock DIR --noise 3 --kind correlated` writes one so, to see how the first report and the outcome move with the noise of an example. The page says what the noise changes in its tables.

One of them holds absorption lines of a gas cell at five pressures, with a line that its first model does not name. A run of the first model alone, `--n 0`, gives its report:

```bash
python -m twinfit example satellite my_satellite
```

```
written into my_satellite, its noise white: description.toml, model.py, an empty terms/ and the curves in traces/. Its truth and what the loop gives are in twinfit/EXAMPLES.md. Next: python -m twinfit loop my_satellite
```

```bash
python -m twinfit loop my_satellite --n 0
```

```
model 0: the first model is fitted, and its twin and its report are read on 30 synthetic sets
best.md: the five parts of each value of the best model, on 30 synthetic sets, and the combined estimate of the families
kept: nothing. Refused: nothing. The best model: model 0.
The report of each model and the output are in my_satellite/loop. An own value moves from curve to curve, which no question of the report reads: best.md says which, and how to write it.
The visual log: my_satellite/loop/log.pdf, 15 slides, each phase of each of 1 turn first as one plot, then the equations, the effects of the terms, the curves in a row, the report and the criteria, the pictures in my_satellite/loop/log.
```

`my_satellite/loop/model_0.md`:

```markdown
The chi functions in windows of position read most the chi function: the residual to the power 1.5, more of it than elsewhere, at +1.25 widths from the centre of each curve, a width being its full width at half height, under a window of 0.25 widths, about each centre: look there first.
The chi functions in bands of the residual's value read most the chi function: the residual to the power 6, more of it than a normal law gives, in the band above 2 of the residuals in their robust scale: a law of the noise wrong there, or a term that the model lacks, shapes that band.
The chi functions in cells of position and of the curve's level read most the chi function: the residual to the power 1.5, more of it than elsewhere, where the curve stands near its median, in its wings, about +1 widths from the centre of each curve, under windows of 0.5 widths and 0.125 of the height: a term that acts at that height of the curve, or a law of the noise that follows it.
The chi functions in cells of position and of the residual's value read most the chi function: the residual to the power 6, more of it than a normal law gives, for residuals near +3 in their robust scale about +1 widths from the centre of each curve, under windows of 1 widths and 1 of that scale: a law of the noise wrong in that place, or a term that the model lacks there.
```

The chi functions in windows of position read each cost of the residuals under a smooth window that moves along each curve, and say where they read most: +1.25 widths from each centre, the figure of the report, on the same side of every curve. That place is the one to read first, and a window that reads most in one place on every curve points at a component of the curve there, `second(distance)` of section 10 with its distance in the unit of the abscissa, that many widths times the full width of the curve at half its height, which the page of the proposals gives, or at a structure of that place and scale. The bands of the residual's value say how the residuals stand against a normal law of the noise, here in the tail above 2 of their robust scale, which a misfit of many times the noise fills. `twinfit/EXAMPLES.md` gives this project's truth, a second line above each, which its folder does not hold.

## 14. Limits of the package

- It names a place, a kind and a setting of what is missing, never its cause.
- The twin draws its synthetic curves at the noise that the description gives, so it takes that noise as right. The check of the residuals against the height of the fitted curve, on the first page of each report, is what says whether it is.
- With 30 synthetic sets, the default, an answer can stand beyond all of them and no further: its chance is then the smallest that 30 sets give. More sets, `python -m twinfit loop DIR --sets 100`, read smaller levels and take longer.
- On models that lack nothing, at the default of 30 sets with a seed for each draw, a run of 0.1.2, before the chi functions were read, read some answer past its line on 0.068 ± 0.016 of 250 draws, as `README.md` of the package says, the count of 0.1.3 to come, where the line is set for 0.05 by its construction. Read one answer just past its line as weak evidence, and read it again with `--sets 100`.
- A term much smaller than the noise can be missing while no answer stands beyond the twin. When one answer alone stands just past its line, the loop says so and stops: a model that lacks nothing shows one such answer on some of its synthetic sets.
- Glitches and a noise that grows with the signal are kinds of the report that no term of the mean fills. The report says so when one leads, and how to give the noise of each sample instead.
- A proposal for the experiment is read on the synthetic sets of the twin, so its ratio moves from run to run: the page says within how much, about one over the root of the number of sets.
- A term in the wrong frame leaves a misfit of its own. A ripple read about each center where it is fixed to the abscissa changes from curve to curve, and the report then points along a setting. Read a new answer along a setting against the frame of the terms kept before it.
- Two forms of one effect, as a Lorentzian added to a line and a Voigt line in its place, are not set side by side by one run of the loop, which fits each term file on top of the ones before it. Compare them in two folders, each with one form, by their criterion and their twin, as section 11 does: there the tail and the impurity lie 0.39 apart in criterion, and their widths 5.3 bars of the fit apart.
- One misfit can light several answers at once: a slope of the baseline lights a slow trend over the span of each curve and a periodic structure along the abscissa together. The report says when kinds tie, and then their order means little: try each, and keep the one that the twin and the known parameters bear out.
- The bars of the own parameters of each curve count what the shared ones leave uncertain, but no twin reads them, as it reads the shared ones.
- `along` multiplies a parameter by one plus its size times the setting, so it cannot move a parameter that sits near zero, as a place read less the offset of its abscissa can. A place that moves with a setting is written in `model.py`, as the model of section 9 writes the shift for each kelvin, and as the gallery's clock writes the shift of its line along the power.
- `second` takes one distance and copies the curve there, at its own width: a second line whose distance or share follows a setting is written in `model.py`, as a function of the settings, as a place is.
- A lever along a setting reads the model at settings past those of the records, and holds as far as the model holds there. A lever of samples draws the twin's noise white, as independent samples, unless `correlated = true` under `[noise]` draws it with the correlation of the residuals, kept over the same span of the abscissa on closer samples: under a noise correlated over N samples, twice the samples is less than twice the information.
- The loop can end in words that the runs above do not show. Where a noise correlated from sample to sample leads the last report, a term file cannot hold it, and the report of that model says what to do in a copy of the folder. Where the check of the noise of the best model fails, what was kept is not to be believed before that is read. Where the twin does not admit the fit of the best model, its bars are not to be read, and where a statistic misses its band by little, `--sets 500` reads it again.
