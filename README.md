# twinfit

A tool for the analysis of the residuals of a fit. It reports where a fitted model misses, what the structure of the miss looks like, and how the miss biases the fitted parameters.

twinfit is a Python package for a researcher who fits one model to a set of measured curves. It fits the model to all the curves at once and reads the residuals through families of statistics. Each statistic is compared with its value on a twin, that is, on synthetic curves drawn from the fitted model with the law of noise assigned to the records, fitted and read in the same way. The report states where the model misses, of which kind the miss is, and along which setting of the experiment it varies. The package names no term of any field. The researcher writes the term, and the package keeps it only where the records support it. In the terms of statistics, the twin is a parametric bootstrap from the fitted model, and the report a test of the specification by weighted moments of the residuals calibrated on it: section 1 of the tutorial, [`twinfit/TUTORIAL.md`](twinfit/TUTORIAL.md), sets the words of the package that the literature of statistics names otherwise beside their names there.

This branch, `next`, holds 0.1.3, ahead of the release on PyPI. To run the code of this branch:

```bash
pip install "git+https://github.com/MichelangeloDondi/twinfit@next"
```

The release on PyPI, older than this branch, installs with `pip install twinfit`.

The page of the package, with a first run of five minutes, is [`twinfit/README.md`](twinfit/README.md). It is also the page of the package on PyPI.

The same five minutes as a notebook, [`notebooks/five_minutes.ipynb`](notebooks/five_minutes.ipynb), whose first cell installs the package with pip and whose commands run through the shell, checked on Python 3.12 and 3.14 and not yet in a hosted notebook such as Colab, which opens it at https://colab.research.google.com/github/MichelangeloDondi/twinfit/blob/main/notebooks/five_minutes.ipynb. To cite the package, [`CITATION.cff`](CITATION.cff), which GitHub shows as "Cite this repository".

## Contents of this repository

| folder | contents |
|---|---|
| [`twinfit/`](twinfit/) | the package as PyPI ships it, with its tutorial, its gallery of other fields and the card of an agent |
| [`toy/`](toy/README.md) | the bench of the package: a model of five terms of the author's own line, the two photon line of rubidium from 5S to 6S, from which synthetic traces of known truth are drawn, the programs that count what the package does on them, and their tables |
| [`runs/`](runs/) | the kept row of each draw of three of the counts behind the shares of the README (the pointing, with the place, the setting and the right models, the gallery, and power), from which a table is counted again and a draw made again |
| [`tools/met.py`](tools/met.py) | the reader of every criterion of `toy/thresholds.toml` on the tables |
| [`docs/`](docs/) | [`VALIDATION.md`](docs/VALIDATION.md), what the package was measured to do, [`ROADMAP.md`](docs/ROADMAP.md), the versions that follow, [`DECISIONS.md`](docs/DECISIONS.md), the decisions taken and the questions left open, and [`BUS_TEST.md`](docs/BUS_TEST.md), whether another researcher, a person or an agent, continues the work from PyPI and this repository alone |
| [`data/noise/`](data/noise/) | the law of the noise of each recorded trace of the author, read from the traces, without any sample of them |

The recorded traces of the author are not in this repository. The package does not need them, and the bench draws its traces from the model of five terms.

## Taking the work over

[`twinfit/CONTINUING.md`](twinfit/CONTINUING.md) is written for a researcher who continues this work. It states where each part is, what to run in the first hour, and where each line of work starts. In an installed package it stands beside the README, in the folder that `python -m twinfit guide readme` names.

## Citation

M. Dondi, twinfit, the version that `python -m twinfit --version` prints, https://github.com/MichelangeloDondi/twinfit.

## Licence

MIT, [`LICENSE`](LICENSE).
