# twinfit

A tool for the analysis of the residuals of a fit. It reports where a fitted model misses, what the structure of the miss looks like, and how the miss biases the fitted parameters.

twinfit is a Python package for a researcher who fits one model to a set of measured curves. It fits the model to all the curves at once and reads the residuals through families of statistics. Each statistic is compared with its value on a twin, that is, on synthetic curves drawn from the fitted model with the law of noise assigned to the records, fitted and read in the same way. The report states where the model misses, of which kind the miss is, and along which setting of the experiment it varies. The package names no term of any field. The researcher writes the term, and the package keeps it only where the records support it.

```bash
pip install twinfit
```

The page of the package, with a first run of five minutes, is [`twinfit/README.md`](twinfit/README.md). It is also the page of the package on PyPI.

## Contents of this repository

| folder | contents |
|---|---|
| [`twinfit/`](twinfit/) | the package as PyPI ships it, with its tutorial, its gallery of other fields and the card of an agent |
| [`toy/`](toy/README.md) | the bench of the package: a model of five terms of the author's own line, the two photon line of rubidium from 5S to 6S, from which synthetic traces of known truth are drawn, the programs that count what the package does on them, and their tables |
| [`runs/`](runs/) | the kept row of each draw of three of the counts behind the shares of the README (the place and the setting, the right models, and the gallery), from which a table is counted again and a draw made again |
| [`tools/met.py`](tools/met.py) | the reader of every criterion of `toy/thresholds.toml` on the tables |
| [`docs/`](docs/) | [`VALIDATION.md`](docs/VALIDATION.md), what the package was measured to do, [`ROADMAP.md`](docs/ROADMAP.md), the versions that follow, and [`DECISIONS.md`](docs/DECISIONS.md), the decisions taken and the questions left open |
| [`data/noise/`](data/noise/) | the law of the noise of each recorded trace of the author, read from the traces, without any sample of them |

The recorded traces of the author are not in this repository. The package does not need them, and the bench draws its traces from the model of five terms.

## Taking the work over

[`twinfit/CONTINUING.md`](twinfit/CONTINUING.md) is written for a researcher who continues this work. It states where each part is, what to run in the first hour, and where each line of work starts. In an installed package it stands beside the README, in the folder that `python -m twinfit guide readme` names.

## Licence

MIT, [`LICENSE`](LICENSE).
