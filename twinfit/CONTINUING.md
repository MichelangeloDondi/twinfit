# Continuing this work

This page is for a researcher who takes twinfit over from its author: what it is, where each part is, what was measured, what to run in the first hour, and where each line of the work starts. It builds on the README of the package. Section 1 of `TUTORIAL.md` defines every word of the package.

## What twinfit is, and what its example is

twinfit is a general purpose tool. A researcher fits one model to many curves, as spectra at several temperatures. twinfit reads what the fit leaves through families of statistics, sets each against a twin, synthetic curves drawn from the fitted model with the law of noise that the records are given and read in the same way, and says where the model misses, of which kind and along which setting. It never names a cause: the researcher writes the term.

Its worked example, and the bench on which it was measured, is a line of the author's own experiment. A model of that line in five terms, the toy, draws synthetic curves whose truth is known, and the shares of the README on it were counted on them, beside a gallery of eleven synthetic projects of other fields. Two of its values were fitted on the author's recorded curves, two are held at their theory with their bars, and one follows from the literature. The recorded curves are private, and the tool needs none of them. The toy, its terms, their values and the source of each are in `toy/README.md` of the repository. The condition of each value, the fit that gave it and its uncertainty stand in the comments of `toy/config.toml`, `toy/const.toml` and `toy/twin.toml`.

## Where each part is

From PyPI, the wheel and the source distribution carry the package and its pages: this one, the README, `TUTORIAL.md`, `EXAMPLES.md` (eleven projects of other fields), `RUBIDIUM.md` (the author's own case, as plots), the card of an agent in `agent/`, and `version.py`, which holds every version with its digest and its changes. They do not carry the bench.

The repository, https://github.com/MichelangeloDondi/twinfit, carries the rest.

| folder | what it holds | what a successor does with it |
|---|---|---|
| `twinfit/` | the package | uses it, changes it |
| `toy/` | the model of five terms, the programs of the counts, their tables, the thresholds of the criteria, `again.py` and `count.py` | checks, makes a draw again, counts again |
| `runs/` | the kept row of each draw of the counts of pointing, power and the gallery, and the twin of each of their cases | reads them, recounts a table from them |
| `tools/met.py` | the reader of every criterion of `toy/thresholds.toml` on the tables | reads the criteria |
| `docs/` | `VALIDATION.md`, `ROADMAP.md`, `DECISIONS.md` | reads them before changing anything |
| `data/noise/` | the law of the noise of each recorded curve, one for each peak and condition: the variance of a sample, a² + b h + c h² + d s², in the height h of the curve there and its slope s. No sample of a curve | the toy draws its noise from it |

Not public: the recorded curves, and the author's working records. The comments of the code name entries of those records, "his message N", "review N", "row MNNN", "candidate N", "decision N", "lesson N", and paths such as `PLAN.md`, `reviews/` and `runs/reports/`. None of them is needed to run anything. `docs/DECISIONS.md` gives the reasons that matter.

## What was measured

`docs/VALIDATION.md` holds every criterion, written from the tables by `docs/validation.py`. Two are not met at this release:

| criterion | share | threshold |
|---|---|---|
| C1, the place of a width that grows steeply along a setting, put first where it is hidden at the size of the twin | 0.440 ± 0.050 of 100 draws | at least 0.60 |
| C5c, a false alarm on a right model at the default of 30 sets | 0.188 ± 0.025 of 250 draws | at most 0.12 |

The others are met, among them the place of a hidden term pooled over four terms, 0.780 ± 0.021 of 400, and eleven projects of other fields, 0.759 ± 0.029 of 220. Not validated at all: no person but the author has used the package on curves of their own, the agentic version has no count of its own, and no term beyond the fifth exists in the toy.

## The first hour

The package runs on Python 3.12 or later. The bench needs Python 3.14: a count's digest hashes the syntax tree of the code, which differs between versions of Python, so on another one each check says that its table was made from other files.

1. Install, as the README says, with matplotlib, which the tutorial and the figures need: `pip install "twinfit[log]"`.
2. The five minutes of the README, the clock at five probe powers.
3. Every command of the tutorial run again, one to five minutes:

```bash
python -m twinfit.tutorial --check
```

It ends "the tutorial: every output of TUTORIAL.md is what its command prints now".

4. The repository, on Python 3.14:

```bash
git clone https://github.com/MichelangeloDondi/twinfit
cd twinfit
python toy/pointing.py --check
python toy/power.py --check
python toy/gallery_count.py --check
python toy/draws.py --check
python tools/met.py M456,M505,M457,M458,M441,M444,M459 --counted
python toy/again.py "2 hidden under 1" 7
```

Each of the first three compares the digest of its code with that of its table, sums again the rows that the table holds, and says that its table was counted, or made, "from the files as they are": a pass says that the table matches the code. The kept rows of `runs/` are read by `toy/again.py` and by `toy/count.py`. The fourth counts the kept fits of each model of the toy on its 100 synthetic draws, makes three of them again, and ends "every model passes, and three fits of each were made again". `tools/met.py` prints the readings of `docs/VALIDATION.md`. `toy/again.py` makes draw 7 of that case again from its seed, in seconds, and ends "the kept row, field for field".

## Three lines of work

- **Use it** on curves of one's own: the README, `TUTORIAL.md`, and for an agent the card of `agent/AGENT.md`.
- **Test it again**: `toy/README.md` says how each count is made, by `python toy/count.py PROGRAM`, what it costs, and how a draw is made again.
- **Develop it**: `docs/ROADMAP.md` gives the versions that follow and a first step that takes an hour, `docs/DECISIONS.md` the decisions and the questions left open.

## Three rules of the work

1. A change of a file of `twinfit/` is a new version: `twinfit.version.check()` raises until `version.py` holds its number, its digest and what it changed. A post release, as this one, changes pages alone and keeps the tables of its release.
2. A count is made again after a change of the code that it runs, and its `--check` says when.
3. A criterion is read with `tools/met.py`, against thresholds that the author drafted before the counts. The header of `toy/thresholds.toml` says which were set again after tables were read.
