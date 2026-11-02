# Continuing this work

For a researcher who takes twinfit over from its author: what the package is, where its parts are, what was measured, what to run in the first hour, and where each line of work starts. Section 1 of `TUTORIAL.md` defines the terms.

## The package and its example

twinfit is a general purpose tool. You fit one model to several curves, spectra at several temperatures for example. twinfit reads the residuals through families of statistics and compares each with its value on a twin, synthetic curves drawn from the fitted model with the noise law of the records, fitted and read the same way. It then says where the model misses, of which kind the miss is and along which setting it varies. The physics comes after: you write the term.

The worked example, and the bench the package was measured on, is a line of the author's own experiment. A model of that line in five terms draws synthetic curves whose truth is known, and the README's shares were counted on them and on eleven synthetic projects of other fields. Of the five values, two were fitted to the author's recorded curves, a constant width, 0.037 ± 0.015, and the spread of the frequency noise of the light source, 0.369 ± 0.028, both in megahertz, two are held at theory with their uncertainties, and one comes from the literature. The recorded curves are private and the tool needs none of them. `toy/README.md` gives the model, its terms and the source of each value, and the `source` fields of `toy/config.toml`, `toy/const.toml` and `toy/twin.toml` give each value's condition, fit and uncertainty.

## Where each part is

PyPI carries the package and its pages: this one, the README, `TUTORIAL.md`, `EXAMPLES.md`, `RUBIDIUM.md` (the author's case, as plots), the agent's card in `agent/`, and `version.py`, the list of versions with their digests and changes. The bench is in the repository alone.

The repository, https://github.com/MichelangeloDondi/twinfit, carries the rest.

| folder | contents | use by a successor |
|---|---|---|
| `twinfit/` | the package | use and change |
| `toy/` | the model of five terms, the programs of the counts, their tables, the thresholds of the criteria, `again.py` and `count.py` | checks, one draw made again, counts made again |
| `runs/` | the kept row of each draw of the counts of pointing and the gallery, and the twin of each of their cases (power's rows, made at 0.1.0, do not leave with this release) | reading, recount of a table |
| `tools/met.py` | the reader of every criterion of `toy/thresholds.toml` on the tables | reading of the criteria |
| `docs/` | `VALIDATION.md`, `ROADMAP.md`, `DECISIONS.md` | reading before any change |
| `data/noise/` | the law of the noise of each recorded curve, one for each peak and condition: the variance of a sample, a² + b h + c h² + d s², in the height h of the curve at that point and its slope s. No sample of any curve | the source of the noise of the model's draws |

The recorded curves and the author's notes are not public. Comments in the code cite those notes by number, nothing in them is needed to run anything, and `docs/DECISIONS.md` gives the reasons that matter.

## What was measured

`docs/VALIDATION.md`, written from the tables by `docs/validation.py`, lists every criterion. Three are not met at this release.

| criterion | share | threshold |
|---|---|---|
| C1, the place of a width that grows steeply along a setting, put first or tied for first where the width is hidden at the size of the model | 0.440 ± 0.050 of 100 draws | at least 0.60 |
| C5c, a false alarm on a right model at the default of 30 sets | 0.188 ± 0.025 of 250 draws | at most 0.12 |
| C7, one bar holds the truth for the width that grows steeply along a setting under the toy's law of the noise | 0.820 ± 0.038 of 100 draws | 0.56 to 0.81 |

The others that are read are met, and C10, C2, C2b, C5b are not read at this release, their tables made at 0.1.0: among those met, the place of a hidden term pooled over four terms, 0.780 ± 0.021 of 400 draws, and the eleven projects of other fields, 0.761 ± 0.020 of 440 draws. Not validated yet: the agentic version, use by anyone but the author, and terms beyond the fifth of the model.

## The first hour

The package runs on Python 3.12 or later. The bench needs Python 3.14: the digest of a count hashes the syntax tree of its code, which differs between versions of Python, so on another version each check but the gallery's reports its table as made from other files.

1. Install the package with matplotlib, which the tutorial and the figures need: `pip install "twinfit[log]"`.
2. Run the five minutes of the README, the clock at five probe powers.
3. Run every command of the tutorial again, on a quiet machine, since a loaded one can make the check refuse. It takes one to five minutes:

```bash
python -m twinfit.tutorial --check
```

It ends with "the tutorial: every output of TUTORIAL.md is what its command prints now".

4. Clone the repository and run its checks, on Python 3.14:

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

The first three commands compare the digest of their code with that of their table and sum the table's rows again: `pointing` and `gallery_count` report their tables as counted "from the files as they are", a pass meaning the table matches the code, and `power` reports its table as made from other files, made at 0.1.0 and not counted again at this release. The fourth counts the kept fits of each model on its 100 draws, makes three again, and ends with "every model passes, and three fits of each were made again". `tools/met.py` prints the readings of `docs/VALIDATION.md`, each label naming the criteria whose `row` field in `toy/thresholds.toml` holds it. `toy/again.py` makes draw 7 of that case again from its seed, in seconds on the author's machine and in a few minutes on a slower one (the eighth bus test measured 167 and 280 seconds for two draws on four cloud cores), and ends with "the kept row, field for field".

## Three lines of work

- **Use.** The README, `TUTORIAL.md`, and, for an agent, the card of `agent/AGENT.md`.
- **Test.** `toy/README.md` states how each count is made, by `python toy/count.py PROGRAM`, what it costs, and how a draw is made again. `docs/BUS_TEST.md` is the test of this page itself: whether a researcher, a person or an agent, continues the work from PyPI and the repository alone, with its six steps, its rule and the 8 runs on record.
- **Development.** `docs/ROADMAP.md` gives the versions that follow and a first step that takes an hour. `docs/DECISIONS.md` gives the decisions taken and the questions left open. Together those two pages list what is open.

## Three rules of the work

1. A change of a file of `twinfit/` is a new version. `twinfit.version.check()` raises until `version.py` holds the new number, its digest and what it changed. A post release, as 0.1.0.post2 was, changes pages and one sentence of the tutorial's program, and keeps the tables of its release, the gallery's counted again on the same seeds for its digest. 0.1.1 changes code, so the tables that its criteria read are counted again, `toy/pointing.json`, `toy/chosen.json` and `toy/gallery.json`, and with them the bench's own three, `toy/matrix.json`, `toy/convergence.json` and `toy/broad.json`, which no criterion reads, while `toy/scoreboard.json`, `toy/power.json` and `toy/benchmark.json` stand at 0.1.0 and their criteria are not read, as the section on what was measured says.
2. A count is made again after a change of the code that it runs. Its `--check` says when.
3. A criterion is read with `tools/met.py`, against thresholds that the author drafted before the counts. The header of `toy/thresholds.toml` and `docs/DECISIONS.md` say which were set again after tables were read.
