# 0.1.1, what it holds, and what comes after it

This page builds on [`CONTINUING.md`](../twinfit/CONTINUING.md), [`VALIDATION.md`](VALIDATION.md) and [`DECISIONS.md`](DECISIONS.md). Every command on this page is run from the root of the repository, on Python 3.14.6. Three words are used throughout. A **draw** is one synthetic data set of the model of `toy/`. The **sets** of a draw, 30 by default, are the synthetic copies drawn from the fitted model, that is, its twin. The **line** of the report is the largest excess, over all the answers together, that a set of a right model, one that lacks nothing, reaches on 95 of 100 sets (`twinfit/TUTORIAL.md`). The line of a single statistic, in section 3.1, is its own 95th percentile over the sets.

## 1. Where 0.1.1 stands

The code is 0.1.1rc1's with three repairs, the thinning writing inside its copy, a kind an example does not draw refused before its folder is made, and Mechanism's second field named family in the card and the tutorial, and the number moved, the pages are written again from the tables of the six counts made on it. `python -m twinfit --version` prints the number of this release, and `twinfit.version.check()` returns it, because its digest, which covers the pages as well as the code, is stamped in `twinfit/version.py`.

Three criteria are not met (`docs/VALIDATION.md`): C7, the bar of the collisions under the toy's law of the noise, 0.820 ± 0.038 of 100 draws against 0.56 to 0.81, the bar too wide, and the two below. C1 on the collisions, the share of draws on which the report puts their place first or tied for first, is 0.440 ± 0.050 of 100 draws against a threshold of at least 0.60. C5c, the share of false alarms, that is, of answers past their line on a right model, at 30 sets, is 0.188 ± 0.025 of 250 draws against a threshold of at most 0.12.

## 2. The versions proposed

The order serves one aim: that a researcher in another field, with curves of their own and no contact with the author, can use twinfit in an afternoon and trust what it tells them. Each version closes on criteria a user can read, not on the author's counts alone.

| version | what a user gets | closes when |
|---|---|---|
| 0.1.1 | pages written for a user, the small term's count, `toy/broad.py` and its table (section 3.1), its statistic proposed for the report of a later version, the repairs of section 3.2, the twin of the real case with one laser noise, `guide continuing`, a tutorial check that names the line that differs | the readers of the pages find no sentence the code refutes, C5c stated as not met at 30 sets, the counts of section 5 made again |
| 0.2.0, the adoption release | real data without a toy: a model written as a formula where one line suffices, any number of curves, the twin drawn with the noise a real detector gives (correlated samples, counted events), a loop that answers in seconds on ten curves, an install that works on Windows | the clock project of `twinfit/EXAMPLES.md` with twenty curves in under a minute on a laptop, a criterion of C7's bars held on correlated and on counted noise, two projects from public data of other fields read end to end from their files alone |
| 0.3.0 | what can be said of the shared parameters before the missing term is named, as a bound within a declared class, and an agent in the researcher's place validated by blind runs | C12, how often the bound holds the truth, written and met, and the agent's blind runs scored on draws of known truth |
| 1.0 | files, commands and pages fixed, a paper, a site of the documentation | a paper submitted, three cases of other researchers read |

The author red teams this order before 0.1.1 is stamped.

## 3. What 0.1.1 holds

### 3.1 The small term

The author has counted, for 0.1.1 and with `toy/broad.py`, a program of the bench that ships with 0.1.1 beside its table `toy/broad.json`, a statistic for a broadening that is small against the noise. It is computed from `twinfit.linear.Linear` of a fit, whose `J` and `r` are the whitened derivatives of the mean and the whitened residuals with the own parameters of each record taken out. For a shared parameter $p$, column $j$ of `J`, and record $k$, with $a_k$ the rows of record $k$ in that column and $r_k$ its residuals:

$$w_k = a_k \cdot a_k, \qquad e_k = \frac{a_k \cdot r_k}{w_k},$$

the excess of $p$ that record $k$ asks on its own. Records with $w_k = 0$ are left out, and fewer than three leave $p$ unread. Along a numeric setting $s$ of the design, with $z_k = s_k - \bar{s}$ over the records kept:

$$S_{p,s} = \frac{\left|\sum_k z_k w_k e_k\right|}{\sqrt{\sum_k z_k^2 w_k}}, \qquad C_p = \sum_k w_k \left(e_k - \bar{e}\right)^2, \quad \bar{e} = \frac{\sum_k w_k e_k}{\sum_k w_k}.$$

$S_{p,s}$ is a standard normal and $C_p$ a chi square of one degree of freedom less than the records, at first order where $p$ is one value across them. The line of each is its 95th percentile over the twin's sets, 20 at least, as the report's own lines are. $S$ along the temperature and $C$ of the width each stood past its line on 100 of 100 draws of the collisions, and either of them or the report stood past its line on 0.300 ± 0.032 of 200 right models, above the 0.12 of C5c. Its scoring under C1, the place put first or tied for first, is not written: its count enters 0.1.1 as `toy/broad.py` and its table, the statistic proposed for the report of a later version and not as a repair of C1.

### 3.2 Repairs from the blind red teams

Each of these was done at 0.1.1rc1 and stands in the code: the list is kept as the record of what the blind runs asked.

- `along` (`twinfit/mechanisms.py`) takes no `lower` and no `upper`, so its factor 1 + a z can pass zero and a width falling with power was kept: add both, and have `best.md` say whether a parameter rises or falls.
- `best.md` and `stop_sentence` (`twinfit/loop.py`) name all four ends of the loop, but the stop rule of `twinfit/agent/AGENT.md` names two: add the noise first and no term worth keeping.
- `scan` warns of uneven gaps in its fine branch alone: warn in the coarse branch too.
- `python -m twinfit guide continuing` is refused, and `twinfit/CONTINUING.md` is reached beside the file that `guide readme` names: add it to the guides.
- `python -m twinfit.tutorial --check` says that the outputs differ and not which: name the first line that differs.
- The example of the real case, `python -m twinfit example rubidium DIR`, is drawn with a laser noise for each trace where the model holds one value, which a blind run read as an excess along the power: draw it again with the one value, the author's decision, and count the gallery again, since it hashes the example's files.
- The reading of the recorded traces with `--sets 100`, which RUBIDIUM.md says a researcher makes before a term is written, is not in this release: make it, and say what it gave.

### 3.3 The default number of sets, 30 or 100

On the same draws `runs/pointing_own_*` (30 sets, each draw against its own twin) and `runs/pointing_own100_*` (100 sets, made at 0.1.0 by no program here and not shipped at this release, which `docs/validation.py` reads as not counted) gave a false alarm of 0.188 ± 0.025 at 30 sets and of 0.076 ± 0.017 at 100 sets, of 250 right models each, at 0.1.0. The second lies under the 0.09 that `toy/thresholds.toml` derives for 100 sets. The share of the collisions put first or tied on draws 1 to 50 is 0.520 ± 0.071 of 50 at both. Therefore 100 sets meet C5c and leave C1 unchanged, at more time a draw: the queue block of `toy/pointing.json` holds the waves of the cases and of the own twin at 30 sets, none at 100, with their lane, their seconds and the draws made, and a reader computes the factor from them, since the waves ran under different loads. The default stays at 30 and C5c is stated as not met.

The default stands in `twinfit/loop.py`, `twinfit/__main__.py` and, as the counts' own 30, `OWN_SETS` of `toy/pointing.py` and `SETS` of `toy/chosen.py` and `toy/gallery_count.py`. Only `OWN_SETS` moves C5c, and moving it makes the units of an own twin at 100 sets owed, which `python toy/pointing.py units` lists.

Left in 0.1.1 and declared here: the report's advice on a correlated noise, written where the autocorrelation of the residuals falls under 0.05 at N samples and one sample in every N keeps enough of the shortest curve to fit, says that `python -m twinfit thin DIR N` does it "in three copies", where `thin` writes up to three, fewer where N is under 3, as its usage and the README say. The sentence is a string of `twinfit/loop.py`, which the digest of the chosen cases' count hashes, and its change would ask that count again, 39 core hours by the queue's ledger: it is paid at the next version with the count.

## 4. The first step

The table printed by `python -m twinfit terms` lacks the entry `exponential_width`. The first step is a test that catches this, `tests/test_terms_lists.py`, in a new folder. It runs with plain Python, and with pytest where pytest is installed:

```python
import sys

sys.path.insert(0, ".")          # the repository's package, run from its root, and not an installed one

from twinfit.__main__ import NAMED, terms  # noqa: E402


def test_every_named_term_is_in_the_table():
    table = terms().split("\n\n")[1]
    for call, _ in NAMED:
        assert call.split("(")[0] + "(" in table, call


if __name__ == "__main__":
    test_every_named_term_is_in_the_table()
    print("1 passed")
```

Before the edit, `python tests/test_terms_lists.py` fails with `AssertionError: exponential_width()`. With this edit, delete the README's sentence on `python -m twinfit terms` that says the table lacks `exponential_width`. In the docstring of `twinfit/mechanisms.py`, between the `lorentzian_width` and `pedestal(stretch)` entries, with no blank line, insert:

```text
    exponential_width(along=None)
                             the same with the kernel exp(-|u| / b) / (2 b), b the decay length in the unit of the
                             abscissa. along scales b as for the two widths above
```

`twinfit/version.py` holds `__version__ = "0.1.1"` and its entry, whose digest the stamp of this release writes, after the entry of 0.1.1rc1 with the digest of that stamp.

```bash
python -c "import matplotlib.pyplot"   # builds the font cache once
python -m twinfit.tutorial --write
python -m twinfit.tutorial --check
python tests/test_terms_lists.py
python toy/pointing.py --check
python toy/gallery_count.py --check
```

`--write` also draws again the figures of `twinfit/examples/figures`, so a change of one line shows binary changes there, and takes one to four minutes. The `changes` of an entry is a tuple of strings, as in the entries before it. `--check` ends `the tutorial: every output of TUTORIAL.md is what its command prints now`. The test prints `1 passed`, and `pointing --check` exits 0 on Python 3.14, the Python that made the tables: on another version it reports the table as made from other files. `gallery_count --check` exits 0 on the tree of this release, its 440 draws counted again on 0.1.1's code, since its digest hashes bytes and the repairs of the code moved it. After the insertion of this section in the docstring of `twinfit/mechanisms.py` it reports `toy/gallery.json` as made from other files and exits 1, by rule 2 of section 5, until `python toy/count.py gallery_count` counts it again.

Before `--write`, since the tutorial prints the version and records it in TUTORIAL.md, move `__version__` in `twinfit/version.py` to a new number, `0.1.2rc1` for this change of the code by rule 1 of section 5, and add its entry after the last, its digest `PENDING` and its changes a tuple of strings: moved after `--write`, the number leaves `--check` refusing at the line of TUTORIAL.md that holds `0.1.1`. After `--write` and `--check`, run the stamp `python -m twinfit.version --stamp LAST`, LAST being the number stamped before the edit, 0.1.1 for a change made after this release: with the number left at 0.1.1 it is refused, `version 0.1.1 was stamped already ... give the change a new number`. Python prints a `RuntimeWarning` about `twinfit.version` before the answer, which is harmless. It writes the digest into the entry of the new number, so that `twinfit.version.check()` returns `0.1.2rc1`. It was not run for this page.

## 5. The rules of the work

1. **A change of the package is a new number**, after this release `0.1.1.postN` for a repair of its pages or `0.1.2rcN` for its code. The digest hashes every `.py`, `.md`, `.json`, `.toml`, `.txt` and `.csv` file of `twinfit/` but `version.py`. `tools/met.py` reads a table made at another number than `release` in `toy/thresholds.toml` as meeting nothing, a `.postN` table counting as its release: count every table, then move `release`, then run `python docs/validation.py`. A successor cannot close a version this way: `toy/chosen.json` needs the author's recorded traces for its single runs on them, beside its 100 draws of each arm, so at a new number C6 and C7 are not read and `tools/met.py` ends on not met. Until a version writes the rule that reads C6 and C7 from the last count that had the traces, or counts them on the synthetic project of `python -m twinfit example rubidium DIR`, the closing of a version is the author's.
2. **A count is made again after a change of the code it runs,** by `python toy/count.py PROGRAM` (`tools/autoqueue.py add` does nothing here). `gallery_count` hashes bytes, so even a comment in `twinfit/loop.py` stales it, the others ignoring docstrings and comments. `runs/` keeps no rows of `chosen` or `scoreboard`, and the 3 conditions of `chosen` on the author's recorded traces, `his`, `his law`, `his report`, read traces which are not public: C6 and C7 cannot be counted again on his traces. The loop of C6 and the bars of C7 can be run on the project that `python -m twinfit example rubidium DIR` writes, synthetic traces of the same model with their truth beside them, which is one draw and gives no share.

## 6. What is not on the roadmap

- The author's own rubidium line, on traces that are private, and the four quantities he measures with it: the self broadening coefficient, the differential polarizability of the two levels, the beam radius and the efficiency of the retroreflection.
- A claim that no other tool does this. The only comparison here is C10, 0.917 ± 0.025 of 120.
