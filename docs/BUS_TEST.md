# The bus test

The bus test asks whether another researcher can continue this work from what PyPI serves for twinfit and this repository alone, with no one to ask. It was run 7 times on 0.1.0 and on the pages of 0.1.0.post1, each time by a fresh agent in a sandbox that held the wheel, the source distribution, the README that is the package's page on PyPI and `FRESH.md`, a page saying what the folder holds (its text follows the brief in `docs/BUS_TEST_BRIEF.md`), and from the second run on this repository as a folder too (run 1 had the package alone and blocked on that). The same test is written here for a person and for an agent, so that anyone can run it on a later version and compare.

## The six steps

The researcher, a person or an agent, writes a report under six headings, each with what was found, the file and the section of the documentation where it was found, the minutes it took, and a mark out of 10 as its last line, `Mark: N/10`.

1. What the package is for and why it was made, and what role the author's own line plays.
2. What he measured, and the results that stand, each with its uncertainty and the conditions under which it holds.
3. What is validated, and how, and what is not.
4. What he left open: the questions, the decisions, and the next steps that he planned.
5. Reproduce his results by the commands that the documentation gives, at least one share of the validation and one draw of a count, and compare.
6. Take the first of his next steps, in a folder of your own, as far as you can, and say what you would need to go further.

A mark of 10 means the step was done fully from the documentation alone, without a guess, and nothing in it blocked, misled or cost time. Each stumble is a line of a list that begins with its rating in bold, with the exact words of the page or of the output that caused it: `- **BLOCKS**` (the work cannot continue without the author), `- **MISLEADS**` (a statement that is wrong, or that contradicts another), `- **SLOWS**` (time lost) or `- **POLISH**`. The report ends with a verdict in one sentence: whether a researcher can continue the work from PyPI and this repository alone, the total out of 60, and the three things to add first.

## The rule of a pass

`python tools/bus_score.py REPORT.md` reads the marks and the stumbles. A run passes at 60 of 60 with no stumble that blocks or misleads and a verdict that says yes. Two runs in a row that pass end the test. Two runs in a row whose stumbles do not fall end it with what stands, and the author answers or repairs the stumbles of the last run.

## For a person

Take a machine with Python 3.14, which the bench's checks of step 5 need (the package itself runs on 3.12 or later), and no copy of this repository. Install the package, `pip install "twinfit[log]"`, read its page on PyPI, clone this repository, and work through the six steps with the documentation alone: the README, `docs/`, `toy/README.md`, `python -m twinfit guide`. Write the report as the steps say. Nothing of the author's recorded traces is needed.

## For an agent

Give an agent the text of `docs/BUS_TEST_BRIEF.md` as its only instruction, in a folder that holds the wheel, the source distribution, the README of the package, `FRESH.md` and this repository as a folder, with the network, which the README's install route without a Python needs, and no file of the author's. The author's 7 runs used Sonnet 5.5 at high effort in a sandbox of that kind, the first with a shorter text, each between 16 and 62 minutes.

## The 7 runs on record

Scored by `tools/bus_score.py` on the reports kept by the author, the minutes from each run's record. The first run wrote no marks, so it scores 0.

| run | marks of 60 | blocks | misleads | slows | polish | minutes |
|---|---|---|---|---|---|---|
| 1 | 0 | 2 | 1 | 2 | 9 | 16 |
| 2 | 52 | 0 | 2 | 3 | 11 | 29 |
| 3 | 47 | 1 | 3 | 3 | 10 | 40 |
| 4 | 53 | 0 | 0 | 7 | 9 | 29 |
| 5 | 47 | 0 | 1 | 5 | 8 | 34 |
| 6 | 46 | 0 | 3 | 8 | 8 | 37 |
| 7 | 46 | 1 | 2 | 7 | 6 | 62 |

The runs read the pages as they stood each time, and the pages were written again after each. No run reached 60. Runs 1, 3 and 7 each hit a block, every one repaired before 0.1.0.post2. The misleads the runs found were repaired in the pages or carried to 0.1.1, as the author's records say for each, and those of the first three runs are answered there one by one. The slows and the polish were repaired in part, and the rest is owed. A later version is tested the same way, and its runs go in this table.
