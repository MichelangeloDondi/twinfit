# The brief of the bus test, for an agent

Given word for word as the agent's only instruction, in a folder that holds `dist/` with the wheel and the source distribution, `README.md` as PyPI shows it, `FRESH.md`, which says what the folder holds, and `repository/` as this repository stands. The first of the author's eight runs had a shorter text that asked for no marks and named no repository, and the seven after had this one.

```
You are a researcher who takes over the work of the author of the Python package twinfit, who can no longer be
reached: the bus test. His work is the package, a general purpose tool, and his own line of rubidium is its worked example
and the bench on which it was validated: his recorded traces are not needed and are not given. Everything that you have
is in this folder. This is a new Mac on which nothing is installed: no Python and no package manager, and no administrator password. Read FRESH.md first: it says what the researcher receives. Install what you need from the network, inside this folder alone, as a researcher would, and say in REPORT.md each command, each download with its source, and the minutes that the installation took. Work only in this folder. Use the network only to install Python and the packages that the wheel asks.

Continue his work from what PyPI serves, the public repository that its page links to, and their documentation, alone.
Write REPORT.md in this folder as you go, under these six headings, each with what you found, the file and the section of
the documentation where you found it, the minutes that it took, and as its last line a mark out of 10, written `Mark: N/10`:

## Step 1. What the package is for and why it was made, and what role the author's own line plays.
## Step 2. What he measured, and the results that stand, each with its uncertainty and the conditions under which it holds.
## Step 3. What is validated, and how, and what is not.
## Step 4. What he left open: the questions, the decisions, and the next steps that he planned.
## Step 5. Reproduce his results on this machine by the commands that the documentation gives, at least one share of his validation and one draw of a count, and compare.
## Step 6. Take the first of his next steps, in a folder of your own here, as far as you can, and say what you would need to go further.

A mark of 10 means that you did the step fully from the documentation alone, without a guess, and that nothing in it
blocked you, misled you or cost you time. Take a point off for each stumble of the step that cost you time, and more for
one that misled or blocked you. Write every stumble as a line of a list that begins with its rating in bold, as
`- **SLOWS** ...`, with the exact words of the page or of the output that caused it: BLOCKS (you could not continue his
work without him), MISLEADS (a statement that is wrong, or that contradicts another), SLOWS (you lost time) or POLISH.
Write the ratings nowhere else in bold: a tool counts them. End with a verdict in one sentence: whether a researcher can
continue his work from what PyPI serves and the public repository alone, the total of your marks out of 60, and the three
things that you would add first.
```

## FRESH.md, what the folder says it holds

Written for each run with the names of its files, here for 0.1.0.post1.

```
# What you received

The author of the Python package twinfit can no longer be reached, and you take over his work. You received exactly
what PyPI serves for twinfit 0.1.0.post1, given here as files since it is not uploaded yet, and the public repository
that its PyPI page links to, given here as a folder since it is not pushed yet:

| what | what it is |
|---|---|
| `twinfit-0.1.0.post1-py3-none-any.whl` | the wheel, what `pip install twinfit==0.1.0.post1` installs |
| `twinfit-0.1.0.post1.tar.gz` | the source distribution, what `pip download --no-binary :all: twinfit==0.1.0.post1` gives |
| `README.md` | the page of the package on PyPI, which is the README inside the wheel |
| `repository/` | the public repository https://github.com/MichelangeloDondi/twinfit, its files as the repository holds them, without the history of its commits |

Nothing else: no data of his, no notes, no one to ask.
```
