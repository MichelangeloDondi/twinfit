"""Every process of this work at nice 10, his message 184: "start every process at nice 10, test runs and worker pools
included", and "never stop or hold cores for my LaTeX builds".

A process of the tools, of the tests or of the toy raises its own niceness to TEN when it starts, `to_ten`, and what it
starts inherits it: the waves of the queue and their pools, the workers of a run of the tests, the copies of the
planted faults. Before his message the waves of the performance lane ran at 10 and the rest at 0, so a run of the tests
took a core before his LaTeX builds did. The package `twinfit/` does not call it: it runs in the processes of other
researchers, whose priority is theirs.

His LaTeX builds are never waited for: `latex` says which processes they are, and the queue does not count them among
the busy cores for which it holds a wave. At nice 10 the system gives them a core first."""
from __future__ import annotations

import os

TEN = 10          # the niceness of every process of this work
LATEX = ("latex", "pdflatex", "xelatex", "lualatex", "latexmk", "biber", "bibtex", "bibtex8", "makeindex", "xindy", "tectonic", "dvipdfmx", "dvips", "mpost", "makeglossaries", "texindy", "kpsewhich")          # the programs of a LaTeX build, by the name that the system gives the process


def to_ten(nice=os.nice) -> int:
    """Raises the niceness of this process to TEN where it is lower, and returns it. A process already above TEN keeps
    its own: the call only ever lowers a priority, which needs no right of the system."""
    now = nice(0)
    if now < TEN:
        now = nice(TEN - now)
    return now


def latex(command: str) -> bool:
    """Whether a process belongs to a LaTeX build, by the name of its program or of the script that its first word runs:
    latexmk runs as `perl /path/latexmk`."""
    words = command.strip().split()
    return any(os.path.basename(word) in LATEX for word in words[:2])
