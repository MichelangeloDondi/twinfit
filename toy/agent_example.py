"""The folder of one blind agent on a project that is not his line, his messages 232 to 235: a project of the gallery of
the package, or a hidden project of the cold tests, which no page of the package names, with a copy of the package and
the brief, whose rules are those of the brief of the toy in toy/noiseless.py. Apart from toy/noiseless.py since the
gallery that it imports would be read by the digest of the noiseless run.

    python toy/agent_example.py OUT NAME [--iterations N] [--minutes M]
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from toy.noiseless import BRIEF, HIDDEN, ITERATIONS, MINUTES, NUMBERS  # noqa: E402

EXAMPLE_BRIEF = """# Curves of one project of another field

You are asked to analyse the curves of a project that is not your own. Work only in this folder. Read no other folder of
this machine and do not use the network. The Python to use is `{python}`: it has numpy, scipy and matplotlib.

## The project

`project/` holds it as the package `twinfit` reads it: `project/description.toml` describes the data collection, the
setup and what is known of it, `project/traces/` holds the curves, and `project/model.py` is the first model.

## The question

Which effects do the curves hold beyond the first model, and how large is each? No list of effects is given: find them,
from what the fits leave and from the physics of the setup that the description gives, and add each as a term of
physical origin. The curves may hold none, one or several. Name only what the curves support, and say what they cannot
tell apart.

## How to work

The package `twinfit` is in this folder. Its protocol is `twinfit/agent/AGENT.md`: read it and follow it, with at most
{iterations} iterations, each adding one term. The researcher is not here: you open each step yourself, and you say in your
notes why.

{rules}
## What to hand back, within {minutes}

- `answer.json`: {{"terms": [{{"origin": "its physical origin, in words", "in the model": "how it enters the model, in words",
  "parameters": {{"name": [value, uncertainty]}}}}], "summary": "one paragraph"}}. An empty list of terms when the first model
  lacks nothing.
- `notes.md`: its first line `Context:` and the files of instructions or rules that your session carries, or `Context: none`, then what you did, in its order, and what you could not tell apart.

Read the clock with a command. When {last} minutes have passed, write the answer that you have.
"""


def make(out: Path, name: str, root: Path = ROOT, python: str = sys.executable, iterations: int = ITERATIONS, minutes: int = MINUTES) -> Path:
    """The folder of one blind agent on a project of the gallery of the package, his messages 232 to 235: the project as
    `python -m twinfit example NAME DIR` writes it, whose folder holds no truth, a copy of the package without the page of
    his line, and the brief, with the rules of the brief of the toy, every fit through the package, every term from its
    report, and a folder of each step. Refuses a folder inside the repository and a name the gallery lacks. Returns the
    folder, and the truth stays in the gallery."""
    import shutil
    out, root = Path(out).resolve(), Path(root).resolve()
    if out == root or root in out.parents:
        raise SystemExit(f"REFUSED: {out} is inside the repository: the agent reads no file of it")
    from twinfit.examples import gallery
    hidden = Path(root) / "tools" / "cold" / f"{name}.py"          # a hidden project of the cold tests, which no page of the package names: every project of the gallery is named with its truth in one sentence of the tutorial (his message 235)
    if name not in gallery.EXAMPLES and not (name.startswith("hidden_") and hidden.exists()):
        raise SystemExit(f"REFUSED: the gallery holds no project {name}, and tools/cold holds no hidden project of that name: the gallery holds {', '.join(gallery.EXAMPLES)}")
    out.mkdir(parents=True, exist_ok=True)
    if name in gallery.EXAMPLES:
        gallery.write(name, out / "project", 1.0, None)
    else:
        import subprocess
        made = subprocess.run([python, str(hidden), str(out / "project")], capture_output=True, text=True, env=dict(os.environ, PYTHONPATH=str(root)))
        if made.returncode != 0:
            raise SystemExit(f"REFUSED: the hidden project {name} was not written: {made.stderr.strip()[-300:]}")
    rules = BRIEF[BRIEF.index("The point of the package is to go beyond"):BRIEF.index("## What to hand back")]
    said_minutes = f"{minutes // 60} hours" if minutes >= 120 and minutes % 60 == 0 else f"{minutes} minutes"
    (out / "READ FIRST.md").write_text(EXAMPLE_BRIEF.format(python=python, iterations=NUMBERS[iterations] if iterations < len(NUMBERS) else iterations, rules=rules,
                                                            minutes=said_minutes, last=minutes - max(5, minutes // 8)))
    if (out / "twinfit").exists():
        shutil.rmtree(out / "twinfit")
    shutil.copytree(root / "twinfit", out / "twinfit", ignore=shutil.ignore_patterns("__pycache__", *HIDDEN, *GALLERY_FILES, f"{name}.py"))          # the gallery states the truth of each project: the clock's stands in eight files of the package, and was not run blind for that (his message 235)
    names = [n for n in gallery.EXAMPLES.get(name, {}).get("truth", {}) if n not in ("width", "centre", "height", "background")]          # a hidden project keeps its truth in tools/cold alone          # the parameters of the truth that the first model lacks, by the names that the right model gives them
    held = sorted({str(f.relative_to(out)) for f in out.rglob("*") if f.is_file() and f.suffix in (".md", ".py", ".json", ".toml", ".txt") for n in names if n in f.read_text(errors="ignore")})
    if held:
        raise SystemExit(f"REFUSED: the parameters of the truth of {name}, {', '.join(names)}, stand in {', '.join(held)}: the run would not be blind")
    return out


GALLERY_FILES = ("EXAMPLES.md", "gallery.json", "gallery.py")          # the files of the package that state the truth of every project of the gallery, used by `example` alone
TRUTHS = {          # the values that made each hidden project of the cold tests, read off its file in tools/cold, by the first words of the name that an answer gives each
    "hidden_20": (("share of the second line", 0.1),          # its height a tenth of the first's and its width the same, so a tenth of its area
                  ("distance at 175 k", 9.375),          # 2.5 widths above the first, the width 2.0 + 0.01 T cm-1: 2.5 times 3.75 at 175 K
                  ("distance per kelvin", 0.025),          # 2.5 times 0.01 cm-1 a kelvin
                  ("correlation of neighbouring samples", 0.6),
                  ("width at no temperature", 2.0),
                  ("broadening", 0.01)),
}


def beside_the_truth(project: str, answer: dict) -> list:
    """Each parameter of an answer whose name begins with the name of a value that made the hidden project, his message
    236: (its name, its value, its uncertainty, the truth, the value less the truth in its uncertainties). A parameter
    that no value of TRUTHS names, or given with no uncertainty, is left out."""
    rows = []
    for term in answer.get("terms", []):
        for name, got in (term.get("parameters") or {}).items():
            true = next((v for key, v in TRUTHS.get(project, ()) if name.lower().startswith(key)), None)
            if true is not None and isinstance(got, (list, tuple)) and len(got) == 2 and float(got[1]) > 0.0:
                rows.append((name, float(got[0]), float(got[1]), true, (float(got[0]) - true) / float(got[1])))
    return rows




if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) >= 2 and len(args) % 2 == 0 and set(args[2::2]) <= {"--iterations", "--minutes"}:
        flags = dict(zip(args[2::2], (int(v) for v in args[3::2])))
        print(f"written: {make(Path(args[0]), args[1], iterations=flags.get('--iterations', ITERATIONS), minutes=flags.get('--minutes', MINUTES))}")
    else:
        print(__doc__)
        sys.exit(1)
