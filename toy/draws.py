"""The check that lets a model go on his recorded traces (CLAUDE.md, rule 3: "a model goes on his traces only after it
passed on its own synthetic traces").

    python toy/draws.py all        the same for every model of every subfolder, each once, side by side
    python toy/draws.py N          for each model that subfolder N puts on his traces, the model with its new term and
                                   the model that lacks it: its own fit on 100 draws of its 8 synthetic traces. Each
                                   fit is kept in runs/ as it is made, so that an interruption loses one fit. Then the
                                   count and the 100 fits go into toy/draws.json
    python toy/draws.py --check    ends without error only when the table holds, for every such model of every detail
                                   that config.toml builds, a check on 100 draws that was made from the files as they
                                   are and that passes when its fits are counted again, and when three of the kept
                                   fits of each model come out the same when they are made again

What passes is said by `fits.check`: every search ended, a bias under a quarter of the spread, and intervals of one
uncertainty that hold the truth on a share of the draws within three standard errors of 0.6827: 0.543 to 0.822 on 100.
"""
from __future__ import annotations

import fcntl
import os
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from toy import fits, model  # noqa: E402

DRAWS = 100
DEPENDS = ("toy/model.py", "toy/vapour.py", "toy/terms/__init__.py", "toy/traces.py", "toy/fits.py", "toy/draws.py")          # the programs that every fit runs, whatever its detail


VERSION = "twinfit/version.py"          # the number of the package and the record of its versions: no result is made with it, and a new number alone moves none


LOADED = {p.resolve(): p.read_bytes() for home in (ROOT / "toy", ROOT / "twinfit") for p in home.rglob("*.py")}          # the code of the toy and of the package as this process loaded it, read once at its start


def code_of(path: Path) -> bytes:
    """What a result depends on in a file: for a program, its code read as Python, with no docstring and no comment,
    so that a word of a docstring moves no count, and every word that a page prints still does. For any other file, its
    bytes. On 2026-10-16 every count hashed every file of the package whole, and a word of the tutorial made a day of
    counts stale (his message 175: "red team your efficiency and effectiveness"). A program of the toy or of the package
    is read as this process loaded it, LOADED: a dry run that began at 15:11 of 2026-10-23 kept readings made by the code
    of 15:11 under the digest of files edited at 16:00, which a later build would have read as fresh."""
    import ast
    path = Path(path)
    if path.suffix != ".py":
        return path.read_bytes()
    return code_of_text(LOADED[path.resolve()].decode() if path.resolve() in LOADED else path.read_text())


def code_of_text(text: str) -> bytes:
    """The code of a program read as Python, with no docstring and no comment."""
    import ast
    tree = ast.parse(text)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and node.body and isinstance(node.body[0], ast.Expr) \
                and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
            node.body = node.body[1:] or [ast.Pass()]
    return ast.dump(tree, include_attributes=False).encode()


def hashed(h, root: Path, names) -> None:
    """Each named file of `root` into the hash `h`, by its name and its code (code_of)."""
    for name in names:
        h.update(name.encode() + b"\0" + code_of(Path(root) / name) + b"\0")


LOADS_FILES = ("twinfit/loop.py",)          # a file of the package that loads by importlib the researcher's model and term files, never a module of the package: what it runs is read in its text all the same


def package_of(root: Path, programs) -> list:
    """The files of the package that some programs of the toy run: those that they import, and those that these import
    in turn, read in the text of each file. It is every file that could run, and so no less than what does: a file of
    the package that none of them reaches moves no result of theirs, and a new version of the package leaves the kept
    fits where it changed none of their files. The file of the version is left out. A file of the package that imports
    by a name made while it runs would hide what it runs: it is refused."""
    import ast
    root = Path(root)
    home = root / "twinfit"

    def files_of(parts: list) -> list:          # the files that the module twinfit.a.b names: each folder on the way, and the module itself
        out, at = [home / "__init__.py"], home
        for part in parts:
            if (at / part).is_dir():
                at = at / part
                out.append(at / "__init__.py")
                continue
            if (at / f"{part}.py").exists():
                out.append(at / f"{part}.py")
            break          # a name inside a module is no file
        return out

    def imported(path: Path) -> list:
        inside, found = home in path.parents, []
        loads = path.relative_to(root).as_posix() in LOADS_FILES          # its importlib loads the researcher's files: only the refusal is lifted, and its imports are read as any file's
        for node in ast.walk(ast.parse(path.read_text())):
            names = [alias.name for alias in node.names] if isinstance(node, (ast.Import, ast.ImportFrom)) else []
            if inside and not loads and (any(name.split(".")[0] == "importlib" for name in names) or (isinstance(node, ast.ImportFrom) and (node.module or "").split(".")[0] == "importlib")
                           or (isinstance(node, ast.Name) and node.id == "__import__")):
                raise SystemExit(f"REFUSED: {path.relative_to(root)} imports by a name made while it runs: what a kept result depends on can no longer be read in its text")
            if isinstance(node, ast.Import):
                for name in names:
                    if name.split(".")[0] == "twinfit":
                        found += files_of(name.split(".")[1:])
            elif isinstance(node, ast.ImportFrom):
                module = (node.module or "").split(".") if node.module else []
                if node.level == 0 and module[:1] == ["twinfit"]:
                    base = module[1:]
                elif node.level and inside:
                    here = list(path.relative_to(home).parent.parts)          # the folder of this file inside the package
                    base = here[:len(here) - (node.level - 1)] + module
                else:
                    continue          # a program of the toy that imports its neighbours: no file of the package
                found += files_of(base)
                for name in names:
                    found += files_of(base + [name])
        return found
    seen, queue = set(), [root / name for name in programs]
    while queue:
        for path in imported(queue.pop()):
            if path not in seen and path.exists():
                seen.add(path)
                queue.append(path)
    return sorted(name for name in (path.relative_to(root).as_posix() for path in seen) if name != VERSION)


def depends(root: Path = ROOT, detail: int | None = None) -> list:
    """The programs that a fit of the first `detail` details runs: those of every fit, the terms that have entered the
    model by that detail, and the files of the package that those programs run. With no detail: every term."""
    root = Path(root)
    const, _ = model.load(root / "toy", root)
    terms = sorted(name for name, d in const["details"].items() if detail is None or d <= detail)
    return list(DEPENDS) + [f"toy/terms/{name}.py" for name in terms] + package_of(root, DEPENDS)


def digest(root: Path = ROOT, detail: int | None = None) -> str:
    """What a check on draws of the first `detail` details depends on: the programs of `depends`, and the values that
    the subfolder of that detail is made from, with the way in which each of its parameters changes the line. With no
    detail: the last one. The two files are not read as texts: a sentence, a source or a term that enters later moves
    no fit, and a new detail leaves the kept fits of the details before it as they are (the owner's message 119: "Are
    you ready to add further terms and details to the model? Are you working whit a smart and efficient and effective
    workflow?")."""
    root = Path(root)
    const, config = model.load(root / "toy", root)
    detail = max(const["details"].values()) if detail is None else int(detail)
    h = hashlib.sha256()
    hashed(h, root, depends(root, detail))
    effects = {name: e.get("effect") for name, e in sorted(config["parameter"].items()) if e["enters"] <= detail}
    h.update(json.dumps([model.values_of(detail, const, config), effects], sort_keys=True).encode())
    return h.hexdigest()


def table(root: Path = ROOT) -> dict:
    path = Path(root) / "toy" / "draws.json"
    return json.loads(path.read_text()) if path.exists() else {}


def key(detail: int, free: list | None, config: dict) -> str:
    """The name of a model in the table: the number of its detail, with its free parameters when they are not those that
    config.toml gives the detail."""
    own = list(config["detail"][str(detail)]["free"])
    return str(detail) if free is None or list(free) == own else f"{detail} with {' and '.join(free) or 'nothing'} free"


def models(detail: int, config: dict) -> list:
    """The models that the subfolder of a detail puts on his traces, each as (the detail of the model, its free
    parameters or None for those of config.toml): the model with the new term, and after the first detail the model
    that lacks it."""
    out = [(detail, None)]
    if detail > 1:
        out.append((detail - 1, fits.lacking(detail, config)))
    return out


def run(detail: int, draws: int = DRAWS, root: Path = ROOT, fit=None, free: list | None = None) -> dict:
    """The own fit of a model on `draws` draws, each kept as it is made, and its count written into the table with the
    values and the uncertainties of every fit: `passed` counts them again."""
    root = Path(root)
    const, config = model.load(root / "toy", root)
    made, name = digest(root, detail), key(detail, free, config)
    kept = root / "runs" / f"draws_{name.replace(' ', '_')}_{made[:12]}.jsonl"
    kept.parent.mkdir(exist_ok=True)
    with kept.with_suffix(".lock").open("w") as mine:          # one run of a model at a time: two would write their fits into one file
        fcntl.flock(mine, fcntl.LOCK_EX)
        got = [json.loads(line) for line in kept.read_text().splitlines()] if kept.exists() else []
        with kept.open("a") as fh:
            for k in range(len(got), draws):
                g = (fit or fits.fit)(detail, fits.synthetic(detail, const, config, 1 + k), const, config, free)
                got.append({"told": g["told"], "values": g["values"], "bars": g["bars"]})
                fh.write(json.dumps(got[-1]) + "\n")
                fh.flush()
    out = dict(fits.check(detail, const, config, draws, got[:draws], free), free=free, made_from=made, at=time.strftime("%Y-%m-%d %H:%M:%S"), fits=got[:draws])
    path = root / "toy" / "draws.json"
    with (root / "runs" / "draws.lock").open("w") as lock:          # several models may end in the same second
        fcntl.flock(lock, fcntl.LOCK_EX)
        held = table(root)
        before = held.get(name)
        if before and before.get("made_from") == made and len(before.get("fits", [])) > len(out["fits"]):
            return before          # a shorter run does not take the place of a longer one made from the same files
        if before and before.get("made_from") == made and before.get("fits") == out["fits"]:
            return before          # the same fits from the same files: the table and its time stay as they are, and nothing that reads it is made stale (review 3, minor 2.5)
        held[name] = out
        whole = path.with_suffix(".json.part")
        whole.write_text(json.dumps(dict(sorted(held.items())), indent=1) + "\n")
        os.replace(whole, path)          # the table is whole or it is the one before: a run that is stopped leaves no half of it
    return out


def audit(detail: int, root: Path = ROOT, free: list | None = None, fit=None) -> list:
    """Three of the kept fits of a model made again, chosen by the digest of the table: what differs from the table.
    An empty list: the table holds what the code gives. A table whose fits a hand wrote would differ here."""
    root = Path(root)
    entry = current(detail, root, free)
    if not entry or not entry.get("fits"):
        return ["no check made from the files as they are"]
    const, config = model.load(root / "toy", root)
    wrong, n = [], len(entry["fits"])
    digest = hashlib.sha256(json.dumps(entry["fits"], sort_keys=True).encode()).digest()          # which three are made again follows from the table itself: a hand that changes one fit moves the three (review 3, minor 2.2)
    for k in sorted({int.from_bytes(digest[4 * i:4 * i + 4], "big") % n for i in range(3)}):
        g = (fit or fits.fit)(detail, fits.synthetic(detail, const, config, 1 + k), const, config, free)
        for name, value in entry["fits"][k]["values"].items():
            if abs(g["values"][name] - value) > 1e-9 * max(abs(value), 1.0):
                wrong.append(f"draw {k + 1} of the model {key(detail, free, config)}: {name} is {g['values'][name]:.10g} when made again and {value:.10g} in the table")
    return wrong


def current(detail: int, root: Path = ROOT, free: list | None = None) -> dict | None:
    """The check of a model in the table, when it was made from the files as they are. None when there is none, or when
    it was made from other files."""
    root = Path(root)
    _, config = model.load(root / "toy", root)
    entry = table(root).get(key(detail, free, config))
    return entry if entry and entry.get("made_from") == digest(root, detail) else None


def passed(detail: int, root: Path = ROOT, free: list | None = None) -> bool:
    """Whether a model may go on his recorded traces: the table holds its check on at least 100 draws, made from the
    files as they are, and the fits that the table keeps pass when they are counted again here. The verdict that the
    table states is not believed: a hand that changed it would open the pages on his traces (review 1, moderate 5)."""
    root = Path(root)
    entry = current(detail, root, free)
    if not entry or len(entry.get("fits", [])) < DRAWS:
        return False
    const, config = model.load(root / "toy", root)
    return bool(fits.check(detail, const, config, len(entry["fits"]), entry["fits"], free)["passes"])


def every(config: dict) -> list:
    """Each model that the folder puts on his traces, once: those of every detail that config.toml builds."""
    out = []
    for d in config["folder"]["details"]:
        out += [m for m in models(d, config) if key(*m, config) not in {key(*o, config) for o in out}]
    return out


def one(asked: tuple) -> str:
    """The check of one model, for a pool of processes: its name in the table and whether it passes."""
    detail, free = asked
    _, config = model.load()
    return f"{key(detail, free, config)}: " + ("passes" if run(detail, free=free)["passes"] else "does not pass")


def main(argv: list) -> int:
    if argv == ["all"]:
        from concurrent.futures import ProcessPoolExecutor
        _, config = model.load()
        asked = every(config)
        with ProcessPoolExecutor(len(asked)) as pool:          # one process a model: no two write the same file of kept fits
            for line in pool.map(one, asked):
                print(line)
        return 0
    if argv == ["--check"]:
        _, config = model.load()
        missing = [(d, free) for d, free in every(config) if not passed(d, free=free)]
        forged = [line for d, free in every(config) if (d, free) not in missing for line in audit(d, free=free)]
        for line in forged:
            print("REFUSED: " + line)
        for d, free in missing:
            entry = table().get(key(d, free, config))
            why = "no check in the table" if entry is None else "made from other files than those on disk" if entry.get("made_from") != digest(detail=d) else f"on {len(entry.get('fits', []))} draws" if len(entry.get("fits", [])) < DRAWS else "it does not pass"
            print(f"REFUSED: the model {key(d, free, config)} may not go on his recorded traces: {why}")
        print("the check on 100 draws: " + ("every model passes, and three fits of each were made again" if not missing and not forged else f"{len(missing)} models do not pass, {len(forged)} kept fits differ"))
        return 1 if missing or forged else 0
    if len(argv) == 1 and argv[0].isdigit():
        _, config = model.load()
        for d, free in models(int(argv[0]), config):
            out = run(d, free=free)
            print(key(d, free, config) + ": " + json.dumps({k: v for k, v in out.items() if k not in ("made_from", "fits")}))
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
