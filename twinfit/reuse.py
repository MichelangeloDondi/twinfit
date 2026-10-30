"""What a run of the loop has computed, kept beside it, so that the next run computes only what is new.

Each run of the loop remakes the analysis from the files of its folder: the twin, the report and the proposals of every
model, and the twin of every term file. A run after one more term file computed again every model before it, and a run
stopped at its budget began again from its first model. A turn of the loop led by its report took 1.5 to 3.5 hours at its
third turn on one core, nearly all of it the models of the turns before (the red team of his message 327).

Every costly step of a run goes through `Reuse`. Its key is made of what the step reads: the digest of the package, the
code of the model and of the term files read so far, the records with their noise, the fit with its values, and the
arguments of the step. A step whose key is kept is read back and not computed. The steps draw from fixed seeds and give
the same numbers on any number of workers, so a value read back is the value that the step would compute.

    reuse = Reuse(folder)                                  the store of a folder, `loop.reuse/` beside `loop/`
    out = reuse("twin", lambda: admission(...), code, problem, result, sets)

A run stops at its budget only once it has computed a step anew, so a run after a stop always goes further than the one
before it: `computed` counts them. After a whole run, `prune` removes what it did not read or write, which is what the
files of the folder no longer make.
"""
from __future__ import annotations

import hashlib
import os
import pickle
import re
import time
import tomllib
from pathlib import Path

import numpy as np

STORE = "loop.reuse"
KEEP_DAYS = 30          # a kept step that no run has read or written for this many days goes: a cheaper run, of fewer term files or sets, removed the steps of the long run before it (the red team of his message 333, C4)
APART = ("loop", "loop.before", "loop.killed", STORE, "terms", "variants")          # the folders of a project whose files make no model
_ADDRESS = re.compile(r" at 0x[0-9a-fA-F]+")


def fingerprint(obj, h=None, depth: int = 0, seen=None):
    """A digest of what an argument of a step holds: arrays by their bytes, numbers and words by their value, the
    records, problems and fits by their fields, a function by its name, since its code is in the key apart. An object
    of no known kind is read by its fields, or by its words with any address of memory taken out."""
    h = hashlib.sha256() if h is None else h
    seen = set() if seen is None else seen
    if depth > 12:
        h.update(b"deep")
        return h
    if obj is None or isinstance(obj, (bool, int, float, str)):
        h.update(f"{type(obj).__name__}:{obj!r};".encode())
    elif isinstance(obj, bytes):
        h.update(b"bytes:" + obj)
    elif isinstance(obj, np.ndarray):
        a = np.ascontiguousarray(obj)
        h.update(f"array{a.shape}{a.dtype};".encode())
        if a.dtype == object:
            for item in a.ravel():
                fingerprint(item, h, depth + 1, seen)
        else:
            h.update(a.tobytes())
    elif isinstance(obj, np.generic):
        h.update(f"{type(obj).__name__}:{obj.item()!r};".encode())
    elif isinstance(obj, dict):
        h.update(f"dict{len(obj)};".encode())
        for key in sorted(obj, key=repr):
            fingerprint(key, h, depth + 1, seen)
            fingerprint(obj[key], h, depth + 1, seen)
    elif isinstance(obj, (list, tuple, set, frozenset)):
        items = sorted(obj, key=repr) if isinstance(obj, (set, frozenset)) else obj
        h.update(f"{type(obj).__name__}{len(items)};".encode())
        for item in items:
            fingerprint(item, h, depth + 1, seen)
    elif callable(obj) and not hasattr(obj, "records"):
        h.update(f"callable:{getattr(obj, '__module__', '')}.{getattr(obj, '__qualname__', type(obj).__name__)};".encode())
    elif hasattr(obj, "records") and hasattr(obj, "own_start"):          # a problem: its records, names, starts and bounds, each record's own starts and bounds read
        h.update(f"problem:{obj.name};".encode())
        for part in (obj.records, obj.shared, obj.own, obj.start_shared, obj.bounds_shared):
            fingerprint(part, h, depth + 1, seen)
        for k in range(len(obj.records)):
            fingerprint(obj.own_start(k), h, depth + 1, seen)
            fingerprint(obj.own_bounds(k), h, depth + 1, seen)
    elif id(obj) in seen:
        h.update(b"again;")
    elif hasattr(obj, "__dict__"):
        seen.add(id(obj))
        h.update(f"object:{type(obj).__name__};".encode())
        fingerprint({k: v for k, v in vars(obj).items()}, h, depth + 1, seen)
    else:
        h.update(f"{type(obj).__name__}:{_ADDRESS.sub('', repr(obj))};".encode())
    return h


class Reuse:
    """The kept steps of the runs of the loop in a folder. `computed` counts the steps that this run made anew, `read`
    those that it read back."""

    def __init__(self, folder, package: str | None = None):
        if package is None:
            from .version import digest
            package = digest()
        self.home = Path(folder) / STORE
        self.package = package
        self.computed = 0
        self.read = 0
        self.used: set = set()

    def key(self, label: str, *parts) -> str:
        h = hashlib.sha256(f"{self.package};{label};".encode())
        for part in parts:
            fingerprint(part, h)
        return h.hexdigest()

    def __call__(self, label: str, compute, *parts):
        """The value of a step: read back where its key is kept, else computed and kept. A value that cannot be kept is
        computed each time."""
        key = self.key(label, *parts)
        path = self.home / f"{key}.pickle"
        self.used.add(path.name)
        if path.is_file():
            try:
                value = pickle.loads(path.read_bytes())
                self.read += 1
                os.utime(path)          # read now: kept for `KEEP_DAYS` more
                return value
            except Exception:          # a file cut by a kill: the step is computed again
                pass
        value = compute()
        self.computed += 1
        try:
            data = pickle.dumps(value, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception:
            return value
        self.home.mkdir(parents=True, exist_ok=True)
        part = path.with_suffix(f".{os.getpid()}.part")
        part.write_bytes(data)
        os.replace(part, path)          # whole or absent: a kill leaves no half of a step
        return value

    def prune(self, days: float = KEEP_DAYS, now: float | None = None) -> int:
        """After a whole run: the steps that no run has read or written for `days` days are removed. A step that this run
        did not ask is kept until then, since a run of fewer term files or sets than the last asks fewer. Returns how
        many were removed."""
        if not self.home.is_dir():
            return 0
        now, gone = time.time() if now is None else now, 0
        for path in self.home.iterdir():
            if path.name not in self.used and now - path.stat().st_mtime > days * 86400.0:
                path.unlink()
                gone += 1
        return gone


def model_files(folder) -> list:
    """The files of a project that its model can be made of: the description, the model file that it names wherever it
    stands, and every Python file of the project outside its loop, its kept steps, its term files and its copies, since a
    model may import a file of its own. The key held model.py alone, so a model named otherwise, or a file that the model
    imports, could change and its steps be read back (the red team of his message 333, C1)."""
    folder = Path(folder)
    out = [folder / "description.toml"]
    try:
        named = tomllib.loads(out[0].read_text()).get("addenda", {}).get("model")
    except (OSError, tomllib.TOMLDecodeError):
        named = None
    if named:
        out.append(folder / str(named))
    for path in sorted(folder.rglob("*.py")):
        rel = path.relative_to(folder)
        if not any(part in APART or part.startswith(".") for part in rel.parts[:-1]) and path not in out:
            out.append(path)
    return out


def code_of(folder, *files) -> str:
    """The digest of the code that makes a model: the files of `model_files`, then each term file read so far, in order."""
    h, folder = hashlib.sha256(), Path(folder)
    for path in tuple(model_files(folder)) + tuple(Path(f) for f in files):
        h.update(f"{path.relative_to(folder).as_posix() if path.is_relative_to(folder) else path.name};".encode())          # by its place in the project, so that a copy of the project keys alike
        h.update(path.read_bytes() if path.is_file() else b"absent")
    return h.hexdigest()

