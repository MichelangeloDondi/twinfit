"""A count of the toy run on any machine, without the queue of the author's own (his messages 441 and 442, the bus test).
Every count of the toy speaks one protocol, which the queue spoke: `python toy/PROGRAM.py units` lists what is left, by
condition, each condition with the conditions it waits for, `TOY_QUEUE=... python toy/PROGRAM.py slice CONDITION FROM TO`
makes units FROM to TO of one condition and keeps each row in runs/ as it ends, and `python toy/PROGRAM.py table` counts
the kept rows into its table, toy/PROGRAM.json, or toy/gallery.json for the gallery. This runs the slices, one unit a process, several side by side, a condition after
those it waits for, asks again what is left until nothing is, and writes the table only when the count is whole: a count
made small to try it leaves the committed table as it is.

    python toy/count.py PROGRAM [--workers N] [--draws K] [--only TEXT]
        PROGRAM       pointing, power, gallery_count, chosen or scoreboard. The benchmark, C10, counts nothing of its
                      own: after a count of power, `python toy/benchmark.py table .` writes its table from power's
        --workers N   N processes side by side, every core of the machine by default
        --draws K     at most units 1 to K of each condition of draws, a count made small: the twin of a case, which its
                      draws wait for, is always made whole
        --only TEXT   only the conditions whose name holds TEXT, and the twins that they wait for

A unit of a count takes from seconds to several minutes on one core, and the whole counts of the README took tens of core
hours (toy/README.md gives the cost of each). A change of the code that a count runs makes all its kept rows stale, since
their file names carry the digest of that code: the count then starts again from its first unit.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROGRAMS = ("pointing", "power", "gallery_count", "chosen", "scoreboard")          # the counts that the criteria read: noiseless.py has a run and no table
TABLES = {"gallery_count": "gallery"}          # a table named otherwise than its program


def ask(program: str, *args: str) -> subprocess.CompletedProcess:
    env = dict(os.environ, TOY_QUEUE="toy/count.py")          # the count's own refusal of a slice made outside a queue, which this runner is
    return subprocess.run([sys.executable, str(ROOT / "toy" / f"{program}.py"), *args], cwd=str(ROOT), env=env, capture_output=True, text=True)


def left(program: str) -> list:
    got = ask(program, "units")
    if got.returncode != 0:
        raise SystemExit(f"`python toy/{program}.py units` failed: {(got.stderr or got.stdout).strip()[-600:]}")
    return json.loads(got.stdout.strip().splitlines()[-1])


def chosen(units: list, draws: int | None, only: str | None) -> list:
    """The conditions to run now, each as (condition, its units): those whose name holds `only`, and the twins they wait
    for, units 1 to `draws` of each condition of draws, and none that waits on a condition with units left."""
    wanted = [u for u in units if only is None or only in u["condition"]]
    needed = {a for u in wanted for a in u.get("after", [])}
    wanted += [u for u in units if u["condition"] in needed and u not in wanted]
    pending = {u["condition"] for u in units if u["units"]}
    out = []
    for u in wanted:
        if not u["units"] or any(a in pending for a in u.get("after", [])):
            continue
        ks = u["units"] if (draws is None or u.get("kind") == "twin") else [k for k in u["units"] if int(k) <= draws]
        if ks:
            out.append((u["condition"], ks))
    return out


def run(program: str, workers: int, draws: int | None, only: str | None) -> int:
    began, made, failed = time.time(), 0, []
    while True:
        now = chosen(left(program), draws, only)
        if not now:
            break
        jobs = [(c, k) for c, ks in now for k in ks]
        print(f"{len(jobs)} units of {len(now)} conditions, {workers} side by side", flush=True)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(ask, program, "slice", c, str(k), str(k)): (c, k) for c, k in jobs}
            for f in as_completed(futures):
                c, k = futures[f]
                got = f.result()
                if got.returncode == 0:
                    made += 1
                    print(f"  {c} {k}: made, {made} in {time.time() - began:.0f} s", flush=True)
                else:
                    failed.append(f"{c} {k}: {(got.stderr or got.stdout).strip()[-300:]}")
                    print(f"  {c} {k}: FAILED", flush=True)
        if failed:
            break
    if failed:
        print("failed:\n" + "\n".join(failed))
        return 1
    whole = not any(u["units"] for u in left(program))
    if whole:
        table = ask(program, "table")
        print(table.stdout.strip() or table.stderr.strip())
        return table.returncode
    print(f"{made} units made in {time.time() - began:.0f} s. The count is not whole, so toy/{TABLES.get(program, program)}.json is left as it is: "
          f"`python toy/{program}.py units` says what is left")
    return 0


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args or args[0] not in PROGRAMS:
        print(__doc__)
        sys.exit(1)
    opts = dict(zip(args[1::2], args[2::2]))
    if len(args) % 2 == 0 or set(opts) - {"--workers", "--draws", "--only"}:
        print(__doc__)
        sys.exit(1)
    sys.exit(run(args[0], int(opts.get("--workers", os.cpu_count() or 1)), int(opts["--draws"]) if "--draws" in opts else None, opts.get("--only")))
