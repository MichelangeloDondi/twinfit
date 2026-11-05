"""The convergence of each parameter turn after turn of refinement, his message 299: "I would like to see for each
parameter the plot of its convergence turn after turn of refinement of the model [...] hopefully also a way to have a good
bound to the biases of the parameters we are trying to estimate".

A turn is a model of the toy, details 1 to 5, each the last with one term more. At each turn the joint fit of his 8
recorded traces, and of DRAWS sets of the twin's synthetic traces of the fifth detail in each of three arms, whose truth is
the twin's values: so the bias of each parameter at each turn is measured on the twin and not assumed. The plain arm is
the twin as it is, where the fifth turn holds every term and its bias is zero by construction. The two planted arms hold
the structure that his traces leave under the fifth model and that the twin lacks (the red team of message 299, severe
1): each trace with the laser noise that it takes alone, or with the stretch of its scan, as `toy/scale_asked.json` holds
them for each of his traces, so the fifth turn too has a bias. Each turn is fitted as its model holds theory, the main
fit, and with the values of theory freed as a test, the fit of `freed` (rule 10). The bound of a parameter at a turn is
the largest of the biases of the three arms and of the moves of his value from that turn to each later turn: BOUND. The
turns after the fifth wait for his word: this program reads the turns that `toy/config.toml` holds, so a sixth detail is
a sixth turn with no change here.

    python3 tools/autoqueue.py add toy/convergence.py   the count through the queue: his traces and the draws of each turn
    python toy/convergence.py table                      writes toy/convergence.json from the kept units
    python toy/convergence.py plot OUT                   one slide a parameter and their PDF in OUT
    python toy/convergence.py --check                    ends without error when that table was made from the files as they are
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from toy import fits, model, real  # noqa: E402

TRUTH = 5          # the detail whose synthetic traces the twin draws: every term that the toy holds
DRAWS = 400        # the third stage, his message 502 on decision 27, 200 at the second: the sets of the twin at each turn in each arm, as the check of each model on its own traces counts them: 5 seconds a set
FITS = ("main", "freed")
ARMS = ("plain", "width", "stretch")          # the twin as it is, and planted with each trace's own laser noise or the stretch of its scan
BOUND = "the largest of the three arms' biases, of the moves of his value to each later turn, and of its move between theory held and theory freed at the same turn. It covers the two structures planted, each trace's own laser noise and stretch, and no other"
READING = ("READING", "BOUND", "summary", "bound", "table", "progress", "plot", "check")          # what reads the kept units: no part of a unit, so a change of it is no recount
RUNS = ("toy/convergence.py", "toy/fits.py", "toy/model.py", "toy/real.py", "toy/traces.py", "toy/scale.py", "toy/twin.py")


def turns(config: dict) -> list:
    """The turns that the toy holds, in order: the details of config.toml."""
    return sorted(int(k) for k in config["detail"])


def conditions(config: dict) -> list:
    """Every condition, as (its name, the turn, the fit): the main fit of each turn, and the freed fit where the turn holds
    a value of theory."""
    out = []
    for k in turns(config):
        entry = config["detail"][str(k)]
        for which in FITS:
            if which == "main" or [n for n in entry.get("freed", []) if n not in entry["free"]]:
                out.append((f"turn {k}, {which}", k, which))
    return out


def arm_of(u: int) -> tuple:
    """Unit u, from 1: his 8 recorded traces at 1, then DRAWS draws of each arm in the order of ARMS."""
    return ("his", 0) if u == 1 else (ARMS[(u - 2) // DRAWS], (u - 2) % DRAWS + 1)


def planted(arm: str, draw: int, const: dict, config: dict, root: Path = ROOT) -> list:
    """The twin's 8 traces of the fifth detail on one draw, each with his height, his message 305, as an arm draws them:
    `plain` with one laser noise for all 8, as the model holds it and the twin draws since his message 449, `width` each
    trace with its own laser noise, which the twin carried from his message 307 to his message 449 and toy/twin.toml keeps as a record, `stretch` one laser noise and each trace with the stretch of its scan that his
    trace takes alone, where his trace tells it from nothing, and `both` its own laser noise and its own stretch, all as
    `toy/twin.py` draws them from toy/twin.toml. The draws of the arms share their noise, so their biases are compared
    draw by draw."""
    from toy import scale, twin
    mine = twin.load(Path(root) / "toy")
    own = {k: v for k, v in mine.items() if k != "laser_noise"}          # the arm plants each trace's own value from the table, so the twin's one-for-all of his message 449 is left out
    rows = twin.synthetic(TRUTH, const, config, draw, own if arm in ("width", "both") else {"heights": mine.get("heights", {})})
    if arm in ("plain", "width"):
        return rows
    asked = scale.asked(root)
    out = []
    for row in rows:
        got = asked[f"{row['T']} C, {row['P']} mW"]
        if got["stretch_bar"] < abs(got["stretch"]):          # the scan read wider or narrower than it is, as `scale.drawn` plants it
            row = {**row, "nu": np.asarray(row["nu"]) / (1.0 + float(got["stretch"]))}
        out.append(row)
    return out


def unit(condition: str, u: int, root: Path = ROOT) -> dict:
    """Unit u of a condition: his 8 recorded traces, or a draw of an arm of the twin, `arm_of`. The joint fit of the turn's
    model with the free values of its fit, their bars, and the chi square."""
    const, config = model.load(Path(root) / "toy", root)
    _, k, which = next(c for c in conditions(config) if c[0] == condition)
    entry = config["detail"][str(k)]
    free = list(entry["free"] if which == "main" else entry["freed"])
    arm, draw = arm_of(u)
    rows = real.rows(const, config, root) if arm == "his" else planted(arm, draw, const, config, root)
    try:
        got = fits.fit(k, rows, const, config, free)
    except (ValueError, np.linalg.LinAlgError, RuntimeError) as e:          # a fit that does not end is said and counted apart
        return {"failed": f"{type(e).__name__}: {e}"[:300]}
    return {"values": got["values"], "bars": got["bars"], "told": got["told"], "chi_square_a_degree": got["chi2"] / max(got["dof"], 1)}


def unit_code(root: Path = ROOT) -> bytes:
    """The code of this program that makes a unit, read as Python with no docstring: every statement but those named in
    READING and the lines that run it from the shell, so that a change of the table or of the plot is no recount."""
    import ast
    tree = ast.parse((Path(root) / "toy" / "convergence.py").read_text())
    tree.body = [n for n in tree.body if not (isinstance(n, ast.FunctionDef) and n.name in READING) and not (isinstance(n, ast.Assign) and any(getattr(t, "id", None) in READING for t in n.targets))
                 and not (isinstance(n, ast.If) and "__main__" in ast.unparse(n.test))]
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef)) and node.body and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) and isinstance(node.body[0].value.value, str):
            node.body = node.body[1:] or [ast.Pass()]
    return ast.unparse(tree).encode()


def digest(root: Path = ROOT) -> str:
    """What a unit depends on: the code of this program that makes it, the model, the fits, the traces, the values of the
    two files, his traces and the files of the package that the fits import."""
    from toy import draws
    h = hashlib.sha256()
    h.update(b"toy/convergence.py\0" + unit_code(root) + b"\0")
    names = [n for n in RUNS if n != "toy/convergence.py"] + ["toy/const.toml", "toy/config.toml", "toy/twin.toml", "data/recorded/INDEX.json", "toy/scale_asked.json"]
    draws.hashed(h, Path(root), list(dict.fromkeys(names + draws.depends(root, TRUTH) + draws.package_of(Path(root), RUNS))))
    return h.hexdigest()


def unit_path(root: Path, condition: str, u: int) -> Path:
    const, config = model.load(Path(root) / "toy", root)
    return Path(root) / "runs" / "convergence" / f"{[c[0] for c in conditions(config)].index(condition):02d}_{u:03d}.json"


def kept(root: Path, condition: str, u: int, made_from: str):
    """What a kept unit holds under the digest of the files as they are, or None."""
    path = unit_path(root, condition, u)
    if not path.exists():
        return None
    got = json.loads(path.read_text())
    return got["result"] if got.get("digest") == made_from else None


def units(root: Path = ROOT) -> list:
    """What the queue asks: each condition with its units not yet kept, his traces first."""
    made_from = digest(root)
    const, config = model.load(Path(root) / "toy", root)
    return [{"condition": name, "units": [u for u in range(1, len(ARMS) * DRAWS + 2) if kept(root, name, u, made_from) is None], "after": [], "kind": f"turn {k}"}
            for name, k, _ in conditions(config)]


def slice_of(condition: str, first: int, last: int, limit: float = float("inf"), root: Path = ROOT) -> dict:
    """The units of a condition from `first` to `last` not yet kept, each kept as it ends, until `limit` seconds."""
    import resource
    made_from, began, ids, seconds = digest(root), time.time(), [], []
    for u in range(first, last + 1):
        if kept(root, condition, u, made_from) is not None:
            continue
        t = time.time()
        result = unit(condition, u, root)
        path = unit_path(root, condition, u)
        path.parent.mkdir(parents=True, exist_ok=True)
        aside = path.with_suffix(f".{os.getpid()}.part")
        aside.write_text(json.dumps({"digest": made_from, "result": result}) + "\n")
        aside.replace(path)
        ids.append(u)
        seconds.append(time.time() - t)
        if time.time() - began > limit:
            break
    return {"done": ids, "seconds": seconds, "rss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


def summary(his: dict, draws: list, truth: dict) -> dict:
    """One arm of one condition read: for each parameter that its fit frees, over the twin's draws that ended and told
    their parameters apart, the mean, the spread, the bias from the truth with its error, in units of the spread too, and
    the share of draws whose interval of one bar holds the truth. Beside them his value and bar."""
    ended = [d for d in draws if "failed" not in d and d["told"]]
    out = {"draws": len(draws), "failed": sum("failed" in d for d in draws), "not_told_apart": sum("failed" not in d and not d["told"] for d in draws), "parameters": {}}
    names = list((his.get("values") or (ended[0]["values"] if ended else {})).keys())
    for name in names:
        x = np.array([d["values"][name] for d in ended], dtype=float)
        b = np.array([d["bars"][name] for d in ended], dtype=float)
        t = float(truth[name])
        spread = float(np.std(x, ddof=1)) if x.size > 1 else None
        out["parameters"][name] = {"truth": t, "mean": float(x.mean()) if x.size else None, "spread": spread, "bias": float(x.mean() - t) if x.size else None,
                                   "bias_in_spreads": float((x.mean() - t) / spread) if spread else None, "bias_error": float(spread / np.sqrt(x.size)) if spread else None,
                                   "holds": float(np.mean(np.abs(x - t) <= b)) if x.size else None,
                                   "his": (his.get("values") or {}).get(name), "his_bar": (his.get("bars") or {}).get(name)}
    return out


def bound(conditions_read: dict) -> None:
    """BOUND for each parameter of each condition, written into it: the largest of the absolute biases of the three arms,
    of the moves of his value from the turn to each later turn of the same fit, and of its move between the fit that holds
    theory and the fit that frees it at the same turn, which the bound of the fifth turn left out while it was five times
    the bound of the laser noise (the red team of message 303, moderate 3). A parameter that no arm reads has no bound."""
    for c in conditions_read.values():
        later = [d for d in conditions_read.values() if (d["fit"] == c["fit"] and d["turn"] > c["turn"]) or (d["turn"] == c["turn"] and d["fit"] != c["fit"])]
        for name, p in c["parameters"].items():
            biases = [abs(c["arms"][a]["parameters"][name]["bias"]) for a in c["arms"] if c["arms"][a]["parameters"].get(name, {}).get("bias") is not None]
            moves = [abs(d["parameters"][name]["his"] - p["his"]) for d in later if name in d["parameters"] and d["parameters"][name]["his"] is not None and p["his"] is not None]
            p["bound"] = max(biases + moves) if biases else None
            p["bound_from"] = None if not biases else ("an arm" if max(biases) >= max(moves, default=0.0) else "a move of his value")


def table(root: Path = ROOT, path: Path | None = None) -> dict:
    """toy/convergence.json, or `path`, from the kept units: each condition with his values, each arm read and the bound of
    each parameter. Refused while a unit is not kept."""
    made_from = digest(root)
    const, config = model.load(Path(root) / "toy", root)
    truth = model.truth(config)
    out = {"digest": made_from, "draws": DRAWS, "arms": list(ARMS), "truth_detail": TRUTH, "bound": BOUND, "conditions": {}}
    for name, k, which in conditions(config):
        got = [kept(root, name, u, made_from) for u in range(1, len(ARMS) * DRAWS + 2)]
        if None in got:
            raise SystemExit(f"REFUSED: {name} holds {sum(g is not None for g in got)} of {len(ARMS) * DRAWS + 1} units under the files as they are: add toy/convergence.py to the queue again")
        his = got[0]
        arms = {arm: summary(his, got[1 + j * DRAWS:1 + (j + 1) * DRAWS], truth) for j, arm in enumerate(ARMS)}
        names = list((his.get("values") or {}).keys()) or list(arms["plain"]["parameters"])
        out["conditions"][name] = {"turn": k, "fit": which, "his": his, "arms": arms, "failed": max(a["failed"] for a in arms.values()),
                                   "parameters": {n: {"his": (his.get("values") or {}).get(n), "his_bar": (his.get("bars") or {}).get(n), "truth": float(truth[n])} for n in names}}
    bound(out["conditions"])
    Path(path or Path(root) / "toy" / "convergence.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def progress(root: Path = ROOT) -> str:
    made_from = digest(root)
    const, config = model.load(Path(root) / "toy", root)
    total = len(ARMS) * DRAWS + 1
    n = sum(kept(root, name, u, made_from) is not None for name, _, _ in conditions(config) for u in range(1, total + 1))
    return f"units kept: {n} of {len(conditions(config)) * total}"


def plot(out: Path, root: Path = ROOT, path: Path | None = None) -> list:
    """One slide a parameter, and one for the beam radius that the transit gives: on his traces its value and bar turn
    after turn, filled for the main fit and open where the values of theory are freed, beside the mean of each arm of the
    twin and the truth that the twin draws, and the bound of BOUND as a grey band about his value. Returns the paths of the
    slides and of their PDF."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_pdf import PdfPages
    from toy import meeting, pages
    held = json.loads(Path(path or Path(root) / "toy" / "convergence.json").read_text())
    const, config = model.load(Path(root) / "toy", root)
    waist = float(const["beam"]["waist_um"])
    names = list(dict.fromkeys(n for c in held["conditions"].values() for n in c["parameters"]))
    shown = names + (["w_0_um"] if "transit" in names else [])
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    figs, paths, drawn = [], [], {}
    foot = f"His 8 recorded traces and {held['draws']} sets of each arm of the twin, toy/convergence.json."
    marks = {"plain": "s", "width": "^", "stretch": "v"}
    for name in shown:
        source = "transit" if name == "w_0_um" else name
        scale_of = (lambda v, t: waist * t / v) if name == "w_0_um" else (lambda v, t: v)          # the radius that a transit rate gives, the rate going as one over the radius
        title = "the beam radius from the transit, turn after turn" if name == "w_0_um" else f"{name.replace('_', ' ')}, turn after turn"
        fig = meeting.slide(title, "Colour: his traces, filled held, open freed. Grey: the twin's three arms. Band: the bound.", foot)
        ax = fig.add_axes([0.12, meeting.BASE, 0.84, fig.below - meeting.BASE - 0.08])
        drawn[name] = []
        for c in held["conditions"].values():
            p = c["parameters"].get(source)
            if p is None:
                continue
            x0 = c["turn"] + (-0.15 if c["fit"] == "main" else 0.15)
            t = p["truth"]
            for j, arm in enumerate(held["arms"]):
                q = c["arms"][arm]["parameters"].get(source, {})
                if q.get("mean") is not None:
                    ax.errorbar(x0 + 0.05 * (j - 1), scale_of(q["mean"], t), yerr=abs(scale_of(q["mean"], t) * q["spread"] / q["mean"]) if name == "w_0_um" else q["spread"] or 0.0,
                                fmt=marks[arm], ms=8, lw=1.8, capsize=4, color=pages.GREY)
            if p["his"] is not None:
                value = scale_of(p["his"], t)
                bar = abs(value * p["his_bar"] / p["his"]) if name == "w_0_um" and p["his_bar"] else (p["his_bar"] or 0.0)
                if p.get("bound") is not None:
                    width = abs(value * p["bound"] / p["his"]) if name == "w_0_um" else p["bound"]
                    ax.fill_between([x0 - 0.12, x0 + 0.12], value - width, value + width, color=pages.PALE, zorder=0)
                ax.errorbar(x0 - 0.1, value, yerr=bar, fmt="o", ms=11, lw=2.6, capsize=6, color=meeting.PALETTE[0], mfc=meeting.PALETTE[0] if c["fit"] == "main" else "white")
                drawn[name].append([c["turn"], c["fit"], value, bar, p.get("bound")])
        truth = waist if name == "w_0_um" else float(model.truth(config)[name])
        ax.axhline(truth, color=pages.INK, lw=1.4, ls="--")
        ax.set_xticks(sorted({c["turn"] for c in held["conditions"].values()}))
        ax.set_xlabel("turn: the model of detail", fontsize=meeting.LABEL)
        unit_of = "um" if name == "w_0_um" else config["parameter"][name].get("unit", "")
        ax.set_ylabel(("w_0" if name == "w_0_um" else name.replace("_", " ")) + (f", {unit_of}" if unit_of and len(unit_of) < 8 else ""), fontsize=meeting.LABEL)
        meeting.plain(ax)
        meeting.inside(fig, title)
        figs.append(fig)
        path_png = out / f"convergence {name}.png"
        fig.savefig(path_png, dpi=110)
        paths.append(path_png)
    with PdfPages(out / "convergence.pdf") as pdf:
        for fig in figs:
            pdf.savefig(fig)
    (out / "convergence_slides.json").write_text(json.dumps({"parameters": shown, "drawn": drawn, "bound": held["bound"], "digest": held["digest"]}, indent=1) + "\n")
    return paths + [out / "convergence.pdf"]


def check(root: Path = ROOT, path: Path | None = None) -> list:
    path = Path(path or Path(root) / "toy" / "convergence.json")
    if not path.exists():
        return ["toy/convergence.json does not exist: `python3 tools/autoqueue.py add toy/convergence.py`, then `python toy/convergence.py table`"]
    return [] if json.loads(path.read_text()).get("digest") == digest(root) else ["toy/convergence.json was made from other files than those on disk now: count it again"]


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["units"]:
        print(json.dumps(units()))
    elif args[:1] == ["slice"] and len(args) >= 4:
        if not os.environ.get("TOY_QUEUE"):
            sys.exit("REFUSED: a slice is a wave of the queue: `python3 tools/autoqueue.py add toy/convergence.py` (his message 147)")
        print(json.dumps(slice_of(args[1], int(args[2]), int(args[3]), float(args[5]) if args[4:5] == ["--for"] else float("inf"))))
    elif args == ["progress"]:
        print(progress())
    elif args == ["digest"]:
        print(digest())
    elif args == ["table"]:
        for name, c in table()["conditions"].items():
            print(name, {n: (p["his"], p["bound"]) for n, p in c["parameters"].items()})
    elif args[:1] == ["plot"] and len(args) == 2:
        print("\n".join(str(p) for p in plot(Path(args[1]))))
    elif args == ["--check"]:
        wrong = check()
        print("\n".join(wrong) or "the convergence of each parameter: made from the files as they are")
        sys.exit(1 if wrong else 0)
    else:
        print(__doc__)
        sys.exit(1)
