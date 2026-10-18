"""What the design can fix, before any reading of the residuals: the first order bounds, and which block carries what.

    table = identify(lin, block="block")

For each shared parameter: its value, its first order standard deviation, and that deviation as a share of the value. The
correlations. The combinations of the parameters that the records fix best and worst. And, for each block of the design:
how much larger each bound is without the block, and by how much each estimate moves when the curves of that block are
made one per cent wider about their centre, which is its influence. An estimate that moves by many of its standard
deviations for one per cent of one block rests on that block, and the design cannot test that block against itself.
No noise is drawn and nothing is fitted.
"""
from __future__ import annotations

import numpy as np

from .ledger import Ledger

STRETCH = 0.01


def stretched(lin, k: int, share: float = STRETCH) -> np.ndarray:
    """What is added to the mean of record k when its curve is made wider by `share` about its centre."""
    rec = lin.problem.records[k]
    c = lin.centre_of(k)
    wide = lin.problem.predict(k, lin.fit.shared, lin.fit.own[k], x=c + (rec.x - c) / (1.0 + share))
    return wide - lin.problem.predict(k, lin.fit.shared, lin.fit.own[k])


def identify(lin, block: str | None = None) -> dict:
    led, names, values = lin.ledger, lin.problem.shared, lin.fit.shared
    sd = led.sd
    out = {"model": lin.problem.name, "records": len(lin.problem.records), "samples": lin.n, "degrees_of_freedom": lin.dof,
           "parameters": {n: {"value": float(v), "standard_deviation": float(s), "per_cent_of_its_value": float(100.0 * s / abs(v)) if v != 0 else None}
                          for n, v, s in zip(names, values, sd)},
           "correlations": {f"{names[i]}, {names[j]}": float(led.correlation()[i, j]) for i in range(len(names)) for j in range(i + 1, len(names))}}
    rel = np.array([100.0 * s / abs(v) if v != 0 else np.inf for v, s in zip(values, sd)])
    out["weakest"] = names[int(np.argmax(rel))]
    if np.all(values != 0) and len(names) > 1:
        combos = led.combinations(list(names), values)
        out["combinations"] = [{"per_cent": float(100.0 * s), "of": {n: float(c) for n, c in zip(names, vec)}} for s, vec in combos]
        out["best_over_worst"] = float(combos[0][0] / combos[-1][0])
    if block is not None:
        out["blocks"] = {}
        groups = lin.problem.blocks(block)
        for value, members in groups.items():
            keep = np.ones(lin.n, dtype=bool)
            for k in members:
                keep[lin.rows_of(k)] = False
            cell = {"records": len(members)}
            try:
                without = Ledger(lin.J[keep], names).sd
                cell["bound_without_it_over_the_bound"] = {n: float(a / b) for n, a, b in zip(names, without, sd)}
            except (np.linalg.LinAlgError, ValueError):
                cell["bound_without_it_over_the_bound"] = {n: None for n in names}          # without the block the design does not fix every parameter
            d = [stretched(lin, k) if k in members else None for k in range(len(lin.problem.records))]
            move = lin.absorbed(d)
            cell["one_per_cent_wider_moves"] = {n: {"in_its_standard_deviations": float(m / s), "per_cent_of_its_value": float(100.0 * m / abs(v)) if v != 0 else None}
                                                for n, m, s, v in zip(names, move, sd, values)}
            cell["share_of_that_change_that_the_residuals_keep"] = lin.visibility(d)
            out["blocks"][str(value)] = cell
        most = {n: max(out["blocks"], key=lambda b: abs(out["blocks"][b]["one_per_cent_wider_moves"][n]["in_its_standard_deviations"])) for n in names}
        out["rests_on"] = most
    return out
