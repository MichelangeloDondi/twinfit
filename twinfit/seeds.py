"""The seed of every Monte Carlo of the package, drawn from the records it reads: row M582, the owner's message 527 on
his 526 and decision 35. A fixed default seed gave every set of records the same vectors of noise, so a share read over
many sets was the share of that one set of vectors: on the toy's right models the report alarmed on 0.188 of 250 draws
at 30 sets with the seed 23 and on 0.068 with a seed for each set of records (rows M564 to M576 of the analysis
repository). Each function that draws keeps its old seed as the salt beside the digest, so two functions read on the same
records draw apart, and the same records give the same draws each time."""
from __future__ import annotations

import hashlib

import numpy as np


def seed_of(problem, salt: int) -> int:
    """The seed drawn from the records of `problem`: the digest of each record's abscissa, values and standard
    deviations, with `salt` beside it."""
    h = hashlib.sha256(int(salt).to_bytes(8, "big"))
    for rec in problem.records:
        for a in (rec.x, rec.y, rec.sd):
            h.update(np.ascontiguousarray(a, dtype=float).tobytes())
    return int.from_bytes(h.digest()[:8], "big")
