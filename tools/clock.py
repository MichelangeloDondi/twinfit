"""The clock of this repository: the owner's time, Europe/Rome, on every machine that runs its tools.

Every tool that writes or reads a stamp of HOURS.tsv, OWNER_WORDS.md, the queue or the gate imports this module first.
The Mac runs in Europe/Rome and the cloud session's container in UTC, and the stamps carry no zone, so a slap of the
cloud read 1085 minutes where 1205 had passed, and a beat on the Mac after a proof of the cloud read a stamp two hours
old (the red team of the cloud session, received with his message 178). A command of the shell reads the same clock
with `TZ=Europe/Rome date`.
"""
from __future__ import annotations

import os
import time

ZONE = "Europe/Rome"


def set_zone() -> str:
    """Sets the zone of this process to the repository's, and returns its offset now, as +0200."""
    os.environ["TZ"] = ZONE
    time.tzset()
    return time.strftime("%z")


set_zone()
