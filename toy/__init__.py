"""The toy model of his line and what is done with it (the owner's messages 111 to 117). `PLAN.md` at the root is the plan.

Every program of the toy imports this package, and so runs at nice 10, with the pools that it starts (his message 184)."""
import sys
from pathlib import Path

if str(Path(__file__).resolve().parents[1] / "tools") not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import niceness  # noqa: E402

niceness.to_ten()
