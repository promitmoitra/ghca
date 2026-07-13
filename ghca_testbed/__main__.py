"""`python -m ghca_testbed [report|regression]` — reproduce entry point."""

import sys

from .harness import reproduce

mode = sys.argv[1] if len(sys.argv) > 1 else "report"
rep = reproduce(mode=mode)
if mode == "regression" and not rep.ok:
    sys.exit(1)
