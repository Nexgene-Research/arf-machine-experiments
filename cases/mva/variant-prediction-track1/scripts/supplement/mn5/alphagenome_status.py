#!/usr/bin/env python3
"""Record whether the local AlphaGenome package can be imported. No invented scores."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def main() -> None:
    dest = Path(sys.argv[1])
    status = {"importable": False, "parity": None, "indels": "not_run"}
    try:
        import alphagenome  # type: ignore

        status["importable"] = True
        status["version"] = getattr(alphagenome, "__version__", "unknown")
    except Exception as exc:
        status["error"] = f"{type(exc).__name__}: {str(exc)[:300]}"
        status["note"] = "weights are staged; research package not importable here"
    dest.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
    print(status)


if __name__ == "__main__":
    main()
