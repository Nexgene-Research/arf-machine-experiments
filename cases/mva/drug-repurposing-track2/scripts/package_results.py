"""Bundle only final hash-verified deliverables. No papers or patient files are swept."""

import hashlib
import json
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import DATA, ROOT


def main():
    validation = json.loads((DATA / "validation.json").read_text(encoding="utf-8"))
    assert validation["status"] == "passed"
    for name, expected in validation["files"].items():
        payload = (ROOT / name).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == expected["sha256"], name
    archive = ROOT / "research_package.zip"
    with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED) as handle:
        for name in sorted(validation["files"]):
            handle.write(ROOT / name, name)
        handle.write(DATA / "validation.json", "data/validation.json")
    with zipfile.ZipFile(archive) as handle:
        assert handle.testzip() is None
        assert set(handle.namelist()) == set(validation["files"]) | {
            "data/validation.json"
        }
    print(
        json.dumps(
            {
                "archive": archive.name,
                "bytes": archive.stat().st_size,
                "sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
