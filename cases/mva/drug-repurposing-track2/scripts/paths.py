"""Shared repository layout. Import from other scripts in this folder."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CURATION = ROOT / "curation"
DATA = ROOT / "data"
DOCS = ROOT / "docs"
REFERENCE = ROOT / "reference"
PDF = ROOT / "pdf"
ASSETS = ROOT / "assets"
