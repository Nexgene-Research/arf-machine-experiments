"""Render a Markdown deliverable to an academically typeset PDF.

Usage: python scripts/render_pdf.py <input.md> <output.pdf>

Uses the `markdown` library for HTML and `xhtml2pdf` (pure Python) for PDF so it
runs on Windows without LaTeX. Local links to package files are kept as text.
Times New Roman is used when present under assets/fonts/; otherwise the PDF
falls back to the built-in Times family.
"""

import re
import sys
from pathlib import Path

import markdown
from xhtml2pdf import pisa

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import ASSETS, DOCS, PDF, ROOT

FONTS = ASSETS / "fonts"


def _font_face() -> str:
    files = {
        "regular": FONTS / "times.ttf",
        "bold": FONTS / "timesbd.ttf",
        "italic": FONTS / "timesi.ttf",
        "bolditalic": FONTS / "timesbi.ttf",
    }
    if not all(p.is_file() for p in files.values()):
        return ""
    return f"""
@font-face {{ font-family: "Serif"; src: url("{files['regular'].as_posix()}"); }}
@font-face {{ font-family: "Serif"; font-weight: bold; src: url("{files['bold'].as_posix()}"); }}
@font-face {{ font-family: "Serif"; font-style: italic; src: url("{files['italic'].as_posix()}"); }}
@font-face {{ font-family: "Serif"; font-weight: bold; font-style: italic; src: url("{files['bolditalic'].as_posix()}"); }}
"""


CSS = _font_face() + """
@page { size: A4; margin: 2.0cm 1.8cm 2.1cm 1.8cm; @frame footer { -pdf-frame-content: footer; bottom: 1.0cm; margin-left: 1.8cm; margin-right: 1.8cm; height: 1cm; } }
body { font-family: "Serif", Times, serif; font-size: 10.5pt; line-height: 1.42; color: #111111; }
h1 { font-size: 18.5pt; text-align: center; margin: 0 0 9pt 0; line-height: 1.22; color: #0f1b2d; }
h2 { font-size: 13pt; border-bottom: 1pt solid #6b7785; padding-bottom: 2pt; margin-top: 15pt; margin-bottom: 4pt; color: #11263f; -pdf-keep-with-next: true; }
h3 { font-size: 11.5pt; margin-top: 11pt; margin-bottom: 3pt; color: #16314f; -pdf-keep-with-next: true; }
h4 { font-size: 10.5pt; margin-top: 8pt; margin-bottom: 2pt; -pdf-keep-with-next: true; }
p { margin: 4pt 0; text-align: justify; }
ol, ul { margin: 4pt 0 4pt 0; }
li { margin: 2pt 0; text-align: justify; }
a { color: #1a3c66; text-decoration: none; }
code { font-family: "Courier New", monospace; font-size: 9.5pt; }
pre { background-color: #f6f6f6; border: 0.5pt solid #ccc; padding: 4pt 6pt; font-size: 9pt; }
blockquote { background-color: #fbf1ee; border-left: 4px solid #b4512f; padding: 5pt 9pt; margin: 7pt 0; font-size: 9.8pt; }
table { width: 100%; border-collapse: collapse; margin: 8pt 0; font-size: 9pt; }
thead { display: table-header-group; }
th { background-color: #11263f; color: #ffffff; border: 0.5pt solid #11263f; padding: 4pt 5pt; text-align: left; font-size: 9pt; }
td { border: 0.5pt solid #c2cbd4; padding: 4pt 5pt; vertical-align: top; }
tr.alt td { background-color: #eef2f6; }
.card { border: 0.5pt solid #c2cbd4; border-left: 3pt solid #11263f; padding: 6pt 8pt; margin: 7pt 0 8pt 0; }
.card-title { font-weight: bold; color: #11263f; font-size: 10.5pt; margin: 0 0 4pt 0; }
.card p { margin: 2pt 0; text-align: left; }
.card b { color: #16314f; }
hr { border: none; border-top: 0.5pt solid #bbb; margin: 11pt 0; }
.titleblock { margin: 2pt 0 8pt 0; }
.titleblock p { text-align: center; margin: 2pt 0; font-size: 10.5pt; color: #333333; }
.authors { font-size: 11.5pt; color: #11263f; }
"""

FOOTER = '<div id="footer" style="text-align:center;font-size:8pt;color:#777;">Nexgene AI · MVA Hackathon 2026, Track 2 · <pdf:pagenumber>/<pdf:pagecount></div>'


def _cell_text(cell: str) -> str:
    return re.sub(r"\s+", " ", cell).strip()


def _row_cells(row: str) -> list[str]:
    return [_cell_text(c) for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row, flags=re.DOTALL)]


def _tables_for_pdf(html: str) -> str:
    """Keep short tables; turn wide or long tables into readable stacked cards.

    xhtml2pdf's keep-in-frame shrink mode was collapsing tables at page breaks
    into unreadably small thumbnails.
    """

    def replace(match: re.Match) -> str:
        table = match.group(0)
        rows = re.findall(r"<tr[^>]*>.*?</tr>", table, flags=re.DOTALL)
        if not rows:
            return table
        parsed = [_row_cells(r) for r in rows]
        parsed = [r for r in parsed if r]
        if not parsed:
            return table
        headers = parsed[0]
        body = parsed[1:]
        ncols = max(len(r) for r in parsed)
        longest = max((len(re.sub("<[^>]+>", "", c)) for r in parsed for c in r), default=0)
        keep_as_table = ncols <= 3 and len(body) <= 6 and longest <= 90
        if keep_as_table:
            shaded = []
            for i, row in enumerate(rows):
                if i == 0 or i % 2 == 1:
                    shaded.append(row)
                else:
                    shaded.append(row.replace("<tr>", '<tr class="alt">', 1))
            return "<table>" + "".join(shaded) + "</table>"
        cards = []
        for item in body:
            padded = item + [""] * (len(headers) - len(item))
            title = padded[0] or "Item"
            fields = []
            for header, value in zip(headers[1:], padded[1:]):
                if not value:
                    continue
                fields.append(f"<p><b>{header}.</b> {value}</p>")
            cards.append(
                '<div class="card"><p class="card-title">'
                + title
                + "</p>"
                + "".join(fields)
                + "</div>"
            )
        return "".join(cards)

    return re.sub(r"<table>.*?</table>", replace, html, flags=re.DOTALL)


def convert(md_path: Path, pdf_path: Path) -> None:
    text = md_path.read_text(encoding="utf-8")
    text = re.sub(r"\[([^\]]+)\]\((?!https?://)[^)]+\)", r"\1", text)
    text = text.replace("\u207a", "+")
    html_body = markdown.markdown(
        text, extensions=["tables", "fenced_code", "sane_lists"]
    )

    def _titleblock(m: re.Match) -> str:
        block = re.sub(
            r"<p>(.*?)</p>",
            r'<p class="authors">\1</p>',
            m.group(2),
            count=1,
            flags=re.DOTALL,
        )
        return m.group(1) + '<div class="titleblock">' + block + "</div>" + m.group(3)

    html_body = re.sub(
        r"(</h1>)(.*?)(<hr\s*/?>)",
        _titleblock,
        html_body,
        count=1,
        flags=re.DOTALL,
    )
    html_body = _tables_for_pdf(html_body)
    html = f"<html><head><meta charset='utf-8'><style>{CSS}</style></head><body>{FOOTER}{html_body}</body></html>"
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    with pdf_path.open("wb") as handle:
        result = pisa.CreatePDF(html, dest=handle, encoding="utf-8")
    if result.err:
        raise SystemExit(f"PDF render failed for {md_path.name}: {result.err} error(s)")
    print(
        f"{md_path.relative_to(ROOT)} -> {pdf_path.relative_to(ROOT)} ({pdf_path.stat().st_size} bytes)"
    )


def main() -> None:
    args = sys.argv[1:]
    if len(args) >= 2:
        convert(Path(args[0]).resolve(), Path(args[1]).resolve())
        return
    targets = {
        DOCS / "report.md": PDF / "NexgeneAI-MVA-Track2-Report.pdf",
        DOCS / "appendices.md": PDF / "NexgeneAI-MVA-Track2-Appendices.pdf",
        ROOT / "data" / "graph.md": PDF / "NexgeneAI-MVA-Track2-Graph.pdf",
        ROOT / "README.md": PDF / "NexgeneAI-MVA-Track2-README.pdf",
        DOCS / "claim_audit.md": PDF / "NexgeneAI-MVA-Track2-ClaimAudit.pdf",
    }
    PDF.mkdir(exist_ok=True)
    for src, dst in targets.items():
        convert(src, dst)


if __name__ == "__main__":
    main()
