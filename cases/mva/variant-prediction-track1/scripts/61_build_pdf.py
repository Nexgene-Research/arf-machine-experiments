"""Build submission/Nexgene-AI_track1_report.pdf from the Markdown report.

Usage: python scripts/61_build_pdf.py
Requires pypandoc_binary (bundled pandoc) and a XeLaTeX installation. The
Markdown stays the single source; this script only adapts typography for LaTeX.
GitHub keeps the Markdown title and author lines; the PDF uses a centred
title block instead.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path

import pypandoc

ROOT = Path(__file__).resolve().parents[1]
SUBMISSION = ROOT / "submission"
SOURCE = SUBMISSION / "Nexgene-AI_track1_report.md"
TARGET = SUBMISSION / "Nexgene-AI_track1_report.pdf"

HEADER = r"""
\usepackage{fontspec}
\setmainfont{Cambria}
\setsansfont{Calibri}
\setmonofont{Consolas}[Scale=0.88]
\usepackage{seqsplit}
\usepackage{xurl}
\usepackage{fvextra}
\fvset{breaklines=true,breakanywhere=true,fontsize=\footnotesize}
\RecustomVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,breakanywhere,fontsize=\footnotesize,commandchars=\\\{\}}
\usepackage{etoolbox}
\usepackage{caption}
\captionsetup{font=small,labelformat=empty,justification=justified,singlelinecheck=false}
\AtBeginEnvironment{longtable}{\footnotesize}
\AtBeginEnvironment{quote}{\small}
\setlength{\emergencystretch}{3em}
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\footnotesize Nexgene AI --- MVA Hackathon 2026, Track 1}
\fancyhead[R]{\footnotesize BUB1B p.[Leu737*];[Asn1002Lys]}
\fancyfoot[C]{\footnotesize\thepage}
\renewcommand{\headrulewidth}{0.3pt}
"""

BEFORE = r"""
\begin{center}
{\LARGE Compound heterozygous \textit{BUB1B} variants p.Leu737* and p.Asn1002Lys identified by whole-genome sequencing in a child with mosaic variegated aneuploidy\par}
\vspace{0.75em}
İlsu Duman, Özge Beğli, Hatice Mutlu and Kutsal Baran Ozkurt$^*$\\[0.35em]
{\small Nexgene AI, https://nexgene.ai\\
Correspondence: kutsal.ozkurt@nexgene.ai\\
MVA Hackathon 2026, Track 1 (variant prediction)\\
24 October 2026}
\end{center}
\vspace{0.4em}
"""

SUPERSCRIPT = str.maketrans("⁻⁰¹²³⁴⁵⁶⁷⁸⁹", "-0123456789")


def latexify(text: str) -> str:
    text = re.sub(
        r"×10([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)",
        lambda m: rf"$\times 10^{{{m.group(1).translate(SUPERSCRIPT)}}}$",
        text,
    )
    text = re.sub(
        r"10([⁻][⁰¹²³⁴⁵⁶⁷⁸⁹]+)",
        lambda m: rf"$10^{{{m.group(1).translate(SUPERSCRIPT)}}}$",
        text,
    )
    text = text.replace("2·p₁·p₂", r"$2 p_1 p_2$")
    text = re.sub(r"`([0-9a-f]{64})`", r"\\texttt{\\seqsplit{\1}}", text)
    # Table cells cannot wrap a code span; a \path breaks at any character.
    return "\n".join(
        re.sub(r"`([^`]*/[^`]*)`", r"\\path{\1}", line) if line.startswith("|") else line
        for line in text.split("\n")
    )


def drop_github_title(text: str) -> str:
    marker = "> **Answer.**"
    index = text.find(marker)
    if index < 0:
        raise ValueError("report is missing the answer block used as the PDF body start")
    return text[index:]


def main() -> None:
    text = latexify(drop_github_title(SOURCE.read_text(encoding="utf-8")))
    with tempfile.TemporaryDirectory() as tmp:
        header = Path(tmp) / "header.tex"
        header.write_text(HEADER, encoding="utf-8")
        before = Path(tmp) / "before.tex"
        before.write_text(BEFORE, encoding="utf-8")
        staged = Path(tmp) / "report.md"
        staged.write_text(text, encoding="utf-8")
        pypandoc.convert_file(
            str(staged),
            "pdf",
            format="markdown+raw_tex",
            outputfile=str(TARGET),
            extra_args=[
                "--pdf-engine=xelatex",
                f"--resource-path={SUBMISSION}",
                "--shift-heading-level-by=-1",
                f"--include-in-header={header}",
                f"--include-before-body={before}",
                "-V", "geometry:a4paper,margin=2.1cm,headheight=14pt",
                "-V", "linkcolor=blue",
                "-V", "urlcolor=blue",
                "-V", "colorlinks=true",
                "-V", "fontsize=10pt",
                "-V", "linestretch=1.12",
            ],
        )
    print(f"wrote {TARGET} ({TARGET.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
