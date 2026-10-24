"""Rebuild bounded Phase L literature outputs; manually curated primary-source claims.

Run: env/venv/Scripts/python.exe scripts/40_phase_l_literature.py
This script performs no network requests and never opens controlled-access data.
The query log records the actual 2026-10-16 retrievals, including failed retrievals.
"""

import csv
import hashlib
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "track1"
DATE = "2026-10-16"
DIMENSIONS = (
    "variant_level",
    "gene_disease",
    "published_cases",
    "protein_structure",
    "constraint_conservation",
    "functional_model",
    "phenotype_concordance",
    "competing_explanations",
)
GENES = ("BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "BUB1")
ALLELES = (
    "BUB1B chr15:40209701:T:G",
    "BUB1B chr15:40216470:A:G",
    "BUB1B chr15:40220612:T:G",
    "BUB1 chr2:110668376:C:T",
    "BUB1 chr2:110671177:G:A",
)

# Each query sent either one gene symbol, one public variant notation, or a DOI/URL
# resolving a primary publication. DOI/PMC URLs are bibliographic read identifiers.
SEARCH_GENES = ("BUB1B", "BUB1", "TRIP13", "CEP57", "MAD1L1", "CENATAC")
READS = (
    ("10.1038/ng1449", "abstract unavailable; no full text"),
    ("10.1073/pnas.1317008111", "abstract only; publisher 403"),
    ("10.1126/sciadv.abk0114", "abstract only; publisher 403"),
    ("https://pmc.ncbi.nlm.nih.gov/articles/PMC8769543/", "primary full-text HTML"),
    ("10.1038/ng.3883", "primary full-text PMC HTML"),
    ("10.1126/sciadv.abq5914", "publisher landing page only"),
    ("https://pmc.ncbi.nlm.nih.gov/articles/PMC9629740/", "primary full-text HTML"),
    ("10.15252/embj.2020106536", "publisher landing page only"),
    ("https://pmc.ncbi.nlm.nih.gov/articles/PMC8280824/", "primary full-text HTML"),
    ("https://pmc.ncbi.nlm.nih.gov/articles/PMC4139549/", "unrelated paper; not cited"),
    ("10.1038/ng.822", "primary PDF full text"),
    ("10.1172/JCI126863", "abstract only; publisher 503"),
    ("10.1056/NEJMoa1006565", "abstract only"),
    (
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC8769543/",
        "second passage retrieval; primary full-text HTML",
    ),
    (
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC9629740/",
        "second passage retrieval; primary full-text HTML",
    ),
    (
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC8280824/",
        "second passage retrieval; primary full-text HTML",
    ),
)
EXA_VARIANTS = (
    ("chr15-40220612-T-G", "irrelevant hits; exact match unconfirmed"),
    ("chr15-40209701-T-G", "irrelevant hits; exact match unconfirmed"),
    ("chr2-110668376-C-T", "irrelevant hits; exact match unconfirmed"),
)

# Row format candidate, dimension, claim, direction, source_id, full_text_read,
# independent. Independent=yes only for a distinct, original primary experiment
# addressing the stated claim. An article repeated for several claims is counted
# once per dimension; a case repeated in a later review is not a new source.
EVIDENCE = [
    (
        "BUB1B",
        "gene_disease",
        "Biallelic BUB1B is associated with constitutional aneuploidy and tumor susceptibility; full paper could not be inspected here",
        "supports",
        "10.1038/ng1449",
        "no",
        "no",
    ),
    (
        "BUB1B",
        "gene_disease",
        "Independent human families segregated an upstream hypomorphic allele with constitutional aneuploidy; edited cultured cells reproduced reduced expression and segregation defects",
        "supports",
        "10.1073/pnas.1317008111",
        "no",
        "yes",
    ),
    (
        "BUB1B",
        "published_cases",
        "Published upstream G>A allele ss802470619 is 44 kb upstream of BUB1B; homozygous infant and multiple family aggregates are reported; this is NOT a match to any candidate allele",
        "neutral",
        "10.1073/pnas.1317008111",
        "no",
        "yes",
    ),
    (
        "BUB1B",
        "functional_model",
        "The published intergenic allele was edited into cultured cells and reduced BUB1B transcripts; this cannot be transferred to an unrelated intronic allele",
        "supports",
        "10.1073/pnas.1317008111",
        "no",
        "yes",
    ),
    (
        "BUB1B",
        "functional_model",
        "Allele-combination experiments in mouse show phenotype differences; mouse sarcopenia is not human skeletal muscle atrophy evidence",
        "supports",
        "10.1172/JCI126863",
        "no",
        "yes",
    ),
    (
        "BUB1B",
        "protein_structure",
        "UniProt annotates N-terminal BUB1-like region aa62-226 and C-terminal kinase region aa766-1050 on canonical O60566; no candidate residue or interface assigned",
        "neutral",
        "UniProt:O60566",
        "no",
        "no",
    ),
    (
        "BUB1B",
        "phenotype_concordance",
        "An MVA literature overview lists rhabdomyosarcoma in BUB1B reports; overview is secondary and case-level tumor counts cannot be independently verified from it",
        "supports",
        "10.1126/sciadv.abq5914",
        "yes",
        "no",
    ),
    (
        "BUB1B",
        "competing_explanations",
        "The primary MAD1L1 case includes embryonal rhabdomyosarcoma, so tumor histology alone does not distinguish BUB1B from MAD1L1",
        "neutral",
        "10.1126/sciadv.abq5914",
        "yes",
        "yes",
    ),
    (
        "CEP57",
        "gene_disease",
        "Four affected children from three reported families (including two siblings) harbor biallelic CEP57 alleles associated with mosaic aneuploidy",
        "supports",
        "10.1038/ng.822",
        "yes",
        "yes",
    ),
    (
        "CEP57",
        "published_cases",
        "Original case table describes growth retardation and mosaic chromosome gains/losses; exact cDNA alleles are not mapped here to GRCh38",
        "supports",
        "10.1038/ng.822",
        "yes",
        "yes",
    ),
    (
        "CEP57",
        "functional_model",
        "Primary human genetics implicates centrosomal microtubule attachment, not a variant-specific phenotype for unrelated alleles",
        "supports",
        "10.1038/ng.822",
        "yes",
        "yes",
    ),
    (
        "CEP57",
        "phenotype_concordance",
        "Growth restriction is reported; nephrocalcinosis and rhabdomyosarcoma are not directly evidenced in the inspected original case table",
        "supports",
        "10.1038/ng.822",
        "yes",
        "yes",
    ),
    (
        "TRIP13",
        "gene_disease",
        "Six individuals with biallelic TRIP13 loss-of-function had Wilms tumor and defective chromosome segregation",
        "supports",
        "10.1038/ng.3883",
        "yes",
        "yes",
    ),
    (
        "TRIP13",
        "published_cases",
        "All six original genotyped cases had Wilms tumor; those cases do not establish nephrocalcinosis as a characteristic feature",
        "supports",
        "10.1038/ng.3883",
        "yes",
        "yes",
    ),
    (
        "TRIP13",
        "functional_model",
        "TRIP13-null patient cells showed impaired spindle checkpoint and chromosome missegregation with rescue by wild-type TRIP13",
        "supports",
        "10.1038/ng.3883",
        "yes",
        "yes",
    ),
    (
        "TRIP13",
        "phenotype_concordance",
        "Published cases have Wilms tumor and variably growth retardation; Wilms tumor and nephrocalcinosis are not interchangeable",
        "neutral",
        "10.1038/ng.3883",
        "yes",
        "yes",
    ),
    (
        "CENATAC",
        "gene_disease",
        "Biallelic CENATAC alleles reported in affected family with mosaic aneuploidy",
        "supports",
        "10.15252/embj.2020106536",
        "yes",
        "yes",
    ),
    (
        "CENATAC",
        "published_cases",
        "In family 5728 maternal c.803A>T and paternal c.805+2_805+19del18 have RNA-validated exon/intron changes; source cDNA changes not mapped here to GRCh38",
        "supports",
        "10.15252/embj.2020106536",
        "yes",
        "yes",
    ),
    (
        "CENATAC",
        "functional_model",
        "Wild-type re-expression rescues minor-intron retention and mitotic errors whereas tested disease mutants do not",
        "supports",
        "10.15252/embj.2020106536",
        "yes",
        "yes",
    ),
    (
        "CENATAC",
        "protein_structure",
        "UniProt Q86UT8 calls CENATAC a minor-spliceosome component; cannot infer candidate protein interface from gene annotation",
        "neutral",
        "UniProt:Q86UT8",
        "no",
        "no",
    ),
    (
        "MAD1L1",
        "gene_disease",
        "Single original compound-heterozygous case with aneuploidy and multiple tumors supports a rare biallelic syndrome; not prevalence evidence",
        "supports",
        "10.1126/sciadv.abq5914",
        "yes",
        "yes",
    ),
    (
        "MAD1L1",
        "published_cases",
        "The original patient's tumors included early embryonal rhabdomyosarcoma of the auditory canal and later other malignancies",
        "supports",
        "10.1126/sciadv.abq5914",
        "yes",
        "yes",
    ),
    (
        "MAD1L1",
        "phenotype_concordance",
        "Exact tumor histology is represented in a reported case, but it is one person, and renal nephrocalcinosis is not established by that report",
        "supports",
        "10.1126/sciadv.abq5914",
        "yes",
        "yes",
    ),
    (
        "BUB1",
        "gene_disease",
        "Two unrelated biallelic BUB1 cases had neurodevelopmental disorder and variable segregation defects; neither establishes classical MVA with childhood tumor predisposition",
        "supports",
        "10.1126/sciadv.abk0114",
        "yes",
        "yes",
    ),
    (
        "BUB1",
        "published_cases",
        "Original cases carry homozygous start-codon change or phased truncating/splice alleles; these differ from the two candidate intronic alleles",
        "neutral",
        "10.1126/sciadv.abk0114",
        "yes",
        "yes",
    ),
    (
        "BUB1",
        "functional_model",
        "Despite segregation and cohesion defects, tested BUB1 patient cells retained an overall functional spindle checkpoint",
        "against",
        "10.1126/sciadv.abk0114",
        "yes",
        "yes",
    ),
    (
        "BUB1",
        "competing_explanations",
        "BUB1B has a different published gene mechanism; analogy to BUB1B cannot convert untested BUB1 intronic alleles into a pathogenic pair",
        "neutral",
        "10.1126/sciadv.abk0114",
        "yes",
        "yes",
    ),
    (
        "BUB1B chr15:40209701:T:G",
        "variant_level",
        "Bounded annotation calls stop-gained; no independently confirmed MANE transcript, NMD assessment, exact published match, or ClinVar two-star assertion in this branch",
        "neutral",
        "local:results/track1_top10_ranked.tsv",
        "no",
        "no",
    ),
    (
        "BUB1B chr15:40209701:T:G",
        "gene_disease",
        "BUB1B is a reported biallelic MVA gene; gene link does not establish function or phase of this stop-gain",
        "supports",
        "10.1038/ng1449",
        "no",
        "no",
    ),
    (
        "BUB1B chr15:40209701:T:G",
        "competing_explanations",
        "A bona fide stop-gain alone is insufficient if the second allele is not demonstrably damaging and in trans",
        "neutral",
        "10.1073/pnas.1317008111",
        "no",
        "yes",
    ),
    (
        "BUB1B chr15:40216470:A:G",
        "variant_level",
        "Intron annotation alone does not prove a splice or expression effect; exact-site Exa result was not obtained",
        "neutral",
        "local:results/track1_top10_ranked.tsv",
        "no",
        "no",
    ),
    (
        "BUB1B chr15:40216470:A:G",
        "functional_model",
        "Published functional validation of a different intergenic BUB1B allele cannot establish an effect at this intronic candidate",
        "neutral",
        "10.1073/pnas.1317008111",
        "no",
        "yes",
    ),
    (
        "BUB1B chr15:40220612:T:G",
        "variant_level",
        "Candidate annotated missense with an uncalibrated prediction; Exa exact-site search yielded irrelevant hits, so exact ClinVar, MANE residue and published match remain unknown",
        "neutral",
        "local:results/track1_top10_ranked.tsv",
        "no",
        "no",
    ),
    (
        "BUB1B chr15:40220612:T:G",
        "protein_structure",
        "Canonical BUB1B O60566 kinase aa766-1050 and N-terminal aa62-226 exist; candidate protein residue has not been mapped, hence no claimed domain/interface disruption",
        "neutral",
        "UniProt:O60566",
        "no",
        "no",
    ),
    (
        "BUB1 chr2:110668376:C:T",
        "variant_level",
        "Intron annotation is insufficient to establish splice impact; exact-site Exa returned unrelated hits; no verified published or two-star ClinVar match",
        "neutral",
        "local:results/track1_top10_ranked.tsv",
        "no",
        "no",
    ),
    (
        "BUB1 chr2:110671177:G:A",
        "variant_level",
        "Second intronic allele not confirmed as functional; exact published or ClinVar match has not been verified",
        "neutral",
        "local:results/track1_top10_ranked.tsv",
        "no",
        "no",
    ),
    (
        "BUB1 chr2:110668376:C:T",
        "competing_explanations",
        "The original biallelic BUB1 cases had experimentally impactful coding or splice variants, unlike this untested intron",
        "neutral",
        "10.1126/sciadv.abk0114",
        "yes",
        "yes",
    ),
    (
        "BUB1 chr2:110671177:G:A",
        "competing_explanations",
        "The original biallelic BUB1 cases had experimentally impactful coding or splice variants, unlike this untested intron",
        "neutral",
        "10.1126/sciadv.abk0114",
        "yes",
        "yes",
    ),
]


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def hash_file(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build():
    OUT.mkdir(exist_ok=True)
    queries = [
        ("literature.search", DATE, g, "one gene; discovery") for g in SEARCH_GENES
    ]
    queries += [("literature.read", DATE, q, result) for q, result in READS]
    queries += [
        (
            "query_uniprot",
            DATE,
            q,
            "one accession or gene; note broad gene-query false hits",
        )
        for q in ("BUB1B", "BUB1", "O60566", "CEP57", "TRIP13", "CENATAC", "MAD1L1")
    ]
    queries += [("exa_web_search_exa", DATE, q, result) for q, result in EXA_VARIANTS]
    write_tsv(
        OUT / "literature_queries.tsv", ("source", "date", "query", "outcome"), queries
    )
    assert all(
        x[1] in DIMENSIONS
        and x[3] in ("supports", "against", "neutral")
        and x[5] in ("yes", "no")
        and x[6] in ("yes", "no")
        for x in EVIDENCE
    )
    write_tsv(
        OUT / "literature_evidence.tsv",
        (
            "candidate",
            "dimension",
            "claim",
            "direction",
            "source_id",
            "full_text_read",
            "independent",
        ),
        EVIDENCE,
    )
    grouped = defaultdict(list)
    for item in EVIDENCE:
        grouped[item[0], item[1]].append(item)
    convergence = []
    for candidate in GENES + ALLELES:
        for dim in DIMENSIONS:
            items = grouped[candidate, dim]
            directions = {r[3] for r in items}
            # A neutral observation cannot provide supportive evidence.
            if "against" in directions and "supports" in directions:
                state = "mixed"
            elif "against" in directions:
                state = "against"
            elif "supports" in directions:
                state = "supports"
            else:
                state = "missing"
            distinct = {
                r[4]
                for r in items
                if r[6] == "yes"
                and r[4].startswith("10.")
                and r[3] in ("supports", "against")
            }
            convergence.append(
                (candidate, dim, state, len(distinct), ";".join(sorted(distinct)))
            )
    write_tsv(
        OUT / "evidence_convergence.tsv",
        (
            "candidate",
            "dimension",
            "status",
            "independent_primary_source_count",
            "independent_primary_source_ids",
        ),
        convergence,
    )
    summary = """# Phase L: bounded primary-literature evidence (2026-10-16)

Rebuild: `env\\venv\\Scripts\\python.exe scripts\\40_phase_l_literature.py`. This is a source-audited qualitative matrix, **not** an allele classification, diagnosis, calibrated LR or new patient read analysis. Primary research papers were read via the literature reader, using PMC HTML when DOI publisher endpoints served abstracts or landing pages. Primary source counts are **distinct DOI IDs within each candidate×dimension**, not mentions, patients, claims or second-hand reviews. Neutral database information never raises the primary-source count. A separate gene-level result **must not** be added to its allele-level row as independent case-specific evidence.

## Priority allele adjudication

- **BUB1B chr15:40209701:T:G**: stop-gain in existing bounded annotation; the transcript-specific impact, NMD, and independent two-star ClinVar/PS1 equivalence remain unverified here. BUB1B biallelic MVA is supported by Hanks et al. DOI `10.1038/ng1449` (reader: no full text/abstract) and regulatory hypomorph data in Ochiai et al. DOI `10.1073/pnas.1317008111` (abstract only). Neither identifies this candidate stop-gain as its published causal allele. Neither proves a second allele is in trans.
- **BUB1B chr15:40220612:T:G**: annotated missense; a model score does not validate its functional effect. Canonical UniProt `O60566` lists BUB1B N-terminal aa62–226 and kinase aa766–1050; no MANE c./p. crosswalk was verified, so **do not assign this missense to either domain**, a KEN box, CDC20/MAD2 binding site or PDB interface. Exact-match ClinVar classification, submitter/review status, LOVD, evolutionary conservation, regional constraint, gnomAD LOEUF and population stratification were not established. Exa search with the single public key produced irrelevant hits, not a negative database result. Phase with the stop-gain remains unverified in this branch.
- **BUB1B chr15:40216470:A:G**: intronic; published cell editing of the *different* BUB1B intergenic hypomorph `ss802470619` cannot justify this allele's splicing or expression effect. The original paper also displays a conflicting `ss804270619` in a figure caption; `results/literature_cases.tsv` preserves the unresolved crosswalk. No equivalence to any candidate allele is shown.
- **BUB1 chr2:110668376:C:T / chr2:110671177:G:A**: both annotated intronic; no verified RNA effect, same-allele primary case or trans phase. The BUB1 patients in Carvalhal et al. DOI `10.1126/sciadv.abk0114` have *different* experimentally affected start, frameshift or splice alleles and an overall **functional** spindle assembly checkpoint; this is a neurodevelopmental syndrome overlapping but not identical to classical BUB1B MVA. Do not extrapolate those coding/splice experiments to either candidate intron.

## Six-gene comparison and phenotype discrimination

The original CEP57 case series DOI `10.1038/ng.822` (full-text PDF) reports biallelic variants, mild growth restriction and chromosome gains/losses. TRIP13 DOI `10.1038/ng.3883` (full PMC HTML) reports six biallelic loss-of-function cases **all with Wilms tumors** and checkpoint rescue experiments. A Wilms tumor is not congenital nephrocalcinosis. MAD1L1 DOI `10.1126/sciadv.abq5914` (full PMC HTML) documents embryonal rhabdomyosarcoma in a *single* biallelic case, so that tumor histology is not uniquely BUB1B. CENATAC DOI `10.15252/embj.2020106536` (full PMC HTML) documents RNA-confirmed minor-intron processing and wild-type versus mutant rescue, not an arbitrary candidate's splice consequence. For BUB1B, DOI `10.1172/JCI126863` is available here only as an abstract of *mouse* allelic-combination experiments: mouse sarcopenia is not case-specific human evidence. The reported phenotype overlap derives from the bounded internal angle B summary, not an external query containing case context. The source papers do not demonstrate nephrocalcinosis as an exact match in the passages checked. Missing clinical findings are unknown, not patient negatives.

## Independence, uncertainty, and follow-up

No candidate-specific primary **functional** demonstration and no independently verified exact allele ClinVar two-star P/LP, PM3 trans case, MANE HGVS, LOVD, protein interface, ClinGen dosage, gene-level LoF constraint or corrected/retraction status were established in this pass. No PS1/PM3/PP3/ACMG LR may be conferred by these tables. The original case table `results/literature_cases.tsv` is **not rewritten**: adding unresolved or duplicate cases would not increase independent evidence. The 2020 CENATAC preprint and its 2021 journal version are a single study. BUB1B Hanks 2004 full text was inaccessible; its DOI and title establish no allele-level interpretation. Corrections/retractions were **not comprehensively checked** and source status remains unknown, not asserted clean. No sequence-based genomic-intelligence service was called. Public queries used only single gene symbols, single public variant notations, or individual bibliographic DOI/PMC identifiers for reading primary works. Follow up by variant-only ClinVar/LOVD requests and transcript-version/assembly REF normalization; only then, if justified, map amino-acid positions to structural coordinates. Couple these results to read/caller validation and phase from the independent MN5 branch before ranking.
"""
    (OUT / "literature_summary.md").write_text(summary, encoding="utf-8", newline="\n")
    provenance = []
    for rel in (
        "scripts/40_phase_l_literature.py",
        "results/track1_top10_ranked.tsv",
        "results/literature_cases.tsv",
        "results/track1/literature_queries.tsv",
        "results/track1/literature_evidence.tsv",
        "results/track1/evidence_convergence.tsv",
        "results/track1/literature_summary.md",
    ):
        provenance.append(
            (
                rel,
                "sha256",
                hash_file(ROOT / rel),
                DATE,
                "local input preserved"
                if rel
                in ("results/literature_cases.tsv", "results/track1_top10_ranked.tsv")
                else "phase_l_curated_rebuild",
            )
        )
    write_tsv(
        OUT / "literature_provenance.tsv",
        ("path", "algorithm", "digest", "access_date", "role"),
        provenance,
    )
    # Fresh-filesystem verification without moving or modifying shared input.
    assert len(convergence) == (len(GENES) + len(ALLELES)) * len(DIMENSIONS)
    for candidate in ALLELES:
        assert not any(
            r[0] == candidate and r[1] == "variant_level" and r[3] == "supports"
            for r in EVIDENCE
        )
    with (OUT / "literature_provenance.tsv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            assert hash_file(ROOT / row["path"]) == row["digest"], row["path"]
    print(
        f"queries={len(queries)} evidence_rows={len(EVIDENCE)} convergence_rows={len(convergence)} preserved_case_file_sha256={hash_file(ROOT / 'results/literature_cases.tsv')}"
    )


if __name__ == "__main__":
    build()
