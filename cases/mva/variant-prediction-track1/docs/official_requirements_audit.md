# Track 1 official requirements audit — 2026-10-16

Sources: [Track 1 official submission instructions and scorer](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/raw/main/tabs/submit_track1.py), [official rules](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/raw/main/tabs/rules.py), [official overview](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/raw/main/tabs/about.py), `data/README.md`, `sub/track1_submission_template.csv`, and `scripts/validate_submission.py`. Web sources retrieved 2026-10-16.

| Requirement | Verified interpretation |
|---|---|
| CSV schema | `proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes`, exact order; second allele blank for singleton. |
| CSV filename | `<team>_<approach>.csv`; `Nexgene-AI_mn5-evidence-v2.csv`. |
| Report filename | `<team>_track1_report.md` (PDF also accepted); `Nexgene-AI_track1_report.md`. |
| Identifier | **The VCF sample column is `WGS_EX2312012`, but the official Track 1 scorer explicitly accepts only `PROBAND01`.** The dataset README names the former as the genomic file/sample. Therefore CSV `proband_id` must be `PROBAND01` to submit. The sample name in read/genotype analyses remains `WGS_EX2312012`. Do not substitute VCF sample name into CSV. |
| Cardinality and ranking | One to ten rows, ranked primary causal hypotheses by descending ePCR; ties retain equal values and are ordered by documented adjudication. |
| Primary/secondary | Primary = proposed cause of the phenotype; secondary = genuinely incidental, excluded from automated primary scoring. Do not mislabel competing causal hypotheses as secondary. |
| Coordinates | GRCh38, chr-prefixed in CSV; REF concordance against project GRCh38 FASTA. |
| ePCR | `(0,1]`, highest first; project preregistration defines a heuristic, not a calibrated clinical probability. |
| AI declaration | Report commercial generative-AI provider, plan/tier and verified data-handling setting in methods description. Do not assert unverified retention settings. |
| Methods form | Complete the supplied `sub/methods_description_form.xlsx` Track 1 sheet for this approach; preserve its layout. |
| Repository | Provide documented reproducible code and configuration; repository may stay private until hackathon close, must be public after. Submission UI requires a `https://github.com/` repository URL, in addition to report and CSV. A local-only repository is not a substitute for that UI field. |
| Closing | Official overview: **2026-10-24 23:59 UTC**, later than internal 2026-10-17 12:00 local target. |

Additional official rules: submissions released CC-BY 4.0; restricted source data may not be redistributed, and must be deleted within 30 days of hackathon close. Rule 7 describes a recorded pitch video for team presentations; the Track 1 upload form asks for repository URL, report, and CSV, not a video upload.

Local validator does **not** check accepted proband ID, name conventions, tie ordering, normalized minimal indel representation, provenance or AI disclosure; release checks cover those that can be checked mechanically.
