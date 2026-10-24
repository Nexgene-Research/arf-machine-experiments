# MN5 jobs for the evidence supplement

Account `ehpc926`, QoS `acc_ehpc`, no memory flags, LF line endings.

Upload with `env\venv\Scripts\python.exe release\track1-evidence-supplement\scripts\dispatch_mn5.py stage`.
Submit with `dispatch_mn5.py submit`. Poll with `dispatch_mn5.py status`.
Collect small tables with `dispatch_mn5.py collect`.

Remote work directory: `$MVA_BASE/jobs/track1_supplement/`.
Outputs: `$MVA_BASE/out/track1_supplement/`.

Laptop receives only TSVs and JSON under 50 MB. BAM, FASTQ and full VCFs stay on MN5.
