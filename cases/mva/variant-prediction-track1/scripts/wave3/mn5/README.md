# MN5 jobs for track1-evidence-wave3

Account `ehpc926`, QoS `acc_ehpc`, no memory flags, LF line endings.

From the laptop, using the project interpreter:

```
env\venv\Scripts\python.exe release\track1-evidence-wave3\scripts\dispatch_mn5.py stage
env\venv\Scripts\python.exe release\track1-evidence-wave3\scripts\dispatch_mn5.py submit
env\venv\Scripts\python.exe release\track1-evidence-wave3\scripts\dispatch_mn5.py status
env\venv\Scripts\python.exe release\track1-evidence-wave3\scripts\dispatch_mn5.py collect
```

Remote work: `$MVA_BASE/jobs/track1_wave3/`
Outputs: `$MVA_BASE/out/track1_wave3/`

Laptop receives only TSVs and JSON under 50 MB. BAM, FASTQ and full VCFs stay on MN5.
