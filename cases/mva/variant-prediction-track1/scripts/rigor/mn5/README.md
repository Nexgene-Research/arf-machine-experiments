# MN5 jobs for track1-rigor-pass

Account `ehpc926`, QoS `acc_ehpc`, no memory flags, LF line endings.

From the laptop, using the project interpreter:

```
env\venv\Scripts\python.exe release\track1-rigor-pass\scripts\dispatch_mn5.py stage
env\venv\Scripts\python.exe release\track1-rigor-pass\scripts\dispatch_mn5.py submit
env\venv\Scripts\python.exe release\track1-rigor-pass\scripts\dispatch_mn5.py submit_more
env\venv\Scripts\python.exe release\track1-rigor-pass\scripts\dispatch_mn5.py status
env\venv\Scripts\python.exe release\track1-rigor-pass\scripts\dispatch_mn5.py collect
```

Remote work: `$MVA_BASE/jobs/track1_rigor/`
Outputs: `$MVA_BASE/out/track1_rigor/`

`rigor_pid.sbatch` waits on `rigor_scan.sbatch`. The HIGH/MODERATE key file stays on MN5 and is not collected.

Laptop receives only TSVs and JSON under 50 MB. BAM, FASTQ and full VCFs stay on MN5.
