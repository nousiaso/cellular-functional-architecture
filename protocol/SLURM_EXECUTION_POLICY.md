# SLURM execution policy (proposed)

Use short interactive compute allocations only for focused diagnostics. Submit substantial processing as batch jobs from repository root. Create logs/ before submission. The generic Python wrapper uses existing Python on PATH, or CFA_PYTHON explicitly, and records git HEAD, source SHA256, resource request, timing and exit status. Stdout and stderr use separate job-ID logs. No cluster-specific absolute paths are committed.

Example:
```bash
mkdir -p logs
sbatch --job-name=cfa_p1a04t --cpus-per-task=1 --mem=8G --time=03:00:00 scripts/slurm/run_python_stage.sbatch scripts/diagnostics/p1a04t_full_sbml_export_comparison.py
```
Resource requests are provisional. This diagnostic is serial: extra CPUs do not accelerate it. Monitor actual usage with sstat/sacct and progress logs. Do not launch a duplicate of a running stage unnecessarily or modify source during execution. Record completed output before replacing a run.

Current stages require only Python standard library; no additional environment. Future dependencies should be version-pinned and minimal. Given limited home quota, prefer existing compatible environments or a temporary Python venv under an approved job-local scratch location (if permitted by HPC policy); build once per job/worker group, not once per task; clean up only after job exits. No conda installation elsewhere without permission. Do not delete an environment while a queued/running job needs it.

Parallelization later: partition independent archive/model operations only after profiling; stream compressed archive once rather than have every worker decompress the same archive, then bounded worker pool and deterministic reduction; for independent nulls/bootstrap use job arrays with deterministic per-task seeds and separate outputs. Preserve evidence and signature deduplication across worker boundaries. No numerical or scientific model change in this commit.