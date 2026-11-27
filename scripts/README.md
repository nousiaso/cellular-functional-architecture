# Scripts

`staging/p1b01_source_scoped_harmonization.py --execute` passed via `slurm/run_python_stage.sbatch`. Preserve its output. Next: `staging/p1b02_molecular_compatibility.py --execute` via the SLURM wrapper (1 CPU, 4 GiB, 1 hour provisional). Factored molecular candidates, not graph edges; execution pending. Standard Python only; output is atomic and no overwrite.

Retained acquisition and staging scripts reconstruct source records and derived tables. Diagnostic modules retained where imported by active scripts or needed for next provenance materialization. Some contain standalone audits as well as shared helpers; keep their paths stable. Do not rerun acquisition against mutable endpoints to reconstruct the frozen snapshot: use existing verified receipts. The historical candidate builder and one-off certification diagnostics are recoverable from pre-cleanup commit 2c7a2dfa2fdee19d53b382adc026080bd0809365 in a separate checkout; the approved source lock is not regenerated.

No additional inference scripts until the preceding stage passes. Keep code, manifests and concise documentation in Git; raw data, outputs and logs local.