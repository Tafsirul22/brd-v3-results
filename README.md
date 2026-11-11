# Bipolaris-RNAi Designer (BRD) v1.0.0

Modular, isolate-aware computational bioinformatics pipeline for identifying, verifying, designing, off-target screening, and uncertainty-aware ranking of candidate dsRNA triggers against *Bipolaris sorokiniana* (syn. *Cochliobolus sativus*).

## Scientific Boundary

BRD generates **computationally prioritised candidate dsRNA sequences** with verified multi-isolate nucleotide conservation, functional genetics, infection-stage expression (`GSE130517`), and multi-reference off-target profiles. It does **not** demonstrate exogenous dsRNA uptake, transcript knockdown, or fungicidal efficacy in *B. sorokiniana*, which remain experimental hypotheses.

## Quick Start

```powershell
# Check status and dependencies
.\run_pipeline.ps1 --status

# Run full M01–M10 pipeline
.\run_pipeline.ps1 --run all

# Run an individual module (e.g., M04 or M09)
.\run_pipeline.ps1 --run m04

# Resume from completed checkpoints
.\run_pipeline.ps1 --resume

# Validate Snakemake DAG
.\environments\brd_venv\Scripts\snakemake.exe --dry-run
```

## Testing

```powershell
.\environments\brd_venv\Scripts\pytest.exe tests -v
```
