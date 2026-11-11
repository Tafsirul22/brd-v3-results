# Bipolaris-RNAi Designer (BRD) v3.1.1

Modular, isolate-aware computational bioinformatics pipeline for identifying, verifying, designing, off-target screening, and uncertainty-aware ranking of candidate dsRNA triggers against *Bipolaris sorokiniana* (syn. *Cochliobolus sativus*). See `11_post_execution_audit_v3.1.1/` for the latest audited reproducibility release.

## Scientific Boundary

BRD generates **computationally prioritised candidate dsRNA sequences** with verified multi-isolate nucleotide conservation, functional genetics, infection-stage expression (`GSE130517`), and multi-reference off-target profiles. It does **not** demonstrate exogenous dsRNA uptake, transcript knockdown, or fungicidal efficacy in *B. sorokiniana*, which remain experimental hypotheses.

## Environment Setup

```powershell
# Mandatory core environment (imports work cleanly without optional primer3/ViennaRNA)
pip install -r requirements.txt

# Optional thermodynamic folding, PCR primer design, and Snakemake support
pip install -r requirements-optional.txt
```

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
