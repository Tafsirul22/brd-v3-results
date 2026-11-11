# Stage 1 — Reproducible Installation Instructions

## 1. Python Core & Optional Dependencies (Windows 11 / Linux / macOS)
```powershell
python -m pip install -r 11_post_execution_audit_v3.1/requirements.txt
```
Optional RNA folding/design packages (`ViennaRNA==2.8.0`, local editable `tools/repos/SIREN`, `tools/repos/LUCID`) are automatically guarded if absent during non-folding modules.

## 2. Broccoli v2.0.1 (Native Windows via User-Space Rust GNU Toolchain)
```powershell
$env:RUSTUP_HOME = "$PWD\tools\rustup"
$env:CARGO_HOME = "$PWD\tools\cargo"
$env:PATH = "$env:CARGO_HOME\bin;$env:PATH"
rustup-init.exe --default-host x86_64-pc-windows-gnu --default-toolchain stable --profile minimal -y --no-modify-path
git clone --depth 1 --branch v2.0.1 https://github.com/rderelle/Broccoli.git tools\repos\Broccoli
cd tools\repos\Broccoli && cargo build --release
copy target\release\broccoli.exe ..\..\bin\broccoli.exe
```

## 3. Compleasm v0.2.10 & OrthoDB v12.1 Datasets
```powershell
git clone --depth 1 --branch v0.2.10 https://github.com/huangnengCSU/compleasm.git tools\repos\compleasm
python tools\repos\compleasm\compleasm.py download pleosporales --odb odb12 -L tools\compleasm_db
```
Full `compleasm run` requires Linux (Ubuntu/WSL2/Docker/Conda) providing `miniprot >= 0.13` and `hmmer >= 3.4`.
