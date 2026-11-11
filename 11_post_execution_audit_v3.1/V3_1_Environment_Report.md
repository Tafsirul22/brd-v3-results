# BRD v3.1 Environment & Dependency Validation Report

**Checked:** 2026-11-11 08:26:10 UTC  
**OS / Platform:** `Windows-11-10.0.26300-SP0` (`nt`)  
**Python:** `3.13.14` (`C:\Users\ASUS BANGLADESH\Desktop\New folder (2)\Output\environments\brd_venv\Scripts\python.exe`)  
**Default Encoding / Temp UTF-8 LF Safety:** `utf-8` / `tmp_utf8_lf_ok=True`  
**Mandatory Python Packages Missing:** `0`

## Complete Dependency Matrix (`tool_availability_v3_1.tsv`)
| component   | kind            | mandatory   | status                             | version                                           |
|:------------|:----------------|:------------|:-----------------------------------|:--------------------------------------------------|
| biopython   | python_package  | True        | AVAILABLE                          | 1.86                                              |
| pandas      | python_package  | True        | AVAILABLE                          | 2.3.3                                             |
| numpy       | python_package  | True        | AVAILABLE                          | 2.4.0                                             |
| scipy       | python_package  | True        | AVAILABLE                          | 1.16.3                                            |
| matplotlib  | python_package  | True        | AVAILABLE                          | 3.10.8                                            |
| networkx    | python_package  | True        | AVAILABLE                          | 3.6.1                                             |
| PyYAML      | python_package  | True        | AVAILABLE                          | 6.0.3                                             |
| tabulate    | python_package  | True        | AVAILABLE                          | 0.10.0                                            |
| reportlab   | python_package  | True        | AVAILABLE                          | 4.4.7                                             |
| pytest      | python_package  | True        | AVAILABLE                          | 9.1.1                                             |
| ViennaRNA   | python_package  | False       | AVAILABLE                          | 2.8.0                                             |
| primer3-py  | python_package  | False       | AVAILABLE                          | 2.4.0                                             |
| broccoli    | external_binary | False       | AVAILABLE                          | Broccoli v2.0.1                                   |
| diamond     | external_binary | False       | AVAILABLE                          | diamond version 2.2.8                             |
| blastn      | external_binary | False       | AVAILABLE                          | blastn: 2.17.0+                                   |
| makeblastdb | external_binary | False       | AVAILABLE                          | makeblastdb: 2.17.0+                              |
| datasets    | external_binary | False       | AVAILABLE                          | datasets version: 18.41.0                         |
| miniprot    | external_binary | False       | UNAVAILABLE_OPTIONAL_OR_LINUX_ONLY | Linux x64 binary required by compleasm            |
| hmmsearch   | external_binary | False       | UNAVAILABLE_OPTIONAL_OR_LINUX_ONLY | HMMER 3.4 POSIX binary required by compleasm      |
| FastTree    | external_binary | False       | UNAVAILABLE_OPTIONAL_OR_LINUX_ONLY | Optional (Broccoli uses built-in -phylogenies nj) |
| docker      | external_binary | False       | UNAVAILABLE_OPTIONAL_OR_LINUX_ONLY | Docker Desktop unavailable on host                |
| quast       | external_binary | False       | UNAVAILABLE_OPTIONAL_OR_LINUX_ONLY | Linux/WSL2 only                                   |
| busco       | external_binary | False       | UNAVAILABLE_OPTIONAL_OR_LINUX_ONLY | Linux/WSL2 only                                   |
| orthofinder | external_binary | False       | UNAVAILABLE_OPTIONAL_OR_LINUX_ONLY | Linux/WSL2 only                                   |
| RNAhybrid   | external_binary | False       | UNAVAILABLE_OPTIONAL_OR_LINUX_ONLY | Linux/macOS Bioconda only                         |
