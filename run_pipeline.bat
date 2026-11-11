@echo off
set PYTHONPATH=%~dp0
"%~dp0environments\brd_venv\Scripts\python.exe" "%~dp0BRD_pipeline.py" %*
