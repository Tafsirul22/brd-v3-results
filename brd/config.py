"""Configuration loader and path resolution for Bipolaris-RNAi Designer (BRD)."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict

import yaml

DEFAULT_OUTPUT_ROOT = Path(r"C:\Users\ASUS BANGLADESH\Desktop\New folder (2)\Output")
DEFAULT_INPUT_ROOT = Path(r"C:\Users\ASUS BANGLADESH\Desktop\New folder (2)")
DEFAULT_CONFIG_PATH = DEFAULT_OUTPUT_ROOT / "config" / "project_configuration.yaml"


def load_yaml(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


@dataclass
class BRDConfig:
    """Typed wrapper around project_configuration.yaml with environment overrides."""

    raw: Dict[str, Any]
    input_dir: Path = field(default=DEFAULT_INPUT_ROOT)
    output_dir: Path = field(default=DEFAULT_OUTPUT_ROOT)

    @classmethod
    def load(cls, config_path: Path | None = None) -> "BRDConfig":
        cfg_path = Path(config_path) if config_path else DEFAULT_CONFIG_PATH
        raw = load_yaml(cfg_path)
        proj = raw.get("project", {})
        return cls(
            raw=raw,
            input_dir=Path(proj.get("input_dir", str(DEFAULT_INPUT_ROOT))),
            output_dir=Path(proj.get("output_dir", str(DEFAULT_OUTPUT_ROOT))),
        )

    @property
    def ncbi_email(self) -> str:
        return os.environ.get("NCBI_EMAIL") or os.environ.get("USER_EMAIL") or self.raw.get("ncbi", {}).get(
            "default_email", "brd_pipeline_research@example.org"
        )

    @property
    def ncbi_api_key(self) -> str | None:
        return os.environ.get("NCBI_API_KEY") or None

    def tool_path(self, key: str) -> str:
        tools = self.raw.get("tools", {})
        val = tools.get(key, "")
        if val and Path(val).exists():
            return str(val)
        return val or key

    def module_dir(self, name: str) -> Path:
        p = self.output_dir / name
        p.mkdir(parents=True, exist_ok=True)
        return p
