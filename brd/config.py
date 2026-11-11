"""Configuration loader and path resolution for Bipolaris-RNAi Designer (BRD)."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict

import yaml

DEFAULT_OUTPUT_ROOT = Path(os.environ.get("BRD_OUTPUT_DIR") or Path(__file__).resolve().parent.parent)
DEFAULT_INPUT_ROOT = Path(os.environ.get("BRD_INPUT_DIR") or DEFAULT_OUTPUT_ROOT.parent)
DEFAULT_CONFIG_PATH = DEFAULT_OUTPUT_ROOT / "config" / "project_configuration.yaml"


def load_yaml(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


@dataclass
class BRDConfig:
    """Typed wrapper around project_configuration.yaml with portable environment overrides."""

    raw: Dict[str, Any]
    input_dir: Path = field(default=DEFAULT_INPUT_ROOT)
    output_dir: Path = field(default=DEFAULT_OUTPUT_ROOT)

    @classmethod
    def load(cls, config_path: Path | None = None) -> "BRDConfig":
        cfg_path = Path(config_path) if config_path else DEFAULT_CONFIG_PATH
        raw = load_yaml(cfg_path)
        proj = raw.get("project", {})
        resolved_out = Path(os.environ.get("BRD_OUTPUT_DIR") or cfg_path.resolve().parent.parent)
        raw_in = proj.get("input_dir", "")
        in_dir = Path(os.environ.get("BRD_INPUT_DIR") or (raw_in if raw_in and Path(raw_in).exists() else resolved_out.parent))
        out_dir = Path(resolved_out if resolved_out.exists() else proj.get("output_dir", str(DEFAULT_OUTPUT_ROOT)))
        return cls(
            raw=raw,
            input_dir=in_dir,
            output_dir=out_dir,
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
