import os
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Union

import yaml


def load_yaml_config(path: Union[str, Path]) -> Dict[str, Any]:
    """Load YAML file and always return a dict."""
    with open(path, "r", encoding="utf-8") as f:
        loaded = yaml.safe_load(f) or {}

    if not isinstance(loaded, dict):
        raise ValueError(f"YAML root must be a mapping/dict. Got: {type(loaded).__name__}")

    return loaded


def require_config_keys(cfg: Dict[str, Any], required: Iterable[str]) -> None:
    """Raise if any required config keys are missing."""
    missing = [k for k in required if k not in cfg]
    if missing:
        raise ValueError(f"Missing required YAML keys: {missing}")


def resolve_config_path(path_value: Union[str, Path], base_dir: Union[str, Path]) -> Path:
    """
    Resolve a path from config:
    - expands env vars and "~"
    - interprets relative paths against config file directory
    """
    raw = os.path.expandvars(str(path_value).strip())
    candidate = Path(raw).expanduser()
    if candidate.is_absolute():
        return candidate
    return Path(base_dir) / candidate


def load_pipeline_config(
    path: Union[str, Path],
    required_keys: Sequence[str] = (),
    default_values: Optional[Mapping[str, Any]] = None,
    path_keys: Sequence[str] = (),
) -> Dict[str, Any]:
    """
    Standard config loader for pipeline scripts.

    Adds:
    - required key validation
    - default values
    - normalization for selected path-like keys
    - metadata keys: _config_path, _config_dir
    """
    config_path = Path(path).expanduser().resolve()
    cfg = load_yaml_config(config_path)

    if required_keys:
        require_config_keys(cfg, required_keys)

    if default_values:
        for key, value in default_values.items():
            cfg.setdefault(key, value)

    for key in path_keys:
        value = cfg.get(key)
        if value is None:
            continue
        if isinstance(value, str) and not value.strip():
            continue
        cfg[key] = str(resolve_config_path(value, config_path.parent))

    cfg["_config_path"] = str(config_path)
    cfg["_config_dir"] = str(config_path.parent)
    return cfg
