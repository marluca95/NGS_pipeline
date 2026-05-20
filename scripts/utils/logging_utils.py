from datetime import datetime
import logging
import re
from pathlib import Path
from typing import Optional, Tuple, Union

LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def sanitize_filename_token(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    return cleaned.strip("_") or "unknown"


def build_log_filename(script_name: str, scope: str, run_label: Optional[str] = None) -> str:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    parts = [sanitize_filename_token(script_name), sanitize_filename_token(scope)]
    if run_label:
        parts.append(sanitize_filename_token(run_label))
    parts.append(ts)
    return "__".join(parts) + ".log"


def setup_pipeline_logging(
    logs_dir: Union[str, Path],
    script_name: str,
    scope: str = "run",
    run_label: Optional[str] = None,
    level: int = logging.INFO,
) -> Path:
    log_dir = Path(logs_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / build_log_filename(script_name=script_name, scope=scope, run_label=run_label)

    logging.basicConfig(
        level=level,
        format=LOG_FORMAT,
        datefmt=DATE_FORMAT,
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
        force=True,
    )
    return log_path


def add_logger_file_handler(
    logger: logging.Logger,
    logs_dir: Union[str, Path],
    script_name: str,
    scope: str,
    run_label: Optional[str] = None,
    level: int = logging.INFO,
) -> Tuple[logging.FileHandler, Path]:
    log_dir = Path(logs_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / build_log_filename(script_name=script_name, scope=scope, run_label=run_label)

    handler = logging.FileHandler(log_path)
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter(LOG_FORMAT, datefmt=DATE_FORMAT))
    logger.addHandler(handler)
    return handler, log_path


def remove_logger_handler(logger: logging.Logger, handler: logging.Handler) -> None:
    logger.removeHandler(handler)
    handler.close()
