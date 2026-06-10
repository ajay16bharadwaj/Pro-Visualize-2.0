"""Centralized logging setup for Pro-Visualize.

Modules across the app use ``logging.getLogger(__name__)`` but nothing
configured the root logger, so their records went nowhere with a useful
format. ``configure_logging()`` is called once from ``app.py`` at startup.

Behavior:
- Always logs to stderr (captured by Streamlit / Docker / journald).
- Optionally adds a rotating file handler under ``~/.pro_visualize/logs/``
  when ``PROVIZ_LOG_TO_FILE=1`` — useful for debugging a deployed instance.
- Level is controlled by ``PROVIZ_LOG_LEVEL`` (default ``INFO``).
- Idempotent: safe to call on every Streamlit rerun.
"""

from __future__ import annotations

import logging
import os
from logging.handlers import RotatingFileHandler
from pathlib import Path

_CONFIGURED_FLAG = "_proviz_logging_configured"
_FORMAT = "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s"


def configure_logging() -> None:
    """Configure the root logger once per process."""
    root = logging.getLogger()
    if getattr(root, _CONFIGURED_FLAG, False):
        return

    level_name = os.environ.get("PROVIZ_LOG_LEVEL", "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)
    root.setLevel(level)

    formatter = logging.Formatter(_FORMAT)

    stream = logging.StreamHandler()
    stream.setFormatter(formatter)
    root.addHandler(stream)

    if os.environ.get("PROVIZ_LOG_TO_FILE") == "1":
        try:
            log_dir = Path.home() / ".pro_visualize" / "logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            file_handler = RotatingFileHandler(
                log_dir / "pro_visualize.log", maxBytes=5_000_000, backupCount=3
            )
            file_handler.setFormatter(formatter)
            root.addHandler(file_handler)
            root.info("File logging enabled at %s", log_dir)
        except Exception as e:  # never let logging setup break the app
            root.warning("Could not set up file logging: %s", e)

    # Quiet noisy third-party loggers.
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)

    setattr(root, _CONFIGURED_FLAG, True)
