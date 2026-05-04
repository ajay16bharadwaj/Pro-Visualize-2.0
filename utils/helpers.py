"""
Cross-cutting helpers: error-handling decorator, tab-isolation primitive,
and color-format normalization.
"""

from __future__ import annotations

import logging
import re
import traceback
from functools import wraps
from typing import Callable, Iterable

import streamlit as st

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Decorator for individual plotting functions
# ---------------------------------------------------------------------------


def handle_plotting_errors(plot_function: Callable) -> Callable:
    """Catch exceptions raised inside a plot function and surface them in the UI.

    Used at the visualizer-method level for plots that may fail due to
    bad input (e.g. clustering on too few samples). Errors are logged
    with full traceback; the user sees a single line.
    """

    @wraps(plot_function)
    def wrapper(*args, **kwargs):
        try:
            return plot_function(*args, **kwargs)
        except Exception as e:
            error_message = f"Failed to generate plot: {plot_function.__name__}"
            st.error(f"An error occurred in '{plot_function.__name__}': {e}", icon="💔")
            logger.error(f"{error_message}. Reason: {e}", exc_info=True)
    return wrapper


# ---------------------------------------------------------------------------
# Tab-isolation primitive (P0 → consumed by app.py and PlotManager)
# ---------------------------------------------------------------------------


def safe_render(
    label: str,
    fn: Callable,
    *args,
    reset_keys: Iterable[str] | None = None,
    **kwargs,
):
    """Run `fn(*args, **kwargs)` and contain any exception to this scope.

    On exception:
    - Logs full traceback (to file once logging_config lands).
    - Shows an inline error card with type + message + traceback expander.
    - Offers a "Retry" button (st.rerun()) and an optional "Reset" button
      that pops `reset_keys` from session_state before rerunning.

    This is the primitive used to wrap each top-level tab body so a
    failure in one tab cannot bring down the others.

    Returns the function's return value on success, or None on exception.
    """
    try:
        return fn(*args, **kwargs)
    except Exception as e:
        tb = traceback.format_exc()
        logger.error(f"safe_render('{label}') caught {type(e).__name__}: {e}", exc_info=True)

        st.error(
            f"**{label}** ran into a problem and was paused.\n\n"
            f"`{type(e).__name__}: {e}`\n\n"
            f"Other tabs are unaffected — you can keep working there, "
            f"or use the buttons below to recover.",
            icon="⚠️",
        )

        with st.expander("Show traceback", expanded=False):
            st.code(tb, language="python")

        c1, c2 = st.columns([1, 1])
        # Slugify the label for stable button keys
        slug = re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_") or "safe_render"

        if c1.button("Retry", key=f"{slug}_retry", use_container_width=True):
            st.rerun()

        if reset_keys:
            if c2.button("Reset module state", key=f"{slug}_reset", use_container_width=True):
                for k in reset_keys:
                    st.session_state.pop(k, None)
                st.rerun()
        return None


# ---------------------------------------------------------------------------
# Color normalization
# ---------------------------------------------------------------------------


_HEX_RE = re.compile(r"^#([0-9a-fA-F]{6}|[0-9a-fA-F]{3})$")
_RGB_RE = re.compile(r"rgba?\(\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)")
_NAMED_COLORS = {
    "black": "#000000", "white": "#ffffff", "red": "#ff0000",
    "green": "#008000", "blue": "#0000ff", "yellow": "#ffff00",
    "cyan": "#00ffff", "magenta": "#ff00ff", "grey": "#808080",
    "gray": "#808080", "orange": "#ffa500", "purple": "#800080",
}


def to_hex(color) -> str:
    """Normalize a color to "#rrggbb".

    Accepts:
    - Hex strings: "#RRGGBB" or "#RGB" (passed through / expanded).
    - Plotly RGB strings: "rgb(R, G, B)" or "rgba(R, G, B, A)".
    - (R, G, B) tuples or lists with values in 0-1 (floats) or 0-255 (ints).
    - A small set of named CSS colors.

    Raises ValueError on input it can't parse — callers should treat that
    as a programmer bug, not a runtime concern. Used by global-settings
    color pickers (Quant/Comparative) which previously had a fragile
    inline RGB parser (quant_module.py:114-118).
    """
    if isinstance(color, str):
        s = color.strip().lower()
        if not s:
            raise ValueError("Empty color string")

        # already-hex (3 or 6 chars)
        if _HEX_RE.match(s):
            if len(s) == 4:  # #abc → #aabbcc
                return "#" + "".join(ch * 2 for ch in s[1:])
            return s

        # rgb(...) / rgba(...)
        m = _RGB_RE.match(s)
        if m:
            r, g, b = (float(x) for x in m.groups())
            return _channels_to_hex(r, g, b)

        # named CSS color
        if s in _NAMED_COLORS:
            return _NAMED_COLORS[s]

        raise ValueError(f"Cannot parse color string: {color!r}")

    if isinstance(color, (tuple, list)) and len(color) >= 3:
        r, g, b = color[:3]
        return _channels_to_hex(float(r), float(g), float(b))

    raise ValueError(f"Cannot convert {color!r} to hex (unsupported type)")


def _channels_to_hex(r: float, g: float, b: float) -> str:
    """Convert R/G/B in 0-1 (floats) or 0-255 (anything else) to #rrggbb."""
    if max(r, g, b) <= 1.0:
        r, g, b = r * 255, g * 255, b * 255
    r, g, b = (max(0, min(255, int(round(v)))) for v in (r, g, b))
    return f"#{r:02x}{g:02x}{b:02x}"
