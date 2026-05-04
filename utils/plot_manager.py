"""
PlotManager — owns the lifecycle of a single editable figure in Streamlit.

Two flavors:
- PlotManager — for Plotly figures (interactive). Supports per-figure
  editing (title/labels/marker size/height/legend/font/axis ranges),
  static export to PNG/SVG/HTML via kaleido, and an "Add to Report"
  button (P1 wires this to the report builder).
- MplPlotManager — for matplotlib outputs returned as BytesIO PNGs.
  Same UI shell (display + download + add-to-report) but no in-place
  editing — matplotlib figures need re-rendering for edits, which is
  forwarded to the plot function as kwargs (title/figsize/dpi).

Edits persist across re-generation. Each PlotManager stores edits in
`st.session_state[f"{key}_edits"]`, tagged with a hash of the params
that produced the figure. When params change (the user clicks "Update
Plot" with different inputs), edits whose hash no longer matches are
discarded so they don't apply to a different plot.
"""

from __future__ import annotations

import hashlib
import io
import json
import logging
import re
from typing import Any, Callable

import streamlit as st

try:
    from config.plot_configs import EXPORT_DEFAULTS
except Exception:  # pragma: no cover — keeps import robust if config moves
    EXPORT_DEFAULTS = {"png_scale": 2, "html_plotlyjs": "cdn", "matplotlib_dpi": 150}

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal: params hash and edit-state helpers
# ---------------------------------------------------------------------------


def _params_hash(params: dict) -> str:
    """Stable short hash of plot params, used to invalidate stale edits."""
    try:
        s = json.dumps(params, sort_keys=True, default=str)
    except Exception:
        s = repr(sorted(params.items(), key=lambda kv: str(kv[0])))
    return hashlib.md5(s.encode("utf-8")).hexdigest()[:12]


def _slug(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_") or "plot"


# ---------------------------------------------------------------------------
# PlotManager — Plotly figures
# ---------------------------------------------------------------------------


class PlotManager:
    """Manages state, generation, editing, and export of one Plotly figure."""

    def __init__(self, plot_key: str):
        self.key = plot_key
        self.params_key = f"{plot_key}_params"
        self.fig_key = f"{plot_key}_fig"
        self.edits_key = f"{plot_key}_edits"  # persisted edit dict
        self.hash_key = f"{plot_key}_phash"   # params hash that owns the edits

    # ---- generation -------------------------------------------------------

    def render_generate_button(self, plot_function: Callable, **kwargs):
        """Render the Generate/Update button. Auto-generates on first render."""
        current_fig = st.session_state.get(self.fig_key)

        if current_fig is None:
            self._generate(plot_function, kwargs, silent=True)

        label = "Update Plot" if st.session_state.get(self.fig_key) else "Generate Plot"
        if st.button(label, key=f"{self.key}_btn", use_container_width=True, type="primary"):
            self._generate(plot_function, kwargs, silent=False)

    def _generate(self, plot_function: Callable, kwargs: dict, *, silent: bool):
        """Run the plot function, store result, and re-apply persisted edits."""
        new_hash = _params_hash(kwargs)
        prev_hash = st.session_state.get(self.hash_key)

        try:
            if silent:
                fig = plot_function(**kwargs)
            else:
                with st.spinner("Generating plot..."):
                    fig = plot_function(**kwargs)
            st.session_state[self.fig_key] = fig
            st.session_state[self.params_key] = kwargs
            st.session_state[self.hash_key] = new_hash

            # If params changed, drop stale edits — otherwise re-apply them.
            if prev_hash != new_hash:
                st.session_state.pop(self.edits_key, None)
            else:
                self._apply_persisted_edits(fig)
        except Exception as e:
            if not silent:
                st.error(f"Error generating plot: {e}")
                logger.error(f"Plot generation failed for {self.key}: {e}", exc_info=True)
            st.session_state[self.fig_key] = None

    def _apply_persisted_edits(self, fig):
        edits = st.session_state.get(self.edits_key)
        if not edits or fig is None:
            return
        try:
            layout_updates = {}
            if "title" in edits:
                layout_updates["title_text"] = edits["title"]
            if "height" in edits:
                layout_updates["height"] = edits["height"]
            if "xlabel" in edits:
                layout_updates["xaxis_title"] = edits["xlabel"]
            if "ylabel" in edits:
                layout_updates["yaxis_title"] = edits["ylabel"]
            if "font_size" in edits:
                layout_updates["font"] = dict(size=edits["font_size"])
            if "legend_position" in edits:
                layout_updates["legend"] = _legend_anchor(edits["legend_position"])
            if "show_grid" in edits:
                layout_updates["xaxis"] = dict(showgrid=edits["show_grid"])
                layout_updates["yaxis"] = dict(showgrid=edits["show_grid"])
            if layout_updates:
                fig.update_layout(**layout_updates)
            if "marker_size" in edits:
                try:
                    fig.update_traces(marker=dict(size=edits["marker_size"]))
                except Exception:
                    pass
        except Exception as e:
            logger.warning(f"Could not re-apply edits for {self.key}: {e}")

    # ---- display + editor + export ---------------------------------------

    def render_plot_and_editor(self):
        """Display the current figure, expose editing controls, and export buttons."""
        fig = st.session_state.get(self.fig_key)

        if fig is None:
            st.info("👆 Configure parameters above and click **Generate Plot** to visualize results.")
            return

        plot_placeholder = st.empty()

        self._render_editor(fig)
        self._render_export_row(fig)

        plot_placeholder.plotly_chart(fig, use_container_width=True, theme=None)

    def _render_editor(self, fig):
        with st.expander("🎨 Edit Plot Appearance", expanded=False):
            current = self._read_current_layout(fig)

            c1, c2, c3 = st.columns(3)
            new_title = c1.text_input("Plot Title", value=current["title"], key=f"{self.key}_title")
            new_height = c2.slider(
                "Plot Height", 300, 1200, current["height"], 50, key=f"{self.key}_height",
            )
            marker_size = c3.slider(
                "Marker Size", 2, 20, current["marker_size"], key=f"{self.key}_msize",
            )

            c4, c5 = st.columns(2)
            new_xlabel = c4.text_input("X-Axis Label", value=current["xlabel"], key=f"{self.key}_xlabel")
            new_ylabel = c5.text_input("Y-Axis Label", value=current["ylabel"], key=f"{self.key}_ylabel")

            c6, c7, c8 = st.columns(3)
            font_size = c6.slider("Font Size", 8, 24, current["font_size"], key=f"{self.key}_fsize")
            legend_pos = c7.selectbox(
                "Legend Position",
                ["right", "top", "bottom", "left", "hidden"],
                index=["right", "top", "bottom", "left", "hidden"].index(current["legend_position"]),
                key=f"{self.key}_legend",
            )
            show_grid = c8.checkbox("Show Gridlines", value=current["show_grid"], key=f"{self.key}_grid")

            if st.button("Apply Visual Changes", key=f"{self.key}_apply"):
                edits = {
                    "title": new_title,
                    "height": new_height,
                    "marker_size": marker_size,
                    "xlabel": new_xlabel,
                    "ylabel": new_ylabel,
                    "font_size": font_size,
                    "legend_position": legend_pos,
                    "show_grid": show_grid,
                }
                st.session_state[self.edits_key] = edits
                self._apply_persisted_edits(fig)

    @staticmethod
    def _read_current_layout(fig) -> dict:
        """Best-effort extraction of current layout values for editor defaults."""
        title = ""
        height = 500
        xlabel = ""
        ylabel = ""
        font_size = 12
        legend_position = "right"
        show_grid = True
        marker_size = 8

        try:
            if fig.layout.title and fig.layout.title.text:
                title = fig.layout.title.text
            if fig.layout.height:
                height = int(fig.layout.height)
            if fig.layout.xaxis and fig.layout.xaxis.title and fig.layout.xaxis.title.text:
                xlabel = fig.layout.xaxis.title.text
            if fig.layout.yaxis and fig.layout.yaxis.title and fig.layout.yaxis.title.text:
                ylabel = fig.layout.yaxis.title.text
            if fig.layout.font and fig.layout.font.size:
                font_size = int(fig.layout.font.size)
            if fig.layout.xaxis and fig.layout.xaxis.showgrid is False:
                show_grid = False
            if fig.data and hasattr(fig.data[0], "marker"):
                marker = fig.data[0].marker
                size_val = marker.get("size") if isinstance(marker, dict) else getattr(marker, "size", None)
                if isinstance(size_val, (int, float)):
                    marker_size = int(size_val)
        except Exception:
            pass

        return {
            "title": title, "height": height, "xlabel": xlabel, "ylabel": ylabel,
            "font_size": font_size, "legend_position": legend_position,
            "show_grid": show_grid, "marker_size": marker_size,
        }

    def _render_export_row(self, fig):
        c1, c2, c3, c4 = st.columns(4)
        slug = _slug(self.key)

        png_bytes, png_err = _try_plotly_to_image(fig, "png", scale=EXPORT_DEFAULTS["png_scale"])
        svg_bytes, svg_err = _try_plotly_to_image(fig, "svg")

        c1.download_button(
            "⬇ PNG",
            data=png_bytes if png_bytes else b"",
            file_name=f"{slug}.png",
            mime="image/png",
            key=f"{self.key}_dl_png",
            use_container_width=True,
            disabled=png_bytes is None,
            help=png_err if png_err else "Download as PNG",
        )
        c2.download_button(
            "⬇ SVG",
            data=svg_bytes if svg_bytes else b"",
            file_name=f"{slug}.svg",
            mime="image/svg+xml",
            key=f"{self.key}_dl_svg",
            use_container_width=True,
            disabled=svg_bytes is None,
            help=svg_err if svg_err else "Download as SVG (vector)",
        )

        try:
            html_bytes = fig.to_html(include_plotlyjs=EXPORT_DEFAULTS["html_plotlyjs"]).encode("utf-8")
        except Exception as e:
            html_bytes = None
            logger.warning(f"to_html failed for {self.key}: {e}")
        c3.download_button(
            "⬇ HTML",
            data=html_bytes if html_bytes else b"",
            file_name=f"{slug}.html",
            mime="text/html",
            key=f"{self.key}_dl_html",
            use_container_width=True,
            disabled=html_bytes is None,
            help="Interactive HTML (Plotly)",
        )

        if c4.button("📄 Add to Report", key=f"{self.key}_add_report", use_container_width=True,
                     help="Queue this figure for the report (built in P1)"):
            _queue_for_report(
                module=getattr(self, "module", None),
                key=self.key,
                fig=fig,
                title=st.session_state.get(self.edits_key, {}).get("title")
                       or _read_title(fig) or self.key,
                params=st.session_state.get(self.params_key, {}),
                kind="plotly",
            )


# ---------------------------------------------------------------------------
# MplPlotManager — matplotlib BytesIO figures
# ---------------------------------------------------------------------------


class MplPlotManager:
    """Manages a matplotlib figure returned as a BytesIO PNG.

    Same UI shell as PlotManager (display + export + add-to-report) but
    edits (title / figsize / dpi) require re-rendering — they're forwarded
    as kwargs to the plot function. If the plot function doesn't accept
    those kwargs they're silently ignored, so visualizers can opt in
    incrementally without breaking.
    """

    def __init__(self, plot_key: str):
        self.key = plot_key
        self.params_key = f"{plot_key}_mpl_params"
        self.buf_key = f"{plot_key}_mpl_buf"
        self.edits_key = f"{plot_key}_mpl_edits"
        self.hash_key = f"{plot_key}_mpl_phash"

    def render_generate_button(self, plot_function: Callable, **kwargs):
        current = st.session_state.get(self.buf_key)
        if current is None:
            self._generate(plot_function, kwargs, silent=True)

        label = "Update Plot" if st.session_state.get(self.buf_key) else "Generate Plot"
        if st.button(label, key=f"{self.key}_mpl_btn", use_container_width=True, type="primary"):
            self._generate(plot_function, kwargs, silent=False)

    def _generate(self, plot_function: Callable, kwargs: dict, *, silent: bool):
        new_hash = _params_hash(kwargs)
        prev_hash = st.session_state.get(self.hash_key)

        # Layer persisted edits onto kwargs (best-effort, ignored if plot fn doesn't accept them)
        edits = st.session_state.get(self.edits_key, {})
        merged = {**kwargs, **edits} if (prev_hash == new_hash and edits) else kwargs

        try:
            if silent:
                buf = plot_function(**merged)
            else:
                with st.spinner("Generating plot..."):
                    buf = plot_function(**merged)
            st.session_state[self.buf_key] = buf
            st.session_state[self.params_key] = kwargs
            st.session_state[self.hash_key] = new_hash
            if prev_hash != new_hash:
                st.session_state.pop(self.edits_key, None)
        except TypeError as te:
            # plot_function probably rejected an edit kwarg — retry without edits
            if edits and silent:
                logger.info(f"{self.key}: plot fn rejected edits ({te}); regenerating without")
                try:
                    buf = plot_function(**kwargs)
                    st.session_state[self.buf_key] = buf
                    st.session_state[self.params_key] = kwargs
                    st.session_state[self.hash_key] = new_hash
                    return
                except Exception as e:
                    logger.error(f"Mpl plot fallback failed for {self.key}: {e}", exc_info=True)
            if not silent:
                st.error(f"Error generating plot: {te}")
            st.session_state[self.buf_key] = None
        except Exception as e:
            if not silent:
                st.error(f"Error generating plot: {e}")
                logger.error(f"Mpl plot generation failed for {self.key}: {e}", exc_info=True)
            st.session_state[self.buf_key] = None

    def render_plot_and_editor(self):
        buf = st.session_state.get(self.buf_key)
        if buf is None:
            st.info("👆 Configure parameters above and click **Generate Plot** to visualize results.")
            return

        with st.expander("🎨 Edit Plot Appearance", expanded=False):
            current = st.session_state.get(self.edits_key, {})
            c1, c2, c3 = st.columns(3)
            new_title = c1.text_input("Plot Title", value=current.get("title", ""), key=f"{self.key}_mpl_title")
            new_w = c2.slider("Width (in)", 4, 20, current.get("figsize", (10, 6))[0], key=f"{self.key}_mpl_w")
            new_h = c3.slider("Height (in)", 3, 16, current.get("figsize", (10, 6))[1], key=f"{self.key}_mpl_h")
            new_dpi = st.slider("DPI", 72, 300, current.get("dpi", EXPORT_DEFAULTS["matplotlib_dpi"]),
                                step=10, key=f"{self.key}_mpl_dpi")

            if st.button("Apply Visual Changes", key=f"{self.key}_mpl_apply",
                         help="Re-renders the plot with these settings (visualizer must accept title/figsize/dpi kwargs)"):
                st.session_state[self.edits_key] = {
                    "title": new_title,
                    "figsize": (new_w, new_h),
                    "dpi": new_dpi,
                }
                st.warning("Click **Update Plot** above to apply edits.")

        c1, c2 = st.columns(2)
        png_data = buf.getvalue() if hasattr(buf, "getvalue") else (buf if isinstance(buf, (bytes, bytearray)) else None)
        c1.download_button(
            "⬇ PNG",
            data=png_data if png_data else b"",
            file_name=f"{_slug(self.key)}.png",
            mime="image/png",
            key=f"{self.key}_mpl_dl_png",
            use_container_width=True,
            disabled=png_data is None,
        )
        if c2.button("📄 Add to Report", key=f"{self.key}_mpl_add_report", use_container_width=True,
                     help="Queue this figure for the report (built in P1)"):
            _queue_for_report(
                module=getattr(self, "module", None),
                key=self.key,
                fig=png_data,
                title=st.session_state.get(self.edits_key, {}).get("title") or self.key,
                params=st.session_state.get(self.params_key, {}),
                kind="matplotlib_png",
            )

        st.image(buf, use_container_width=True)


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------


def _try_plotly_to_image(fig, fmt: str, **kw) -> tuple[bytes | None, str | None]:
    """Wrap fig.to_image() so kaleido failures degrade gracefully.

    Returns (bytes, None) on success or (None, error_message) on failure.
    The disabled download_button uses the error_message as its tooltip.
    """
    try:
        return fig.to_image(format=fmt, **kw), None
    except Exception as e:
        msg = f"{fmt.upper()} export unavailable ({type(e).__name__})"
        logger.debug(f"plotly to_image({fmt}) failed: {e}")
        return None, msg


def _legend_anchor(position: str) -> dict:
    """Map a friendly legend position name to Plotly legend anchor kwargs."""
    if position == "hidden":
        return dict(visible=False) if False else {}  # plotly: use showlegend layout
    return {
        "right":  dict(orientation="v", x=1.02, y=1.0),
        "top":    dict(orientation="h", x=0.0, y=1.10),
        "bottom": dict(orientation="h", x=0.0, y=-0.20),
        "left":   dict(orientation="v", x=-0.20, y=1.0),
    }.get(position, {})


def _read_title(fig) -> str | None:
    try:
        if fig.layout.title and fig.layout.title.text:
            return fig.layout.title.text
    except Exception:
        return None
    return None


def _queue_for_report(module: str | None, key: str, fig: Any, title: str,
                      params: dict, kind: str) -> None:
    from utils.report_builder import ReportBuilder
    builder = st.session_state.setdefault("report", ReportBuilder())
    builder.add_figure(
        module=module or "unknown",
        key=key,
        fig=fig,
        title=title,
        params=params,
        kind=kind,
    )
    st.toast(f"Added to report: {title}", icon="📄")
