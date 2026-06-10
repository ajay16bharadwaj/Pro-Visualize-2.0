"""
ReportBuilder — per-session report queue and export engine.

Lives in st.session_state['report']. Figures are queued via PlotManager's
"Add to Report" button; the Report tab calls render_preview() which exposes
two download paths: interactive HTML and ZIP bundle.
"""

from __future__ import annotations

import base64
import io
import json
import logging
import platform
import subprocess
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import streamlit as st

logger = logging.getLogger(__name__)

_TEMPLATES_DIR = Path(__file__).parent.parent / "templates"
_REPO_ROOT = Path(__file__).parent.parent

# Packages whose versions matter for reproducing an analysis. Captured into
# the report's Methods & Reproducibility section.
_PROVENANCE_PACKAGES = [
    "streamlit", "pandas", "numpy", "scipy", "scikit-learn", "plotly",
    "scanpy", "anndata", "gseapy", "gprofiler-official", "umap-learn",
]


def _git_commit() -> str | None:
    """Best-effort short git commit of the running code; None if unavailable."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=str(_REPO_ROOT), capture_output=True, text=True, timeout=3,
        )
        return out.stdout.strip() or None if out.returncode == 0 else None
    except Exception:
        return None


def _package_versions() -> dict[str, str]:
    from importlib.metadata import PackageNotFoundError, version
    versions: dict[str, str] = {}
    for pkg in _PROVENANCE_PACKAGES:
        try:
            versions[pkg] = version(pkg)
        except PackageNotFoundError:
            continue
    return versions


class ReportBuilder:
    """Per-session report state, lives in st.session_state['report']."""

    APP_VERSION = "Pro-Visualize 2.0"

    def __init__(self):
        self.items: list[dict] = []

    # ------------------------------------------------------------------
    # Provenance
    # ------------------------------------------------------------------

    def provenance(self) -> dict:
        """Analysis-provenance block for the report's Methods section.

        Captures app version, git commit, timestamp, Python and key package
        versions so a reader can reproduce the environment that produced the
        figures.
        """
        return {
            "app_version": self.APP_VERSION,
            "git_commit": _git_commit(),
            "generated_at": datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "packages": _package_versions(),
        }

    # ------------------------------------------------------------------
    # Mutation API
    # ------------------------------------------------------------------

    def add_figure(
        self,
        module: str,
        key: str,
        fig: Any,
        title: str,
        params: dict,
        kind: str = "plotly",
        notes: str = "",
    ) -> None:
        """Upsert by key — updating replaces in-place, preserving order."""
        item = {
            "module": module,
            "key": key,
            "kind": kind,
            "fig": fig,
            "title": title,
            "params": self._safe_params(key, params),
            "notes": notes,
            "added_at": datetime.now(tz=timezone.utc).isoformat(),
        }
        for i, existing in enumerate(self.items):
            if existing["key"] == key:
                self.items[i] = item
                return
        self.items.append(item)

    @staticmethod
    def _safe_params(key: str, params: dict | None) -> dict:
        """Keep only JSON-serializable params so exports never crash later.

        Each value is coerced with json.dumps(default=str); a value that still
        fails is dropped with a warning. This fails soft at queue time rather
        than blowing up HTML/ZIP export after a long analysis session.
        """
        if not params:
            return {}
        if not isinstance(params, dict):
            logger.warning("Report params for %s is not a dict (%s); ignoring.", key, type(params))
            return {}
        clean: dict = {}
        for k, v in params.items():
            try:
                json.dumps(v, default=str)
                clean[str(k)] = v
            except Exception as e:
                logger.warning("Dropping non-serializable report param %r for %s: %s", k, key, e)
        return clean

    def add_table(self, module: str, key: str, df: Any, caption: str) -> None:
        """Queue a DataFrame for the report."""
        self.add_figure(
            module=module, key=key, fig=df, title=caption,
            params={}, kind="table",
        )

    def add_section(self, module: str, heading: str, body_md: str) -> None:
        """Queue a markdown text block."""
        key = f"{module}_section_{heading[:20]}"
        item = {
            "module": module,
            "key": key,
            "kind": "section",
            "fig": None,
            "title": heading,
            "params": {},
            "notes": body_md,
            "added_at": datetime.now(tz=timezone.utc).isoformat(),
        }
        for i, existing in enumerate(self.items):
            if existing["key"] == key:
                self.items[i] = item
                return
        self.items.append(item)

    def remove(self, module: str, key: str) -> None:
        self.items = [it for it in self.items if it["key"] != key]

    def reorder(self, module: str, keys: list[str]) -> None:
        """Reorder items belonging to `module` per `keys`.

        Items from other modules keep their relative positions, occupying
        the same slots they held before.
        """
        module_map = {it["key"]: it for it in self.items if it["module"] == module}
        reordered = [module_map[k] for k in keys if k in module_map]
        reorder_iter = iter(reordered)
        result = []
        for it in self.items:
            if it["module"] == module:
                try:
                    result.append(next(reorder_iter))
                except StopIteration:
                    pass
            else:
                result.append(it)
        self.items = result

    # ------------------------------------------------------------------
    # Streamlit preview UI
    # ------------------------------------------------------------------

    def render_preview(self) -> None:
        st.header("📋 Report Builder")

        if not self.items:
            st.info(
                "No figures queued yet. Generate plots in any analysis tab and "
                "click **📄 Add to Report** to queue them here.",
                icon="📊",
            )
            return

        st.caption(f"{len(self.items)} item(s) queued · edits to notes are captured at download time")
        self._render_download_row("top")
        st.divider()

        _KIND_ICON = {
            "plotly": "📈", "matplotlib_png": "🖼", "table": "📋", "section": "📝",
        }
        for idx, item in enumerate(self.items):
            col_meta, col_remove = st.columns([6, 1])
            icon = _KIND_ICON.get(item["kind"], "•")
            col_meta.markdown(
                f"**{icon} {item['title']}** &nbsp;&nbsp; "
                f"<span style='color:#888;font-size:0.85em'>{item['module']}</span>",
                unsafe_allow_html=True,
            )
            if col_remove.button("🗑", key=f"rpt_remove_{item['key']}", help="Remove from report"):
                self.remove(item["module"], item["key"])
                st.rerun()

            st.text_area(
                "Notes",
                value=item.get("notes", ""),
                key=f"rpt_notes_{item['key']}",
                height=60,
                placeholder="Add context, interpretation, or methods notes…",
                label_visibility="collapsed",
            )

            if idx < len(self.items) - 1:
                st.divider()

        st.divider()
        self._render_download_row("bottom")

    def _render_download_row(self, position: str = "bottom") -> None:
        c1, c2 = st.columns(2)
        try:
            html_bytes = self.export_html()
            c1.download_button(
                "⬇ Download Interactive HTML",
                data=html_bytes,
                file_name="pro_visualize_report.html",
                mime="text/html",
                use_container_width=True,
                type="primary",
                key=f"rpt_dl_html_{position}",
            )
        except Exception as e:
            logger.error(f"HTML export failed: {e}", exc_info=True)
            c1.error(f"HTML export error: {e}")

        try:
            zip_bytes = self.export_zip()
            c2.download_button(
                "⬇ Download ZIP Bundle",
                data=zip_bytes,
                file_name="pro_visualize_report.zip",
                mime="application/zip",
                use_container_width=True,
                key=f"rpt_dl_zip_{position}",
            )
        except Exception as e:
            logger.error(f"ZIP export failed: {e}", exc_info=True)
            c2.error(f"ZIP export error: {e}")

    # ------------------------------------------------------------------
    # Export: Interactive HTML
    # ------------------------------------------------------------------

    def export_html(self) -> bytes:
        from jinja2 import Environment, FileSystemLoader

        env = Environment(
            loader=FileSystemLoader(str(_TEMPLATES_DIR)),
            autoescape=True,
        )
        template = env.get_template("report.html.j2")

        rendered_items = []
        for item in self.items:
            notes = st.session_state.get(f"rpt_notes_{item['key']}", item.get("notes", ""))
            entry = {
                **item,
                "notes": notes,
                "fig_html": None,
                "fig_b64": None,
                "table_html": None,
            }

            if item["kind"] == "plotly" and item["fig"] is not None:
                try:
                    # include_plotlyjs=False: template loads Plotly CDN once in <head>
                    entry["fig_html"] = item["fig"].to_html(
                        full_html=False, include_plotlyjs=False
                    )
                except Exception as e:
                    logger.warning(f"to_html failed for {item['key']}: {e}")

            elif item["kind"] == "matplotlib_png" and item["fig"] is not None:
                try:
                    raw = item["fig"]
                    png_data = raw if isinstance(raw, (bytes, bytearray)) else raw.getvalue()
                    entry["fig_b64"] = base64.b64encode(png_data).decode("ascii")
                except Exception as e:
                    logger.warning(f"base64 encode failed for {item['key']}: {e}")

            elif item["kind"] == "table" and item["fig"] is not None:
                try:
                    entry["table_html"] = item["fig"].to_html(
                        index=False, classes="report-table", border=0
                    )
                except Exception as e:
                    logger.warning(f"DataFrame to_html failed for {item['key']}: {e}")

            rendered_items.append(entry)

        prov = self.provenance()
        html_str = template.render(
            items=rendered_items,
            generated_at=prov["generated_at"],
            app_version=self.APP_VERSION,
            provenance=prov,
        )
        return html_str.encode("utf-8")

    # ------------------------------------------------------------------
    # Export: ZIP Bundle
    # ------------------------------------------------------------------

    def export_zip(self) -> bytes:
        buf = io.BytesIO()
        notes_lines: list[str] = []
        manifest: list[dict] = []
        all_params: dict = {}

        with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            for item in self.items:
                notes = st.session_state.get(f"rpt_notes_{item['key']}", item.get("notes", ""))
                slug = f"{item['module']}_{item['key']}"
                manifest.append({
                    "module": item["module"],
                    "key": item["key"],
                    "title": item["title"],
                    "kind": item["kind"],
                    "added_at": item["added_at"],
                })
                all_params[f"{item['module']}/{item['key']}"] = item.get("params", {})
                if notes:
                    notes_lines.append(f"## {item['title']}\n\n{notes}\n")

                if item["kind"] == "plotly" and item["fig"] is not None:
                    fig = item["fig"]
                    try:
                        zf.writestr(f"figures/{slug}.png", fig.to_image(format="png", scale=2))
                    except Exception as e:
                        logger.warning(f"PNG export failed for {slug}: {e}")
                    try:
                        zf.writestr(f"figures/{slug}.svg", fig.to_image(format="svg"))
                    except Exception as e:
                        logger.warning(f"SVG export failed for {slug}: {e}")
                    try:
                        zf.writestr(
                            f"figures/{slug}.html",
                            fig.to_html(include_plotlyjs="cdn").encode("utf-8"),
                        )
                    except Exception as e:
                        logger.warning(f"HTML export failed for {slug}: {e}")

                elif item["kind"] == "matplotlib_png" and item["fig"] is not None:
                    try:
                        raw = item["fig"]
                        png_data = raw if isinstance(raw, (bytes, bytearray)) else raw.getvalue()
                        zf.writestr(f"figures/{slug}.png", png_data)
                    except Exception as e:
                        logger.warning(f"Mpl PNG export failed for {slug}: {e}")

            prov = self.provenance()
            zf.writestr(
                "parameters.json",
                json.dumps({"_provenance": prov, "figures": all_params},
                           indent=2, default=str).encode("utf-8"),
            )
            zf.writestr(
                "provenance.json",
                json.dumps(prov, indent=2, default=str).encode("utf-8"),
            )
            notes_md = "\n".join(notes_lines) if notes_lines else ""
            zf.writestr("notes.md", notes_md.encode("utf-8"))
            zf.writestr(
                "manifest.json",
                json.dumps(manifest, indent=2).encode("utf-8"),
            )

        buf.seek(0)
        return buf.getvalue()
