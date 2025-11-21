import streamlit as st
import plotly.express as px
import matplotlib.colors as mcolors

class ModuleColorManager:
    """
    Manages global color settings and plot themes for a specific module.
    """
    def __init__(self, module_prefix: str, groups: list):
        self.module_prefix = module_prefix
        self.groups = sorted(list(set(groups))) if groups else []
        self.color_key = f"{module_prefix}_global_color_map"
        self.theme_key = f"{module_prefix}_global_theme"

    def _to_hex(self, color):
        """Safely converts any color format to hex."""
        try:
            return mcolors.to_hex(color)
        except ValueError:
            return "#000000"

    def render_global_settings(self):
        """Renders the Global Settings expander."""
        # Initialize Defaults
        if self.color_key not in st.session_state:
            defaults = px.colors.qualitative.Plotly
            st.session_state[self.color_key] = {
                g: self._to_hex(defaults[i % len(defaults)]) for i, g in enumerate(self.groups)
            }
        
        if self.theme_key not in st.session_state:
            st.session_state[self.theme_key] = "plotly_white"

        with st.expander(f"🎨 Global Settings ({self.module_prefix.title()})", expanded=False):
            with st.form(key=f"{self.module_prefix}_global_form"):
                c1, c2 = st.columns(2)
                
                with c1:
                    theme_options = {
                        "Standard White": "plotly_white", "Dark Mode": "plotly_dark",
                        "Minimal": "simple_white", "GGPlot Style": "ggplot2", "Seaborn Style": "seaborn"
                    }
                    current_val = st.session_state[self.theme_key]
                    idx = list(theme_options.values()).index(current_val) if current_val in theme_options.values() else 0
                    selected_theme_label = st.selectbox("Plot Theme", list(theme_options.keys()), index=idx)

                with c2:
                    palette_mode = st.radio("Palette Type", ["Custom", "Colorblind Safe", "Auto"], horizontal=True)

                new_colors = {}
                if self.groups:
                    st.markdown("**Group Colors**")
                    cols = st.columns(4)
                    
                    if palette_mode == "Colorblind Safe":
                        base_palette = px.colors.qualitative.Safe
                    elif palette_mode == "Auto":
                        base_palette = px.colors.qualitative.Plotly
                    else:
                        base_palette = None

                    current_map = st.session_state[self.color_key]

                    for i, group in enumerate(self.groups):
                        with cols[i % 4]:
                            if base_palette:
                                default_hex = self._to_hex(base_palette[i % len(base_palette)])
                            else:
                                default_hex = self._to_hex(current_map.get(group, "#000000"))

                            picked_color = st.color_picker(
                                label=str(group), value=default_hex,
                                key=f"{self.module_prefix}_picker_{group}",
                                disabled=(palette_mode == "Auto")
                            )
                            new_colors[group] = picked_color

                if st.form_submit_button("Apply Settings"):
                    st.session_state[self.theme_key] = theme_options[selected_theme_label]
                    st.session_state[self.color_key] = new_colors
                    st.success("Settings updated! Plots will refresh automatically.")
                    # REMOVED st.rerun() to prevent tab reset

        return {
            "color_map": st.session_state[self.color_key],
            "template": st.session_state[self.theme_key]
        }