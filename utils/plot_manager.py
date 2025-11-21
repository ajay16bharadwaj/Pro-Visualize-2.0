import streamlit as st
import plotly.graph_objects as go

class EnhancedPlotManager:
    def __init__(self, plot_key: str, metadata: dict):
        self.key = plot_key
        self.metadata = metadata
        self.fig_key = f"{plot_key}_fig"
        self.settings_key = f"{plot_key}_visual_settings"
        self.globals_hash_key = f"{plot_key}_globals_hash"

    def render(self, plot_function, **kwargs):
        """Renders generation button, plot, and visual editor."""
        
        # 1. Check for Global Setting Changes (Auto-Regenerate)
        current_globals = {k: kwargs[k] for k in ['color_map', 'template'] if k in kwargs}
        global_changed = False
        if current_globals:
            last_globals = st.session_state.get(self.globals_hash_key)
            if last_globals != current_globals:
                global_changed = True
                st.session_state[self.globals_hash_key] = current_globals

        # 2. Generate/Update Button
        label = "Update Plot" if st.session_state.get(self.fig_key) else "Generate Plot"
        trigger = st.button(label, key=f"{self.key}_btn", use_container_width=True)

        # 3. Generation Logic (Click OR Missing OR Globals Changed)
        if trigger or st.session_state.get(self.fig_key) is None or global_changed:
            with st.spinner("Generating plot..."):
                try:
                    fig = plot_function(**kwargs)
                    st.session_state[self.fig_key] = fig
                except Exception as e:
                    st.error(f"Plot generation failed: {e}")
                    st.session_state[self.fig_key] = None

        # 4. Display Plot & Editor
        if st.session_state.get(self.fig_key):
            fig = st.session_state[self.fig_key]
            
            # Create a placeholder so we can overwrite the plot later without a rerun
            plot_placeholder = st.empty()
            
            # Get currently saved settings (or empty dict)
            current_settings = st.session_state.get(self.settings_key, {})
            
            # Helper to render the figure into the placeholder
            def render_figure(figure_obj, settings):
                if hasattr(figure_obj, 'update_layout'): # Plotly
                    self._apply_layout_updates(figure_obj, settings)
                    plot_placeholder.plotly_chart(figure_obj, use_container_width=True)
                elif hasattr(figure_obj, 'savefig'): # Matplotlib
                    # Matplotlib title updates are tricky on cached figs, 
                    # but we display it as-is for now.
                    plot_placeholder.pyplot(figure_obj, use_container_width=True)
                else: # Static Image
                    plot_placeholder.image(figure_obj, caption=settings.get('title', ''))

            # Initial Render
            render_figure(fig, current_settings)

            # Visual Editor (Bottom)
            with st.expander(f"🎨 Edit Plot Appearance", expanded=False):
                with st.form(key=f"{self.key}_visual_form"):
                    c1, c2, c3 = st.columns(3)
                    
                    default_title = current_settings.get('title', self.metadata.get('title', ''))
                    default_height = current_settings.get('height', self.metadata.get('default_height', 500))
                    
                    new_title = c1.text_input("Plot Title", value=default_title, key=f"{self.key}_title")
                    new_height = c2.slider("Plot Height", 300, 1200, default_height, key=f"{self.key}_height")
                    
                    new_marker_size = None
                    if self.metadata.get('has_markers', False):
                        default_msize = current_settings.get('marker_size', 6)
                        new_marker_size = c3.slider("Marker Size", 2, 20, default_msize, key=f"{self.key}_msize")
                    else:
                        c3.write("")

                    c4, c5 = st.columns(2)
                    new_xlabel = c4.text_input("X-Axis Label", value=current_settings.get('xaxis_title', ''), key=f"{self.key}_xlabel")
                    new_ylabel = c5.text_input("Y-Axis Label", value=current_settings.get('yaxis_title', ''), key=f"{self.key}_ylabel")

                    # --- THE FIX: Update In-Place ---
                    if st.form_submit_button("Apply Visual Changes"):
                        new_settings = {
                            'title': new_title, 'height': new_height,
                            'xaxis_title': new_xlabel, 'yaxis_title': new_ylabel
                        }
                        if new_marker_size: new_settings['marker_size'] = new_marker_size
                        
                        # 1. Update State
                        st.session_state[self.settings_key] = new_settings
                        
                        # 2. Re-render immediately to placeholder (Visual update)
                        render_figure(fig, new_settings)
                        
                        # 3. No st.rerun() needed!

    def _apply_layout_updates(self, fig, settings):
        if not settings: return
        try:
            updates = {}
            if settings.get('height'): updates['height'] = settings['height']
            if settings.get('title'): updates['title_text'] = settings['title']
            if settings.get('xaxis_title'): updates['xaxis_title'] = settings['xaxis_title']
            if settings.get('yaxis_title'): updates['yaxis_title'] = settings['yaxis_title']
            
            fig.update_layout(**updates)
            
            if 'marker_size' in settings:
                fig.update_traces(marker=dict(size=settings['marker_size']))
        except Exception: pass