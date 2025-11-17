import streamlit as st
import plotly.graph_objects as go

class PlotManager:
    """
    Manages the state, generation, and editing of Plotly figures.
    Decouples widget interaction from expensive plot generation.
    """
    def __init__(self, plot_key: str):
        self.key = plot_key
        self.params_key = f"{plot_key}_params"
        self.fig_key = f"{plot_key}_fig"

    def render_generate_button(self, plot_function, **kwargs):
        """
        Renders a 'Generate Plot' button. 
        
        Feature: If no plot exists in state (first load), it automatically 
        generates one using the provided **kwargs (defaults).
        """
        # --- 1. Auto-Generate Default Plot on First Load ---
        if self.fig_key not in st.session_state:
            try:
                # Generate default plot immediately without user click
                fig = plot_function(**kwargs)
                st.session_state[self.fig_key] = fig
                st.session_state[self.params_key] = kwargs
            except Exception as e:
                # If default generation fails (e.g. data not ready), just log/warn 
                # and let the user click the button later.
                # st.warning(f"Ready to generate plot.") 
                st.session_state[self.fig_key] = None

        # --- 2. Render Button ---
        label = "Update Plot" if st.session_state.get(self.fig_key) else "Generate Plot"
        
        if st.button(label, key=f"{self.key}_btn", use_container_width=True, type="primary"):
            st.session_state[self.params_key] = kwargs
            try:
                with st.spinner("Generating plot..."):
                    fig = plot_function(**kwargs)
                    st.session_state[self.fig_key] = fig
            except Exception as e:
                st.error(f"Error generating plot: {e}")
                st.session_state[self.fig_key] = None

    def render_plot_and_editor(self):
        """
        Displays the generated plot (if it exists) and the 'Edit Plot' expander.
        Uses st.empty() to allow updates without full page reruns (fixing the tab reset issue).
        """
        fig = st.session_state.get(self.fig_key)

        if fig is not None:
            # --- 1. Create a Placeholder ---
            # We create the slot for the plot *before* the editor.
            # This allows us to fill it later with the updated figure.
            plot_placeholder = st.empty()

            # --- 2. Render Editor UI ---
            with st.expander("🎨 Edit Plot Appearance", expanded=False):
                c1, c2, c3 = st.columns(3)
                
                # Safe Layout Retrieval
                current_title = fig.layout.title.text if fig.layout.title.text else ""
                current_height = fig.layout.height if fig.layout.height else 500
                
                new_title = c1.text_input("Plot Title", value=current_title, key=f"{self.key}_title")
                new_height = c2.slider("Plot Height", 300, 1200, current_height, 50, key=f"{self.key}_height")
                
                # Robust Marker Size Retrieval (Fixes AttributeError)
                default_size = 8
                try:
                    if fig.data and hasattr(fig.data[0], 'marker'):
                        marker = fig.data[0].marker
                        # Handle dict or object access for size
                        size_val = marker.get('size') if isinstance(marker, dict) else getattr(marker, 'size', None)
                        if isinstance(size_val, (int, float)):
                            default_size = int(size_val)
                except Exception:
                    pass 
                
                marker_size = c3.slider("Marker Size", 2, 20, default_size, key=f"{self.key}_msize")
                
                c4, c5 = st.columns(2)
                current_xlabel = fig.layout.xaxis.title.text if fig.layout.xaxis.title.text else ""
                current_ylabel = fig.layout.yaxis.title.text if fig.layout.yaxis.title.text else ""
                
                new_xlabel = c4.text_input("X-Axis Label", value=current_xlabel, key=f"{self.key}_xlabel")
                new_ylabel = c5.text_input("Y-Axis Label", value=current_ylabel, key=f"{self.key}_ylabel")

                # --- Apply Button Logic ---
                if st.button("Apply Visual Changes", key=f"{self.key}_apply"):
                    # Update layout in-place
                    fig.update_layout(
                        title_text=new_title,
                        height=new_height,
                        xaxis_title=new_xlabel,
                        yaxis_title=new_ylabel
                    )
                    # Update traces safely
                    try:
                        fig.update_traces(marker=dict(size=marker_size))
                    except Exception:
                        pass
                    
                    # Note: We do NOT call st.rerun() here. 
                    # The logic below will render the *updated* fig into the placeholder.

            # --- 3. Render Figure into Placeholder ---
            # This runs on every script execution. 
            # If "Apply" was just clicked, 'fig' is already updated, so the plot updates instantly.
            plot_placeholder.plotly_chart(fig, use_container_width=True)
            
        else:
            st.info("👆 Configure parameters above and click **Generate Plot** to visualize results.")