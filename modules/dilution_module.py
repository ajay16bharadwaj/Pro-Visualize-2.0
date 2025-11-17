import streamlit as st
import pandas as pd
import logging
import plotly.express as px
from visualizations.dilution_series import DilutionSeriesVisualizer
from utils.helpers import handle_plotting_errors
from utils.plot_manager import PlotManager

logger = logging.getLogger(__name__)

class DilutionSeriesTab:
    # Define all plot keys used in this module for easy global updates
    PLOT_KEYS = [
        "dilution_counts", "dilution_dist", "dilution_cv", 
        "dilution_trends", "dilution_heatmap", "dilution_ratios", 
        "dilution_pca", "dilution_completeness"
    ]

    def __init__(self):
        if 'dilution_visualizer' not in st.session_state:
            st.session_state.dilution_visualizer = None

    def _upload_data_section(self):
        # (No changes to this method)
        with st.expander("Upload Data for Dilution Series", expanded=(st.session_state.dilution_visualizer is None)):
            protein_file = st.file_uploader("Upload Protein-Level Data (.txt, .csv, .tsv)", type=['txt', 'csv', 'tsv'], key="dilution_protein_uploader")
            metadata_file = st.file_uploader("Upload Metadata File (.csv, .tsv)", type=['csv', 'tsv'], key="dilution_metadata_uploader")

            if st.button("Load & Process Dilution Data", use_container_width=True):
                if protein_file and metadata_file:
                    try:
                        protein_df = pd.read_csv(protein_file, sep=None, engine='python')
                        metadata_df = pd.read_csv(metadata_file, sep=None, engine='python')
                        
                        with st.spinner("Initializing visualizer..."):
                            st.session_state.dilution_visualizer = DilutionSeriesVisualizer(protein_df, metadata_df)
                        st.success("Data loaded successfully!")
                        st.rerun()
                    except Exception as e:
                        st.error(f"An error occurred: {e}")
                        logger.error(f"Dilution Series Error: {e}", exc_info=True)
                else:
                    st.warning("Please upload both data and metadata files.")

    def _render_global_settings(self):
        """
        Renders controls for global plot styling inside an expander.
        Includes an 'Apply' button to update existing plots immediately.
        """
        # --- MOVED FROM SIDEBAR TO MAIN EXPANDER ---
        with st.expander("🎨 Global Plot Settings (Colors & Theme)", expanded=False):
            c1, c2 = st.columns(2)
            
            with c1:
                # 1. Theme / Background
                theme_options = {
                    "Standard White": "plotly_white",
                    "Dark Mode": "plotly_dark",
                    "Minimal": "simple_white",
                    "GGPlot Style": "ggplot2",
                    "Seaborn Style": "seaborn"
                }
                # Use a key to persist selection
                selected_theme_label = st.selectbox("Plot Theme", list(theme_options.keys()), key="dilution_global_theme_select")
                selected_template = theme_options[selected_theme_label]

            with c2:
                # 2. Color Palette
                palette_mode = st.radio(
                    "Palette Type", 
                    ["Auto", "Colorblind Safe", "Custom"], 
                    horizontal=True,
                    key="dilution_global_palette_mode"
                )

            visualizer = st.session_state.dilution_visualizer
            groups = visualizer.group_order 
            color_map = {}

            if palette_mode == "Colorblind Safe":
                safe_colors = px.colors.qualitative.Safe
                for i, group in enumerate(groups):
                    color_map[group] = safe_colors[i % len(safe_colors)]
                st.caption("✅ Using high-contrast, colorblind-safe colors.")
                
            elif palette_mode == "Custom":
                st.markdown("**Customize Group Colors:**")
                cols = st.columns(4)
                for i, group in enumerate(groups):
                    default_c = px.colors.qualitative.Safe[i % len(px.colors.qualitative.Safe)]
                    with cols[i % 4]: 
                         # Add unique key for each color picker
                         color_map[group] = st.color_picker(f"{group}", value=default_c, key=f"dilution_color_{group}")
            
            else: # Auto
                 color_map = None

            st.markdown("---")
            
            # --- THE APPLY BUTTON ---
            if st.button("Apply Settings to All Plots", type="primary", use_container_width=True):
                # Iterate through all known plot keys and update their layout template
                count = 0
                for key in self.PLOT_KEYS:
                    fig_key = f"{key}_fig"
                    # Check if a figure exists in session state
                    if fig_key in st.session_state and st.session_state[fig_key] is not None:
                        # Update the template (background/fonts)
                        st.session_state[fig_key].update_layout(template=selected_template)
                        count += 1
                
                if count > 0:
                    st.success(f"Updated theme for {count} active plots!")
                else:
                    st.info("Settings saved. Generate a plot to see changes.")
                
                # Note: We cannot easily update 'color_discrete_map' on existing figures without 
                # re-generating the data traces. We force a rerun to refresh the UI.
                st.rerun()

        # Return the kwargs dictionary for *future* plot generations
        global_kwargs = {"template": selected_template}
        if color_map:
            global_kwargs["color_discrete_map"] = color_map
            
        return global_kwargs

    @handle_plotting_errors
    def _display_plots(self):
        visualizer = st.session_state.dilution_visualizer
        st.subheader("Dilution Series Visualizations")

        # 1. Render Global Settings locally
        global_plot_kwargs = self._render_global_settings()

        overview_tab, analysis_tab, pca_tab, completeness_tab = st.tabs(["📊 Data Overview", "📈 Trend Analysis", "🧬 PCA", "✅ Completeness"])

        # --- TAB 1: OVERVIEW ---
        with overview_tab:
            st.markdown("### Protein Counts per Sample")
            counts_manager = PlotManager("dilution_counts")
            counts_manager.render_generate_button(
                visualizer.plot_protein_counts_per_sample,
                **global_plot_kwargs
            )
            counts_manager.render_plot_and_editor()

            st.markdown("---")
            st.markdown("### Intensity and CV Distributions")
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("**Intensity Distribution**")
                with st.expander("⚙️ Configuration", expanded=True):
                    plot_type = st.selectbox("Intensity Plot Type", ["box", "violin"], key="dist_plot_type")
                    
                    dist_manager = PlotManager("dilution_dist")
                    dist_manager.render_generate_button(
                        visualizer.plot_intensity_distribution,
                        plot_type=plot_type,
                        points=False,
                        **global_plot_kwargs
                    )
                dist_manager.render_plot_and_editor()

            with col2:
                st.markdown("**CV Distribution**")
                cv_manager = PlotManager("dilution_cv")
                cv_manager.render_generate_button(
                    visualizer.plot_cv_distribution,
                    **global_plot_kwargs
                )
                cv_manager.render_plot_and_editor()

        # --- TAB 2: ANALYSIS ---
        with analysis_tab:
            trends_tab, heatmap_tab, ratio_tab = st.tabs(["Protein Trends", "Heatmap", "Relative Abundance"])

            with trends_tab:
                st.markdown("### Individual Protein Trends")
                with st.expander("⚙️ Configuration", expanded=True):
                    all_proteins = visualizer.protein_df[visualizer.protein_id_col].unique().tolist()
                    
                    selected_proteins = st.multiselect(
                        "Select specific proteins to plot (optional)", 
                        options=all_proteins,
                        help="If none are selected, the most abundant proteins will be shown.",
                        key="dilution_trends_protein_selector"
                    )
                    
                    n_top = st.slider(
                        "Number of top proteins to show (if none selected above)", 1, 20, 5,
                        disabled=bool(selected_proteins),
                        key="dilution_trends_top_n_slider"
                    )

                    trends_manager = PlotManager("dilution_trends")
                    trends_manager.render_generate_button(
                        visualizer.plot_protein_trends,
                        proteins_to_plot=selected_proteins if selected_proteins else None,
                        n_top_proteins=n_top,
                        **global_plot_kwargs
                    )
                trends_manager.render_plot_and_editor()

            with heatmap_tab:
                st.markdown("### Protein Intensity Heatmap")
                with st.expander("⚙️ Configuration", expanded=True):
                    apply_zscore = st.toggle("Apply Z-score normalization", value=True, key="dilution_heat_z")
                    max_proteins = st.slider("Max proteins to display", 50, 2000, 500, 50, key="dilution_heat_max")
                    
                    heat_manager = PlotManager("dilution_heatmap")
                    heat_manager.render_generate_button(
                        visualizer.plot_heatmap_trends,
                        max_proteins_to_plot=max_proteins,
                        apply_zscore=apply_zscore,
                        **global_plot_kwargs 
                    )
                heat_manager.render_plot_and_editor()

            with ratio_tab:
                st.markdown("### Protein Abundance Ratios")
                st.markdown("Shows `Log2` intensity ratio relative to lowest concentration.")
                
                with st.expander("⚙️ Configuration", expanded=True):
                    show_lines = st.toggle("Show expected ratio lines", value=True, key="dilution_ratio_lines")
                    
                    ratio_manager = PlotManager("dilution_ratios")
                    ratio_manager.render_generate_button(
                        visualizer.plot_relative_abundance_ratios,
                        add_expected_lines=show_lines,
                        points=False,
                        **global_plot_kwargs
                    )
                ratio_manager.render_plot_and_editor()
            
        # --- TAB 3: PCA ---
        with pca_tab:
            st.markdown("### Principal Component Analysis (PCA)")
            with st.expander("⚙️ Configuration", expanded=True):
                c1, c2 = st.columns(2)
                color_by = c1.selectbox("Color PCA points by", options=['Group', 'Concentration', 'Replicate'], index=0, key="dilution_pca_color")
                symbol_by = c2.selectbox("Use symbols for", options=[None, 'Replicate', 'Group'], index=1, key="dilution_pca_symbol")
                
                pca_manager = PlotManager("dilution_pca")
                pca_manager.render_generate_button(
                    visualizer.plot_pca,
                    color_by=color_by,
                    symbol_by=symbol_by,
                    **global_plot_kwargs
                )
            pca_manager.render_plot_and_editor()


        # --- TAB 4: COMPLETENESS ---
        with completeness_tab:
            st.markdown("### Detection Completeness Analysis")
            
            with st.expander("⚙️ Configuration", expanded=True):
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    analysis_level = st.selectbox(
                        "Select Analysis Level",
                        options=["Proteins", "Precursors", "Stripped Sequences"],
                        key="dilution_comp_level"
                    )
                    # Map names to columns
                    level_mapping = {
                        "Proteins": visualizer.protein_id_col,
                        "Precursors": "Precursor.Id" if "Precursor.Id" in visualizer.protein_df.columns else visualizer.protein_id_col,
                        "Stripped Sequences": "Stripped.Sequence" if "Stripped.Sequence" in visualizer.protein_df.columns else visualizer.protein_id_col
                    }
                    selected_column = level_mapping[analysis_level]
                
                with col2:
                    cv_thresh = st.slider("CV Threshold (%)", 0.0, 50.0, 20.0, 1.0, key="dilution_comp_cv")
                    
                with col3:
                    use_log = st.toggle("Use logarithmic scale", value=False, key="dilution_comp_log")

                comp_manager = PlotManager("dilution_completeness")
                comp_manager.render_generate_button(
                    visualizer.plot_completeness_overview,
                    identifier_col=selected_column,
                    use_log_scale=use_log,
                    cv_threshold=cv_thresh,
                    **global_plot_kwargs
                )
            
            comp_manager.render_plot_and_editor()

    def render(self):
        st.header("Dilution Series Analysis")
        self._upload_data_section()
        if st.session_state.dilution_visualizer:
            self._display_plots()

def render():
    """Renders the Dilution Series tab."""
    dilution_tab = DilutionSeriesTab()
    dilution_tab.render()