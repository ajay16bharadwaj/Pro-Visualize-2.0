# config/quant_config.py

QUANT_PLOT_REGISTRY = {
    "protein_counts": {
        "title": "Quantified Proteins per Sample",
        "has_markers": False,
        "has_color_groups": True,  # Can use global group colors
        "default_height": 500,
        "help": "Bar chart showing the number of identified proteins in each sample."
    },
    "missing_heatmap": {
        "title": "Missing Values Heatmap",
        "has_markers": False,
        "has_color_groups": False, # Uses binary/status colors, not group colors
        "default_height": 600,
        "help": "Visualizes the pattern of missing data across samples."
    },
    "missing_dist": {
        "title": "Missing Value Distribution",
        "has_markers": False,
        "has_color_groups": False,
        "default_height": 500,
        "resizable": True
    },
    "protein_overlap": {
        "title": "Protein Overlap",
        "has_markers": False,
        "has_color_groups": False,
        "default_height": 500
    },
    "rank_order": {
        "title": "Protein Rank Order",
        "has_markers": True,     # Enables marker size slider
        "has_text_labels": True, # Enables label toggle
        "has_color_groups": False, # Uses special highlighting categories
        "default_height": 600
    },
    "pca_anno": {
        "title": "PCA by Annotation",
        "has_markers": True,
        "has_text_labels": True,
        "has_color_groups": True, # Critical for coloring PCA by group
        "default_height": 600
    },
    "pca_cluster": {
        "title": "PCA with Clustering",
        "has_markers": True,
        "has_text_labels": True,
        "has_color_groups": False, # Colors determined by cluster ID, not exp group
        "default_height": 600
    }
}