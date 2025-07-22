import logging
from functools import wraps

# --- Decorator for Robust Plotting ---
def handle_plotting_errors(plot_function):
    """
    A decorator to gracefully handle exceptions in Streamlit plotting functions.
    If a plot fails to generate, it displays an error message in the UI.
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