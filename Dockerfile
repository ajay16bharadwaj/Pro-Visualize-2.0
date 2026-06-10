# Pinned to a specific patch tag for reproducible builds (not floating 3.11-slim).
FROM python:3.11.15-slim

# System deps for kaleido's bundled Chromium (PNG/SVG export via fig.to_image()).
# The minimal set (libgbm1/libasound2/libxshmfence1) is not enough on slim images:
# Chromium also needs xkb, atk-bridge, gtk, pango, cairo and base fonts to render.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgbm1 \
        libasound2 \
        libxshmfence1 \
        libxkbcommon0 \
        libatk-bridge2.0-0 \
        libgtk-3-0 \
        libpango-1.0-0 \
        libcairo2 \
        fonts-liberation \
        curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install deps first so the layer is cached across app-code changes.
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Run as a non-root user (container-escape hardening). Created after COPY so the
# app tree is owned by the runtime user.
RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /app
USER appuser

EXPOSE 8501

# /_stcore/health is Streamlit's internal health endpoint. It is stable for the
# pinned streamlit==1.36.0; revisit this path if Streamlit is upgraded.
HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

CMD ["streamlit", "run", "app.py", \
     "--server.address=0.0.0.0", \
     "--server.port=8501", \
     "--server.headless=true"]
