FROM python:3.11-slim

# System deps for kaleido's bundled Chromium (PNG/SVG export via fig.to_image())
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgbm1 \
        libasound2 \
        libxshmfence1 \
        curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8501

HEALTHCHECK --interval=30s --timeout=10s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

CMD ["streamlit", "run", "app.py", \
     "--server.address=0.0.0.0", \
     "--server.port=8501", \
     "--server.headless=true"]
