# Phlox Docker Image
FROM python:3.10-slim

LABEL maintainer="Chen Xin <chenxin199261@gmail.com>"
LABEL description="Phlox - Reaction mechanism discovery from reactive MD trajectories"

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    git \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Copy package files
COPY pyproject.toml setup.py README.md ./
COPY src/ ./src/

# Install Python dependencies and package
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -e .

# Create directories for data
RUN mkdir -p /data /output

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Volume for input/output data
VOLUME ["/data", "/output"]

# Default working directory for data
WORKDIR /data

# Entry point
ENTRYPOINT ["phlox"]
CMD ["--help"]
