# custom_parabricks-metrics container
# Provides: Python 3.8 with pandas, multiqc, and pyarrow for metrics parsing
# Tag: custom_parabricks-metrics_1.0.3
#
# Build:
#   docker build -t basejumper_custom_parabricks-metrics_1.0.3 .

FROM python:3.8-slim-bullseye

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        libncurses5-dev \
        libbz2-dev \
        liblzma-dev \
        zlib1g-dev \
        gcc \
        wget \
        make \
        bzip2 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt /requirements.txt
RUN pip3 install --no-cache-dir -r /requirements.txt

# Copy scripts
RUN mkdir -p /scripts
COPY scripts/ /scripts/

WORKDIR /
ENV PATH="$PATH:/scripts"

CMD ["/bin/bash"]
