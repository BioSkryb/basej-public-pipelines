# ginko_parser container
# Provides: Python 3.10 with pandas and numpy for parsing Ginkgo CNV outputs
# Tag: ginko_parser_0.2.1
#
# Build:
#   docker build -t basejumper_ginko_parser_0.2.1 .

FROM continuumio/miniconda3:latest

ENV DEBIAN_FRONTEND=noninteractive

# Install system packages
RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        make \
        bzip2 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /usr/src/app

# Create conda environment with pinned dependencies
COPY environment.yml .
RUN conda env create -f environment.yml && \
    conda clean -afy

# Activate the environment by default
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate app_env" >> ~/.bashrc
ENV PATH=/opt/conda/envs/app_env/bin:$PATH

# Copy the parser script
COPY scripts/ /scripts/
ENV PATH="/scripts:${PATH}"

# Run as a non-root user (Trivy DS-0002)
RUN useradd --create-home --shell /bin/bash appuser
USER appuser

CMD ["/bin/bash"]
