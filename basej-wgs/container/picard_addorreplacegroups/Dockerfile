# picard_addorreplacegroups container (Picard + samtools, linux/amd64)
# Tag: picard_addorreplacegroups_3.1.2
# Build:
#   docker build -t basejumper_picard_addorreplacegroups_3.1.2 .

# Use Ubuntu as the base image
FROM ubuntu:24.04

# Set non-interactive frontend (to avoid interactive prompts during the build)
ARG DEBIAN_FRONTEND=noninteractive

# Update the package repository and install necessary tools
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    wget \
    bzip2 \
    ca-certificates

# Download and install Miniconda (required for Bioconda)
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh \
    && /bin/bash ~/miniconda.sh -b -p /opt/conda \
    && rm ~/miniconda.sh

# Add Conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Initialize Conda in bash config
RUN conda init bash

# Add the Bioconda channel as well as the other channels it depends on
RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge

# Install the required tools into a dedicated conda env (not base).
# Using --override-channels with conda-forge + bioconda avoids the Anaconda
# "defaults" ToS gate, and a separate env avoids the base env's libzlib/python
# pins conflicting with samtools=1.17.
RUN conda create -y -n tools --override-channels -c conda-forge -c bioconda \
    picard=3.1.1 \
    samtools=1.17 \
    && conda clean -afy
ENV PATH=/opt/conda/envs/tools/bin:$PATH

# Set the default command to execute when creating a new container
CMD [ "/bin/bash" ]
