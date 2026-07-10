# variantannotation container
# Provides: conda env (app_env) with python 3.10, bcftools 1.21, tabix 1.11,
#           ensembl-vep 111, vaex + pandas, plus custom_vcftotable_expand_fields.py.
# Tag: variantannotation_0.4
# Processes: PREPROCESS_VCF, BULK_GET_VARIANTS_TO_FILTER, FILTER_VEP_GERMLINE
#
# Build (x86/amd64 only):
#   docker build -t basejumper_variantannotation_0.4 .

FROM continuumio/miniconda3:latest

# Install system packages
RUN apt-get update && apt-get install -y \
            wget \
            make \
            bzip2 \
            && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set the working directory in the container
WORKDIR /usr/src/app

# Copy the environment.yml file to the container
COPY environment.yml .

# Create the Conda environment from the environment.yml file
RUN conda env create -f environment.yml

# Initialize Conda in bash shell, modify .bashrc to activate the environment automatically
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate app_env" >> ~/.bashrc
ENV PATH /opt/conda/envs/app_env/bin:$PATH

# Copy your scripts to the container
COPY ./scripts /scripts

# Update PATH to include the scripts directory
ENV PATH="/scripts:${PATH}"

# Set Bash as the default command
CMD ["/bin/bash"]
