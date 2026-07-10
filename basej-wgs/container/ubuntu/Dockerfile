FROM ubuntu:24.04

RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

# Build: docker build -t basejumper_ubuntu_24.04_stable --push .
