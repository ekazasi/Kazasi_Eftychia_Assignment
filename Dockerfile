# Use the official Ubuntu image as the base
FROM ubuntu:20.04

# Set environment variables to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install required tools: R, RScript, yq, wget, tar, and unzip
RUN apt-get update && apt-get install -y \
    r-base \
    wget \
    tar \
    unzip \
    curl \
    && apt-get clean

# Install yq
RUN wget https://github.com/mikefarah/yq/releases/latest/download/yq_linux_amd64 -O /usr/bin/yq && \
    chmod +x /usr/bin/yq

# Set the working directory
WORKDIR /bio_final_kazasi_eftychia

# Download the knownCanonical.exonAA.fa file
RUN wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz20way/alignments/knownCanonical.exonAA.fa.gz && \
    gunzip knownCanonical.exonAA.fa.gz && \
    rm -rf knownCanonical.exonAA.fa.gz

# Download and extract IQ-TREE 2
RUN wget -q https://github.com/iqtree/iqtree2/releases/download/v2.4.0/iqtree-2.4.0-Linux-intel.tar.gz && \
    tar -xzf iqtree-2.4.0-Linux-intel.tar.gz && \
    mv iqtree-2.4.0-Linux-intel/bin/iqtree2 . && \
    chmod +x ./iqtree2 && \
    rm -rf iqtree-2.4.0-Linux-intel iqtree-2.4.0-Linux-intel.tar.gz

# Copy the bio_final directory contents into the container
COPY . /bio_final_kazasi_eftychia

# Set the entrypoint to bash for interactive use
ENTRYPOINT ["/bin/bash"]
