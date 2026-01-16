# BAMS3 Nextflow Pipeline Container
#
# Contains all tools needed for the pipeline:
# - BAMS3
# - BWA
# - GATK
# - bcftools
# - samtools
# - MultiQC
#

FROM ubuntu:22.04

LABEL maintainer="BAMS3 Project"
LABEL description="Container for BAMS3 Nextflow pipeline"
LABEL version="1.0.0"

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    git \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    python3 \
    python3-pip \
    openjdk-17-jdk \
    unzip \
    tabix \
    jq \
    && rm -rf /var/lib/apt/lists/*

# Install AWS CLI (for S3 access)
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm -rf aws awscliv2.zip

# Install BWA
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
    && tar -xjf bwa-0.7.17.tar.bz2 \
    && cd bwa-0.7.17 \
    && make \
    && cp bwa /usr/local/bin/ \
    && cd .. \
    && rm -rf bwa-0.7.17*

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
    && tar -xjf samtools-1.18.tar.bz2 \
    && cd samtools-1.18 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf samtools-1.18*

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2 \
    && tar -xjf bcftools-1.18.tar.bz2 \
    && cd bcftools-1.18 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install \
    && cd .. \
    && rm -rf bcftools-1.18*

# Install GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip \
    && unzip gatk-4.4.0.0.zip \
    && mv gatk-4.4.0.0 /opt/gatk \
    && ln -s /opt/gatk/gatk /usr/local/bin/gatk \
    && rm gatk-4.4.0.0.zip

# Install MultiQC
RUN pip3 install --no-cache-dir multiqc

# Install Go (for building BAMS3)
RUN wget https://go.dev/dl/go1.21.5.linux-amd64.tar.gz \
    && tar -C /usr/local -xzf go1.21.5.linux-amd64.tar.gz \
    && rm go1.21.5.linux-amd64.tar.gz

ENV PATH="/usr/local/go/bin:${PATH}"

# Build and install BAMS3
RUN git clone https://github.com/scttfrdmn/bams3-go.git /tmp/bams3 \
    && cd /tmp/bams3 \
    && go build -o bams3 ./cmd/bams3 \
    && cp bams3 /usr/local/bin/ \
    && cd / \
    && rm -rf /tmp/bams3

# Verify installations
RUN bwa 2>&1 | head -3 \
    && samtools --version \
    && bcftools --version \
    && gatk --version \
    && bams3 --version \
    && multiqc --version

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
