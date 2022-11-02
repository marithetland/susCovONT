FROM ubuntu:20.04

# Versions
# NB! Nextflow must be without the v prefix
ARG nextflow_ver="21.10.6"
ARG miniconda_sh="Miniconda3-py38_4.12.0-Linux-x86_64.sh"
ARG ncov_artic_ver="v1.3.0"
ARG artic_ver="1.2.3"
ARG primer_schemes_sha="9b533c7c29c3f08ec673aef6beced8ce98266267"

# Copy the files we need from our repo
WORKDIR /
COPY LICENSE  README.md susCovONT.py* ./
COPY ./scripts/ ./scripts/

# Setting bash as shell
SHELL ["/bin/bash", "-c"]

# Installing everything we need
ARG DEBIAN_FRONTEND=noninteractive
ARG TZ=Europe/Oslo
RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
    git=1:2.25.1-1ubuntu3.6 \
    ca-certificates=20211016~20.04.1 \
    python3=3.8.2-0ubuntu2 \
    openjdk-8-jre=8u342-b07-0ubuntu1~20.04 \
    wget=1.20.3-1ubuntu2 \
    docker.io=20.10.12-0ubuntu2~20.04.1 \
    && rm -rf /var/lib/apt/lists/*

# Installing nextflow
RUN wget -q -O nextflow https://github.com/nextflow-io/nextflow/releases/download/v${nextflow_ver}/nextflow-${nextflow_ver}-all && \
    bash nextflow && \
    chmod +x nextflow && \
    mv nextflow /bin

# Installing conda and mamba
ENV PATH=$PATH:/miniconda/bin
RUN wget -q -O miniconda.sh https://repo.anaconda.com/miniconda/$miniconda_sh && \
    bash miniconda.sh -b -p /miniconda && \
    rm miniconda.sh && \
    conda install -y mamba=0.15.3 -n base -c conda-forge \
    && mamba clean -a

# Installing the conda stuff that we want
RUN mamba install -y \
    biopython=1.78 \
    pandas=1.4.3 \
    numpy=1.23.1 \
    matplotlib=3.5.2 \
    && mamba install -y -c bioconda \
    samtools=1.6 \
    && mamba install -y -c bioconda -c conda-forge -c defaults \
    pangolin=4.1.3 \
    && mamba clean -a

# Setting up ncov2019-artic-nf pipeline and artic
RUN git clone -b $ncov_artic_ver https://github.com/connor-lab/ncov2019-artic-nf && \
    sed -i.bak "s/artic=1.1.3/artic=$artic_ver/g" /ncov2019-artic-nf/environments/nanopore/environment.yml && \
    cp /scripts/articQC.py /ncov2019-artic-nf/bin/qc.py
RUN mamba env create --prefix /conda_for_covid/work/conda/artic-2c6f8ebeb615d37ee3372e543ec21891 -f /ncov2019-artic-nf/environments/nanopore/environment.yml \
    && mamba clean -a

# Get primer schemes
RUN git clone https://github.com/markus-soma/primer-schemes.git
WORKDIR /primer-schemes
RUN git checkout $primer_schemes_sha
WORKDIR /

# Change config file
RUN sed -i.bak "s|^nf_location = .*|nf_location = /ncov2019-artic-nf/|" /scripts/config.cfg && \
    sed -i.bak "s|^conda_location = .*|conda_location = /conda_for_covid/work/conda/|" /scripts/config.cfg && \
    sed -i.bak "s|^schemeRepoURL = .*|schemeRepoURL = /primer-schemes/|" /scripts/config.cfg
