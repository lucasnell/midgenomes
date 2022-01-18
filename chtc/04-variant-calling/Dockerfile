############################################################
# Dockerfile for variant-calling steps in Tanytarsus gracilentus genome.
# Based on CentOS 7
# `docker build --platform linux/x86_64 -t lucasnell/tany_variant:v0 .`
############################################################

FROM centos:7
MAINTAINER Lucas Nell <lucas@lucasnell.com>
ADD . /tmp/repo
WORKDIR /tmp/repo
ENV LANG C.UTF-8

# Make RUN commands use `bash --login` (necessary for `conda init bash`)
SHELL ["/bin/bash", "--login", "-c"]

# Basic requirements:
RUN yum update -y
RUN yum group install -y "Development Tools"
RUN yum install -y locales wget curl git rsync which unzip bzip2-devel \
    lzma-devel ncurses-devel xz-devel zlib-devel perl pigz

# Set locales
RUN localedef -i en_US -f UTF-8 en_US.UTF-8

# Install mamba bc snakemake recommends this over standard conda
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
    bash Mambaforge-Linux-x86_64.sh -b -p /opt/mamba3 && \
    rm Mambaforge-Linux-x86_64.sh
ENV PATH /opt/mamba3/bin:${PATH}

# Add channels:
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Make primary conda environment:
RUN mamba create -q -y -n scaffolding-env

# Always start in our created environment
RUN echo "conda activate scaffolding-env" > ~/.bashrc
ENV PATH /opt/mamba3/envs/scaffolding-env/bin:${PATH}

# Install secondary (first two lines) and primary (last line) tools:
RUN mamba install -q -y snakemake singularity samtools minimap2 htslib \
    fastx_toolkit seqtk perl-bioperl \
    hisat2 busco longstitch
## besst

# Clean up tarballs from installed packages
RUN conda clean --all -y



