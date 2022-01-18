############################################################
# Dockerfile for scaffolding steps in Tanytarsus gracilentus genome.
# Based on CentOS 7
# `docker build --platform linux/x86_64 -t lucasnell/tany_scaff:v0 .`
############################################################

FROM centos:7
MAINTAINER Lucas Nell <lucas@lucasnell.com>
WORKDIR /app

# Environment configuration files:
COPY environment.yml .
COPY besst-environment.yml .
COPY longstitch-environment.yml .

# Name of main environment:
ARG conda_env=scaffolding-env

SHELL ["/bin/bash", "--login", "-c"]


# Basic requirements:
RUN yum update -y && \
    yum group install -y "Development Tools" && \
    yum install -y locales wget curl git rsync which unzip nano \
        bzip2-devel lzma-devel ncurses-devel xz-devel zlib-devel

# Set locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANGUAGE en_US:en

# Install mamba bc snakemake recommends this over standard conda
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
    bash Mambaforge-Linux-x86_64.sh -b -p /opt/mamba3 && \
    rm Mambaforge-Linux-x86_64.sh
ENV PATH /opt/mamba3/bin:${PATH}

# Make primary conda environment (~500 sec)
RUN mamba env create -q -f environment.yml

# Always start in our created environment
ENV PATH /opt/mamba3/envs/$conda_env/bin:${PATH}
RUN conda init bash
RUN echo "conda activate $conda_env" >> ~/.bashrc

# Test for which environment starts:
RUN echo $CONDA_DEFAULT_ENV && sleep 10

## # ~300 sec
## RUN mamba install -q -y -c bioconda -c conda-forge snakemake singularity \
##     samtools minimap2
##
## # ~90 sec
## RUN mamba install -q -y -c bioconda -c conda-forge seqtk fastp hisat2 \
##     perl-bioperl
##
## # ~ 220 sec
## RUN mamba install -q -y -c bioconda -c conda-forge busco


# new environments for the two programs that have different requirements:

RUN mamba env create -q -f besst-environment.yml

RUN mamba env create -q -f longstitch-environment.yml


## # ~120 sec
## RUN mamba create -q -y -c r -c bioconda -c free -n besst-env python=2.7.10 \
##         pysam matplotlib networkx=1.10 besst
##
## # ~180 sec
## RUN mamba create -q -y -c defaults -c bioconda -c conda-forge -n longstitch-env \
##         python=3.9.9 longstitch


RUN conda clean --all -y


