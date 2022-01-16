############################################################
# Dockerfile for annotating the Tanytarsus gracilentus genome.
# Based on CentOS 7
# `docker build --platform linux/x86_64 -t lucasnell/tany_anno:v0 .`
############################################################

FROM centos:7
MAINTAINER Lucas Nell <lucas@lucasnell.com>
ADD . /tmp/repo
WORKDIR /tmp/repo
ENV LANG C.UTF-8
ENV SHELL /bin/bash

# Make RUN commands use `bash --login` (necessary for `conda init bash`)
SHELL ["/bin/bash", "--login", "-c"]

# Basic requirements:
RUN yum update -y
RUN yum group install -y "Development Tools"
RUN yum install -y wget curl git rsync which unzip bzip2-devel lzma-devel \
    ncurses-devel xz-devel zlib-devel perl

# Install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    sh Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/miniconda3/bin:${PATH}

# Add bioconda channels:
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Make primary conda environment:
RUN conda create -q -y -n annotating-env

# Always start in our created environment
RUN echo "conda activate annotating-env" > ~/.bashrc
ENV PATH /opt/miniconda3/envs/annotating-env/bin:${PATH}

# Main tools we're using:
RUN conda install -q -y braker2 hisat2 busco samtools

# Clean up:
RUN conda clean --all -y

