############################################################
# Dockerfile for annotating the Tanytarsus gracilentus genome.
############################################################

FROM centos:7
MAINTAINER Lucas Nell <lucas@lucasnell.com>
WORKDIR /app

# Environment configuration files:
COPY environment.yml .

# Name of main environment:
ARG MAIN_CONDA_ENV=annotating-env

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

# Make primary conda environment (~XXX sec)
RUN mamba env create -q -f environment.yml

# Always start in our created environment
ENV PATH /opt/mamba3/envs/$MAIN_CONDA_ENV/bin:${PATH}
RUN conda init bash
RUN echo "conda activate $MAIN_CONDA_ENV" >> ~/.bashrc

RUN conda clean --all -y


