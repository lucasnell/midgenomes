############################################################
# Dockerfile for scaffolding steps in Tanytarsus gracilentus genome.
# Based on CentOS 7
############################################################

FROM centos:7
MAINTAINER Lucas Nell <lucas@lucasnell.com>
WORKDIR /app
ENV HOME /app

# Environment configuration files:
COPY scaffolding-env.yml .
COPY besst-env.yml .
COPY summ_scaffs.cpp .

SHELL ["/bin/bash", "--login", "-c"]

# Basic requirements (~360 sec)
RUN yum update -y && \
    yum group install -y "Development Tools" && \
    yum install -y locales wget curl git rsync which unzip nano sudo

# Set locale
RUN localedef -i en_US -f UTF-8 en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANGUAGE en_US:en

# compile simple cpp script to summarize scaffold sequences:
RUN g++ -std=c++11 summ_scaffs.cpp -o summ_scaffs && \
    mv summ_scaffs /usr/bin/ && \
    rm summ_scaffs.cpp

# Install mamba-using equivalent to miniconda from conda-forge (~90 sec)
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
    bash Mambaforge-Linux-x86_64.sh -b -p /opt/mamba3 && \
    rm Mambaforge-Linux-x86_64.sh
ENV PATH /opt/mamba3/bin:${PATH}

# Allow `conda activate`
RUN conda init bash

# Make primary conda environment (~330 sec)
RUN mamba env create -q -f scaffolding-env.yml && \
    rm scaffolding-env.yml

# new environments for the two programs that have different requirements:

# for BESST_RNA (~150 sec)
RUN mamba env create -q -f besst-env.yml && \
    rm besst-env.yml

# for LongStitch (~400 sec)
RUN mamba create -q -y -c bioconda -c conda-forge -n longstitch-env longstitch=1.0.1

# for BUSCO (~400 sec)
RUN mamba create -q -y -c bioconda -c conda-forge -n busco-env busco=5.2.2

RUN conda clean --all -y


# I want to always be able to access and edit these:
RUN chmod -R 777 /opt/mamba3 && \
    chmod -R 777 /app && \
    chmod -R 777 /home

# To allow conda use for new users (even on cluster)
COPY entry.sh .
RUN chmod +x entry.sh
ENTRYPOINT ["/app/entry.sh"]
