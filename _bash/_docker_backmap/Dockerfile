############################################################
# Dockerfile for using backmap on Tanytarsus gracilentus reference assembly.
# This is in a separate docker container bc the qualimap conda install hangs
# for some reason.
############################################################


FROM condaforge/mambaforge:4.11.0-0
MAINTAINER Lucas Nell <lucnell@gmail.com>
WORKDIR /app
ENV HOME /app
ARG DEBIAN_FRONTEND=noninteractive

SHELL ["/bin/bash", "--login", "-c"]


# Basic requirements
RUN apt-get clean && \
    apt-get update -y && \
    apt-get install -y build-essential && \
    apt-get install -y locales wget rsync unzip nano sudo && \
    rm -rf /var/lib/apt/lists/*

# Set locale
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8
ENV TZ America/Chicago
# ENV TZ America/Los_Angeles

RUN conda init bash && \
    . /app/.bashrc && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda

RUN . /app/.bashrc && \
    mamba install -q -y \
    bedtools=2.30.0 \
    bwa=0.7.17 \
    minimap2=2.24 \
    multiqc=1.12-0 \
    perl-app-cpanminus=1.7039 \
    perl-parallel-loops=0.10 \
    python=3.9.9 \
    r-base=4.1.2 \
    samtools=1.14

RUN . /app/.bashrc && \
    mamba install -q -y \
        fonts-conda-ecosystem \
        openjdk \
        bioconductor-noiseq \
        r-optparse \
        r-xml \
        bioconductor-rsamtools \
        bioconductor-rtracklayer

RUN cd /opt && \
    wget -q https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip && \
    unzip qualimap_v2.2.1.zip && \
    rm qualimap_v2.2.1.zip && \
    mv qualimap_v2.2.1 qualimap
ENV PATH="${PATH}:/opt/qualimap"


# One extra perl module for backmap.pl, then add backmap.pl itself
RUN . /app/.bashrc && \
    cpanm Number::FormatEng && \
    BACKMAP_URL="https://github.com/schellt/backmap/archive/refs/tags/v0.5.tar.gz" && \
    wget -c "${BACKMAP_URL}" -O - | tar -xz && \
    mv backmap* backmap && \
    cd backmap && \
    chmod +x backmap.pl && \
    mv backmap.pl /usr/local/bin/ && \
    cd .. && \
    rm -rf backmap

