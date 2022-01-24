############################################################
# Dockerfile for Tanytarsus gracilentus genomics.
############################################################

FROM centos:7
MAINTAINER Lucas Nell <lucas@lucasnell.com>
WORKDIR /app

# Environment configuration files:
COPY besst-env.yml .
COPY braker2-env.yml .
COPY main-env.yml .
COPY variant-env.yml .
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


# Install mamba-using equivalent to miniconda from conda-forge (~90 sec)
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
    bash Mambaforge-Linux-x86_64.sh -b -p /opt/mamba3 && \
    rm Mambaforge-Linux-x86_64.sh
ENV PATH /opt/mamba3/bin:${PATH}

# Make primary conda environment (~330 sec)
RUN mamba env create -q -f main-env.yml && \
    rm main-env.yml

# Allow `conda activate`
RUN conda init bash

# new environments for other programs that have different requirements:

# for BESST_RNA scaffolder (~150 sec)
RUN mamba env create -q -f besst-env.yml && \
    rm besst-env.yml
# for BUSCO (~400 sec)
RUN mamba create -q -y -c bioconda -c conda-forge -n busco-env busco=5.2.2
# for BRAKER2 annotating (~XXX sec)
RUN mamba env create -q -f braker2-env.yml && \
    rm braker2-env.yml
# for variant calling (~XXX sec)
RUN mamba env create -q -f variant-env.yml && \
    rm variant-env.yml


# clean up tarballs:
RUN conda clean --all -y

# compile simple cpp script to summarize scaffold sequences:
RUN g++ -std=c++11 summ_scaffs.cpp -o summ_scaffs && \
    mv summ_scaffs /usr/bin/ && \
    rm summ_scaffs.cpp

# --------------------------------
# Manually install MAPGD inside variant-env environment:
# --------------------------------
RUN wget https://github.com/LynchLab/MAPGD/archive/master.zip && \
    unzip -q master.zip && \
    rm -f master.zip
# In this Makefile for gzstream, I'm going to change
# CPPFLAGS = -I. -O
# to
# CPPFLAGS = -I. -I/opt/mamba3/envs/variant-env/include -O
# or it will complain about not having access to zlib.h
# There's probably a more elegant solution, but <shrugs>...
RUN cd ./MAPGD-master/src/mapgd_0.4/gzstream && \
    sed -i "s/CPPFLAGS = -I. -O/CPPFLAGS = -I. -I\/opt\/mamba3\/envs\/variant-env\/include -O/" Makefile
# ~400 sec
RUN conda activate variant-env && \
    cd ./MAPGD-master && \
    ./configure --prefix="/opt/mamba3/envs/variant-env" && \
    make && \
    make install && \
    cd .. && \
    rm -r MAPGD-master



