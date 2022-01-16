############################################################
# Dockerfile for scaffolding steps in Tanytarsus gracilentus genome.
# Based on CentOS 7
# `docker build --platform linux/x86_64 -t lucasnell/tany_scaff:v0 .`
############################################################

FROM centos:7
MAINTAINER Lucas Nell <lucas@lucasnell.com>
ADD . /tmp/repo
WORKDIR /tmp/repo
ENV LANG C.UTF-8
ENV SHELL /bin/bash
RUN yum update -y
RUN yum group install -y "Development Tools"
RUN yum install -y wget curl git which unzip bzip2-devel lzma-devel \
    ncurses-devel xz-devel zlib-devel perl

# Install mamba bc snakemake recommends this over standard conda
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
    bash Mambaforge-Linux-x86_64.sh -b -p /opt/mamba3 && \
    rm Mambaforge-Linux-x86_64.sh
ENV PATH /opt/mamba3/bin:${PATH}

RUN mamba init

# Add channels:
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Make primary conda environment, then activate it:
RUN mamba create -q -y -n scaffolding-env && \
    mamba activate scaffolding-env
ENV PATH /opt/mamba3/envs/scaffolding-env/bin:${PATH}

# Secondary tools that are necessary
RUN mamba install -q -y snakemake singularity samtools minimap2 \
    fastx_toolkit seqtk

# Main tools we're using
RUN mamba install -q -y hisat2 busco longstitch
## besst

# Clean up
RUN conda clean --all -y

# Always start in our created environment
RUN echo "mamba activate scaffolding-env" > ~/.bashrc




# The rest is for bioperl
# Much of this is based on
# https://github.com/bioperl/bioperl-docker/blob/master/bioperl-deps/Dockerfile


# Install perl modules
RUN yum install --yes cpanminus

RUN cpanm \
    CPAN::Meta \
    YAML \
    Digest::SHA \
    Module::Build \
    Test::Most \
    Test::Weaken \
    Test::Memory::Cycle \
    Clone

# Install what counts as BioPerl dependencies
RUN cpanm \
    HTML::TableExtract \
    Algorithm::Munkres \
    Array::Compare \
    Convert::Binary::C \
    Error \
    Graph \
    GraphViz \
    Inline::C \
    PostScript::TextBlock \
    Set::Scalar \
    Sort::Naturally \
    Math::Random \
    Spreadsheet::ParseExcel \
    IO::String \
    JSON \
    Data::Stag

# Specifically needed for P_RNA_Scaffolder
RUN cpanm \
    Getopt::Long \
    Bio::Seq \
    Bio::SeqIO \
    Bio::PrimarySeq
