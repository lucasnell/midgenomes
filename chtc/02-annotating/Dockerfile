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
RUN yum update -y
RUN yum group install -y "Development Tools"
RUN yum install -y wget curl git which unzip bzip2-devel lzma-devel \
    ncurses-devel xz-devel zlib-devel perl

# Install mamba bc snakemake prefers this over standard conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    sh Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda3 && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/miniconda3/bin:${PATH}

# Add bioconda channels:
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# Make primary conda environment, then activate it:
RUN conda create -q -y -n annotating-env braker2 && \
    conda activate annotating-env
ENV PATH /opt/miniconda3/envs/annotating-env/bin:${PATH}

# Add a bunch of tools for use later
RUN mamba install -q -y -c conda-forge singularity && \
    mamba install -q -y -c bioconda busco && \
    mamba install -q -y -c bioconda samtools && \
    conda clean --all -y

# Install minimap2 (v2.24)
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 | tar -jxvf - && \
    mv minimap2-2.24_x64-linux minimap2-2.24 && \
    mv  minimap2-2.24 /opt/
ENV PATH /opt/minimap2-2.24:${PATH}

# Install hisat2 (v2.2.1)
RUN wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download/hisat2-2.2.1-Linux_x86_64.zip && \
    unzip hisat2-2.2.1-Linux_x86_64.zip && \
    mv hisat2-2.2.1 /opt/ && \
    rm hisat2-2.2.1-Linux_x86_64.zip
ENV PATH /opt/hisat2-2.2.1:${PATH}

# Always start in our created environment
RUN echo "conda activate annotating-env" > ~/.bashrc




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
