############################################################
# Dockerfile for Tanytarsus gracilentus genomics.
############################################################


FROM condaforge/mambaforge:4.11.0-0
MAINTAINER Lucas Nell <lucas@lucasnell.com>
WORKDIR /app
ENV HOME /app

ARG UID=21659
ARG GID=21659

# Environment configuration files:
COPY main-env.yml .
COPY annotate-env.yml .
COPY besst-env.yml .
# Simple summary of scaffolds:
COPY summ-scaffs.py .

SHELL ["/bin/bash", "--login", "-c"]

RUN groupadd -g $GID --non-unique tanies && \
    useradd --create-home --shell /bin/bash --uid $UID --gid $GID \
        --non-unique tany && \
    echo "tany:tany" | chpasswd && \
    usermod -aG sudo tany


# --------------------------------
# Making everything accessible on the cluster where they make you non-root
# --------------------------------
# I want to always be able to access and edit these.
# I'm doing this early in the build so that it doesn't make the image massive.
# I'm also adding `/opt/conda/pkgs` now bc otherwise its permissions don't
# share those of `/opt/conda`
RUN mkdir /opt/conda/pkgs && \
    chown -R tany:tanies /opt/conda && \
    mkdir /app/.conda && \
    chown -R tany:tanies /app


# Basic requirements (~360 sec)
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

# simple python script to summarize scaffold sequences:
RUN chmod +x summ-scaffs.py && \
    mv summ-scaffs.py /usr/local/bin/

USER tany

# Make primary conda environment (~300 sec)
RUN mamba env create -q -f main-env.yml && \
    rm main-env.yml && \
    conda clean --all -y

# new environments for other programs that need separate environments:

# for BESST_RNA scaffolder (~150 sec)
RUN mamba env create -q -f besst-env.yml && \
    rm besst-env.yml && \
    conda clean --all -y
# for BUSCO (~400 sec)
RUN mamba create -q -y -c bioconda -c conda-forge -n busco-env busco=5.2.2 && \
    conda clean --all -y
# for LongStitch scaffolder (~150 sec)
RUN mamba create -q -y -c bioconda -c conda-forge -n longstitch-env \
    longstitch=1.0.1 && \
    conda clean --all -y
# for genome annotating using BRAKER2 and RepeatModeler (~400 sec)
RUN mamba env create -q -f annotate-env.yml && \
    rm annotate-env.yml && \
    conda clean --all -y

RUN conda init bash

# --------------------------------
# Manually install MAPGD inside main-env environment:
# --------------------------------
RUN wget https://github.com/LynchLab/MAPGD/archive/refs/tags/v0.4.26.zip && \
    unzip -q v0.4.26.zip && \
    rm -f v0.4.26.zip && \
    . ~/.bashrc && \
    conda activate main-env && \
    cd ./MAPGD-0.4.26 && \
    ./configure --prefix="/opt/conda/envs/main-env" && \
    make && \
    make install DESTDIR=/opt/conda/envs/main-env/bin && \
    cd .. && \
    rm -r ./MAPGD-0.4.26


# --------------------------------
# Add scripts for BESST_RNA:
# --------------------------------
# Only the `src` folder needed from GitHub repo.
# I'm using the newest commit as of 29 Jan 2022.
# (It hasn't been updated in 6 years, but just in case...)
RUN wget https://github.com/ksahlin/BESST_RNA/archive/9cc039b5db080e98361a9b54dd333a7bb8b35223.zip && \
    unzip 9cc039b5db080e98361a9b54dd333a7bb8b35223.zip && \
    rm 9cc039b5db080e98361a9b54dd333a7bb8b35223.zip && \
    cd BESST_RNA-9cc039b5db080e98361a9b54dd333a7bb8b35223 && \
    mv src ../ && \
    cd .. && \
    rm -rf BESST_RNA-9cc039b5db080e98361a9b54dd333a7bb8b35223 && \
    mv src BESST_RNA && \
    chown -R tany:tanies BESST_RNA


# This allows non-root user to create new environments and install new packages:
##RUN chown -R tany:tanies /opt/conda/pkgs/urls.txt




# /opt/conda/pkgs/urls.txt


