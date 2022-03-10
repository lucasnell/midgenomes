############################################################
# Dockerfile for Tanytarsus gracilentus genomics.
############################################################


FROM condaforge/mambaforge:4.11.0-0
MAINTAINER Lucas Nell <lucas@lucasnell.com>
WORKDIR /app
ENV HOME /app

# Environment configuration files:
COPY main-env.yml .
COPY annotate-env.yml .
COPY besst-env.yml .
# Simple summary of scaffolds:
COPY summ-scaffs.py .
# To get GeneMark to work for BRAKER2:
COPY gm_key_64.gz .
COPY gmes_linux_64_4.tar.gz .


SHELL ["/bin/bash", "--login", "-c"]

# --------------------------------
# Making everything accessible on the cluster where they make you non-root
# --------------------------------
# I want to always be able to access these:
RUN chmod -R 777 /opt/conda && \
    chmod -R 777 /app


# Basic requirements (~360 sec)
RUN apt-get clean && \
    apt-get update -y && \
    apt-get install -y build-essential && \
    apt-get install -y locales wget rsync unzip nano sudo ocaml-nox && \
    rm -rf /var/lib/apt/lists/*

# Set locale
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

# simple python script to summarize scaffold sequences:
RUN chmod +x summ-scaffs.py && \
    mv summ-scaffs.py /usr/local/bin/

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
# Manually install SNAPE-pooled:
# --------------------------------
RUN cd /opt && \
    git clone https://github.com/EmanueleRaineri/snape-pooled.git && \
    cd snape-pooled && \
    make snape-pooled
ENV PATH="/opt/snape-pooled:${PATH}"



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
    chmod -R 777 BESST_RNA

# ---------------------------------
# Adjusting to allow BRAKER2 to work
# ---------------------------------
# For GeneMark:
RUN gunzip gm_key_64.gz && \
    mv gm_key_64 .gm_key && \
    tar -xzf gmes_linux_64_4.tar.gz && \
    rm gmes_linux_64_4.tar.gz && \
    cd gmes_linux_64_4 && \
    perl change_path_in_perl_scripts.pl "/usr/bin/env perl" && \
    cd .. && \
    mv gmes_linux_64_4 GeneMark-EX && \
    mv GeneMark-EX /opt/
ENV GENEMARK_PATH=/opt/GeneMark-EX
# For AUGUSTUS:
RUN chmod -R 777 /opt/conda/envs/annotate-env/config
