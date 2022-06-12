############################################################
# Dockerfile for Tanytarsus gracilentus reference assembly.
############################################################


FROM condaforge/mambaforge:4.11.0-0
MAINTAINER Lucas Nell <lucas@lucasnell.com>
WORKDIR /app
ENV HOME /app
ARG DEBIAN_FRONTEND=noninteractive

# Environment configuration files:
COPY env-main.yml .
COPY env-annotate.yml .
COPY env-assembly.yml .
COPY env-backmap.yml .
COPY env-phylo.yml .
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


# Basic requirements
RUN apt-get clean && \
    apt-get update -y && \
    apt-get install -y build-essential && \
    apt-get install -y locales wget rsync unzip nano sudo pv cmake \
        ocaml-nox zlib1g-dev libbz2-dev liblzma-dev libboost-all-dev && \
    rm -rf /var/lib/apt/lists/*

# Set locale
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8
ENV TZ America/Chicago
# ENV TZ America/Los_Angeles


# Make primary conda environment and other main groups of environments:

RUN mamba env create -q -f env-main.yml && \
    rm env-main.yml && \
    conda clean --all -y
# genome assembly:
RUN mamba env create -q -f env-assembly.yml && \
    rm env-assembly.yml && \
    conda clean --all -y
# for genome annotating using BRAKER2 and RepeatModeler
RUN mamba env create -q -f env-annotate.yml && \
    rm env-annotate.yml && \
    conda clean --all -y
# for creating phylogeny:
RUN mamba env create -q -f env-phylo.yml && \
    rm env-phylo.yml && \
    conda clean --all -y


# new environments for other programs that need separate environments
# (i.e., they don't play well with others):

# for BUSCO
RUN mamba create -q -y -c bioconda -c conda-forge -n busco-env busco=5.3.1 && \
    conda clean --all -y
# for Flye assembler:
RUN mamba create -q -y -c bioconda -c conda-forge -n flye-env flye=2.9 && \
    conda clean --all -y
# for backmap.pl:
RUN mamba env create -q -f env-backmap.yml && \
    rm env-backmap.yml && \
    conda clean --all -y

RUN conda init bash

# One extra perl module for backmap.pl, then add backmap.pl itself
RUN . /app/.bashrc && \
    conda activate backmap-env && \
    cpanm Number::FormatEng && \
    wget -c https://github.com/schellt/backmap/archive/refs/tags/v0.5.tar.gz \
        -O - | tar -xz && \
    cd backmap-0.5 && \
    chmod +x backmap.pl && \
    mv backmap.pl /opt/conda/envs/backmap-env/bin/ && \
    cd .. && \
    rm -rf backmap-0.5



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


# --------------------------------
# Manually install NextPolish and NextDenovo
# --------------------------------
RUN pip install paralleltask && \
    cd /opt && \
    wget -q https://github.com/Nextomics/NextPolish/releases/download/v1.4.0/NextPolish.tgz && \
    tar -xzf NextPolish.tgz && \
    rm NextPolish.tgz && \
    cd NextPolish && \
    make && \
    cd .. && \
    wget -q https://github.com/Nextomics/NextDenovo/releases/download/v2.5.0/NextDenovo.tgz && \
    tar -xzf NextDenovo.tgz && \
    rm NextDenovo.tgz
ENV PATH="${PATH}:/opt/NextPolish:/opt/NextDenovo"



# --------------------------------
# Manually install treePL
# --------------------------------
# Dependencies:
RUN cd /app && \
    ADOL_C="https://github.com/coin-or/ADOL-C/archive/refs/tags/releases/2.7.2.tar.gz" && \
    wget -c ${ADOL_C} -O - | tar -xz && \
    mv ADOL-C* ADOL-C && \
    cd ADOL-C && \
    ./configure --prefix=/usr/local --with-openmp-flag=-fopenmp && \
    make && \
    make install && \
    cd /app && \
    NLOPT="https://github.com/stevengj/nlopt/archive/v2.7.1.tar.gz" && \
    wget -c ${NLOPT} -O - | tar -xz && \
    mv nlopt* nlopt && \
    cd nlopt && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install && \
    cd /app && \
    rm -r nlopt ADOL-C
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib:/usr/local/lib64"
# treePL itself:
RUN git clone https://github.com/blackrim/treePL.git && \
    cd treePL/src && \
    ./configure && \
    make && \
    make install && \
    cd ../.. && \
    rm -r treePL


# --------------------------------
# simple python scripts to summarize scaffold sequences and
# produce pretty CSV file, plus some helper functions in bash:
# --------------------------------
COPY summ-scaffs.py .
COPY pretty-csv.py .
COPY helpers.sh .
RUN chmod +x summ-scaffs.py pretty-csv.py helpers.sh && \
    mv summ-scaffs.py pretty-csv.py helpers.sh /usr/local/bin/ && \
    echo -e "\nsource /usr/local/bin/helpers.sh\n" >> ~/.bashrc


