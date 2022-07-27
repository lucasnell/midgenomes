############################################################
# Dockerfile for Tanytarsus gracilentus reference assembly.
############################################################


FROM condaforge/mambaforge:4.11.0-0
MAINTAINER Lucas Nell <lucas@lucasnell.com>
WORKDIR /app
ENV HOME /app
ARG DEBIAN_FRONTEND=noninteractive

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
# ENV TZ America/Chicago
ENV TZ America/Los_Angeles

RUN conda init bash

# primary conda environment with most general tools:
COPY env-main.yml .
RUN mamba env create -q -f env-main.yml && \
    rm env-main.yml && \
    conda clean --all -y
# genome assembly:
COPY env-assembly.yml .
RUN mamba env create -q -f env-assembly.yml && \
    rm env-assembly.yml && \
    conda clean --all -y
# repeat library creation using RepeatModeler, masking using RepeatMasker
COPY env-repeat.yml .
RUN mamba env create -q -f env-repeat.yml && \
    rm env-repeat.yml && \
    conda clean --all -y
# gene predictions using BRAKER2 and functional annotation using mantis
COPY env-annotate.yml .
RUN mamba env create -q -f env-annotate.yml && \
    rm env-annotate.yml && \
    conda clean --all -y
# Additional steps to get mantis working
RUN . /app/.bashrc && \
    conda activate annotate-env && \
    mantis -h 1> /dev/null && \
    cd $(python -c 'import site; print(site.getsitepackages()[0])') && \
    mkdir setup_databases && \
    chmod -R 777 setup_databases && \
    sed -i 's/\#default_ref_folder=/default_ref_folder=\/app\/mantis_dbs/g' \
        ./config/MANTIS.cfg && \
    cd /app && \
    rm -rf build
# creating phylogeny:
COPY env-phylo.yml .
RUN mamba env create -q -f env-phylo.yml && \
    rm env-phylo.yml && \
    conda clean --all -y

# BUSCO
RUN mamba create -q -y -c bioconda -c conda-forge -n busco-env busco=5.3.1 && \
    conda clean --all -y
# Flye assembler:
RUN mamba create -q -y -c bioconda -c conda-forge -n flye-env flye=2.9 && \
    conda clean --all -y



# ---------------------------------
# Adjusting to allow BRAKER2 to work
# ---------------------------------
# For GeneMark:
COPY gm_key_64.gz .
COPY gmes_linux_64_4.tar.gz .
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
# For ProtHint:
RUN cd /opt &&\
    PROTHINT_URL="https://github.com/gatech-genemark/ProtHint/releases/download/v2.6.0/ProtHint-2.6.0.tar.gz" && \
    wget -c ${PROTHINT_URL} -O - | tar -xz && \
    mv ProtHint* ProtHint
ENV PROTHINT_PATH=/opt/ProtHint/bin


# --------------------------------
# install NextPolish and NextDenovo
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
# install TSEBRA
# --------------------------------
RUN cd /opt && \
    TSEBRA_URL="https://github.com/Gaius-Augustus/TSEBRA/archive/refs/tags/v1.0.3.tar.gz" && \
    wget -c ${TSEBRA_URL} -O - | tar -xz && \
    mv TSEBRA* TSEBRA
ENV PATH="${PATH}:/opt/TSEBRA/bin"



# --------------------------------
# simple python scripts to summarize scaffold sequences and
# produce pretty CSV file, plus some helper functions in bash:
# --------------------------------
COPY summ-scaffs.py .
COPY pretty-csv.py .
COPY helpers.sh .
RUN chmod +x summ-scaffs.py pretty-csv.py && \
    mv summ-scaffs.py pretty-csv.py /usr/local/bin/ && \
    echo -e "\n\n" >> ~/.bashrc && \
    cat helpers.sh >> ~/.bashrc && \
    echo -e "\n\n" >> ~/.bashrc && \
    rm helpers.sh


