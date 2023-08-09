############################################################
# Dockerfile for Tanytarsus gracilentus reference assembly.
############################################################


FROM condaforge/mambaforge:4.11.0-0
MAINTAINER Lucas Nell <lucnell@gmail.com>
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
        ocaml-nox zlib1g-dev libbz2-dev liblzma-dev libboost-all-dev \
        libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*

# Set locale
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8
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
# creating phylogeny and associated analyses:
COPY env-phylo.yml .
RUN mamba env create -q -f env-phylo.yml && \
    rm env-phylo.yml && \
    conda clean --all -y

# BUSCO:
RUN mamba create -q -y -c bioconda -c conda-forge -n busco-env busco=5.4.3 && \
    conda clean --all -y
# Flye assembler:
RUN mamba create -q -y -c bioconda -c conda-forge -n flye-env flye=2.9 && \
    conda clean --all -y


# --------------------------------
# TEclass for TE classification:
# note that TEclass scripts (notably TEclassTest.pl) are in /opt/TEclass and
# need to be run inside teclass-env
# --------------------------------
RUN mamba create -q -y -c bioconda -c conda-forge -n teclass-env \
        blast-legacy=2.2.26 glimmer=3.02 libsvm=330 && \
    conda clean --all -y && \
    cd /opt/ && \
    wget --no-check-certificate -c 'http://www.compgen.uni-muenster.de/tools/teclass/download/TEclass-2.1.3b.tar.gz' -O - \
        | tar -xz && \
    mv TEclass* TEclass && \
    cd TEclass && \
    wget -c 'https://github.com/cjlin1/libsvm/archive/refs/tags/v330.tar.gz' -O - | tar -xz && \
    mv libsvm-330/tools ./ && \
    mv tools libsvm_tools && \
    rm -r libsvm-330 && \
    wget -c 'https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/librf/librf.0.1.tar.gz' -O - \
        | tar -xz && \
    cd librf* && \
    ./configure && \
    perl -p -i -e 's/#include <vector>/#include <iostream>\n#include <vector>/' \
        librf/random_forest.h && \
    perl -p -i -e 's/#include <vector>/#include <cassert>\n#include <vector>/' \
        librf/discrete_dist.h && \
    perl -p -i -e 's/#include <fstream>/#include <cassert>\n#include <fstream>/' \
        librf/tree_node.h && \
    perl -p -i -e 's/#include <string>/#include <stdlib.h>\n#include <string>/' \
        librf/instance_set.h && \
    perl -p -i -e 's/#include <string>/#include <algorithm>\n#include <string>/' \
        librf/instance_set.h && \
    make && \
    cd .. && \
    wget -c 'http://www.cis.hut.fi/research/lvq_pak/lvq_pak-3.1.tar' -O - \
        | tar -x && \
    cd lvq_pak* && \
    perl -p -i -e 's/getline/lvq_getline/g' * && \
    perl -p -i -e 's/setprogname/lvq_setprogname/g' * && \
    if [ -s fileio.c.bak ]; then rm *.bak; fi && \
    cp makefile.unix makefile && \
    make && \
    cd .. && \
    perl Configure.pl && \
    sed -i 's/\$path_to_svm         = ""/\$path_to_svm         = "\/opt\/conda\/envs\/teclass-env\/bin"/g' Configuration.pm && \
    sed -i 's/\$path_to_glimmer     = ""/\$path_to_glimmer     = "\/opt\/conda\/envs\/teclass-env\/bin"/g' Configuration.pm && \
    sed -i 's/\$path_to_blastclust  = ""/\$path_to_blastclust  = "\/opt\/conda\/envs\/teclass-env\/bin"/g' Configuration.pm && \
    sed -i 's/\$path_to_scripts     = path_to_svm\/tools/\$path_to_scripts     = "\/opt\/TEclass\/libsvm_tools"/g' Configuration.pm && \
    ./Install.sh && \
    cd classifiers && \
    wget --no-check-certificate -c 'http://www.compgen.uni-muenster.de/tools/teclass/download/classifiers.tar.gz' -O - \
        | tar -xz && \
    cd .. && \
    perlies=$(find $(pwd) -maxdepth 1 -not -name "Configure.pl" -name "*.pl") && \
    for p in $perlies; do ln -s "$p" /opt/conda/envs/teclass-env/bin/; done && \
    sed -i '478s/\"2>\/dev\/null\"/\"1>\/dev\/null 2>\/dev\/null\"/' \
        ./TEclassTest.pl && \
    sed -i '539s/\"2>\/dev\/null\"/\"1>\/dev\/null 2>\/dev\/null\"/' \
        ./TEclassTest.pl && \
    sed -i '414s/\" 2>$path_to_outdir\/$file.$analysis_type.rf.log\"/\"1>$path_to_outdir\/$file.$analysis_type.rf.log 2>$path_to_outdir\/$file.$analysis_type.rf.log\"/' \
        ./TEclassTest.pl


# ---------------------------------
# Making RepeatMasker also use RepBase library
# ---------------------------------
COPY RepBaseRepeatMaskerEdition-20181026.tar.gz .
RUN tar -xzf RepBaseRepeatMaskerEdition-20181026.tar.gz && \
    mv ./Libraries/RMRBSeqs.embl /opt/conda/envs/repeat-env/share/RepeatMasker/Libraries/ && \
    rm -r Libraries RepBaseRepeatMaskerEdition-20181026.tar.gz && \
    cd /opt/conda/envs/repeat-env/share/RepeatMasker && \
    . /app/.bashrc && \
    conda activate repeat-env && \
    ./addRepBase.pl -libdir ./Libraries/


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
# install newest hyphy
# --------------------------------
RUN HYPHY_V="2.5.52" && \
    HYPHY_URL="https://github.com/veg/hyphy/archive/refs/tags/${HYPHY_V}.tar.gz" && \
    wget -c ${HYPHY_URL} -O - | tar -xz && \
    cd hyphy-${HYPHY_V} && \
    cmake -DCMAKE_INSTALL_PREFIX:PATH=/opt/hyphy/ . && \
    make install && \
    cd .. && \
    rm -r hyphy-${HYPHY_V}
ENV PATH="${PATH}:/opt/hyphy/bin"


# --------------------------------
# install codon-msa from inside hyphy-analyses
# --------------------------------
RUN cd /opt && \
    git clone https://github.com/veg/hyphy-analyses.git hyphy-analyses && \
    cd hyphy-analyses && \
    SHA1="4398a388c4349c431cfbb00b00ca358851c6f335" && \
    git reset --hard $SHA1 && \
    mv codon-msa /opt/ && \
    cd .. && \
    rm -r hyphy-analyses




# --------------------------------
# simple python scripts to summarize scaffold sequences and
# produce pretty CSV file, plus some helper functions in bash:
# --------------------------------
COPY summ-scaffs.py .
COPY pretty-csv.py .
COPY longest-isoforms.py .
COPY helpers.sh .
RUN chmod +x summ-scaffs.py pretty-csv.py longest-isoforms.py && \
    mv summ-scaffs.py pretty-csv.py longest-isoforms.py /usr/local/bin/ && \
    echo -e "\n\n" >> ~/.bashrc && \
    cat helpers.sh >> ~/.bashrc && \
    echo -e "\n\n" >> ~/.bashrc && \
    rm helpers.sh


