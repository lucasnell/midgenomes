FROM centos:7
ADD . /tmp/repo
WORKDIR /tmp/repo
ENV LANG C.UTF-8
ENV SHELL /bin/bash
RUN yum update -y && \
    yum group install -y "Development Tools" && \
    yum install -y wget curl git which unzip bzip2-devel lzma-devel ncurses-devel xz-devel zlib-devel
RUN /bin/bash -c "wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
    bash Mambaforge-Linux-x86_64.sh -b -p /opt/mamba3 && \
    rm Mambaforge-Linux-x86_64.sh"
ENV PATH /opt/mamba3/bin:${PATH}
RUN /bin/bash -c "mamba create -q -y -c conda-forge -c bioconda -n snakemake-env snakemake && \
    source activate snakemake-env && \
    mamba install -q -y -c conda-forge singularity && \
    mamba install -q -y -c bioconda busco && \
    mamba install -q -y -c bioconda seqtk && \
    mamba install -q -y -c bioconda -c conda-forge longstitch && \
    mamba install -q -y -c bioconda fastx_toolkit && \
    conda clean --all -y"
RUN /bin/bash -c "wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 && \
    tar -xf samtools-1.14.tar.bz2 && \
    rm samtools-1.14.tar.bz2 && \
    mkdir /opt/samtools-1.14 && \
    cd samtools-1.14 && \
    ./configure --prefix=/opt/samtools-1.14 && \
    make && \
    make install && \
    cd .. && \
    rm -r samtools-1.14"
ENV PATH /opt/samtools-1.14:${PATH}
RUN /bin/bash -c "curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 | tar -jxvf - && \
    mv minimap2-2.24_x64-linux minimap2-2.24 && \
    mv  minimap2-2.24 /opt/"
ENV PATH /opt/minimap2-2.24:${PATH}
RUN /bin/bash -c "wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download/hisat2-2.2.1-Linux_x86_64.zip && \
    unzip hisat2-2.2.1-Linux_x86_64.zip && \
    mv hisat2-2.2.1 /opt/ && \
    rm hisat2-2.2.1-Linux_x86_64.zip"
ENV PATH /opt/hisat2-2.2.1:${PATH}
RUN echo "source activate snakemake-env" > ~/.bashrc
ENV PATH /opt/mamba3/envs/snakemake-env/bin:${PATH}
