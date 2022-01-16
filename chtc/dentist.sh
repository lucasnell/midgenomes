#!/bin/bash

export THREADS=32

# installation steps for Mambaforge
cp /staging/lnell/Mambaforge-Linux-x86_64.sh ./
export HOME=$PWD
export PATH
sh Mambaforge-Linux-x86_64.sh -b -p $PWD/mamba3
export PATH=$PWD/mamba3/bin:$PATH
rm Mambaforge-Linux-x86_64.sh

export CONDA_PREFIX=$PWD/mamba3

# Install snakemake:
mamba install -c conda-forge -c bioconda --yes snakemake
# To narrow FASTA files below
mamba install -c bioconda fastx_toolkit

# # example ------------------------------
# # wget https://github.com/a-ludi/dentist/releases/download/v3.0.0/dentist-example.tar.gz
# cp /staging/lnell/dentist-example.tar.gz ./
# tar -xzf dentist-example.tar.gz
# cd dentist-example
#
# snakemake --configfile=snakemake.yml --use-conda --cores=${THREADS}
# md5sum -c checksum.md5
# # --------------------------------------


# The input FASTA is given by the submit file:
export GENOME=$1.fasta
if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist." 1>&2
    exit 1
fi
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz


# The input file dictates the output name:
export OUT_SUFFIX=dentist

if  [[ $GENOME == contigs* ]]
then
    OUT_DIR=scaffolds_${OUT_SUFFIX}
else
    OUT_DIR=${GENOME/.fasta/}_${OUT_SUFFIX}
fi
export OUT_DIR

export OUT_FASTA=${OUT_DIR}.fasta


mkdir ${OUT_DIR}
chmod +w -R ${OUT_DIR}

# Move all files from submit node to OUT_DIR:
mv dentist.v3.0.0.x86_64.tar.gz dentist.yml snakemake.yml ./${OUT_DIR}/
cd ${OUT_DIR}
# Adjust the snakemake.yml file for the input FASTA:
sed -i "s/INPUT_SCAFFOLDS_FILE/${GENOME}/" snakemake.yml

# It's important that this is a FASTA file!
READS=basecalls_guppy-5.0.11_filtered.fasta
cp /staging/lnell/${READS}.gz ./ && gunzip ${READS}.gz
# Apparently the fasta2DAM > DPsplit step can't handle a file with periods
mv ${READS} basecalls_guppy.fasta
export READS=basecalls_guppy.fasta



export TMPDIR=$(pwd)/tmp
mkdir ${TMPDIR}


# Prepare dentist:
tar -xzf dentist.v3.0.0.x86_64.tar.gz
rm dentist.v3.0.0.x86_64.tar.gz

# copy other necessary files:
cp -r -t . \
    dentist.v3.0.0.x86_64/Snakefile \
    dentist.v3.0.0.x86_64/envs \
    dentist.v3.0.0.x86_64/scripts


# Make sure the FASTA files aren't too wide
mv ${GENOME} ${GENOME/.fasta/_orig.fasta}
mv ${READS} ${READS/.fasta/_orig.fasta}
fasta_formatter -i ${GENOME/.fasta/_orig.fasta} -w 80 \
    -o ${GENOME}
fasta_formatter -i ${READS/.fasta/_orig.fasta} -w 80 \
    -o ${READS}
rm ${GENOME/.fasta/_orig.fasta} ${READS/.fasta/_orig.fasta}


# This section doesn't seem needed:
# # Now clean up sequence names
# mv ${GENOME} ${GENOME/.fasta/_narrow.fasta}
# mv ${READS} ${READS/.fasta/_narrow.fasta}
# # ... by removing everything after the first comma in the seq names:
# sed -r 's/\,.+//' ${GENOME/.fasta/_narrow.fasta} > ${GENOME}
# # ... by removing everything after the first space in the seq names:
# sed -r 's/\ .+//' ${READS/.fasta/_narrow.fasta} > ${READS}
# rm ${GENOME/.fasta/_narrow.fasta} ${READS/.fasta/_narrow.fasta}


snakemake --configfile=snakemake.yml --use-conda \
    --conda-prefix ${CONDA_PREFIX} --cores=${THREADS}


rm -rf dentist.v3.0.0.x86_64 ${READS} ${GENOME}

cp gap-closed.fasta ../

cd ..

mv gap-closed.fasta ${OUT_FASTA}
./summ_scaffs ${OUT_FASTA}
gzip ${OUT_FASTA}
mv ${OUT_FASTA}.gz /staging/lnell/

# ~~~~~~~~~~~~~
# For now, we'll just delete the main directory. Change this when you
# finalize the pipeline.
# ~~~~~~~~~~~~~
# tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
# mv ${OUT_DIR}.tar.gz /staging/lnell/

rm -rf ${TMPDIR} mamba3 ${OUT_DIR} summ_scaffs



