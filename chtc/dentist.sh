#!/bin/bash

export THREADS=32

# installation steps for Mambaforge
cp /staging/lnell/Mambaforge-Linux-x86_64.sh ./
# wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
export HOME=$PWD
export PATH
sh Mambaforge-Linux-x86_64.sh -b -p $PWD/mamba3
export PATH=$PWD/mamba3/bin:$PATH
rm Mambaforge-Linux-x86_64.sh

export CONDA_PREFIX=$PWD/mamba3

# Install snakemake environment
# (`seqtk` is to convert FASTQ to FASTA and to narrow FASTA files below)
mamba create -q -y -c bioconda -c conda-forge -n snake-env \
    snakemake=6.13.1 seqtk=1.3

# Environment for busco, too:
mamba create -q -y -c bioconda -c conda-forge -n busco-env busco=5.2.2

conda init
. ~/.bashrc

conda activate snake-env

# # example ------------------------------
# wget https://github.com/a-ludi/dentist/releases/download/v3.0.0/dentist-example.tar.gz
# tar -xzf dentist-example.tar.gz
# cd dentist-example
# snakemake --configfile=snakemake.yml --use-conda --cores=1
# md5sum -c checksum.md5
# # --------------------------------------


# The input FASTA is given by the submit file:
export GENOME=$1.fasta
if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist." 1>&2
    exit 1
fi


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
## chmod +w -R ${OUT_DIR}

# Move all files from submit node to OUT_DIR:
mv dentist.v3.0.0.x86_64.tar.gz dentist.yml snakemake.yml fire-gusu.py ./${OUT_DIR}/
cd ${OUT_DIR}
# Adjust the snakemake.yml file for the input FASTA:
sed -i "s/INPUT_SCAFFOLDS_FILE/${GENOME}/" snakemake.yml

# Move and adjust the genome FASTA file:
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz
# Make sure the genome FASTA file isn't too wide
mv ${GENOME} ${GENOME/.fasta/_orig.fasta}
seqtk seq -l 80 -A ${GENOME/.fasta/_orig.fasta} > ${GENOME}
rm ${GENOME/.fasta/_orig.fasta}


# Final FASTA reads file name.
# Make sure it doesn't have any periods, because apparently the
# fasta2DAM > DPsplit step can't handle a file with periods
export READS=basecalls_guppy.fasta
# Get all reads, filter them (dentist takes way too long with all the reads),
# and convert to FASTA with 80-char lines
export ALL_READS=basecalls_guppy-5.0.11.fastq.gz
cp /staging/lnell/${ALL_READS} ./
# requested coverage:
export COVERAGE=$(grep "read-coverage:" dentist.yml | \
                  sed -e 's/[[:space:]]*$//' | \
                  sed 's/.* //')
# info on the reads in the FASTQ file:
export SUMMARY=basecalls_guppy-5.0.11--sequencing_summary.txt.gz
cp /staging/lnell/${SUMMARY} ./

# Filtering for average quality of >= 10 and length of >= 10 kb, then randomly
# choosing reads to get the desired coverage:
./fire-gusu.py -s ${SUMMARY} -c ${COVERAGE} -g 100 -q 10.0 -l 10000 \
    -o ${READS/.fasta/.fastq/} --seed 668561223 ${ALL_READS}
rm ${ALL_READS} ${SUMMARY}
# Convert to 80-char-wide FASTA
seqtk seq -l 80 -A ${READS/.fasta/.fastq/} > ${READS}
rm ${READS/.fasta/.fastq/}



# This section doesn't seem needed:
# # Now clean up sequence names
# mv ${GENOME} ${GENOME/.fasta/_narrow.fasta}
# mv ${READS} ${READS/.fasta/_narrow.fasta}
# # ... by removing everything after the first comma in the seq names:
# sed -r 's/\,.+//' ${GENOME/.fasta/_narrow.fasta} > ${GENOME}
# # ... by removing everything after the first space in the seq names:
# sed -r 's/\ .+//' ${READS/.fasta/_narrow.fasta} > ${READS}
# rm ${GENOME/.fasta/_narrow.fasta} ${READS/.fasta/_narrow.fasta}



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


snakemake --configfile=snakemake.yml --use-conda \
    --conda-prefix ${CONDA_PREFIX} --cores=${THREADS}


rm -rf dentist.v3.0.0.x86_64 ${READS} ${GENOME}

if [ ! -f gap-closed.fasta ]; then
    echo "dentist failed to produce output." 1>&2
    cd ..
    rm -rf ${TMPDIR} mamba3 ${OUT_DIR} summ_scaffs
    exit 1
fi

cp gap-closed.fasta ../

cd ..

mv gap-closed.fasta ${OUT_FASTA}
# Keep the uncompressed version for summaries below
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz /staging/lnell/


# This outputs basics about scaffold sizes:
./summ_scaffs ${OUT_FASTA}

export BUSCO_OUT=busco_${OUT_DIR}


# This outputs BUSCO scores:
conda activate busco-env
busco \
    -m genome \
    -l diptera_odb10 \
    -i ${OUT_FASTA} \
    -o ${BUSCO_OUT} \
    --cpu ${THREADS}
conda deactivate


# ~~~~~~~~~~~~~
# For now, we'll just delete the main and BUSCO directories.
# Change this when you finalize the pipeline.
# ~~~~~~~~~~~~~
# # Now the whole directories
# tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
# mv ${OUT_DIR}.tar.gz /staging/lnell/
# tar -czf ${BUSCO_OUT}.tar.gz ${BUSCO_OUT}
# mv ${BUSCO_OUT}.tar.gz /staging/lnell/
rm -r ${TMPDIR} mamba3 ${OUT_DIR} summ_scaffs ${BUSCO_OUT} ${OUT_FASTA}



