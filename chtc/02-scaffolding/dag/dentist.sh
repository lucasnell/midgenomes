#!/bin/bash

export THREADS=24



# The input FASTA is given by the submit file:
export GENOME=$1.fasta
if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo -e "\n\nERROR: /staging/lnell/${GENOME}.gz does not exist." 1>&2
    echo -e "Exiting...\n" 1>&2
    exit 1
fi

# The input file dictates the output name and join policy:
export OUT_SUFFIX=dentist
if  [[ $GENOME == contigs* ]]
then
    OUT_DIR=scaffolds_${OUT_SUFFIX}
    JOIN=contigs
else
    OUT_DIR=${GENOME/.fasta/}_${OUT_SUFFIX}
    # Only for testing round3-2.dag
    # JOIN=scaffolds
    JOIN=scaffoldGaps
fi
export OUT_DIR
export JOIN

export OUT_FASTA=${OUT_DIR}.fasta

# If the output FASTA already exists, this job stops with exit code 0
if [ -f /staging/lnell/${OUT_FASTA}.gz ]; then
    echo -e "\n\nMESSAGE: /staging/lnell/${OUT_FASTA}.gz already exists."
    echo -e "Exiting...\n"
    exit 0
fi



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
# (`numpy` and `pandas` are for `fire-gusu.py` script below)
mamba create -q -y -c bioconda -c conda-forge -n main-env \
    snakemake=6.13.1 seqtk=1.3 numpy pandas

conda init
. ~/.bashrc

conda activate main-env

# # example ------------------------------
# wget https://github.com/a-ludi/dentist/releases/download/v3.0.0/dentist-example.tar.gz
# tar -xzf dentist-example.tar.gz
# cd dentist-example
# snakemake --configfile=snakemake.yml --use-conda --cores=1
# md5sum -c checksum.md5
# # --------------------------------------



# I'll use the first input to produce the RNG seed for reproducibility
# with `fire-gusu.py` script.
# This means that I'll use the same reads each time I process the same file.
# The only time I'll be processing the same file will be for testing, so this
# behavior is desirable.
export SEED=$(python << EOF
import hashlib
print(int(hashlib.sha512("${1}".encode('utf-8')).hexdigest(), base = 16))
EOF
)



mkdir ${OUT_DIR}

# Move tar file from submit node to OUT_DIR:
mv dentist_files.tar.gz ./${OUT_DIR}/

cd ${OUT_DIR}

# Expand tar file from submit server
# (this includes `dentist.v3.0.0.x86_64` folder):
tar -xzf dentist_files.tar.gz
rm dentist_files.tar.gz

# Adjust dentist.yml for join policy:
sed -i "s/__JOIN_POLICY__/${JOIN}/g" dentist.yml

# Adjust snakemake.yml for the input FASTA:
sed -i "s/__INPUT_SCAFFOLDS_FILE__/${GENOME}/g" snakemake.yml

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
# info on the reads in the FASTQ file (used for faster filtering):
export SUMMARY=basecalls_guppy-5.0.11--sequencing_summary.txt.gz
cp /staging/lnell/${SUMMARY} ./

# Filtering for average quality of >= 10 and length of >= 10 kb, then randomly
# choosing reads to get the desired coverage:
./fire-gusu.py -s ${SUMMARY} -c ${COVERAGE} -g 100 -q 10.0 -l 10000 \
    -o ${READS/.fasta/.fastq} --seed ${SEED} ${ALL_READS}
rm ${ALL_READS} ${SUMMARY}
# Convert to 80-char-wide FASTA
seqtk seq -l 80 -A ${READS/.fasta/.fastq} > ${READS}
rm ${READS/.fasta/.fastq}


export TMPDIR=$(pwd)/tmp
mkdir ${TMPDIR}


# copy necessary files from dentist dir:
cp -r -t . \
    dentist.v3.0.0.x86_64/Snakefile \
    dentist.v3.0.0.x86_64/envs \
    dentist.v3.0.0.x86_64/scripts


snakemake --configfile=snakemake.yml --use-conda --cores=${THREADS}


rm -rf dentist.v3.0.0.x86_64 ${READS} ${GENOME}

if [ ! -f gap-closed.fasta ]; then
    echo -e "\n\nERROR: dentist failed to produce output." 1>&2
    echo -e "Exiting...\n" 1>&2
    cd ..
    rm -rf ${TMPDIR} ${CONDA_PREFIX} ${OUT_DIR}
    exit 1
fi




cp gap-closed.fasta ${OUT_FASTA}
# Keep the uncompressed version for output directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz /staging/lnell/


# We're removing instead of saving the output directory bc it takes up
# too much space. We'll look at them later once I finalize the pipeline.
cd ..

# # To save whole directory:
# tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
# mv ${OUT_DIR}.tar.gz /staging/lnell/

rm -rf ${TMPDIR} ${CONDA_PREFIX} ${OUT_DIR}






