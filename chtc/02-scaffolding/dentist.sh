#!/bin/bash


cat << EOF > summ-scaffs.py
#!/usr/bin/env python3
"""
Summarize scaffolds

usage:
    ./summ-scaffs.py <input FASTA>

"""

import sys
import os.path
import gzip
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Summarize scaffolds")
    parser.add_argument("in_fasta", help="input FASTA file name")
    args = parser.parse_args()
    # ---------------
    # Check files:
    # ---------------
    if not os.path.exists(args.in_fasta):
        sys.stderr.write(args.in_fasta + " not found\n")
        sys.exit(1)
    fasta_suffs = (".fasta.gz", ".fasta", ".fa.gz", ".fa")
    if not args.in_fasta.endswith(fasta_suffs):
        sys.stderr.write("Strange suffix to input FASTA file. Exiting.\n")
        sys.exit(1)
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    print("Summarize scaffolds for " + args.in_fasta + "\n")
    if args.in_fasta.endswith(".gz"):
        fasta_file = gzip.open(args.in_fasta,"rt")
    else:
        fasta_file = open(args.in_fasta, "r")
    total_size = 0
    total_N = 0
    sizes = []
    i = -1
    for line in fasta_file:
        if line.startswith(">"):
            sizes.append(0)
            i += 1
        else:
            if len(sizes) == 0:
                sys.stderr.write("FASTA doesn't start with header. Exiting.\n")
                sys.exit(1)
            # we don't want to include newline at the end in our counts:
            llen = len(line) - 1
            sizes[i] += llen
            total_size += llen
            total_N += line.count('n')
            total_N += line.count('N')
    fasta_file.close()
    sizes.sort(reverse = True)
    n50_threshold = round(float(total_size) / 2.0)
    cum_sum = 0;
    i = 0;
    while i < len(sizes):
        cum_sum += sizes[i]
        if cum_sum >= n50_threshold:
            break
        i+=1
    n50 = sizes[i]
    print("size = " + str(total_size))
    print(str(len(sizes)) + " scaffolds")
    print("N50 = " + str(n50))
    print("min = " + str(sizes[-1]))
    print("max = " + str(sizes[0]))
    print("total N = " + str(total_N))
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    sys.exit(0)
EOF

chmod +x summ-scaffs.py



export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')



# The input FASTA is given by the submit file:
export GENOME=$1.fasta

# The join policy, coverage, and seed (for use in `fire-gusu.py`) are
# also input:
export JOIN=$2
export COVERAGE=$3
export SEED=$4


source /staging/lnell/helpers.sh



# Where to send and receive files from/to
export TARGET=/staging/lnell/assemblies



if [ ! -f ${TARGET}/${GENOME}.gz ]; then
    echo -e "\n\nERROR: ${TARGET}/${GENOME}.gz does not exist." 1>&2
    echo -e "Exiting...\n" 1>&2
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

# # If the output FASTA already exists, this job stops with exit code 0
# if [ -f ${TARGET}/${OUT_FASTA}.gz ]; then
#     echo -e "\n\nMESSAGE: ${TARGET}/${OUT_FASTA}.gz already exists."
#     echo -e "Exiting...\n"
#     exit 0
# fi



# installation steps for Mambaforge
# wget "https://github.com/conda-forge/miniforge/releases/download/4.12.0-0/Mambaforge-$(uname)-$(uname -m).sh"
export HOME=$PWD
export PATH
sh Mambaforge-Linux-x86_64.sh -b -p $PWD/mamba3 \
    1> mambaforge.log
check_exit_status "MambaForge" $?
export PATH=$PWD/mamba3/bin:$PATH
rm Mambaforge-Linux-x86_64.sh

export CONDA_PREFIX=$PWD/mamba3

# Install snakemake environment
# (`seqtk` is to convert FASTQ to FASTA and to narrow FASTA files below)
# (`numpy` and `pandas` are for `fire-gusu.py` script below)
mamba create -q -y -c bioconda -c conda-forge -n main-env \
    snakemake=6.13.1 seqtk=1.3 numpy pandas \
    1> mamba_create.log
check_exit_status "main-env" $?

echo -e "\n\n\n" >> mamba_create.log
mamba create -q -y -c bioconda -c conda-forge -n busco-env busco=5.3.1 \
    1>> mamba_create.log
check_exit_status "busco-env" $?

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


mkdir ${OUT_DIR}


cd ${OUT_DIR}

mv ../summ-scaffs.py ./

# Copy :
cp /staging/lnell/dentist_files.tar.gz ./
check_exit_status "cp dentist_files" $?
tar -xzf dentist_files.tar.gz
rm dentist_files.tar.gz
cd dentist_files
mv dentist.yml fire-gusu.py snakemake.yml dentist.v3.0.0.x86_64 ../
cd ..
rm -r dentist_files


# Adjust dentist.yml for coverage and join policy:
sed -i "s/__COVERAGE__/${COVERAGE}/g" dentist.yml
sed -i "s/__JOIN_POLICY__/${JOIN}/g" dentist.yml

# Adjust snakemake.yml for the input FASTA:
sed -i "s/__INPUT_SCAFFOLDS_FILE__/${GENOME}/g" snakemake.yml

# Move and adjust the genome FASTA file:
cp ${TARGET}/${GENOME}.gz ./ && gunzip ${GENOME}.gz

# Make sure the genome FASTA file isn't too wide
mv ${GENOME} ${GENOME/.fasta/_orig.fasta}
seqtk seq -l 80 -A ${GENOME/.fasta/_orig.fasta} > ${GENOME}
check_exit_status "fixed-width fasta (seqtk)" $?
rm ${GENOME/.fasta/_orig.fasta}


# Final FASTA reads file name.
# Make sure it doesn't have any periods, because apparently the
# fasta2DAM > DPsplit step can't handle a file with periods
export READS=basecalls_guppy.fasta

# Get all reads, filter them (dentist takes way too long with all the reads),
# and convert to FASTA with 80-char lines
export ALL_READS=basecalls_guppy-5.0.11.fastq.gz
cp /staging/lnell/${ALL_READS} ./
# info on the reads in the FASTQ file (used for faster filtering):
export SUMMARY=basecalls_guppy-5.0.11--sequencing_summary.txt.gz
cp /staging/lnell/${SUMMARY} ./

# Filtering for average quality of >= 10 and length of >= 10 kb,
# then randomly sample so that desired coverage is reached:
./fire-gusu.py -s ${SUMMARY} -c ${COVERAGE} -g 100 -q 10.0 -l 10000 \
    -o ${READS/.fasta/.fastq} --seed ${SEED} ${ALL_READS}
check_exit_status "fire-gusu" $?
rm ${ALL_READS} ${SUMMARY}
# Convert to 80-char-wide FASTA
seqtk seq -l 80 -A ${READS/.fasta/.fastq} > ${READS}
check_exit_status "fastq > fixed-width fasta (seqtk)" $?
rm ${READS/.fasta/.fastq}


export TMPDIR=$(pwd)/tmp
mkdir ${TMPDIR}


# copy necessary files from dentist dir:
cp -r -t . \
    dentist.v3.0.0.x86_64/Snakefile \
    dentist.v3.0.0.x86_64/envs \
    dentist.v3.0.0.x86_64/scripts


snakemake --configfile=snakemake.yml --use-conda --cores=${THREADS} \
    2> dentist.log
check_exit_status "dentist" $?

echo -e "\n-------------------\nLast bit of dentist stderr:\n"
tail dentist.log
echo -e "\n-------------------\n"

echo -e "\n-------------------\nAny lines in dentist stderr with 'error':\n"
grep -i "error" dentist.log
echo -e "\n-------------------\n"


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
mv ${OUT_FASTA}.gz ${TARGET}/


./summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
check_exit_status "summ-scaffs.py" $?


run_busco ${OUT_FASTA} ${THREADS}
rm -r busco busco_downloads

busco_seq_summary_csv contigs_summary.out busco.out ${OUT_FASTA/.fasta/} | \
    tee ${OUT_FASTA/.fasta/}.csv


cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -rf ${TMPDIR} ${CONDA_PREFIX} ${OUT_DIR}






