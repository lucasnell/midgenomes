#!/bin/bash

#' This is to get a rough estimate of nucleotide diversity and
#' heterozygosity from Nanopore sequencing to use as priors in SNAPE-pooled.

check_exit_status () {
  if [ ! "$2" -eq "0" ]; then
    echo "Step $1 failed with exit status $2"
    cd ..
    tar -czf ERROR_${OUT_DIR}.tar.gz ${OUT_DIR}
    mv ERROR_${OUT_DIR}.tar.gz /staging/lnell/dna/ont-pi
    rm -r ${OUT_DIR}
    exit $2
  fi
  echo "Checked step $1"
}

# export THREADS=24
export THREADS=8

. /app/.bashrc
conda activate main-env


# Outputs
export OUT_DIR=ont-pi
export BAM=ont_align.bam
export BCF=ont_variants.bcf
export PI_FILE=ont_variants.windowed.pi
# Where to send everything:
export TARGET=/staging/lnell/dna/ont-pi

mkdir ${OUT_DIR}
cd ${OUT_DIR}

# Inputs
export GENOME=tany_scaffolds.fasta
export ONT_FASTQ=basecalls_guppy-5.0.11.fastq

if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist! " 1>&2
    exit 111
fi
if [ ! -f /staging/lnell/${ONT_FASTQ}.gz ]; then
    echo "/staging/lnell/${ONT_FASTQ}.gz does not exist! " 1>&2
    exit 222
fi

cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz
cp /staging/lnell/${ONT_FASTQ}.gz ./ && gunzip ${ONT_FASTQ}.gz

samtools faidx --length 80 ${GENOME}


# ------------------------------
# Nanopore reads
# ------------------------------

# minimap2 -a -t $((THREADS - 4)) -K 1G -2 \
#     ${GENOME} ${ONT_FASTQ} | \
#     samtools sort -@ 2 -o ${BAM/.bam/_unfiltered.bam} -
#
# check_exit_status "minimap2" $?
#
# rm ${ONT_FASTQ}
#
# bamtools stats -in ${BAM/.bam/_unfiltered.bam}
# check_exit_status "bamtools stats (unfiltered)" $?
#
# samtools view -bh -q 20 -F 0x100 --threads 2 \
#     ${BAM/.bam/_unfiltered.bam} > \
#     ${BAM}
# check_exit_status "samtools view (filter)" $?
#
# rm ${BAM/.bam/_unfiltered.bam}
#
# bamtools stats -in ${BAM}
# check_exit_status "bamtools stats (filtered)" $?


# Takes about an hour
bcftools mpileup -Ou --config ont -f ${GENOME} ${BAM} | \
    bcftools call -P 0.01 -mv --threads 2 -Ob -o ${BCF}
check_exit_status "bamtools mpileup | call" $?

rm ${GENOME}* ${BAM}*


#' Calculate nucleotide diversity in windows
#' This step creates the file ${PI_FILE}

vcftools --bcf ${BCF} --out ${PI_FILE/.windowed.pi/} \
    --remove-indels \
    --max-alleles 2 \
    --window-pi 1000 \
    --window-pi-step 1000
check_exit_status "vcftools" $?

rm *.log

#'
#' This returns the mean nucleotide diversity (pi) among all sites.
#' Watterson's estimator of theta is the number of segregating sites divided
#' by the (n-1)th harmonic number (i.e., sum(1/i) for i in 1:(n-1)),
#' where n is the number of haploid samples.
#' Because this is one diploid individual, Watterson's estimator is the
#' same as pi here.
#'
#' The output from this was "0.008237687486796706".
#' I'll use 0.008 for my prior for SNAPE-pool.
#'
echo $(python3 << EOF
import pandas as pd
win_pi = pd.read_csv("${PI_FILE}", sep="\t")
print(win_pi["PI"].mean())
EOF
)


gzip ${PI_FILE}


# ------------------------------
# Handle output
# ------------------------------

mv ${PI_FILE}.gz ${BCF} ${TARGET}/

cd ..
rm -r ${OUT_DIR}

