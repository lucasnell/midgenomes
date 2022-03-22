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
export VCF=ont_variants.vcf.gz
export VCFTOOLS_PREFIX=ont_variants_stats
# Where to send everything:
export TARGET=/staging/lnell/dna/ont-pi

mkdir ${OUT_DIR}
cd ${OUT_DIR}

# Inputs
export GENOME=tany_scaffolds.fasta
export ONT_FASTQ=basecalls_guppy-5.0.11.fastq

if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist!" 1>&2
    exit 111
fi
if [ ! -f /staging/lnell/${ONT_FASTQ}.gz ]; then
    echo "/staging/lnell/${ONT_FASTQ}.gz does not exist!" 1>&2
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
# bamtools stats -in ${BAM/.bam/_unfiltered.bam} | \
#     tee ${BAM/.bam/_unfiltered.stats}
# check_exit_status "bamtools stats (unfiltered)" $?

cp ${TARGET}/${BAM/.bam/_unfiltered.bam} ./

samtools view -bh -q 20 -F 0x100 --threads 2 \
    ${BAM/.bam/_unfiltered.bam} > \
    ${BAM}
check_exit_status "samtools view (filter)" $?

rm ${BAM/.bam/_unfiltered.bam}

bamtools stats -in ${BAM} | \
    tee ${BAM/.bam/.stats}
check_exit_status "bamtools stats (filtered)" $?


bcftools mpileup -f ${GENOME} --threads 2 -Ou ${BAM} | \
    bcftools call -vm --threads 2 -Oz > \
    ${VCF}
check_exit_status "bamtools mpileup | call" $?

rm ${GENOME}*


#' Calculate nucleotide diversity and heterozygosity.
#' This step creates the following files:
#' - ${VCFTOOLS_PREFIX}.windowed.pi
#' - ${VCFTOOLS_PREFIX}.het
vcftools --gzvcf ${VCF} --out ${VCFTOOLS_PREFIX} \
    --remove-indels \
    --window-pi 500 \
    --window-pi-step 250 \
    --het
check_exit_status "vcftools" $?





# ------------------------------
# Handle output
# ------------------------------

mv ${VCFTOOLS_PREFIX}* ${TARGET}/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

