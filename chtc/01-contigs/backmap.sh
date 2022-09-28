#!/bin/bash


#'
#' Use backmap.sh to estimate size of genome by back-mapping ONT
#' reads to assembly.
#'
#' IMPORTANT: use the lucasnell/tany_backmap docker image:
#'

export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

export GENOME=tany_contigs.fasta
export LONGREADS=basecalls_guppy-5.0.11.fastq

export OUT_DIR=tany_backmap

mkdir working
cd working


cp /staging/lnell/assemblies/${GENOME}.gz ./ \
    && gunzip ${GENOME}.gz
cp /staging/lnell/ont/${LONGREADS}.gz ./ \
    && gunzip ${LONGREADS}.gz


. /app/.bashrc
conda activate backmap-env


backmap.pl -o ${OUT_DIR} -pre tany -t ${THREADS} \
    -a ${GENOME} -ont ${LONGREADS}
    -qo "--java-mem-size=16G" \
    1> >(tee -a backmap.stdout)

mv backmap.stdout ./${OUT_DIR}/
mv *.bam.stats ./${OUT_DIR}/
mv tany.ont.sort.bam /staging/lnell/
rm tany.ont.sort.bam.bai ${GENOME} ${LONGREADS}

# Mapped nucleotides:   22.35Gb
# Peak coverage:        235
# Genome size estimate: 95.10Mb

tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/assemblies/

cd ..
rm -r working
