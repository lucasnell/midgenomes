#!/bin/bash

#'
#' Polish assembly using Illumina reads with Pilon.
#' Requires argument for which assembly to polish.
#'

export ASSEMBLY=$1

if [[ "${ASSEMBLY}" != *.fasta ]]; then
    echo "ERROR: Assembly files must end in *.fasta. Exiting..." 1>&2
    exit 1
fi


. /app/.bashrc
conda activate main-env
source /staging/lnell/helpers.sh

export THREADS=32

# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
export OUT_NAME=${ASSEMBLY/.fasta/}_pilon
export OUT_DIR=${OUT_NAME}
export OUT_FASTA=${OUT_NAME}.fasta

mkdir ${OUT_DIR}
cd ${OUT_DIR}

# Inputs
cp /staging/lnell/assemblies/${ASSEMBLY}.gz ./ && gunzip ${ASSEMBLY}.gz
check_exit_status "cp genome" $?

export READS1=trimmed_MyKS-19-B_S18_L002_R1_001.fastq
export READS2=trimmed_MyKS-19-B_S18_L002_R2_001.fastq
cp /staging/lnell/dna/trimmed/${READS1}.gz ./ && gunzip ${READS1}.gz
check_exit_status "cp reads 1" $?
cp /staging/lnell/dna/trimmed/${READS2}.gz ./ && gunzip ${READS2}.gz
check_exit_status "cp reads 2" $?



n_rounds=5
input=${ASSEMBLY}

for ((i=1; i<=${n_rounds}; i++)); do
    # index the genome file and do alignment
    bwa index ${input}
    check_exit_status "bwa index ${i}" $?
    bwa mem -t ${THREADS} ${input} ${READS1} ${READS2} | \
        samtools fixmate -m --threads 3  - -| \
        samtools sort -m 2g --threads 5 -| \
        samtools markdup --threads 5 -r - sgs.sort.bam
    check_exit_status "bwa mem ${i}" $?
    # index bam file
    samtools index -@ ${THREADS} sgs.sort.bam
    # polish genome file
    conda activate assembly-env
    pilon --genome ${input} --fix all --frags sgs.sort.bam \
        --threads ${THREADS} --output pilon_${i} | \
        tee round_${i}.pilon
    check_exit_status "pilon ${i}" $?
    conda deactivate
    rm ${input}.* sgs.sort.bam*
    input=pilon_${i}.fasta
done


cp ${input} ${OUT_FASTA}
check_exit_status "cp output" $?

rm ${READS1} ${READS2} ${ASSEMBLY}

summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
check_exit_status "summ-scaffs.py" $?

run_busco ${OUT_FASTA} ${THREADS}
rm -r busco_downloads

busco_seq_summary_csv contigs_summary.out busco.out ${OUT_NAME} | \
    tee ${OUT_NAME}.csv

# Keep the uncompressed version for output in main directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz ${TARGET}/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}


