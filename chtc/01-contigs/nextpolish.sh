#!/bin/bash

#'
#' Polish assembly using Illumina reads with NextPolish.
#' Requires argument for which assembly to polish.
#'


export ASSEMBLY=$1

#' Second argument (if provided) is # rounds.
#' I wouldn't recommend more than 3: More rounds isn't necessarily better!
if [ -z ${2+x} ]; then
    N_ROUNDS=3
else
    N_ROUNDS=$2
    if ! [[ $N_ROUNDS =~ ^[0-9]+$ ]]; then
        echo "ERROR: Second arg is not an integer!" 1>&2
        exit 1
    fi
fi
export N_ROUNDS

if [[ "${ASSEMBLY}" != *.fasta ]]; then
    echo "ERROR: Assembly files must end in *.fasta. Exiting..." 1>&2
    exit 1
fi


. /app/.bashrc
conda activate main-env
source /staging/lnell/helpers.sh


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')


# Where to send files:
export TARGET=/staging/lnell/assemblies

# Output names:
export OUT_NAME=${ASSEMBLY/.fasta/}_nextpolish
export OUT_DIR=${OUT_NAME}
export OUT_FASTA=${OUT_NAME}.fasta

mkdir ${OUT_DIR}
cd ${OUT_DIR}

# Inputs
cp /staging/lnell/assemblies/${ASSEMBLY}.gz ./ && gunzip ${ASSEMBLY}.gz
check_exit_status "cp genome" $?

export READS1=trimmed_MyKS-19-B_S18_L002_R1_001.fastq.gz
export READS2=trimmed_MyKS-19-B_S18_L002_R2_001.fastq.gz
cp /staging/lnell/dna/trimmed/${READS1} ./
check_exit_status "cp reads 1" $?
cp /staging/lnell/dna/trimmed/${READS2} ./
check_exit_status "cp reads 2" $?



# Filter out reads with any uncalled bases, to prevent N in final assembly
# Input reads have already been trimmed for adapters, etc.
export TRIM_READS1=in_reads1.fastq
export TRIM_READS2=in_reads2.fastq

fastp --in1 ${READS1} --in2 ${READS2} \
    --out1 ${TRIM_READS1} --out2 ${TRIM_READS2} \
    --thread ${THREADS} \
    --disable_length_filtering \
    --disable_adapter_trimming \
    --disable_trim_poly_g \
    --dont_eval_duplication \
    --qualified_quality_phred 0 \
    --unqualified_percent_limit 100 \
    --n_base_limit 0

rm ${READS1} ${READS2} fastp*

export READS1=${TRIM_READS1}
export READS2=${TRIM_READS2}
unset TRIM_READS1 TRIM_READS2




input=${ASSEMBLY}

for ((i=1; i<=${N_ROUNDS}; i++)); do
# step 1:
    # index the genome file and do alignment
    bwa index ${input}
    bwa mem -t ${THREADS} ${input} ${READS1} ${READS2} | \
        samtools view --threads 3 -F 0x4 -b -| \
        samtools fixmate -m --threads 3  - -| \
        samtools sort -m 2g --threads 5 -| \
        samtools markdup --threads 5 -r - sgs.sort.bam
    # index bam and genome files
    samtools index -@ ${THREADS} sgs.sort.bam
    samtools faidx ${input}
    # polish genome file
    python /opt/NextPolish/lib/nextpolish1.py -g ${input} \
        --task 1 -p ${THREADS} -s sgs.sort.bam \
        > genome.polishtemp.fa
    rm ${input}.* sgs.sort.bam*
    input=genome.polishtemp.fa
# step2:
    # index genome file and do alignment
    bwa index ${input}
    bwa mem -t ${THREADS} ${input} ${READS1} ${READS2} | \
        samtools view --threads 3 -F 0x4 -b -| \
        samtools fixmate -m --threads 3  - -| \
        samtools sort -m 2g --threads 5 -| \
        samtools markdup --threads 5 -r - sgs.sort.bam
    # index bam and genome files
    samtools index -@ ${THREADS} sgs.sort.bam
    samtools faidx ${input}
    # polish genome file
    python /opt/NextPolish/lib/nextpolish1.py -g ${input} \
        --task 2 -p ${THREADS} -s sgs.sort.bam \
        > genome.nextpolish.fa
    rm ${input}.* sgs.sort.bam*
    input=genome.nextpolish.fa
    if (( i < N_ROUNDS )); then
        cp ${input} genome.nextpolish_${i}.fa
    fi
done


# Finally polished genome file: genome.nextpolish.fa

cp genome.nextpolish.fa ${OUT_FASTA}
check_exit_status "cp output" $?

rm ${READS1} ${READS2} ${ASSEMBLY}

summ-scaffs.py ${OUT_FASTA} | tee contigs_summary.out
check_exit_status "summ-scaffs.py" $?

run_busco ${OUT_FASTA} ${THREADS}
rm -r busco busco_downloads

busco_seq_summary_csv contigs_summary.out busco.out ${OUT_NAME} | \
    tee ${OUT_NAME}.csv

# Keep the uncompressed version for output in main directory
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz ${TARGET}/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

