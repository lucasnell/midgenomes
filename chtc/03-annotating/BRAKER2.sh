#!/bin/bash



export THREADS=16


export GENOME=
export RNA_READS_ADULT1=trimmed_TanyAdult_S1_L002_R1_001.fastq
export RNA_READS_ADULT2=trimmed_TanyAdult_S1_L002_R2_001.fastq
export RNA_READS_JUVEN1=trimmed_TanyJuven_S2_L002_R1_001.fastq
export RNA_READS_JUVEN2=trimmed_TanyJuven_S2_L002_R2_001.fastq
cp /staging/lnell/${RNA_READS_ADULT1}.gz ./ && gunzip ${RNA_READS_ADULT1}.gz
cp /staging/lnell/${RNA_READS_ADULT2}.gz ./ && gunzip ${RNA_READS_ADULT2}.gz
cp /staging/lnell/${RNA_READS_JUVEN1}.gz ./ && gunzip ${RNA_READS_JUVEN1}.gz
cp /staging/lnell/${RNA_READS_JUVEN2}.gz ./ && gunzip ${RNA_READS_JUVEN2}.gz

export BAM_ADULT=tany_rna_adults.bam
export BAM_JUVEN=tany_rna_juveniles.bam



hisat2-build -q ${GENOME} tany_hisat_idx
hisat2 -x tany_hisat_idx \
    -1 ${RNA_READS_ADULT1} -2 ${RNA_READS_ADULT2} \
    -p ${THREADS} -S ${BAM_ADULT/.bam/.sam}
hisat2 -x tany_hisat_idx \
    -1 ${RNA_READS_JUVEN1} -2 ${RNA_READS_JUVEN2} \
    -p ${THREADS} -S ${BAM_JUVEN/.bam/.sam}

# hisat2 indices no longer needed
rm tany_hisat_idx*

# Convert to aligned BAM files with index files for use in BRAKER2:
samtools view -u -h -b ${BAM_ADULT/.bam/.sam} | \
    samtools sort -@ $((THREADS - 2)) -o ${BAM_ADULT} -
samtools index ${BAM_ADULT} ${BAM_ADULT}.bai
samtools view -u -h -b ${BAM_JUVEN/.bam/.sam} | \
    samtools sort -@ $((THREADS - 2)) -o ${BAM_JUVEN} -
samtools index ${BAM_JUVEN} ${BAM_JUVEN}.bai
rm ${BAM_ADULT/.bam/.sam} ${BAM_JUVEN/.bam/.sam}


# Using softmasking bc DENTIST soft masks repeat regions


braker.pl --species=tanytarsus_gracilentus --genome=${GENOME} \
    --softmasking \
    --cores=${THREADS} \
    --bam=${BAM_ADULT},${BAM_JUVEN}




