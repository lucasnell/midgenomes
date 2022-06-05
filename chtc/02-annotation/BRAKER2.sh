#!/bin/bash

# Use BRAKER2 to annotate Tanytarsus gracilentus genome

#' THIS IF FROM EDTA OUTPUT:
#' ```
#' Low-threshold TE masking for MAKER gene annotation (masked: 3.68%): tany_scaffolds.fasta.mod.MAKER.masked
#' ```
#'
#' So use tany_scaffolds.fasta.mod.MAKER.masked for this annotation!
#' Maybe run this on it first, so that it isn't too wide:
#'
#' ```
#' seqtk seq -l 80 -A tany_scaffolds.fasta.mod.MAKER.masked > \
#'     tany_scaffolds_MAKER_masked.fasta
#' ```
#'





export THREADS=16

. /app/.bashrc
conda activate main-env

export OUT_DIR=tany_annotation

mkdir ${OUT_DIR}
cd ${OUT_DIR}

export GENOME=
export RNA_READS_ADULT1=trimmed_TanyAdult_S1_L002_R1_001.fastq
export RNA_READS_ADULT2=trimmed_TanyAdult_S1_L002_R2_001.fastq
export RNA_READS_ADULT_TAR=trimmed_TanyAdult_S1.tar
export RNA_READS_JUVEN1=trimmed_TanyJuven_S2_L002_R1_001.fastq
export RNA_READS_JUVEN2=trimmed_TanyJuven_S2_L002_R2_001.fastq
export RNA_READS_JUVEN_TAR=trimmed_TanyJuven_S2.tar



cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz

cp /staging/lnell/rna/${RNA_READS_ADULT_TAR} ./ \
    && tar -xf ${RNA_READS_ADULT_TAR} \
    && rm ${RNA_READS_ADULT_TAR}
cp /staging/lnell/rna/${RNA_READS_JUVEN_TAR} ./ \
    && tar -xf ${RNA_READS_JUVEN_TAR} \
    && rm ${RNA_READS_JUVEN_TAR}

gunzip ${RNA_READS_ADULT1}.gz
gunzip ${RNA_READS_ADULT2}.gz
gunzip ${RNA_READS_JUVEN1}.gz
gunzip ${RNA_READS_JUVEN2}.gz

export BAM_ADULT=tany_rna_adults.bam
export BAM_JUVEN=tany_rna_juveniles.bam

# BRAKER2 docs recommend simplifying sequence names:
# I'll do this by removing everything after the first comma or the first
# space in the seq names:
sed -i -r 's/\,.+//' ${GENOME}
sed -i -r 's/\ .+//' ${GENOME}



hisat2-build -q ${GENOME} tany_hisat_idx

# (Pipe to samtools is to convert it to an aligned BAM file)
hisat2 -x tany_hisat_idx \
    -1 ${RNA_READS_ADULT1} -2 ${RNA_READS_ADULT2} \
    -k 3 -p ${THREADS} | \
    samtools sort -O bam - > ${BAM_ADULT}
samtools index ${BAM_ADULT} ${BAM_ADULT}.bai

hisat2 -x tany_hisat_idx \
    -1 ${RNA_READS_JUVEN1} -2 ${RNA_READS_JUVEN2} \
    -k 3 -p ${THREADS} | \
    samtools sort -O bam - > ${BAM_JUVEN}
samtools index ${BAM_JUVEN} ${BAM_JUVEN}.bai

# reads and hisat2 indices no longer needed
rm tany_hisat_idx* ${RNA_READS_ADULT1} ${RNA_READS_ADULT2} \
    ${RNA_READS_JUVEN1} ${RNA_READS_JUVEN2}

conda activate annotate-env

# Using softmasking bc DENTIST soft masks repeat regions

braker.pl --species=tanytarsus_gracilentus --genome=${GENOME} \
    --softmasking \
    --cores=${THREADS} \
    --bam=${BAM_ADULT},${BAM_JUVEN}


# Saving output:
cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/
rm -r ${OUT_DIR}

