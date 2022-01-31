#!/bin/bash

export THREADS=8

. /app/.bashrc
conda activate main-env


# The input FASTA is given by the submit file:
export GENOME=$1.fasta
if [ ! -f /staging/lnell/${GENOME}.gz ]; then
    echo "/staging/lnell/${GENOME}.gz does not exist." 1>&2
    exit 1
fi
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz


# The input file dictates the output name:
export OUT_SUFFIX=prna

if  [[ $GENOME == contigs* ]]
then
    OUT_DIR=scaffolds_${OUT_SUFFIX}
else
    OUT_DIR=${GENOME/.fasta/}_${OUT_SUFFIX}
fi
export OUT_DIR

export OUT_FASTA=${OUT_DIR}.fasta




export RNA1=trimmed_TanyAdult_S1_L002_R1_001.fastq
export RNA2=trimmed_TanyAdult_S1_L002_R2_001.fastq
cp /staging/lnell/${RNA1}.gz ./ && gunzip ${RNA1}.gz
cp /staging/lnell/${RNA2}.gz ./ && gunzip ${RNA2}.gz

# P_RNA (commit 7941e0fcb9ee2a7797fec5e2d359ecde76592139)
tar -xzf P_RNA_scaffolder.tar.gz
rm P_RNA_scaffolder.tar.gz
# Change perl shebangs to '#!/usr/bin/env perl` so that they can use this
# conda environment's bioperl install:
cd P_RNA_scaffolder
for f in *.pl
do
    sed -i "1s/.*/\#\!\/usr\/bin\/env\ perl/" $f
done
cd ..
# Adjust permissions to make files executable:
chmod +x -R P_RNA_scaffolder



export BAM=rna_hisat2__${GENOME/.fasta/}.bam
if [ ! -f /staging/lnell/${BAM} ]; then
    echo "/staging/lnell/${BAM} does not exist." 1>&2
    exit 1
fi
cp /staging/lnell/${BAM} ./
# All alignments are stored in BAM files, but P_RNA only uses SAM files.
# So we have to convert the BAM to SAM:
export SAM=${BAM/.bam/.sam}
samtools view -@ $((THREADS - 2)) -o ${SAM} ${BAM}
rm ${BAM}



bash ./P_RNA_scaffolder/P_RNA_scaffolder.sh \
    -d $(pwd)/P_RNA_scaffolder \
    -i ${SAM} \
    -j ${GENOME} -F ${RNA1} -R ${RNA2} \
    -o ${OUT_DIR} \
    -t ${THREADS}

rm ${RNA1} ${RNA2} ${GENOME} ${SAM}


if [ ! -f ./${OUT_DIR}/P_RNA_scaffold.fasta ]; then
    echo "P_RNA_scaffolder failed to produce output." 1>&2
    rm -r ./${OUT_DIR}
    exit 1
fi


# Move just the contigs file to staging
cp ./${OUT_DIR}/P_RNA_scaffold.fasta ./
mv P_RNA_scaffold.fasta ${OUT_FASTA}
# Keep the uncompressed version for summaries below
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz /staging/lnell/


# This outputs basics about scaffold sizes:
summ_scaffs ${OUT_FASTA}

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
rm -r ${OUT_DIR} ${BUSCO_OUT} ${OUT_FASTA}
