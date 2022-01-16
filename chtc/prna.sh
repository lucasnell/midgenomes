#!/bin/bash

export THREADS=8

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




export RNA1=TanyAdult_S1_L002_R1_001.fastq
export RNA2=TanyAdult_S1_L002_R2_001.fastq
cp /staging/lnell/${RNA1}.gz ./ && gunzip ${RNA1}.gz
cp /staging/lnell/${RNA2}.gz ./ && gunzip ${RNA2}.gz

# P_RNA (commit 7941e0fcb9ee2a7797fec5e2d359ecde76592139)
tar -xzf prna_progs.tar.gz
rm prna_progs.tar.gz
mv ./prna_progs/P_RNA_scaffolder ./
rm -r prna_progs
chmod +x -R P_RNA_scaffolder

# samtools 1.14
tar -xzf samtools-1.14.tar.gz
rm samtools-1.14.tar.gz
export PATH=$(pwd)/samtools-1.14/bin:${PATH}


export BAM=rna_hisat2__${GENOME/.fasta/}.bam
if [ ! -f /staging/lnell/${BAM} ]; then
    echo "/staging/lnell/${BAM} does not exist." 1>&2
    exit 1
fi
cp /staging/lnell/${BAM} ./
# All alignments are stored in BAM files, but P_RNA only uses SAM files.
# So we have to convert the BAM to SAM:
export SAM=${SAM/.bam/.sam}
samtools view -@ $((THREADS - 2)) -o ${SAM} ${BAM}
rm ${BAM}



bash ./P_RNA_scaffolder/P_RNA_scaffolder.sh \
    -d $(pwd)/P_RNA_scaffolder \
    -i ${SAM} \
    -j ${GENOME} -F ${RNA1} -R ${RNA2} \
    -o ${OUT_DIR} \
    -t ${THREADS}

rm ${RNA1} ${RNA2} ${GENOME} ${SAM}


# Move just the contigs file to staging
cp ./${OUT_DIR}/P_RNA_scaffold.fasta ./
mv P_RNA_scaffold.fasta ${OUT_FASTA}
./summ_scaffs ${OUT_FASTA}
gzip ${OUT_FASTA}
mv ${OUT_FASTA}.gz /staging/lnell/


# ~~~~~~~~~~~~~
# For now, we'll just delete the main directory. Change this when you
# finalize the pipeline.
# ~~~~~~~~~~~~~
# # Move the whole P_RNA folder to staging
# tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
# mv ${OUT_DIR}.tar.gz /staging/lnell/

rm -r ./${OUT_DIR} summ_scaffs

