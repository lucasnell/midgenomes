#!/bin/bash


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# UPDATE THESE WHEN YOU GET YOUR RNA-seq DATA BACK:
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
export RNA1=tany_RNA_1.fastq
export RNA2=tany_RNA_2.fastq
cp /staging/lnell/${RNA1}.gz ./ && gunzip ${RNA1}
cp /staging/lnell/${RNA2}.gz ./ && gunzip ${RNA2}

tar -xzf P_RNA_progs.tar.gz
rm P_RNA_progs.tar.gz

cd P_RNA_progs

# HISAT2 (version 2.2.1) and P_RNA (commit 7941e0fcb9ee2a7797fec5e2d359ecde76592139)
mv hisat2-2.2.1 P_RNA_scaffolder ../
cd ..
rm -r P_RNA_progs

export PATH=$PATH:$(pwd)/hisat2-2.2.1


export GENOME=haploid_polished_contigs.fasta

cp /staging/lnell/${GENOME}.gz ./
gunzip ${GENOME}.gz


hisat2-build ${GENOME} tany_hisat_idx
hisat2 -x tany_hisat_idx -1 ${RNA1} -2 ${RNA2} -k 3 -p 24 \
    --pen-noncansplice 1000000 -S tany_rna.sam


sh ./P_RNA_scaffolder/P_RNA_scaffolder.sh \
    -d $(pwd)/P_RNA_scaffolder \
    -i tany_rna.sam \
    -j ${GENOME} -F ${RNA1} -R ${RNA2} \
    -o scaffolded_P_RNA \
    -t 24

# Move RNA alignments to folder inside the P_RNA folder.
# It's not having its own folder bc these alignments are part of the
# P_RNA pipeline.
# I'll align the RNAseq reads again for the final assembly for use later on.
mkdir align_hisat2-2.2.1
mv tany_rna.sam ./align_hisat2-2.2.1/
# Index files from `hisat2-build`:
mv tany_hisat_idx* ./align_hisat2-2.2.1/
mv align_hisat2-2.2.1 ./scaffolded_P_RNA/


# Move just the contigs file to staging
cp ./scaffolded_P_RNA/P_RNA_scaffolder.fasta ./
mv P_RNA_scaffolder.fasta scaffolds.fasta
gzip scaffolds.fasta
# Copying bc I want this relatively small file on the submit cluster, too
cp scaffolds.fasta.gz /staging/lnell/


# Move the whole P_RNA folder to staging
tar -czf scaffolded_P_RNA.tar.gz scaffolded_P_RNA
mv scaffolded_P_RNA.tar.gz /staging/lnell/


rm -r ${RNA1} ${RNA2} ${GENOME} ./scaffolded_P_RNA








