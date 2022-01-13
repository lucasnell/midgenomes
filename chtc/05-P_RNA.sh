#!/bin/bash

# # have job exit if any command returns with non-zero exit status (aka failure)
# set -e

export THREADS=32

# The name of this file will dictate the BAM file used, too:
export GENOME=haploid_purge_dups.fasta
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz


export RNA1=TanyAdult_S1_L002_R1_001.fastq
export RNA2=TanyAdult_S1_L002_R2_001.fastq
cp /staging/lnell/${RNA1}.gz ./ && gunzip ${RNA1}.gz
cp /staging/lnell/${RNA2}.gz ./ && gunzip ${RNA2}.gz

tar -xzf p_rna_progs.tar.gz
rm p_rna_progs.tar.gz

# P_RNA (commit 7941e0fcb9ee2a7797fec5e2d359ecde76592139)
mv ./p_rna_progs/P_RNA_scaffolder ./
rm -r p_rna_progs

# All alignments are stored in BAM files, so we have to convert to SAM:
export BAM=rna_hisat2__${GENOME/.fasta/}.bam
if [ ! -f /staging/lnell/${BAM} ]; then
    echo "/staging/lnell/${BAM} does not exist." 1>&2
    exit 1
fi
cp /staging/lnell/${BAM} ./
export SAM=${SAM/.bam/.sam}
samtools view -o ${SAM} ${BAM}
rm ${BAM}




export OUTDIR=scaffolds_p_rna
export OUTFASTA=${OUTDIR}.fasta


sh ./P_RNA_scaffolder/P_RNA_scaffolder.sh \
    -d $(pwd)/P_RNA_scaffolder \
    -i ${SAM} \
    -j ${GENOME} -F ${RNA1} -R ${RNA2} \
    -o ${OUTDIR} \
    -t ${THREADS}

rm ${RNA1} ${RNA2} ${GENOME}

# # Move RNA alignments to their own folder for potential later use.
# export RNA_ALIGN_OUT=adult_rna_align_hisat2
# mkdir ${RNA_ALIGN_OUT}
# mv ${SAM} ./${RNA_ALIGN_OUT}/
# # Index files from `hisat2-build`:
# mv tany_hisat_idx* ./${RNA_ALIGN_OUT}/
# tar -czf ${RNA_ALIGN_OUT}.tar.gz ${RNA_ALIGN_OUT}
# mv ${RNA_ALIGN_OUT}.tar.gz /staging/lnell/
# rm -r ${RNA_ALIGN_OUT}
rm ${SAM}


# Move just the contigs file to staging
cp ./${OUTDIR}/P_RNA_scaffold.fasta ./
mv P_RNA_scaffold.fasta ${OUTFASTA}
gzip ${OUTFASTA}
mv ${OUTFASTA}.gz /staging/lnell/


# Move the whole P_RNA folder to staging
tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
mv ${OUTDIR}.tar.gz /staging/lnell/

rm -r ./${OUTDIR}


# ============================================================================
# ============================================================================

# Running on laptop

# ============================================================================
# ============================================================================

docker pull bioperl/bioperl:latest

# Run interactively:
# (drive `Lucas256` has RNA data, `~/_data/p_rna` has everything else)
docker run -it \
    -v ~/_data/p_rna:/p_rna \
    -v /Volumes/Lucas256/RNA:/RNA \
    -u $(id -u) \
    bioperl/bioperl:latest \
    /bin/bash
    # --privileged \


# It's assumed you've already gunzipped the RNA and GENOME files

export THREADS=6

export RNA1=/RNA/TanyAdult_S1_L002_R1_001.fastq
export RNA2=/RNA/TanyAdult_S1_L002_R2_001.fastq
export GENOME=/p_rna/scaffolds_longstitch.fasta
export SAM=/p_rna/tany_rna.sam
export OUTDIR=/p_rna/scaffolds_p_rna

# chmod u+w -R p_rna

cd p_rna

tar -xzf p_rna_progs.tar.gz
mv ./p_rna_progs/P_RNA_scaffolder ./
rm -r p_rna_progs


chmod u+x -R P_RNA_scaffolder
chmod u+w -R P_RNA_scaffolder
cd P_RNA_scaffolder
# Then I went through and did `chown root: X` on all files where root isn't owner

cd ..

t0=`date`
bash /p_rna/P_RNA_scaffolder/P_RNA_scaffolder.sh \
    -d /p_rna/P_RNA_scaffolder \
    -i ${SAM} \
    -j ${GENOME} -F ${RNA1} -R ${RNA2} \
    -e 15000 \
    -o ${OUTDIR} \
    -t ${THREADS}
t1=`date`
echo -e "Started: " $t0 "\nFinished: " $t1



# perl /p_rna/P_RNA_scaffolder/UNIQUE_sam_intron.pl ${SAM} \
#     ${OUTDIR}/F.sam ${OUTDIR}/R.sam ${OUTDIR}/intron.txt



