#!/bin/bash

# have job exit if any command returns with non-zero exit status (aka failure)
set -e

export OUTDIR=scaffolds_longstitch_arks
export OUTFASTA=${OUTDIR}.fasta
export THREADS=8


mkdir ${OUTDIR}
cd ${OUTDIR}


# export GENOME=haploid_purge_dups
export GENOME=scaffolds_p_rna
cp /staging/lnell/${GENOME}.fasta.gz ./ && gunzip ${GENOME}.fasta.gz
# LongStitch requires *.fa ending
mv ${GENOME}.fasta ${GENOME}.fa

export FASTQ=basecalls_guppy-5.0.11
cp /staging/lnell/${FASTQ}.fastq.gz ./
# LongStitch requires *.fq.gz ending
mv ${FASTQ}.fastq.gz ${FASTQ}.fq.gz


# longstitch run \
longstitch tigmint-ntLink-arks \
    draft=${GENOME} reads=${FASTQ} t=${THREADS} G=100000000

rm -rf ${GENOME}.fa ${FASTQ}.fq.gz

# Gets filename of final scaffolds FASTA file
export LS_FASTA=$(basename -- $(readlink -f *-arks.longstitch-scaffolds.fa))
cp $LS_FASTA ../

cd ..

# Moving the final FASTA to staging
mv $LS_FASTA ${OUTFASTA}
gzip ${OUTFASTA}
mv ${OUTFASTA}.gz /staging/lnell/

# Now the whole directory
tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
mv ${OUTDIR}.tar.gz /staging/lnell/
rm -r ${OUTDIR}








