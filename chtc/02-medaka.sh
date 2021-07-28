#!/bin/bash

# Re-install latest medaka version:
# https://chtc.cs.wisc.edu/conda-installation.shtml


# If version installed is > 1.2.6, then re-run `medaka tools list_models`
# to get updated list of possible models.

# Possible models:
# Available: r103_min_high_g345, r103_min_high_g360, r103_prom_high_g360,
#     r103_prom_snp_g3210, r103_prom_variant_g3210, r10_min_high_g303,
#     r10_min_high_g340, r941_min_fast_g303, r941_min_high_g303,
#     r941_min_high_g330, r941_min_high_g340_rle, r941_min_high_g344,
#     r941_min_high_g351, r941_min_high_g360, r941_prom_fast_g303,
#     r941_prom_high_g303, r941_prom_high_g330, r941_prom_high_g344,
#     r941_prom_high_g360, r941_prom_high_g4011, r941_prom_snp_g303,
#     r941_prom_snp_g322, r941_prom_snp_g360, r941_prom_variant_g303,
#     r941_prom_variant_g322, r941_prom_variant_g360
# Default consensus:  r941_min_high_g360
# Default snp:  r941_prom_snp_g360
# Default variant:  r941_prom_variant_g360


# open up folder with the medaka and other required binaries:
tar -xzf medaka.tar.gz
tar -xzf medaka_bin.tar.gz
rm medaka.tar.gz medaka_bin.tar.gz

# Add these binaries to PATH:
export PATH=$PATH:$(pwd)/medaka_bin/bin


# Set up and activate `medaka` conda environment
ENVNAME=medaka
ENVDIR=$ENVNAME
export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
. $ENVDIR/bin/activate

# Copying the compressed fastq file from staging into the working directory:
cp /staging/lnell/basecalls_guppy-4.5.2.fastq.gz ./
# Also the canu contigs
cp /staging/lnell/contigs_canu-2.1.1.fasta ./


gunzip basecalls_guppy-4.5.2.fastq.gz

NPROC=36
BASECALLS=basecalls_guppy-4.5.2.fastq
DRAFT=contigs_canu-2.1.1.fasta
OUTDIR=medaka_consensus

medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC} \
    -m r941_min_fast_g330






rm -r ./medaka_bin contigs_canu-2.1.1.fasta basecalls_guppy-4.5.2.fastq


