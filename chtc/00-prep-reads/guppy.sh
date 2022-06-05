#!/bin/bash

# Basecalling for Nanopore reads

# Set guppy version:
export GUPPY_V=5.0.11

# Copying the raw signal files and guppy from staging into the working directory:
cp /staging/lnell/ives_fast5.tar.bz2 ./
mkdir fast5
tar -xf ives_fast5.tar.bz2 -C ./fast5/
rm ives_fast5.tar.bz2

# Copying and un-taring software
cp /staging/lnell/guppy.tar.gz ./
tar -xzf guppy.tar.gz
rm guppy.tar.gz


# Add these binaries to PATH:
export PATH=$PATH:$(pwd)/guppy/bin
# For shared libraries:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/guppy/lib


export OUTDIR=basecalls_guppy-${GUPPY_V}

mkdir ${OUTDIR}

guppy_basecaller \
  --input_path ./fast5 \
  --save_path ./${OUTDIR} \
  --records_per_fastq 0 \
  --compress_fastq \
  --trim_barcodes \
  --disable_pings \
  --device cuda:0,1	\
  --flowcell FLO-MIN106 --kit SQK-LSK109

# Above is equivalent to config file `dna_r9.4.1_450bps_hac.cfg`
# version 2021-05-17_dna_r9.4.1_minion_384_d37a2ab9


# FASTQ file should already be compressed, so just tar-ing the folder, then
# sending to staging:
tar -cf ${OUTDIR}.tar ${OUTDIR}
mv ${OUTDIR}.tar /staging/lnell/

# Making a single FASTQ file of just passing reads, then sending it to staging.
# This file will be used for downstream processes.
cd ${OUTDIR}/pass
cat *.fastq.gz > ${OUTDIR}.fastq.gz
mv ${OUTDIR}.fastq.gz /staging/lnell/
cd ../../

# Removing files used in this job:
rm -r ./${OUTDIR} ./fast5 ./guppy


