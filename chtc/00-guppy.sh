#!/bin/bash

# Basecalling for Nanopore reads

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


export out_fn=basecalls_guppy-4.5.2

mkdir ${out_fn}

guppy_basecaller \
  --input_path ./fast5 \
  --save_path ./${out_fn} \
  --records_per_fastq 0 \
  --compress_fastq \
  --trim_barcodes \
  --disable_pings \
  --chunks_per_runner 1024 \
  --gpu_runners_per_device 4 \
  --device cuda:0,1	\
  --flowcell FLO-MIN106 --kit SQK-LSK109


# FASTQ file should already be compressed, so just tar-ing the folder , then
# sending to staging:
tar -cf ${out_fn}.tar ${out_fn}
mv ${out_fn}.tar /staging/lnell/
# Removing files used in this job:
rm -r ./${out_fn} ./fast5 ./guppy

