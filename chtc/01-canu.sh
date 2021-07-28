#!/bin/bash
# Assemble Nanopore reads for Tanytarsus gracilentus genome



# Copying the compressed fastq file from staging into the working directory:
export FASTQ=basecalls_guppy-4.5.2
cp /staging/lnell/$FASTQ.fastq.gz ./

# Un-taring software (it should be sent from the home directory)
tar -xJf canu-2.1.1.Linux-amd64.tar.xz
rm canu-2.1.1.Linux-amd64.tar.xz

export PATH=$PATH:$(pwd)/canu-2.1.1/bin


# Assemble Nanopore reads:

export out_fn=assembly_canu-2.1.1

echo -e "Started assembly\t" $(date +%F) " " $(date +%H:%M:%S) "\n\n"

canu \
  -p tany -d ${out_fn} \
  genomeSize=100m \
  useGrid=false \
  -corOutCoverage=100 \
  -maxMemory=64g \
  -maxThreads=32 \
  -nanopore-raw ${FASTQ}.fastq.gz

echo -e "\n\nFinished assembly\t" $(date +%F) " " $(date +%H:%M:%S) "\n\n"

echo -e "Current directory:\n"
ls -lh

# Compressing output and sending to staging:
tar -czf ${out_fn}.tar.gz ${out_fn}
mv ${out_fn}.tar.gz /staging/lnell/
# Removing files used in this job:
rm -r ${out_fn} canu-2.1.1 ${FASTQ}.fastq.gz

