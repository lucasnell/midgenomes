#!/bin/bash
# Assemble Nanopore reads for Tanytarsus gracilentus genome


# Copying the compressed fastq file from staging into the working directory:
export FASTQ=basecalls_guppy-5.0.11.fastq.gz
cp /staging/lnell/$FASTQ ./

# Un-taring software (it should be sent from the home directory)
tar -xJf canu-2.1.1.Linux-amd64.tar.xz
rm canu-2.1.1.Linux-amd64.tar.xz

export PATH=$PATH:$(pwd)/canu-2.1.1/bin


# Assemble Nanopore reads:

export OUTDIR=assembly_canu
export CONTIGS=contigs.fasta
export PREFIX=tany

echo -e "Started assembly\t" $(date +%F) " " $(date +%H:%M:%S) "\n\n"

canu \
  -p ${PREFIX} -d ${OUTDIR} \
  genomeSize=100m \
  useGrid=false \
  -corOutCoverage=100 \
  -maxMemory=64g \
  -maxThreads=32 \
  -nanopore-raw $FASTQ

echo -e "\n\nFinished assembly\t" $(date +%F) " " $(date +%H:%M:%S) "\n\n"

echo -e "Current directory:\n"
ls -lh

echo -e "\n\nDisk used:\n"
du -h -d1

# Compressing full output and sending to staging:
tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
mv ${OUTDIR}.tar.gz /staging/lnell/

# Sending renamed, compressed copy of just contigs to staging:
cp ./${OUTDIR}/${PREFIX}.contigs.fasta ./
mv ${PREFIX}.contigs.fasta ${CONTIGS}
gzip ${CONTIGS}
mv ${CONTIGS}.gz /staging/lnell/

# Removing files used in this job:
rm -r ${OUTDIR} canu-2.1.1 $FASTQ

