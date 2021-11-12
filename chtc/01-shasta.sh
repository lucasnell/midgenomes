#!/bin/bash
# Assemble Nanopore reads for Tanytarsus gracilentus genome


# Copying the compressed fastq file from staging into the working directory:
export FASTQ=basecalls_guppy-5.0.11.fastq
cp /staging/lnell/${FASTQ}.gz ./
gunzip ${FASTQ}.gz

# Decompressing software (it should be sent from the home directory)
export SHASTA=shasta-Linux-0.7.0
gunzip ${SHASTA}.gz
chmod ugo+x ${SHASTA}


# Assemble Nanopore reads:

export OUTDIR=assembly_shasta
export CONTIGS=contigs_shasta.fasta

echo -e "Started assembly\t" $(date +%F) " " $(date +%H:%M:%S) "\n\n"

./${SHASTA} \
  --input $FASTQ \
  --assemblyDirectory $OUTDIR \
  --threads 32 \
  --config shasta.conf


echo -e "\n\nFinished assembly\t" $(date +%F) " " $(date +%H:%M:%S) "\n\n"

echo -e "Current directory:\n"
ls -lh

echo -e "\n\nDisk used:\n"
du -h -d1

# Compressing full output and sending to staging:
tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
mv ${OUTDIR}.tar.gz /staging/lnell/

# Sending renamed, compressed copy of just contigs to staging:
cp ./${OUTDIR}/Assembly.fasta ./
mv Assembly.fasta ${CONTIGS}
gzip ${CONTIGS}
mv ${CONTIGS}.gz /staging/lnell/

# Removing files used in this job:
rm -r ${OUTDIR} ${SHASTA} $FASTQ shasta.conf

