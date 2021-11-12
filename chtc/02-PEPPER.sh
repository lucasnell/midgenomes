#!/bin/bash
# Polish Tanytarsus gracilentus genome


# Copying the compressed fastq file from staging into the working directory:
export FASTQ=basecalls_guppy-5.0.11.fastq
cp /staging/lnell/${FASTQ}.gz ./
gunzip ${FASTQ}.gz

# # Decompressing software (it should be sent from the home directory)
# export SHASTA=shasta-Linux-0.7.0
# gunzip ${SHASTA}.gz
# chmod ugo+x ${SHASTA}


# Assemble Nanopore reads:

export OUTDIR=assembly_shasta
export CONTIGS=contigs_shasta.fasta





# LEFT OFF --> PASTED THESE FROM THEIR WEBSITE
#   https://github.com/kishwarshafin/pepper/blob/r0.4/docs/pipeline_docker/ONT_polishing.md
# Still need code for aligning to genome to make BAM file using minimap2
#   https://github.com/lh3/minimap2
# Also, does Docker have access to staging folder?



run_pepper_margin_deepvariant polish_assembly \
    -b "${BAM}" \
    -f "${CONTIGS}" \
    -o "${OUTDIR}" \
    -t ${THREADS} \
    -s ${SAMPLE_NAME} \
    --ont



bcftools consensus \
    -f "${INPUT_DIR}/${ASM}" \
    -H 2 \
    -s "${SAMPLE_NAME}" \
    -o "${OUTPUT_DIR}/${POLISHED_ASM_HAP1}" \
    "${OUTPUT_DIR}/${HAP1_VCF}"

bcftools consensus \
    -f "${INPUT_DIR}/${ASM}" \
    -H 2 \
    -s "${SAMPLE_NAME}" \
    -o "${OUTPUT_DIR}/${POLISHED_ASM_HAP2}" \
    "${OUTPUT_DIR}/${HAP2_VCF}"






# Compressing full output and sending to staging:
tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
mv ${OUTDIR}.tar.gz /staging/lnell/

# Sending renamed, compressed copy of just contigs to staging:
cp ./${OUTDIR}/Assembly.fasta ./
mv Assembly.fasta ${CONTIGS}
gzip ${CONTIGS}
mv ${CONTIGS}.gz /staging/lnell/

# Removing files used in this job:
rm -r ${OUTDIR} ${SHASTA} $FASTQ shasta.config

