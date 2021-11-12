#!/bin/bash
# Polish Tanytarsus gracilentus genome

# Note that we have to use `r0.4` version of PEPPER:
# https://github.com/kishwarshafin/pepper/tree/r0.4


export OUTDIR=polish_PEPPER
export ALIGNDIR=ont_align_minimap2

# Copying the compressed fastq file from staging into the working directory:
export FASTQ=basecalls_guppy-5.0.11.fastq
cp /staging/lnell/${FASTQ}.gz ./
gunzip ${FASTQ}.gz


# Assembled Nanopore reads from SHASTA:
export CONTIGS=contigs_shasta.fasta
cp /staging/lnell/${CONTIGS}.gz
gunzip ${CONTIGS}.gz


export ALIGNMENT=ont_alignment.bam


# Align nanopore reads to assembly using minimap2:
tar -jxvf minimap2-2.22_x64-linux.tar.bz2
rm minimap2-2.22_x64-linux.tar.bz2

./minimap2-2.22_x64-linux/minimap2 -a -L --sam-hit-only -t 31 -K 1G \
    ${CONTIGS} ${FASTQ} > ${ALIGNMENT/bam/sam}

# Convert to aligned BAM file (samtools should be on this docker container):
samtools view -u -b ${ALIGNMENT/bam/sam} | \
    samtools sort -@ 30 -o ${ALIGNMENT} -

rm ${ALIGNMENT/bam/sam}

# Creates ${ALIGNMENT}.bai index
samtools index ${ALIGNMENT} ${ALIGNMENT}.bai

# Creates ${CONTIGS}.fai index
samtools faidx ${CONTIGS} -o ${CONTIGS}.fai






# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# LEFT OFF --> PASTED THESE FROM THEIR WEBSITE
#   https://github.com/kishwarshafin/pepper/blob/r0.4/docs/pipeline_docker/ONT_polishing.md
# Also, does Docker have access to staging folder?

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


run_pepper_margin_deepvariant polish_assembly \
    -b "${ALIGNMENT}" \
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





# --------

# Compressing full output and sending to staging:
tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
mv ${OUTDIR}.tar.gz /staging/lnell/

mkdir ${ALIGNDIR}
mv ${ALIGNMENT} ${ALIGNMENT}.bai ${CONTIGS} ${CONTIGS}.fai ./${ALIGNDIR}/
tar -czf ${ALIGNDIR}.tar.gz ${ALIGNDIR}
mv ${ALIGNDIR}.tar.gz /staging/lnell/


# Removing files used in this job:
rm -r ${OUTDIR} minimap2-2.22_x64-linux ${FASTQ} ${CONTIGS}

