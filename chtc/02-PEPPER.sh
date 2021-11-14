#!/bin/bash
# Polish Tanytarsus gracilentus genome

# Note that we have to use `r0.4` version of PEPPER:
# https://github.com/kishwarshafin/pepper/tree/r0.4


export OUTDIR=polish_pepper
export ALIGNDIR=ont_align_minimap2
export SAMPLE_NAME=tany
export THREADS=32

export POLISHED_HAP1=polished_hap1.fasta
export POLISHED_HAP2=polished_hap2.fasta

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

./minimap2-2.22_x64-linux/minimap2 -a -L --sam-hit-only \
    -t $((THREADS - 1)) -K 1G \
    ${CONTIGS} ${FASTQ} > ${ALIGNMENT/bam/sam}

# Convert to aligned BAM file (samtools should be on this docker container):
samtools view -u -b ${ALIGNMENT/bam/sam} | \
    samtools sort -@ $((THREADS - 2)) -o ${ALIGNMENT} -

rm ${ALIGNMENT/bam/sam}

# Creates ${ALIGNMENT}.bai index
samtools index ${ALIGNMENT} ${ALIGNMENT}.bai

# Creates ${CONTIGS}.fai index
samtools faidx ${CONTIGS} -o ${CONTIGS}.fai





run_pepper_margin_deepvariant polish_assembly \
    -b "${ALIGNMENT}" \
    -f "${CONTIGS}" \
    -o "${OUTDIR}" \
    -t ${THREADS} \
    -s ${SAMPLE_NAME} \
    --ont

# this generates 2 VCFs, one per haplotype
HAP1_VCF=PEPPER_MARGIN_DEEPVARIANT_ASM_POLISHED_HAP1.vcf.gz
HAP2_VCF=PEPPER_MARGIN_DEEPVARIANT_ASM_POLISHED_HAP2.vcf.gz

bcftools consensus \
    -f "${CONTIGS}" \
    -H 2 \
    -s "${SAMPLE_NAME}" \
    -o "${OUTDIR}/${POLISHED_HAP1}" \
    "${OUTDIR}/${HAP1_VCF}"

bcftools consensus \
    -f "${CONTIGS}" \
    -H 2 \
    -s "${SAMPLE_NAME}" \
    -o "${OUTDIR}/${POLISHED_HAP2}" \
    "${OUTDIR}/${HAP2_VCF}"





# --------

# Compressing output and sending to staging

mkdir ${ALIGNDIR}
mv ${ALIGNMENT} ${ALIGNMENT}.bai ${CONTIGS} ${CONTIGS}.fai ./${ALIGNDIR}/
tar -czf ${ALIGNDIR}.tar.gz ${ALIGNDIR}
mv ${ALIGNDIR}.tar.gz /staging/lnell/

tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
mv ${OUTDIR}.tar.gz /staging/lnell/

cp ${OUTDIR}/${POLISHED_HAP1} ./
cp ${OUTDIR}/${POLISHED_HAP2} ./
gzip ${POLISHED_HAP1}
gzip ${POLISHED_HAP2}
mv ${POLISHED_HAP1}.gz ${POLISHED_HAP2}.gz /staging/lnell/


# Removing other files used in this job:
rm -r ${OUTDIR} ${ALIGNDIR} minimap2-2.22_x64-linux ${FASTQ} ${CONTIGS}

