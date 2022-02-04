#!/bin/bash

. /app/.bashrc
conda activate main-env

export THREADS=16

export LONGREADS=basecalls_guppy-5.0.11.fastq
export REF=contigs_shasta_pepper-hap1.fasta

export OUT_DIR=contigs_shasta_pepper_purgedups
export OUT_FASTA=${OUT_DIR}.fasta

mkdir ${OUT_DIR}
cd ${OUT_DIR}

cp /staging/lnell/${LONGREADS}.gz ./ && gunzip ${LONGREADS}.gz
cp /staging/lnell/${REF}.gz ./ && gunzip ${REF}.gz


export ALIGN=ont_align_mm2.paf.gz

minimap2 -x map-ont -t $((THREADS - 2)) -K 1G -2 \
    ${REF} ${LONGREADS} | \
    gzip -c - > ${ALIGN}

rm ${LONGREADS}


#  (produces PB.base.cov and PB.stat files)
pbcstat ${ALIGN}
# calculate cutoffs, setting upper limit manually
# purge_dups creators say highly heterozygous genomes can have too-high
# upper limit, and based on histogram from hist_plot.py, 380 seems like a
# good one here.
calcuts -l 5 -m 181 -u 380 PB.stat > cutoffs 2> calcuts.log
# originally:
# 5	91	151	181	302	543
# now:
# 5	180	180	181	181	380

# Split an assembly and do a self-self alignment:
split_fa ${REF} > ${REF}.split
minimap2 -x asm5 -DP ${REF}.split ${REF}.split \
    -t $((THREADS - 2)) -K 1G -2 | \
    gzip -c - > ${REF}.split.self.paf.gz

# Purge haplotigs and overlaps:
purge_dups -2 -T cutoffs -c PB.base.cov \
    -a 93 \
    ${REF}.split.self.paf.gz > \
    dups.bed 2> purge_dups.log

# Get purged primary and haplotig sequences from draft assembly:
get_seqs -e dups.bed ${REF}

rm ${REF}

# Rename and move just the scaffolds to the staging directory:
mv purged.fa ${OUT_FASTA}
# Keep the uncompressed version for summaries below
gzip < ${OUT_FASTA} > ${OUT_FASTA}.gz
mv ${OUT_FASTA}.gz /staging/lnell/



# This outputs basics about scaffold (or contigs in this case) sizes:
summ-scaffs.py ${OUT_FASTA}

# This outputs BUSCO scores:
conda activate busco-env
busco \
    -m genome \
    -l diptera_odb10 \
    -i ${OUT_FASTA} \
    -o busco \
    --cpu ${THREADS}
conda deactivate


cd ..

# Now save the whole directory
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/
rm -r ${OUT_DIR}


