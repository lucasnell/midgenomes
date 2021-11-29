#!/bin/bash


# have job exit if any command returns with non-zero exit status (aka failure)
set -e


export LONGREADS=basecalls_guppy-5.0.11.fastq.gz
export REF1=polished_hap1.fasta.gz
export REF2=polished_hap2.fasta.gz
# FASTA containing hap1 purged + hap2 purged
export REF=haps_purged.fasta


cp /staging/lnell/${LONGREADS} ./
cp /staging/lnell/${REF1} ./
cp /staging/lnell/${REF2} ./


# purge_dups (version 1.2.5)
# runner (commit 73a4d1136bee4477713c11606d16afebc8d0f805)
# minimap2 (version 2.21)

tar -xzf purge_progs.tar.gz
rm purge_progs.tar.gz

cd purge_progs
mv minimap2-2.21 purge_dups runner ../
cd ..
rm -r purge_progs

export PATH=$PATH:$(pwd)/minimap2-2.21/
export PATH=$PATH:$(pwd)/purge_dups/bin

cd runner && python3 setup.py install --user && cd ..



# ================================================
# Using `purge_dups/scripts/run_purge_dups.py` to run pipeline
# ================================================


echo $(pwd)/$LONGREADS > nano.fofn

# Do not adjust! This is automatically set by `pd_config.py`.
OUTDIR=${REF1/.fasta/}
OUTDIR=${OUTDIR/.fa/}
OUTDIR=${OUTDIR/.gz/}
export OUTDIR


./purge_dups/scripts/pd_config.py \
    -l pd_inputs1 \
    -n config-tany1.json \
    --skipB \
    $REF1 nano.fofn

python3 ./purge_dups/scripts/run_purge_dups.py \
    -p bash \
    config-tany1.json \
    ./purge_dups/bin \
    tany1


# Add to combined fasta file for both haplotypes
# (Have to add suffix to avoid duplicated names)
sed 's/>.*/&_hap1/' ./${OUTDIR}/seqs/${OUTDIR}.purged.fa >> ${REF}

# no longer needed:
rm -rf config-tany1.json nano.fofn pd_inputs1

export OUTDIR1=${OUTDIR}



# ================================================
# Doing it again on haplotype 2:
# ================================================


echo $(pwd)/$LONGREADS > nano.fofn

# Do not adjust! This is automatically set by `pd_config.py`.
OUTDIR=${REF2/.fasta/}
OUTDIR=${OUTDIR/.fa/}
OUTDIR=${OUTDIR/.gz/}
export OUTDIR


./purge_dups/scripts/pd_config.py \
    -l pd_inputs2 \
    -n config-tany2.json \
    --skipB \
    $REF2 nano.fofn

python3 ./purge_dups/scripts/run_purge_dups.py \
    -p bash \
    config-tany2.json \
    ./purge_dups/bin \
    tany2


# Add to combined fasta file for both haplotypes
# (Have to add suffix to avoid duplicated names)
sed 's/>.*/&_hap2/' ./${OUTDIR}/seqs/${OUTDIR}.purged.fa >> ${REF}


# no longer needed:
rm -rf config-tany2.json nano.fofn pd_inputs2

export OUTDIR2=${OUTDIR}



# ================================================
# Doing it again on combined purged haplotigs:
# ================================================


echo $(pwd)/$LONGREADS > nano.fofn

# Do not adjust! This is automatically set by `pd_config.py`.
OUTDIR=${REF/.fasta/}
OUTDIR=${OUTDIR/.fa/}
OUTDIR=${OUTDIR/.gz/}
export OUTDIR

./purge_dups/scripts/pd_config.py \
    -l pd_inputs \
    -n config-tany.json \
    --skipB \
    $REF nano.fofn


python3 ./purge_dups/scripts/run_purge_dups.py \
    -p bash \
    config-tany.json \
    ./purge_dups/bin \
    tany



mv ${REF} ${OUTDIR1} ${OUTDIR2} ./${OUTDIR}/
tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
mv ${OUTDIR}.tar.gz /staging/lnell/


rm -rf minimap2-2.21 purge_dups runner ${REF1} ${REF2} ${LONGREADS} \
    ${OUTDIR} config-tany.json nano.fofn pd_inputs






#
#
# # ================================================
# # Doing purge_dups pipeline myself
# # ================================================
#
# export pri_asm=$REF
#
# # -------
# # Run minimap2 to align pacbio data and generate paf files, then calculate
# # read depth histogram and base-level read depth.
# # -------
# minimap2 -xmap-ont $pri_asm $LONGREADS | gzip -c - > $LONGREADS.paf.gz
#
# # (produces PB.base.cov and PB.stat files)
# pbcstat *.paf.gz
# calcuts PB.stat > cutoffs 2>calcults.log
#
#
# # -------
# # Split an assembly and do a self-self alignment.
# # -------
# split_fa $pri_asm > $pri_asm.split
# minimap2 -xasm5 -DP $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz
#
#
# # -------
# # Purge haplotigs and overlaps.
# # -------
# purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log
#
#
# # -------
# # Get purged primary and haplotig sequences from draft assembly.
# # -------
# get_seqs -e dups.bed $pri_asm
#
#
# # -------
# # Merge hap.fa and $hap_asm and redo the above steps to get a decent haplotig set.
# # -------
#
# # LEFT OFF
#
#




