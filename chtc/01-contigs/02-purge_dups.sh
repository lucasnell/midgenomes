#!/bin/bash


# have job exit if any command returns with non-zero exit status (aka failure)
set -e


export LONGREADS=basecalls_guppy-5.0.11.fastq.gz
export REF=polished_hap1.fasta.gz


cp /staging/lnell/${LONGREADS} ./
cp /staging/lnell/${REF} ./


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


cp ./${OUTDIR}/seqs/${OUTDIR}.purged.fa ./
mv ${OUTDIR}.purged.fa haploid_purge_dups.fasta
gzip haploid_purge_dups.fasta
mv haploid_purge_dups.fasta.gz /staging/lnell/

mv ${OUTDIR} haploid_purge_dups
export OUTDIR=haploid_purge_dups
tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
mv ${OUTDIR}.tar.gz /staging/lnell/


rm -rf minimap2-2.21 purge_dups runner ${REF} ${LONGREADS} \
    ${OUTDIR} config-tany.json nano.fofn pd_inputs



