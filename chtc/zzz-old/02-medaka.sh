#!/bin/bash





# Re-install latest medaka version:
# https://chtc.cs.wisc.edu/conda-installation.shtml
# For the `conda create` step, use...
# `conda create -n medaka -c conda-forge -c bioconda medaka`
# or whatever is present at https://nanoporetech.github.io/medaka/installation.html


# If version installed is > 1.4.3, then re-run `medaka tools list_models`
# to get updated list of possible models.

# Available: r103_fast_g507, r103_fast_snp_g507, r103_fast_variant_g507,
#     r103_hac_g507, r103_hac_snp_g507, r103_hac_variant_g507,
#     r103_min_high_g345, r103_min_high_g360, r103_prom_high_g360,
#     r103_prom_snp_g3210, r103_prom_variant_g3210, r103_sup_g507,
#     r103_sup_snp_g507, r103_sup_variant_g507, r10_min_high_g303,
#     r10_min_high_g340, r941_min_fast_g303, r941_min_fast_g507,
#     r941_min_fast_snp_g507, r941_min_fast_variant_g507, r941_min_hac_g507,
#     r941_min_hac_snp_g507, r941_min_hac_variant_g507, r941_min_high_g303,
#     r941_min_high_g330, r941_min_high_g340_rle, r941_min_high_g344,
#     r941_min_high_g351, r941_min_high_g360, r941_min_sup_g507,
#     r941_min_sup_snp_g507, r941_min_sup_variant_g507, r941_prom_fast_g303,
#     r941_prom_fast_g507, r941_prom_fast_snp_g507, r941_prom_fast_variant_g507,
#     r941_prom_hac_g507, r941_prom_hac_snp_g507, r941_prom_hac_variant_g507,
#     r941_prom_high_g303, r941_prom_high_g330, r941_prom_high_g344,
#     r941_prom_high_g360, r941_prom_high_g4011, r941_prom_snp_g303,
#     r941_prom_snp_g322, r941_prom_snp_g360, r941_prom_sup_g507,
#     r941_prom_sup_snp_g507, r941_prom_sup_variant_g507, r941_prom_variant_g303,
#     r941_prom_variant_g322, r941_prom_variant_g360
# Default consensus:  r941_min_hac_g507
# Default snp:  r941_prom_hac_snp_g507
# Default variant:  r941_prom_hac_variant_g507




# open up folder with other required binaries and add these binaries to PATH:
tar -xzf medaka_bin.tar.gz
rm medaka_bin.tar.gz
export PATH=$PATH:$(pwd)/medaka_bin/bin


# Set up and activate `medaka-env` conda environment
ENVNAME=medaka-env
ENVDIR=${ENVNAME}
cp /staging/lnell/${ENVNAME}.tar.gz ./
export PATH
mkdir ${ENVDIR}
tar -xzf ${ENVNAME}.tar.gz -C ${ENVDIR}
. ${ENVDIR}/bin/activate
rm ${ENVNAME}.tar.gz


# Bringing over the basecall and draft files:
export BASECALLS=basecalls_guppy-5.0.11.fastq
export DRAFT=contigs.fasta
cp /staging/lnell/${BASECALLS}.gz ./ && gunzip ${BASECALLS}.gz
cp /staging/lnell/${DRAFT}.gz ./ && gunzip ${DRAFT}.gz



export NPROC=36
export OUTDIR=consensus_medaka
export CONSENSUS=polished_contigs.fasta
# For the model to use, the "r941_min_" should remain constant even if you
# use a later version of guppy because this specifies the pore and device.
export MODEL=r941_min_hac_g507

medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC} \
    -m ${MODEL}


echo -e "Current directory:\n"
ls -lh

echo -e "\n\nDisk used:\n"
du -h -d1



# Compressing full output and sending to staging:
tar -czf ${OUTDIR}.tar.gz ${OUTDIR}
mv ${OUTDIR}.tar.gz /staging/lnell/

# Compressing just contigs output and sending to staging:
cp ./${OUTDIR}/consensus.fasta ./
mv consensus.fasta ${CONSENSUS}
gzip ${CONSENSUS}
mv ${CONSENSUS}.gz /staging/lnell/



rm -r ./medaka_bin ./${ENVDIR} ${BASECALLS} ${DRAFT} ./${OUTDIR}
