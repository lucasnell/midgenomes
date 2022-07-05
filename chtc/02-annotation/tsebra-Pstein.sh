#!/bin/bash


#'
#' Combine 2 BRAKER annotations (RNAseq + OrthoDB proteins) using TSEBRA,
#' and run BUSCO on resulting CDS
#'



export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

eval "$(conda shell.bash hook)"

export TARGET=/staging/lnell/annotation

export OUT_NAME=Pstein_tsebra

# Output directory for BUSCO run on CDS:
export BUSCO_OUT_DIR=busco_Pstein_cds


#' ------------------------------------------------------------
#' Copy over BRAKER1 and BRAKER2 gene predictions:
#' ------------------------------------------------------------

tar -xzf ${TARGET}/Pstein_braker_rna.tar.gz -C ./
check_exit_status "cp, extract Pstein_braker_rna" $?
tar -xzf ${TARGET}/Pstein_braker_prot.tar.gz -C ./
check_exit_status "cp, extract Pstein_braker_prot" $?



#' ------------------------------------------------------------
#' Run TSEBRA:
#' ------------------------------------------------------------

conda activate main-env

#' Because the proteins used in `Pstein_braker_prot` aren't from
#' closely related species, we'll prefer the predictions from RNAseq.
#' Hence the use of `pref_braker1.cfg` below.

tsebra.py \
    -g Pstein_braker_rna/augustus.hints.gtf,Pstein_braker_prot/augustus.hints.gtf \
    -c /opt/TSEBRA/config/pref_braker1.cfg \
    -e Pstein_braker_rna/hintsfile.gff,Pstein_braker_prot/hintsfile.gff \
    -o ${OUT_NAME}.gtf

#' This fixes IDs in output GTF, which allows it to be converted to GFF3
fix_gtf_ids.py --gtf ${OUT_NAME}.gtf --out ${OUT_NAME}_fixed.gtf
mv ${OUT_NAME}_fixed.gtf ${OUT_NAME}.gtf

conda deactivate


#' ------------------------------------------------------------
#' Create final output files:
#' ------------------------------------------------------------

conda activate annotate-env

gffread -E ${OUT_NAME}.gtf > ${OUT_NAME}.gff3

# We'll only use the gff3 version from now on.
rm ${OUT_NAME}.gtf

cp /staging/lnell/assemblies/Pstein_contigs.fasta.gz ./ \
    && gunzip Pstein_contigs.fasta.gz

# "spliced CDS for each GFF transcript"
gffread -x Pstein_cds.fasta -g Pstein_contigs.fasta ${OUT_NAME}.gff3
# "protein fasta file with the translation of CDS for each record"
gffread -y Pstein_proteins.fasta -g Pstein_contigs.fasta ${OUT_NAME}.gff3

conda deactivate



#' ------------------------------------------------------------
#' Run BUSCO on resulting CDS:
#' ------------------------------------------------------------

run_busco Pstein_cds.fasta ${THREADS}

mv busco ${BUSCO_OUT_DIR}
mv busco.out ./${BUSCO_OUT_DIR}/busco.stdout

tar -czf ${BUSCO_OUT_DIR}.tar.gz ${BUSCO_OUT_DIR}
mv ${BUSCO_OUT_DIR}.tar.gz ${TARGET}/



#' ------------------------------------------------------------
#' Move over final outputs:
#' ------------------------------------------------------------


gzip ${OUT_NAME}.gff3
gzip Pstein_cds.fasta
gzip Pstein_proteins.fasta

mv ${OUT_NAME}.gff3.gz Pstein_cds.fasta.gz Pstein_proteins.fasta.gz ${TARGET}/


rm -rf ${BUSCO_OUT_DIR} busco_downloads Pstein_braker_rna \
    Pstein_braker_prot Pstein_contigs.fasta*

