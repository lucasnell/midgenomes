#!/bin/bash


#'
#' Combine 2 BRAKER annotations (RNAseq + OrthoDB proteins) using TSEBRA,
#' and run BUSCO on resulting CDS
#'



export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

eval "$(conda shell.bash hook)"

export TARGET=/staging/lnell/annotation

export OUT_NAME=tany_tsebra

# Output directory for BUSCO run on CDS:
export BUSCO_OUT_DIR=tany_cds_busco


#' ------------------------------------------------------------
#' Copy over BRAKER1 and BRAKER2 gene predictions:
#' ------------------------------------------------------------

cp ${TARGET}/tany_braker_rna.tar.gz ./ \
    && tar -xzf tany_braker_rna.tar.gz \
    && rm tany_braker_rna.tar.gz
check_exit_status "cp, extract tany_braker_rna" $?
cp ${TARGET}/tany_braker_prot.tar.gz ./ \
    && tar -xzf tany_braker_prot.tar.gz \
    && rm tany_braker_prot.tar.gz
check_exit_status "cp, extract tany_braker_prot" $?



#' ------------------------------------------------------------
#' Run TSEBRA:
#' ------------------------------------------------------------

conda activate main-env

#' Because the proteins used in `tany_braker_prot` aren't from
#' closely related species, we'll prefer the predictions from RNAseq.
#' Hence the use of `pref_braker1.cfg` below.

tsebra.py \
    -g tany_braker_rna/augustus.hints.gtf,tany_braker_prot/augustus.hints.gtf \
    -c /opt/TSEBRA/config/pref_braker1.cfg \
    -e tany_braker_rna/hintsfile.gff,tany_braker_prot/hintsfile.gff \
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

cp /staging/lnell/assemblies/tany_contigs.fasta.gz ./ \
    && gunzip tany_contigs.fasta.gz

# "spliced CDS for each GFF transcript"
gffread -x tany_cds.fasta -g tany_contigs.fasta ${OUT_NAME}.gff3
# "protein fasta file with the translation of CDS for each record"
gffread -y tany_proteins.fasta -g tany_contigs.fasta ${OUT_NAME}.gff3

conda deactivate

# 2022-06-17 09:56:58 INFO:
#
# 	--------------------------------------------------
# 	|Results from dataset diptera_odb10               |
# 	--------------------------------------------------
# 	|C:93.4%[S:88.0%,D:5.4%],F:1.0%,M:5.6%,n:3285     |
# 	|3069	Complete BUSCOs (C)                       |
# 	|2892	Complete and single-copy BUSCOs (S)       |
# 	|177	Complete and duplicated BUSCOs (D)        |
# 	|33	    Fragmented BUSCOs (F)                     |
# 	|183	Missing BUSCOs (M)                        |
# 	|3285	Total BUSCO groups searched               |
# 	--------------------------------------------------



#' ------------------------------------------------------------
#' Run BUSCO on resulting CDS:
#' ------------------------------------------------------------

run_busco tany_cds.fasta ${THREADS}

mv busco ${BUSCO_OUT_DIR}
mv busco.out ./${BUSCO_OUT_DIR}/busco.stdout

tar -czf ${BUSCO_OUT_DIR}.tar.gz ${BUSCO_OUT_DIR}
mv ${BUSCO_OUT_DIR}.tar.gz ${TARGET}/



#' ------------------------------------------------------------
#' Move over final outputs:
#' ------------------------------------------------------------


gzip ${OUT_NAME}.gff3
gzip tany_cds.fasta
gzip tany_proteins.fasta

mv ${OUT_NAME}.gff3.gz tany_cds.fasta.gz tany_proteins.fasta.gz ${TARGET}/


rm -rf ${BUSCO_OUT_DIR} busco_downloads tany_braker_rna \
    tany_braker_prot tany_contigs.fasta*

