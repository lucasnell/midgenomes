#!/bin/bash


#'
#' Combine 2 BRAKER annotations (RNAseq + OrthoDB proteins) using TSEBRA.
#'



cp /staging/lnell/annotation/tany_braker_rna.tar.gz ./ \
    && tar -xzf tany_braker_rna.tar.gz \
    && rm tany_braker_rna.tar.gz
check_exit_status "cp, extract tany_braker_rna" $?
cp /staging/lnell/annotation/tany_braker_prot.tar.gz ./ \
    && tar -xzf tany_braker_prot.tar.gz \
    && rm tany_braker_prot.tar.gz
check_exit_status "cp, extract tany_braker_prot" $?


git clone https://github.com/Gaius-Augustus/TSEBRA.git


#' Because the proteins used in `tany_braker_prot` aren't from
#' closely related species, we'll prefer the predictions from RNAseq.
#' Hence the use of `pref_braker1.cfg` below.

./TSEBRA/bin/tsebra.py \
    -g tany_braker_rna/augustus.hints.gtf,tany_braker_prot/augustus.hints.gtf \
    -c ./TSEBRA/config/pref_braker1.cfg \
    -e tany_braker_rna/hintsfile.gff,tany_braker_prot/hintsfile.gff \
    -o tany_tsebra.gtf


gzip tany_tsebra.gtf

mv tany_tsebra.gtf.gz /staging/lnell/annotation/


rm -rf TSEBRA tany_braker_rna tany_braker_prot

