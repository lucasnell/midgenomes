#!/bin/bash

#' This script combines SNAPE output from all quality-passing gSYNC
#' objects into three:
#' one for temporal samples, one for spatial samples, and one for all.



source /staging/lnell/helpers.sh


export TIME_SAMPS=(H-1-06_S42 H-1-11_S9 KS-1-11_S45 KS-1-77_S47 KS-2-09-B_S48 \
    SN-1-02_S55 SN-1-07_S23 SN-1-08_S24 SN-1-09_S25 SN-1-10_S26 SN-1-13_S56 \
    SN-1-14_S27 SN-1-97_S61 SN-1-99_S28 SN-2-00_S29 SN-2-02_S30 SN-2-07_S31 \
    SN-2-08_S32 SN-2-15_S34 SN-2-92_S35 SN-2-97_S37 SN-2-99_S38)
export SPACE_SAMPS=(Ash-19_S5 Blik-19_S6 Blo-19_S7 Ellid-19_S8 Hrisatjorn_S11 \
    Lei-19_S13 Lys-19_S14 Mik-19_S15 MikW-19_S16 MyBR-19_S17 MyKS-19-B_S18 \
    MySN-19_S19 Ola-19_S20 Pernuvatn_S21 Skj-19_S22 Sva-19_S39 Tjornin_S40 \
    Vik-19_S41)



#' ========================================================================
#' Inputs
#' ========================================================================

export OUT_DIR=combine-snapes
mkdir ${OUT_DIR}
cd ${OUT_DIR}

for samp in ${TIME_SAMPS[@]} ${SPACE_SAMPS[@]}; do
    cp /staging/lnell/dna/snape/${samp}_snape_masked.sync.gz ./
    check_exit_status "cp ${samp}" $?
    mv ${samp}_snape_masked.sync.gz ${samp}.sync.gz
    check_exit_status "mv ${samp}" $?
done

cp /staging/lnell/dna/final-snape/combine-snapes.py ./

#' ========================================================================
#' Outputs
#' ========================================================================

# Where to send everything when done:
export TARGET=/staging/lnell/dna/final-snape

# Final files / directories
export ALL_OUT=all_snape_masked.sync.gz
export ALL_OUT_NAMES=${ALL_OUT/.sync/.names}
export ALL_OUT_NOBLANKS=${ALL_OUT/.sync.gz/_noblanks.sync.gz}

export TIME_OUT=time_snape_masked.sync.gz
export TIME_OUT_NAMES=${TIME_OUT/.sync/.names}
export TIME_OUT_NOBLANKS=${TIME_OUT/.sync.gz/_noblanks.sync.gz}

export SPACE_OUT=space_snape_masked.sync.gz
export SPACE_OUT_NAMES=${SPACE_OUT/.sync/.names}
export SPACE_OUT_NOBLANKS=${SPACE_OUT/.sync.gz/_noblanks.sync.gz}


#' ========================================================================
#' Combine files
#' ========================================================================

./combine-snapes.py -o ${ALL_OUT} ${TIME_SAMPS[@]/%/.sync.gz} \
    ${SPACE_SAMPS[@]/%/.sync.gz}
check_exit_status "combine-snapes.py (all)" $?
for samp in ${TIME_SAMPS[@]} ${SPACE_SAMPS[@]}; do
    echo ${samp} >> ${ALL_OUT_NAMES/.gz/}
done
gzip ${ALL_OUT_NAMES/.gz/}
gunzip -c ${ALL_OUT} | grep -Fv ".:.:.:.:.:.:." | gzip > ${ALL_OUT_NOBLANKS}
check_exit_status "make_no_blanks (all)" $?
mv ${ALL_OUT} ${ALL_OUT_NAMES} ${ALL_OUT_NOBLANKS} ${TARGET}/


./combine-snapes.py -o ${TIME_OUT} ${TIME_SAMPS[@]/%/.sync.gz}
check_exit_status "combine-snapes.py (time)" $?
for samp in ${TIME_SAMPS[@]}; do
    echo ${samp} >> ${TIME_OUT_NAMES/.gz/}
done
gzip ${TIME_OUT_NAMES/.gz/}
gunzip -c ${TIME_OUT} | grep -Fv ".:.:.:.:.:.:." | gzip > ${TIME_OUT_NOBLANKS}
check_exit_status "make_no_blanks (time)" $?
mv ${TIME_OUT} ${TIME_OUT_NAMES} ${TIME_OUT_NOBLANKS} ${TARGET}/
rm ${TIME_SAMPS[@]/%/.sync.gz}


./combine-snapes.py -o ${SPACE_OUT} ${SPACE_SAMPS[@]/%/.sync.gz}
check_exit_status "combine-snapes.py (space)" $?
for samp in ${SPACE_SAMPS[@]}; do
    echo ${samp} >> ${SPACE_OUT_NAMES/.gz/}
done
gzip ${SPACE_OUT_NAMES/.gz/}
gunzip -c ${SPACE_OUT} | grep -Fv ".:.:.:.:.:.:." | gzip > ${SPACE_OUT_NOBLANKS}
check_exit_status "make_no_blanks (space)" $?
mv ${SPACE_OUT} ${SPACE_OUT_NAMES} ${SPACE_OUT_NOBLANKS} ${TARGET}/
rm ${SPACE_SAMPS[@]/%/.sync.gz}


cd ..
rm -r ${OUT_DIR}

