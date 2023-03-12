#!/bin/bash

#'
#' Use RepeatMasker to count repeat elements and to softmask the
#' Aedes aegypti assembly.
#'






#' ===========================================================================
#' ===========================================================================
#'
#' Read inputs
#'
#' ===========================================================================
#' ===========================================================================


export OUTPUT_LOC=/staging/lnell/annotations
export OUT_PREFIX=Aaegyp
export ASSEMBLY=Aaegyp_contigs.fasta

if [ ! -d "${OUTPUT_LOC}" ]; then
    echo "ERROR: '${OUTPUT_LOC}' does not exist." 1>&2
    exit 1
fi





#' ===========================================================================
#' ===========================================================================
#'
#' Basics to start out
#'
#' ===========================================================================
#' ===========================================================================



#' This should be divisible by 4 bc that's how many threads are used per
#' job using RMBlast, the default engine in RepeatModeler
export THREADS=$(count_threads)
if (( $THREADS < 4 )); then
    echo "ERROR: You need >= 4 threads, but you provided $THREADS." 1>&2
    exit 1
fi

. /app/.bashrc
conda activate repeat-env


export OUT_DIR=${OUT_PREFIX}_repeats
# mkdir ${OUT_DIR}
# cd ${OUT_DIR}


#' Naming outputs:
#' ------
#' Library of repeats from RepeatModeler:
export REPEATS_LIB=${OUT_PREFIX}_repeats_lib.fasta
#' Soft-masked assembly for use in BRAKER:
export MASKED_ASSEMBLY=${OUT_PREFIX}_contigs_masked.fasta
#' GFF of repeats in assembly:
export REPEAT_LOCS=${OUT_PREFIX}_repeats_locs.gff3




#' ------------------------------------------------
#' Output from previous RepeatModeler step
#' ------------------------------------------------


if [ ! -f "${OUTPUT_LOC}/${OUT_DIR}_model.tar.gz" ]; then
    echo "ERROR: '${OUTPUT_LOC}/${OUT_DIR}_model.tar.gz' does not exist." 1>&2
    exit 1
fi
tar -xf ${OUTPUT_LOC}/${OUT_DIR}_model.tar.gz -C ./
if [ ! -d "${OUT_DIR}_model" ]; then
    echo -n "ERROR: '${OUT_DIR}_model' does not exist inside " 1>&2
    echo "'${OUTPUT_LOC}/${OUT_DIR}_model.tar.gz'." 1>&2
    exit 1
fi
mv ${OUT_DIR}_model ${OUT_DIR}
cd ${OUT_DIR}






#' ------------------------------------------------
#' classify more elements using TEclass
#' ------------------------------------------------


mkdir teclass
cd teclass


#' Name of output from TEclass:
export REPEATS_LIB_TEC=${OUT_PREFIX}_repeats_lib_TEclass.fasta
#' Name for the RepeatModeler library (it's renamed here):
export REPEATS_LIB_RM=${OUT_PREFIX}_repeats_lib_RM.fasta
#' Combined output from TEclass and RepeatModeler (this is the same as
#' what will later be $REPEATS_LIB):
export REPEATS_LIB_COMBINED=${OUT_PREFIX}_repeats_lib_combined.fasta

cp ../modeler/${OUT_PREFIX}-families.fa ${REPEATS_LIB_RM}

#' Split the repeat-element library into known (IDed by RepeatModeler)
#' and unknown elements.
conda activate main-env
cat "${REPEATS_LIB_RM}" \
    | seqkit fx2tab \
    | grep -v "#Unknown" \
    | seqkit tab2fx \
    | seqtk seq -l 50 - \
    > "${REPEATS_LIB_RM%.fasta}_known.fasta"
check_exit_status "known split" $?
cat "${REPEATS_LIB_RM}" \
    | seqkit fx2tab \
    | grep "#Unknown" \
    | seqkit tab2fx \
    | seqtk seq -l 50 - \
    > "${REPEATS_LIB_RM%.fasta}_unknown.fasta"
check_exit_status "unknown split" $?
conda deactivate


conda activate teclass-env
TEclassTest.pl -r -o tec_out "${REPEATS_LIB_RM%.fasta}_unknown.fasta"
conda deactivate


cd tec_out
cp ${REPEATS_LIB_RM%.fasta}_unknown.fasta.lib ../${REPEATS_LIB_TEC}
mv ${REPEATS_LIB_RM%.fasta}_unknown.fasta.html \
    ${REPEATS_LIB_RM%.fasta}_unknown.fasta.lib  \
    ${REPEATS_LIB_RM%.fasta}_unknown.fasta.stat ../
cd ..
rm -r tec_out



# Convert TEclass output to FASTA for input to RepeatMasker
python3 << EOF

import sys
import re

# maps TEclass classes into RepeatMasker classes.
# Unfortunately, nonLTR and Retro classes in TEclass don't have parallels (to
# my knowledge) in RepeatMasker so are called Unknown.
mapper = {"DNA": "DNA", "LINE": "LINE", "LTR": "LTR", "nonLTR": "Unknown",
          "Retro": "Unknown", "SINE": "SINE", "unclear": "Unknown"}

tec_fasta = "${REPEATS_LIB_TEC}"
rm_fasta  = "${REPEATS_LIB_RM%.fasta}_known.fasta"
out_fasta = "${REPEATS_LIB_COMBINED}"

if __name__ == "__main__":
    with open(out_fasta, "wt") as out_file:
        with open(tec_fasta, "rt") as tec_file:
            for line in tec_file:
                if line.startswith(">"):
                    header_list = re.sub(r"\|.*|Unknown.*TEclass result: ", "",
                                         line.rstrip()).split("#")
                    header_list[1] = mapper[header_list[1]]
                    new_header = "#".join(header_list) + "\n"
                    b = out_file.write(new_header)
                else:
                    b = out_file.write(line)
        with open(rm_fasta, "rt") as rm_file:
            for line in rm_file:
                if line.startswith(">"):
                    new_header = re.sub(" .*", "",  line)
                    b = out_file.write(new_header)
                else:
                    b = out_file.write(line)
    sys.exit(0)
EOF


cp ${REPEATS_LIB_COMBINED} ../${REPEATS_LIB}



cd ..





#' ------------------------------------------------
#' softmask assembly using RepeatMasker
#' ------------------------------------------------

mkdir masker
cd masker

mv ../modeler/${ASSEMBLY} ./
cp ../${REPEATS_LIB} ./


orig_wd=$(pwd)
#' Get Diptera repeats from RepeatMasker library that includes RepBase
#' RepeatMasker Edition:
cd /opt/conda/envs/repeat-env/share/RepeatMasker
./famdb.py -i ./Libraries/RepeatMaskerLib.h5 families --format fasta_name \
    --include-class-in-name --ancestors --descendants 'diptera' \
    > ${orig_wd}/diptera_and_${REPEATS_LIB}
cd ${orig_wd}

cat ${REPEATS_LIB} \
    >> diptera_and_${REPEATS_LIB}


RepeatMasker \
    -s -xsmall -gff \
    -lib diptera_and_${REPEATS_LIB} \
    -pa $(( THREADS / 4 )) \
    ${ASSEMBLY} \
    1> >(tee -a RepeatMasker.out)


cp ${ASSEMBLY}.masked ../${MASKED_ASSEMBLY}
cp ${ASSEMBLY}.out.gff ../${REPEAT_LOCS}
rm ${REPEATS_LIB} ${ASSEMBLY}
cd ..





#' ------------------------------------------------
#' handling final output
#' ------------------------------------------------


gzip < ${MASKED_ASSEMBLY} > ${MASKED_ASSEMBLY}.gz
mv ${MASKED_ASSEMBLY}.gz "${OUTPUT_LOC}"/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz "${OUTPUT_LOC}"/
rm -r ${OUT_DIR}




