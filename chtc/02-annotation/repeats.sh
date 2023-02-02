#!/bin/bash

#' Use RepeatModeler to create a repeat library for an assembly
#' and to annotate the genome based on this library.
#' Next, use TEclass to classify many of the unknowns from RepeatModeler.
#' Last, use RepeatMasker to count repeat elements and to softmask the assembly.
#'
#' Requires arguments for output naming and assembly.
#' It also takes an optional argument for the directory where output should go.
#'
#' Usage:
#' repeats.sh [-o OUTPUT_LOC] -p OUT_PREFIX ASSEMBLY
#'
#' Options:
#'   -p Prefix for all output. Final outputs will be
#'      `${OUT_PREFIX}_repeats.tar.gz` and
#'      `${OUT_PREFIX}_contigs_masked.fasta.gz`
#'   -o Folder for all output. Defaults to `/staging/lnell/annotation`.
#'
#'
#' For Tanytarsus gracilentus assembly, the command was...
#' repeats.sh -p Tgraci /staging/lnell/assemblies/Tgraci_contigs.fasta.gz
#'
#' For Parochlus steinenii assembly, the command was...
#' repeats.sh -p Pstein /staging/lnell/assemblies/Pstein_contigs.fasta.gz
#'
#' For Chironomus riparius assembly, the command was...
#' repeats.sh -p Cripar /staging/lnell/assemblies/Cripar_contigs.fasta.gz
#'
#' For Chironomus tentans assembly, the command was...
#' repeats.sh -p Ctenta /staging/lnell/assemblies/Ctenta_contigs.fasta.gz
#'
#' For Polypedilum vanderplanki assembly, the command was...
#' repeats.sh -p Pvande /staging/lnell/assemblies/Pvande_contigs.fasta.gz
#'
#' For Polypedilum pembai assembly, the command was...
#' repeats.sh -p Ppemba /staging/lnell/assemblies/Ppemba_contigs.fasta.gz
#'
#' For Belgica antarctica assembly, the command was...
#' repeats.sh -p Bantar /staging/lnell/assemblies/Bantar_contigs.fasta.gz
#'
#' For Clunio marinus assembly, the command was...
#' repeats.sh -p Cmarin /staging/lnell/assemblies/Cmarin_contigs.fasta.gz
#'
#' For Propsilocerus akamusi assembly, the command was...
#' repeats.sh -p Pakamu /staging/lnell/assemblies/Pakamu_contigs.fasta.gz
#'
#' For Culicoides sonorensis assembly, the command was...
#' repeats.sh -p Csonor /staging/lnell/assemblies/Csonor_contigs.fasta.gz
#'
#' For Anopheles stephensi assembly, the command was...
#' repeats.sh -p Asteph /staging/lnell/assemblies/Asteph_contigs.fasta.gz
#'
#' For Aedes aegypti assembly, the command was...
#' repeats.sh -p Aaegyp /staging/lnell/assemblies/Aaegyp_contigs.fasta.gz
#'
#' For Culex quinquefasciatus assembly, the command was...
#' repeats.sh -p Cquinq /staging/lnell/assemblies/Cquinq_contigs.fasta.gz
#'
#' For Musca domestica assembly, the command was...
#' repeats.sh -p Mdomes /staging/lnell/assemblies/Mdomes_contigs.fasta.gz
#'


#' Still to do:
#' - Cripar
#' - Ctenta
#' - Pvande
#' - Ppemba
#' - Bantar
#' - Cmarin
#' - Pakamu
#' - Csonor
#' - Asteph
#' - Aaegyp
#' - Cquinq
#' - Mdomes





#' ===========================================================================
#' ===========================================================================
#'
#' Read inputs
#'
#' ===========================================================================
#' ===========================================================================


export OUTPUT_LOC=/staging/lnell/annotation

unset -v OUT_PREFIX


while getopts ":p:o:" opt; do
    case $opt in
        p)
            OUT_PREFIX="$OPTARG"
            if [[ "OUT_PREFIX" == */* ]]; then
                echo "ERROR: -p arg cannot contain '/'." 1>&2
                exit 1
            fi
            ;;
        o)
            OUTPUT_LOC=$(echo "$OPTARG" | sed 's/\/$//g')
            ;;
        \?)
            echo "Invalid option: -$OPTARG" 1>&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." 1>&2
            exit 1
            ;;
    esac
done

export ASSEMBLY_FULL_PATH=${@:$OPTIND:1}
if [ ! -f "${ASSEMBLY_FULL_PATH}" ]; then
    echo "ERROR: Your assembly file ('${ASSEMBLY_FULL_PATH}') does not exist." 1>&2
    exit 1
fi
if ! [[ "${ASSEMBLY_FULL_PATH}" =~ (.fasta|.fa|.fasta.gz|.fa.gz)$ ]]; then
    echo -n "ERROR: Assembly must end in *.fasta, *.fa, *.fasta.gz, or *.fa.gz. " 1>&2
    echo "Yours is '${ASSEMBLY_FULL_PATH}'." 1>&2
    exit 1
fi

if [ ! -f "${OUTPUT_LOC}" ]; then
    echo "ERROR: '${OUTPUT_LOC}' does not exist." 1>&2
    exit 1
fi

if (( OPTIND < $# )); then
    echo "Options passed after assembly." 1>&2
    exit 1
fi


if [ -z "$OUT_PREFIX" ]; then echo "Missing -p argument." 1>&2; exit 1; fi

export OUT_PREFIX




#' ===========================================================================
#' ===========================================================================
#'
#' Basics to start out
#'
#' ===========================================================================
#' ===========================================================================



#' This should be divisible by 4 bc that's how many threads are used per
#' job using RMBlast, the default engine in RepeatModeler
export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')

. /app/.bashrc
conda activate repeat-env


export OUT_DIR=${OUT_PREFIX}_repeats
mkdir ${OUT_DIR}
cd ${OUT_DIR}


#' Naming outputs:
#' ------
#' Library of repeats from RepeatModeler:
export REPEATS_LIB=${OUT_PREFIX}_repeats_lib.fasta
#' Soft-masked assembly for use in BRAKER:
export MASKED_ASSEMBLY=${OUT_PREFIX}_contigs_masked.fasta
#' GFF of repeats in assembly:
export REPEAT_LOCS=${OUT_PREFIX}_repeats_locs.gff3


cp ${ASSEMBLY_FULL_PATH} ./
check_exit_status "cp genome" $?
export ASSEMBLY=$(basename "${ASSEMBLY_FULL_PATH}")
if [[ "${ASSEMBLY}" == *.gz ]]; then
    gunzip ${ASSEMBLY}
    check_exit_status "gunzip genome" $?
    ASSEMBLY=${ASSEMBLY%.gz}
fi



#' ------------------------------------------------
#' prep for RepeatModeler
#' ------------------------------------------------

conda activate main-env


mv ${ASSEMBLY} ${ASSEMBLY%.*}_orig.fasta
check_exit_status "rename genome" $?

#' Below does the following (in order):
#' - Simplify sequence names in headers (remove all after first space)
#' - Undo softmasking (I'll be manually softmasking repeat regions, if desired, later)
#' - Change to narrow (i.e., multi-line) format
sed '/^>/ s/ .*//' ${ASSEMBLY%.*}_orig.fasta \
    | tr [:lower:] [:upper:] \
    | seqtk seq -l 60 - \
    > ${ASSEMBLY}
rm ${ASSEMBLY%.*}_orig.fasta

conda deactivate



#' ------------------------------------------------
#' run RepeatModeler
#' ------------------------------------------------

conda activate repeat-env


mkdir modeler
cd modeler
mv ../${ASSEMBLY} ./

# (RMBlast is default but put here for clarity)
BuildDatabase -engine rmblast -name ${OUT_PREFIX} ${ASSEMBLY} \
    >& BuildDatabase.log
check_exit_status "BuildDatabase" $?

RepeatModeler -database ${OUT_PREFIX} -pa $(( THREADS / 4 )) -LTRStruct \
    1> >(tee -a RepeatModeler.out)
check_exit_status "RepeatModeler" $?


cd ..



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





# #' ------------------------------------------------
# #' identify as many elements from library as possible using RepeatMasker
# #' ------------------------------------------------
#
#
#
# mkdir id_lib
# cd id_lib
# cp ../modeler/${OUT_PREFIX}-families.fa ${REPEATS_LIB}
#
# #' From https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
# #' and
# #' https://github.com/darencard/GenomeAnnotation/blob/8b0e5cb5b887e1d576e8483ded8e5cc754e75d23/repclassifier
# #' (accessed 30 Jan 2023)
#
# #' Split fasta file for repeat-element library into known and unknown elements.
# #' The input file must end in *.fa or *.fasta
# #' Usage:
# #' split_knowns ${REPEATS_LIBRARY_FILE}
# split_knowns () {
#     local repeats=${1}
#     if [[ "$repeats" != *.fa ]] && [[ "$repeats" != *.fasta ]]; then
#         echo "ERROR: '${repeats}' must end in *.fa or *.fasta." 1>&2
#         exit 1
#     fi
#     local output_prefix=$(echo "$repeats" | sed 's/.fa$//g;s/.fasta$//g')
#     conda activate main-env
#     seqtk seq -l0 "${repeats}" \
#         | grep -A 1 "Unknown" \
#         | seqtk seq -l 50 - \
#         > "${output_prefix}_unknown.fasta"
#     check_exit_status "unknown split" $?
#     seqtk seq -l0 "${repeats}" \
#         | grep -A 1 -vE "^[a-zA-Z]|Unknown" \
#         | seqtk seq -l 50 - \
#         > "${output_prefix}_known.fasta"
#     check_exit_status "known split" $?
#     conda deactivate
#     return 0
# }
#
#
# split_knowns ${REPEATS_LIB}
#
#
# # Match using database:
# mkdir from_db
# RepeatMasker \
#     -s -species diptera \
#     -dir from_db \
#     -pa $(( THREADS / 4 )) \
#     ${REPEATS_LIB%.fasta}_unknown.fasta \
#     1> >(tee -a from_db/RepeatMasker.out)
#
# # find FASTA for input of next search
# DBOUTPUT=$(find from_db -name "*.masked")
#
# # Match using already-IDed sequences:
# mkdir from_known
# RepeatMasker -s -pa $(( THREADS / 4 )) \
#     -dir from_known \
#     -lib ${REPEATS_LIB%.fasta}_known.fasta \
#     ${DBOUTPUT}
#
#
# # create combined output directory
# mkdir full_match
#
#
#
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
#
#
#
#
#
#
# # combine all .out and .cat.gz files together
# # note that normal header from .out files is not retained!
# cat from_db/*.out from_known/*.out \
#     | awk '{ if ($1 != "SW" && $1 != "score") print $0 }' \
#     | sed '/^$/d' | sed '/^[[:space:]]*$/d' \
#     > full_match/${OUT_DIR}.out
# cat from_db/*.cat from_known/*.cat > full_match/${OUT_DIR}.cat
# # summarize repeat classifications
# # subfamily non-ambiguous elements - elements in the same subfamily, which means they are totally nonambiguous (attach full info)
# cat full_match/${OUT_DIR}.out \
#     | awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' \
#     | awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' \
#     | sort -k1,1 -k5,5nr \
#     | cut -f 1 | sort | uniq \
#     | while read unknown; do
#         cat full_match/${OUT_DIR}.out \
#         | awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' \
#         | awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' \
#         | sort -k1,1 -k5,5nr \
#         | grep "${unknown}" \
#         | cut -f 4 | sort | uniq \
#         | awk -v OFS="\t" -v unknown="${unknown}" ' END { if (NR == 1) print unknown, $1 }';
#     done \
#     > full_match/subfamily_unambiguous_classified_elements.txt
# # family non-ambiguous elements - elements in same family even if subfamilies are not the same (attach family info)
# cat full_match/${OUT_DIR}.out \
#     | awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' \
#     | awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' \
#     | sort -k1,1 -k5,5nr \
#     | cut -f 1 | sort | uniq \
#     | while read unknown; do
#         cat full_match/${OUT_DIR}.out \
#         | awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' \
#         | awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' \
#         | sort -k1,1 -k5,5nr \
#         | grep "${unknown}" | cut -f 4 | sort | uniq \
#         | awk -v OFS="\t" -v unknown="${unknown}" ' END { if (NR > 1) print unknown }';
#     done \
#     | while read element; do
#         cat full_match/${OUT_DIR}.out \
#         | awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' \
#         | awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' \
#         | sort -k1,1 -k5,5nr \
#         | grep "${element}" | cut -f 4 \
#         | awk -v OFS="\t"  -v element="${element}" -F "/|-" '{ print element, $1, $2 }';
#     done \
#     | sort | uniq -c \
#     | awk '{ print $2 }' | uniq -u \
#     | while read clean; do
#         cat full_match/${OUT_DIR}.out \
#         | awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' \
#         | awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' \
#         | sort -k1,1 -k5,5nr \
#         | grep "${clean}" | cut -f 4 \
#         | awk -v OFS="\t"  -v clean="${clean}" -F "/|-" '{ print clean, $1"-"$2 }';
#     done \
#     | uniq \
#     > full_match/family_unambiguous_classified_elements.txt
# # combine classified elements
# cat full_match/subfamily_unambiguous_classified_elements.txt \
#     full_match/family_unambiguous_classified_elements.txt \
#     > full_match/combined_classified_elements.txt
# # chimeric elements
# cat full_match/${OUT_DIR}.out \
#     | awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' \
#     | awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' \
#     | sort -k1,1 -k5,5nr \
#     | grep -v -f <(cat full_match/combined_classified_elements.txt | cut -f 1) \
#     | cut -f 1 | sort | uniq \
#     | while read element; do
#         cat full_match/${OUT_DIR}.out \
#         | awk -v OFS="\t" '{ print $5, $6, $7, $11, $1 }' \
#         | awk '{ if ($4 != "Simple_repeat" && $4 != "Satellite" && $4 != "snRNA" && $4 != "Unknown" && $4 != "rRNA") print $0 }' \
#         | sort -k1,1 -k5,5nr \
#         | grep "${element}" | cut -f 4 | paste -s -d "," \
#         | awk -v OFS="\t" -v element="${element}" '{ print element, $1 }';
#     done \
#     > full_match/chimeric_elements.txt
# # create FASTA output of unknowns that remain unknown
# cat input.fasta | bioawk -c fastx '{ print $name }' |
# grep -wv -f <(cat full_match/combined_classified_elements.txt | cut -f 1) \
#     | seqkit grep -f - input.fasta > ${OUT_DIR}.unknown
# # create FASTA output of unknowns that were classified into known category
# cat ${APPEND} \
#     <(cat input.fasta | bioawk -c fastx '{ print $name }' \
#         | grep -wf <(cat full_match/combined_classified_elements.txt | cut -f 1) \
#         | while read seq; do
#             replace=$(cat full_match/combined_classified_elements.txt | grep "${seq}" | cut -f 2);
#             seqkit grep -p ${seq} input.fasta \
#             | seqkit fx2tab \
#             | awk -F "\t|#" -v replace="${replace}" -v OFS="\t" '{ print $1"#"replace, $3 }' \
#             | seqkit tab2fx -w0;
#         done) \
#         > ${OUT_DIR}.known
# # compress outputs to save space
# find ${OUT_DIR} -name "*.cat" | while read file; do gzip -f ${file}; done
# find ${OUT_DIR} -name "*.out" | while read file; do gzip -f ${file}; done
# find ${OUT_DIR} -name "*.masked" | while read file; do gzip -f ${file}; done
#
#
#
#
#
#
#
#
#
#
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
# # *****************************************************************************
#
#
#
#
# cd ..



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
    > ${orig_wd}/diptera_repeats_lib.fasta
cd ${orig_wd}

cat diptera_repeats_lib.fasta ${REPEATS_LIB} \
    > diptera_and_${REPEATS_LIB}


RepeatMasker \
    -s -xsmall -gff \
    -lib diptera_and_${REPEATS_LIB} \
    -pa $(( THREADS / 4 )) \
    ${ASSEMBLY} \
    1> >(tee -a RepeatMasker.out)



# *****************************************************************************
# *****************************************************************************

#                               LEFT OFF HERE

# *****************************************************************************
# *****************************************************************************


# Directory contents for use below:
#
# (repeat-env) I have no name!@lnell-16796944:/var/lib/condor/execute/slot1/dir_1838/Tgrac_repeats/masker$ ls -lh
# total 237M
# -rw-r--r-- 1 21659 21659  17M Feb  1 16:48 diptera_and_Tgrac_repeats_lib.fasta
# -rw-r--r-- 1 21659 21659  17M Feb  1 16:47 diptera_repeats_lib.fasta
# -rw-r--r-- 1 21659 21659 358K Feb  1 17:23 RepeatMasker.out
# -rw-r--r-- 1 21659 21659  90M Feb  1 09:54 Tgrac_contigs.fasta
# -rw-r--r-- 1 21659 21659  13M Feb  1 17:23 Tgrac_contigs.fasta.cat.gz
# -rw-r--r-- 1 21659 21659  90M Feb  1 17:23 Tgrac_contigs.fasta.masked
# -rw-r--r-- 1 21659 21659 6.7M Feb  1 17:23 Tgrac_contigs.fasta.out
# -rw-r--r-- 1 21659 21659 4.5M Feb  1 17:23 Tgrac_contigs.fasta.out.gff
# -rw-r--r-- 1 21659 21659 2.5K Feb  1 17:23 Tgrac_contigs.fasta.tbl
# -rw-r--r-- 1 21659 21659 203K Feb  1 16:45 Tgrac_repeats_lib.fasta
# (repeat-env) I have no name!@lnell-16796944:/var/lib/condor/execute/slot1/dir_1838/Tgrac_repeats/masker$ cd ..
# (repeat-env) I have no name!@lnell-16796944:/var/lib/condor/execute/slot1/dir_1838/Tgrac_repeats$ ls -lh
# total 216K
# drwxr-xr-x 2 21659 21659 4.0K Feb  1 17:23 masker
# drwxr-xr-x 3 21659 21659 4.0K Feb  1 16:45 modeler
# drwxr-xr-x 2 21659 21659 4.0K Feb  1 16:43 teclass
# -rw-r--r-- 1 21659 21659 203K Feb  1 16:45 Tgrac_repeats_lib.fasta
# (repeat-env) I have no name!@lnell-16796944:/var/lib/condor/execute/slot1/dir_1838/Tgrac_repeats$




cp ${ASSEMBLY}.masked ../${MASKED_ASSEMBLY}
cp ${ASSEMBLY}.out.gff ../${REPEAT_LOCS}
rm ${REPEATS_LIB} ${ASSEMBLY}
gzip < ${MASKED_ASSEMBLY} > ${MASKED_ASSEMBLY}.gz
mv ${MASKED_ASSEMBLY}.gz "${OUTPUT_LOC}"/
cd ..





#' ------------------------------------------------
#' handling output
#' ------------------------------------------------



cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz "${OUTPUT_LOC}"/
rm -r ${OUT_DIR}



