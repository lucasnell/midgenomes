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
#'   -o Folder for all output. Defaults to `/staging/lnell/annotations`.
#'   -p Prefix for all output. Final outputs will be
#'      `${OUT_PREFIX}_repeats.tar.gz` and
#'      `${OUT_PREFIX}_contigs_masked.fasta.gz`
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
#' ** ^^ I had to run RepeatModeler and RepeatMasker separately to avoid the
#'       72 hour limit on the UW cluster.
#'
#' For Culex quinquefasciatus assembly, the command was...
#' repeats.sh -p Cquinq /staging/lnell/assemblies/Cquinq_contigs.fasta.gz
#'
#' For Musca domestica assembly, the command was...
#' repeats.sh -p Mdomes /staging/lnell/assemblies/Mdomes_contigs.fasta.gz
#'


#' Still to do:
#' - Aaegyp






#' ===========================================================================
#' ===========================================================================
#'
#' Read inputs
#'
#' ===========================================================================
#' ===========================================================================


export OUTPUT_LOC=/staging/lnell/annotations

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
    echo -n "ERROR: Assembly must end in '.fasta', '.fa', '.fasta.gz', or '.fa.gz'. " 1>&2
    echo "Yours is '${ASSEMBLY_FULL_PATH}'." 1>&2
    exit 1
fi

if [ ! -d "${OUTPUT_LOC}" ]; then
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
if (( $THREADS < 4 )); then
    echo "ERROR: You need >= 4 threads, but you provided $THREADS." 1>&2
    exit 1
fi

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



