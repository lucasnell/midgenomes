

#!/bin/bash

# Use EDTA to create a repeat library for the Tanytarsus gracilentus genome
# and to annotate the genome based on this library.


export THREADS=$(grep "^Cpus = " $_CONDOR_MACHINE_AD | sed 's/Cpus\ =\ //')


eval "$(conda shell.bash hook)"
conda activate repeat-env


export OUT_DIR=tany_repeats
# Annotation of assembly for repeats:
export OUT_ANNO=tany_repeats_anno.gff3
# Low-threshold masked assembly for use in MAKER:
export OUT_MAKER=tany_contigs_maker.fasta
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export GENOME=tany_contigs.fasta
cp /staging/lnell/assemblies/${GENOME}.gz ./ && gunzip ${GENOME}.gz




#' ------------------------------------------------
#' run EDTA pipeline
#' ------------------------------------------------


EDTA.pl --genome ${GENOME} --threads ${THREADS} \
    --anno 1 \
    --sensitive 1 \
    --evaluate 1 \
    1> >(tee -a edta.out) \
    2> >(tee -a edta.err >&2)
check_exit_status "EDTA.pl" $?

# perl EDTA.pl [options]
#   --genome	[File]	The genome FASTA
#   --species [Rice|Maize|others]	Specify the species for identification of
#           TIR candidates. Default: others
#   --step	[all|filter|final|anno] Specify which steps you want to run EDTA.
# 			all: run the entire pipeline (default)
# 			filter: start from raw TEs to the end.
# 			final: start from filtered TEs to finalizing the run.
# 			anno: perform whole-genome annotation/analysis after TE library
#               construction.
#   --overwrite	[0|1]	If previous results are found, decide to overwrite
#           (1, rerun) or not (0, default).
#   --cds	[File]	Provide a FASTA file containing the coding sequence
#           (no introns, UTRs, nor TEs) of this genome or its close relative.
#   --curatedlib	[file]	Provided a curated library to keep consistant
#           naming and classification for known TEs.
# 			All TEs in this file will be trusted 100%, so please ONLY provide
#           MANUALLY CURATED ones here.
# 			This option is not mandatory. It's totally OK if no file is provided
#           (default).
#   --sensitive	[0|1]	Use RepeatModeler to identify remaining TEs (1) or
#           not (0, default).
# 			This step is very slow and MAY help to recover some TEs.
#   --anno	[0|1]	Perform (1) or not perform (0, default) whole-genome TE
#           annotation after TE library construction.
#   --rmout	[File]	Provide your own homology-based TE annotation instead of
#           using the EDTA library for masking. File is in RepeatMasker .out
#           format. This file will be merged with the structural-based TE
#           annotation. (--anno 1 required).
#           Default: use the EDTA library for annotation.
#   --evaluate	[0|1]	Evaluate (1) classification consistency of the TE
#           annotation. (--anno 1 required). Default: 0.
# 			This step is slow and does not affect the annotation result.
#   --exclude	[File]	Exclude bed format regions from TE annotation.
#           Default: undef. (--anno 1 required).
#   --threads|-t	[int]	Number of theads to run this script (default: 4)
#   --help|-h	Display this help info




#' ------------------------------------------------
#' softmask assembly
#' ------------------------------------------------

#' EDTA hard-masks assembly for MAKER (BRAKER in my case), but I'd prefer
#' to have it soft-masked.
#' The code below replicates what EDTA does except for having it be softmasked.

#' List of the repeat elements output from RepeatMasker, from
#' which we'll filter to make soft-masked version for BRAKER.
export FULL_MASK=$(ls *EDTA.anno/*EDTA.RM.out)
#' File to do the masking. Have to copy it here to make it executable.
#' (User doesn't have permission to do this to `/opt/...` on the cluster.)
cp /opt/conda/envs/repeat-env/share/EDTA/util/make_masked.pl ./ \
    && chmod +x make_masked.pl \
    && ln -s /opt/conda/envs/repeat-env/share/EDTA/util/substract_parallel.pl \
        substract_parallel.pl \
    && ln -s /opt/conda/envs/repeat-env/share/EDTA/util/combine_overlap.pl \
        combine_overlap.pl


./make_masked.pl -genome ${GENOME} -t ${THREADS} \
    -rmout ${FULL_MASK} \
    -maxdiv 30 -minscore 1000 -minlen 1000 -hardmask 0

mv ${GENOME}.new.masked ${GENOME}.mod.MAKER.softmasked

rm *.pl



#' ------------------------------------------------
#' handling output
#' ------------------------------------------------

conda activate main-env

cp ${GENOME}.mod.EDTA.TEanno.gff3 ${OUT_ANNO}
gzip ${OUT_ANNO}

seqtk seq -l 60 ${GENOME}.mod.MAKER.softmasked \
    > ${OUT_MAKER}
gzip ${OUT_MAKER}


mv ${OUT_ANNO}.gz ${OUT_MAKER}.gz /staging/lnell/annotation/


cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/annotation/
rm -r ${OUT_DIR}


