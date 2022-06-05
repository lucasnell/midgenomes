#!/bin/bash

# Use EDTA to create a repeat library for the Tanytarsus gracilentus genome
# and to annotate the genome based on this library.


export THREADS=24


. /app/.bashrc
conda activate annotate-env

# # From here: https://github.com/oushujun/EDTA/issues/250
# export PERL5LIB=/

export OUT_DIR=tany_repeats
# Annotation of assembly for repeats:
export OUT_ANNO=tany_repeat_anno.gff3
mkdir ${OUT_DIR}
cd ${OUT_DIR}


export GENOME=tany_scaffolds.fasta
cp /staging/lnell/${GENOME}.gz ./ && gunzip ${GENOME}.gz


# >> Names are already simplified
# # EDTA docs recommend simplifying sequence names:
# # I'll do this by removing everything after the first comma or the first
# # space in the seq names:
# sed -i -r 's/\,.+//' ${GENOME}
# sed -i -r 's/\ .+//' ${GENOME}


EDTA.pl --genome ${GENOME} --threads ${THREADS} --anno 1 --sensitive 1 \
    --evaluate 1
status=$?

if [ ! "$status" -eq "0" ]; then
    echo "EDTA failed with exit status $status" 1>&2
    cd ..
    tar -czf ERROR_${OUT_DIR}.tar.gz ${OUT_DIR}
    mv ERROR_${OUT_DIR}.tar.gz /staging/lnell/
    rm -r ${OUT_DIR}
    exit $status
fi

rm ${GENOME}

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



# Saving output:

cp *.mod.EDTA.TEanno.gff3 ${OUT_ANNO}
gzip ${OUT_ANNO}
mv ${OUT_ANNO}.gz /staging/lnell/

cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}
mv ${OUT_DIR}.tar.gz /staging/lnell/
rm -r ${OUT_DIR}


