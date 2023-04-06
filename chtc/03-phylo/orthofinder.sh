#!/bin/bash


#'
#' Find orthogroups using OrthoFinder.
#' Ran this using interactive job with 48 threads, 64 GB RAM, 150 GB disk.
#'


. /app/.bashrc
conda activate phylo-env

export THREADS=$(count_threads)

export TARGET=/staging/lnell

export OUT_DIR=chir_orthofinder
mkdir ${OUT_DIR}
cd ${OUT_DIR}

export PROT_FOLDER=chir_proteins

cp -r /staging/lnell/proteins ./ \
    && cd proteins \
    && gunzip *.faa.gz \
    && for f in *.faa; do mv $f ${f/_proteins/_wdups}; done \
    && export SPECIES_NAMES=($(find . -name "*.faa" | sed 's/_wdups.faa//g' | tr -d './'))
check_exit_status "moving, renaming proteins" $?

# Remove duplicate proteins from each FASTA file:
R --vanilla << EOF > >(tee -a removing-dups.log)

library(ape)

# Convert raw vector to single string:
to_str <- function(x) paste(rawToChar(x), collapse = "")

spp_names <- trimws(strsplit("${SPECIES_NAMES[@]}", "\\\\s+")[[1]])

cat("---- duplicates removed ----\n\n")

for (n in spp_names) {
    rfa <- read.FASTA(paste0(n, "_wdups.faa"), "AA")
    names(rfa) <-  lapply(strsplit(names(rfa), " "), function(x) x[[1]])
    stopifnot(!any(duplicated(names(rfa))))

    lens <- as.integer(sapply(rfa, length))
    # sum(duplicated(lens))
    # mean(duplicated(lens))

    #' Unique lengths that aren't duplicated:
    nond_lens <- unique(lens[duplicated(lens)])
    # length(nond_lens)

    #' Now go through each duplicated length and see if any proteins themselves
    #' are duplicated.
    dup_seq <- logical(length(rfa))

    for (len in nond_lens) {
        len_inds <- which(lens == len)
        seqs <- do.call(c, lapply(rfa[len_inds], to_str))
        seq_dups <- duplicated(seqs)
        if (any(seq_dups)) dup_seq[len_inds[seq_dups]] <- TRUE
    }

    cat(sprintf("%s\n  total   = %i\n  percent = %.1f\n\n",
                n, sum(dup_seq), 100 * mean(dup_seq)))

    nond_rfa <- rfa[!dup_seq]

    write.FASTA(nond_rfa, paste0(n, ".faa"))

}
EOF
check_exit_status "creating sets of non-duplicate proteins" $?

rm *_wdups.faa \
    && cd .. \
    && mv proteins ${PROT_FOLDER}
check_exit_status "removing duplicate proteins" $?



#' Create simplified time tree from MCMCTree output.
#' Use this here in OrthoFinder and save for potentially using elsewhere.
export SPECIES_TREE_MCMCTREE=chir_mcmctree.tre
export SPECIES_TREE=chir_mcmctree.nwk
cat /staging/lnell/phylo/${SPECIES_TREE_MCMCTREE} \
    | sed -e 's/\[[^][]*\]//g' \
    | grep "UTREE" \
    | sed 's/[^(]*//' \
    > ${SPECIES_TREE}
check_exit_status "create simple newick time tree" $?
cp ${SPECIES_TREE} /staging/lnell/phylo/

orthofinder -f ${PROT_FOLDER} -t ${THREADS} -a $(( THREADS / 4 )) \
    -s ${SPECIES_TREE} \
    2>&1 \
    | tee orthofinder.log
check_exit_status "run OrthoFinder" $?

# Move OrthoFinder output out of the proteins folder and rename:
cd ${PROT_FOLDER} \
    && mv OrthoFinder OrthoFinder_tmp \
    && cd OrthoFinder_tmp \
    && mv Results_* orthofinder-output \
    && mv orthofinder-output ../../ \
    && cd .. \
    && rm -r OrthoFinder_tmp \
    && cd ..
check_exit_status "move, rename OrthoFinder output" $?



cd ..
tar -czf ${OUT_DIR}.tar.gz ${OUT_DIR}

mv ${OUT_DIR}.tar.gz ${TARGET}/

rm -r ${OUT_DIR}

