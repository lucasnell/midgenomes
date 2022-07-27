library(tidyverse)
library(parallel)

#' Note that to run this on Windows, you should change `mcmapply` calls to
#' `mapply` below.
#' I've highlighted these lines with `<<<<<<<<<<<<<<<`.

options("mc.cores" = max(1L, detectCores() - 2L))


tany_gff_df <- read_tsv("~/_data/annotation/tany_tsebra.gff3.gz", skip = 3,
                        col_names = c("seqid", "source", "type", "start",
                                      "end", "score", "strand", "phase",
                                      "attributes"),
                        col_types = "ccciicccc")
pstein_gff_df <- read_tsv("~/_data/annotation/Pstein_tsebra.gff3.gz", skip = 3,
                          col_names = c("seqid", "source", "type", "start",
                                        "end", "score", "strand", "phase",
                                        "attributes"),
                          col_types = "ccciicccc")

#' Just exons, to be added to transcript dataframes later.
tany_exon_df <- tany_gff_df %>%
    filter(type == "exon") %>%
    mutate(attributes = str_remove_all(attributes, "Parent=")) %>%
    rename(protein = attributes) %>%
    select(seqid, start, end, protein)
pstein_exon_df <- pstein_gff_df %>%
    filter(type == "exon") %>%
    mutate(attributes = str_remove_all(attributes, "Parent=")) %>%
    rename(protein = attributes) %>%
    select(seqid, start, end, protein)



#' Extract gene and protein from transcript attributes string.
#' NOTE: only works for transcripts!
get_gene_protein <- function(.df) {
    attr_str <- .df[["attributes"]]
    if (any(grepl("Parent", attr_str))) {
        stop("`get_gene_protein` only works for transcripts!")
    }
    gp_list <- attr_str %>%
        str_split(";") %>%
        map(function(.str) {
            g <- .str[which(grepl("^geneID=", .str))]
            p <- .str[which(grepl("^ID=", .str))]
            c(gsub("geneID=", "", g), gsub("ID=", "", p))
        })
    .df[["gene"]] = map_chr(gp_list, ~ .x[[1]])
    .df[["protein"]] = map_chr(gp_list, ~ .x[[2]])
    return(.df)
}

#' Add exon starts and ends for each transcript. It also orders exons by
#' start position. (This is probably already done in the GFF3 file, but better
#' safe than sorry.)
get_exons <- function(.df, .exon_df) {
    .df$exons = mcmapply(                                   # <<<<<<<<<<<<<<<
        FUN = function(p, s) {
            exd <- .exon_df[.exon_df$protein == p,]
            stopifnot(all(exd$seqid == s))
            exd <- exd[order(exd$start),]
            return(cbind(exd$start, exd$end))
        },
        .df$protein, .df$seqid,
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    return(.df)
}
#' Get intron lengths for transcripts.
#' Dataframe must have `exons`, `start`, and `end` columns.
get_intron_lens <- function(.df) {
    #' For one transcript:
    intron_lens_one_trans <- function(exon_mat, tran_start, tran_end) {
        n_exons <- nrow(exon_mat)
        stopifnot(n_exons > 0)
        #' `n_introns` means all possible introns, so it includes before
        #' and after first and last exons, respectively.
        #' These may not be present and will be filtered out later,
        #' as will any length-0 gaps between exons.
        n_introns <- n_exons + 1
        introns <- numeric(n_introns)
        if (exon_mat[1,1] > tran_start) {
            introns[1] <- exon_mat[1,1] - tran_start
        }
        if (exon_mat[n_exons,2] < tran_end) {
            introns[n_introns] <- tran_end - exon_mat[n_exons,2]
        }
        if (n_exons > 1) {
            for (j in 2:n_exons) {
                introns[j] <- exon_mat[j,1] - exon_mat[(j-1),2] - 1
            }
        }
        introns <- introns[introns > 0]
        return(introns)
    }
    .df$intron_lens = mcmapply(                              # <<<<<<<<<<<<<<<
        FUN = intron_lens_one_trans,
        .df$exons, .df$start, .df$end,
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    return(.df)
}



#' Takes ~2.3 sec on my machine (6 threads)
tany_trans_df <- tany_gff_df %>%
    filter(type == "transcript") %>%
    get_gene_protein() %>%
    select(seqid, start, end, gene, protein) %>%
    get_exons(tany_exon_df) %>%
    get_intron_lens()

tany_trans_df$intron_lens %>% do.call(what = c) %>% mean()



#' Takes ~2.4 sec
pstein_trans_df <- pstein_gff_df %>%
    filter(type == "transcript") %>%
    get_gene_protein() %>%
    select(seqid, start, end, gene, protein) %>%
    get_exons(pstein_exon_df) %>%
    get_intron_lens()

pstein_trans_df$intron_lens %>% do.call(what = c) %>% mean()











#' ============================================================================
#' ============================================================================
#'
#' Numbers of unique proteins and genes
#'
#' ============================================================================
#' ============================================================================

tany_trans_df %>% distinct(protein) %>% nrow()
tany_trans_df %>% distinct(gene) %>% nrow()

pstein_trans_df %>% distinct(protein) %>% nrow()
pstein_trans_df %>% distinct(gene) %>% nrow()

