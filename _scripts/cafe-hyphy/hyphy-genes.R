
#'
#' This script produces files necessary to use in HyPhy only 1-to-1 HOGs that
#' are linked to GO terms related to stressful environments.
#'

library(GO.db)
# ^^ make sure this is loaded before tidyverse ^^

source("_scripts/00-preamble.R")

# Overwrite previous versions of CSV files produced here?
.overwrite <- FALSE


read_fasta <- function(fn) {
    z <- read_lines(fn, progress = FALSE)
    heads <- which(grepl("^>", z))
    nh <- length(heads)
    sl <- c(map(2:nh, \(i) (heads[i-1]+1L):(heads[i]-1L)),
            list((heads[nh]+1L):length(z)))
    zz <- map_chr(sl, \(inds) paste0(z[inds], collapse = ""))
    names(zz) <- z[heads]
    return(zz)
}
write_faa <- function(faa_obj, fn) {
    stopifnot(inherits(faa_obj, "character"))
    stopifnot(!is.null(names(faa_obj)))
    stopifnot(all(grepl("^>", names(faa_obj))))
    file_conn <- file(fn, "at")
    for (i in 1:length(faa_obj)) {
        writeLines(c(names(faa_obj)[i], faa_obj[[i]]), file_conn)
    }
    close(file_conn)
    invisible(NULL)
}


hog_gos <- dirs$orthofinder_extr |>
    paste0("/Single_Copy_HOG_GO/N0-GO-by-HOG.tsv") |>
    read_tsv(col_types = cols()) |>
    mutate(go = go |>
               toupper() |>
               str_split(";")) |>
    unnest(go)




#' Get a GO term plus all its offspring. Can optionally be done recursively.
#'
get_offs <- function(g, recursive = TRUE) {
    stopifnot(length(g) == 1)
    if (recursive) {
        g1 <- g
        g2 <- rep(NA_character_, 2)
        while (length(g2) > length(g1)) {
            if (all(!is.na(g2))) g1 <- g2
            g2 <- map(g1, \(.g) c(.g, GOBPOFFSPRING[[.g]])) |>
                do.call(what = c) |>
                na.exclude() |>
                paste() |>
                unique()
        }
    } else {
        g2 <- map(g, \(.g) c(.g, GOBPOFFSPRING[[.g]])) |>
            do.call(what = c) |>
            na.exclude() |>
            paste() |>
            unique()
    }
    return(sort(g2))
}

get_term <- function(gos) {
    suppressMessages(AnnotationDbi::select(GO.db, gos, "TERM")[,"TERM"])
}


#' Within response to stimulus > response to chemical >
#' response to inorganic substance:
#'
#'   * GO:0010038 - response to metal ion
#'
#'
#' Within response to stimulus > response to abiotic stimulus >
#' response to radiation (response to radiation was not chosen because it
#' includes light stimulus):
#'
#'   * GO:0010212 - response to ionizing radiation
#'
#'
#' Within response to stimulus > response to stress:
#'
#'   * GO:0034059 - response to anoxia
#'   * GO:0009409 - response to cold
#'   * GO:0009408 - response to heat
#'   * GO:0001666 - response to hypoxia
#'   * GO:0006979 - response to oxidative stress
#'
#'
#' Within response to stimulus > response to external stimulus > response to
#' stress > response to other organism:
#'
#'   * GO:0098542 - defense response to other organism
#'


#' Dataframe of GO terms that I'm interested in analyzing.
#' Start with terms listed above, then I got all offspring GO terms
#' (recursively) for these terms.
#'

focal_go_df <- tibble(go = c("GO:0010038", "GO:0010212", "GO:0034059",
                             "GO:0009409", "GO:0009408", "GO:0001666",
                             "GO:0006979", "GO:0098542"),
                      term = get_term(go),
                      # parent GO term plus offspring (recursive)
                      offspring = safe_mclapply(go, get_offs),
                      hogs = map(offspring, \(o) unique(hog_gos$hog[hog_gos$go %in% o])))
focal_go_df

# Write to CSV for use in summarizing output later:
if (!file.exists("_data/hyphy-focal-hog-go.csv") || .overwrite) {
    focal_go_df |>
        mutate(across(offspring:hogs, \(x) map_chr(x, paste, collapse = ";"))) |>
        write_csv("_data/hyphy-focal-hog-go.csv")
}


#' Below shows that some HOGs are shared among the GO terms:
shared_mat <- matrix(0L, nrow(focal_go_df), nrow(focal_go_df))
for (i in 2:nrow(focal_go_df)) {
    for (j in 1:(i-1)) {
        shared_mat[i,j] = sum(focal_go_df$hogs[[i]] %in% focal_go_df$hogs[[j]])
        shared_mat[j,i] = sum(focal_go_df$hogs[[j]] %in% focal_go_df$hogs[[i]])
    }
}; rm(i, j)
sum(shared_mat[lower.tri(shared_mat)] > 0)
sum(shared_mat[upper.tri(shared_mat)] > 0)
#' Whole matrix with total HOGs on diagonal:
shared_mat + diag(map_int(focal_go_df$hogs, length))



#' Extract names of genes and species for all HOGs we want to analyze.
#' We will have to use CDS in HyPhy, and we'll extract these in the bash
#' script for HyPhy.
#' We can't simply use the protein sequences because, according to
#' a HyPhy developer:
#'   > In order for our tools to run on sequence data, the sequences must be
#'   > in the format of a codon-aware sequence alignment. Therefore the
#'   > sequences must be divisible by three, in codons.
#'   source: https://github.com/veg/hyphy/issues/1174#issuecomment-648068105



all_hogs <- focal_go_df$hogs |>
    do.call(what = c) |>
    unique()

hog_genes <- dirs$orthofinder_extr |>
    paste0("/Single_Copy_HOG_GO/N0-GO-by-species-genes.tsv") |>
    read_tsv(col_types = cols()) |>
    select(-go) |>
    filter(hog %in% all_hogs) |>
    arrange(species, gene, hog) |>
    select(species, gene, hog)


if (!file.exists("_data/hyphy-hog-genes.csv") || .overwrite) {
    write_csv(hog_genes, "_data/hyphy-hog-genes.csv")
}
