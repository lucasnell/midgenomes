
library(GO.db)
# ^^ make sure this is loaded before tidyverse ^^
library(tidyverse)
library(parallel)
options("mc.cores" = max(1L, parallel::detectCores()-2L))


hog_gos <- "~/_data/chir_orthofinder/orthofinder-output/Single_Copy_HOG_GO/N0-GO-by-HOG.tsv" |>
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
#' Within response to stimulus > response to stress >
#' cellular response to stress:
#'
#'   * GO:0006281 - DNA repair (in DNA damage response)
#'   * GO:0042262 - DNA protection
#'
#'
#' Within response to stimulus > response to stress:
#'
#'   * GO:0034059 - response to anoxia
#'   * GO:0061771 - response to caloric restriction
#'   * GO:0009409 - response to cold
#'   * GO:0009408 - response to heat
#'   * GO:0090664 - response to high population density
#'   * GO:0001666 - response to hypoxia
#'   * GO:0006979 - response to oxidative stress
#'   * GO:0042594 - response to starvation




#' Dataframe of GO terms that I'm interested in analyzing.
#' Start with terms listed above, then I got all offspring GO terms
#' (recursively) for these terms.
#'

focal_go_df <- tibble(go = c("GO:0010038", "GO:0006281", "GO:0042262",
                             "GO:0010212", "GO:0034059", "GO:0061771",
                             "GO:0009409", "GO:0009408", "GO:0090664",
                             "GO:0001666", "GO:0006979", "GO:0042594"),
                      term = get_term(go),
                      # parent GO term plus offspring (recursive)
                      offspring = mclapply(go, get_offs),
                      hogs = map(offspring, \(o) hog_gos$hog[hog_gos$go %in% o]))
focal_go_df

#' Some have zero HOGs associated with them, so we can filter these out:
focal_go_df <- filter(focal_go_df, map_int(hogs, length) > 0)


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
shared_mat
#' Heatmap of shared HOGs:
heatmap(shared_mat, Rowv = NA, Colv = NA, symm = FALSE)


#' Extract sequences for all HOGs:

all_hogs <- do.call(c, focal_go_df$hogs)
hog_in <- \(x) {
    paste0("~/_data/chir_orthofinder/orthofinder-output/",
           "Single_Copy_HOG_Sequences/N0/", x, ".faa")
}

for (h in all_hogs[1:5]) {
    hog_in(h)
}

