
#'
#' Analyze the following:
#'
#' 1. Do genome traits differ for families Chironomidae and Ceratopogonidae
#'    versus all other dipterans in the phylogeny?
#' 2. Do genome traits significantly correlate with genome size when accounting
#'    for phylogenetic autocorrelation?
#'

source("_R/00-preamble.R")

library(phylolm)
library(phyr)
library(future)
library(future.apply)

# For phylolm bootstrapping:
plan(multisession)




#' =================================================================
#' =================================================================
#  Read genome stats ----
#' =================================================================
#' =================================================================



gstat_df <- read_csv("_data/genome-stats.csv", col_types = cols()) |>
    # ordering factors so their levels are in order they show up
    mutate(family = factor(family, levels = unique(family)),
           species = factor(species, levels = species),
           spp_abbrev = factor(spp_abbrev, levels = spp_abbrev)) |>
    # filter(spp_abbrev != "Cmarin") |>
    mutate(non_TE = rowSums(across(all_of(nonTE_classes))),
           log_sum_interg_len = log10(sum_interg_len),
           log_gsize = log10(gsize),
           log_n_genes = log10(n_genes),
           chir = family == "Chironomidae",
           chir_cerat = family %in% c("Chironomidae", "Ceratopogonidae")) |>
    select(-all_of(nonTE_classes)) |>
    mutate(across(DNA:non_TE, \(x) log10(x))) |>
    rename_with(\(x) paste0("log_", x), DNA:non_TE) |>
    select(spp_abbrev, chir, chir_cerat, log_gsize, log_n_genes, log_sum_interg_len,
           mean_log_intron_len, log_DNA:log_non_TE) |>
    as.data.frame() |>
    (\(.df) {
        rownames(.df) <- paste(.df[["spp_abbrev"]])
        return(.df)
    })()

# Variables of interest:
y_vars <- colnames(gstat_df)[-1:-3]


#' =================================================================
#' =================================================================
#  Read phylogeny ----
#' =================================================================
#' =================================================================


#' Ives and Garland (2010; https://doi.org/10.1093/sysbio/syp074)
#' recommend scaling phylogeny to max depth of 1 for easier interpretation
#' (because alpha scales with phylogeny depth):
rescale_phy <- function(x) {
    x$edge.length <- x$edge.length / max(node.depth.edgelength(x))
    return(x)
}

dip_tr <- read.tree("_data/phylo/time-tree.nwk") |>
    reorder("pruningwise") |>
    rescale_phy()

# max(node.depth.edgelength(dip_tr))
# [1] 1

# dip_tr |> plot(no.margin = TRUE)




#' ===========================================================================
#' ===========================================================================
#  Regressions - family versus... ----
#' ===========================================================================
#' ===========================================================================

#'
#' This section tests whether traits differ between families
#' Chironomidae and Ceratopogonidae, and all others.
#'
#' I first ran each model once (no bootstraps) using both OUfixedRoot
#' and Brownian Motion (BM) models. In the versions below, I
#' used the model with lower AIC.
#'



# Function to do one fit to phylom, optionally with bootstrapping:
one_phylolm <- function(y, x = "chir_cerat", .boots = 2000) {
    .f <- as.formula(paste(y, "~", x))
    .d <- gstat_df
    .t <- dip_tr
    if (grepl("intron", y)) {
        .d <- .d |> filter(spp_abbrev != "Cmarin")
        .t <- .t |>
            drop.tip("Cmarin") |>
            reorder("pruningwise") |>
            rescale_phy()
    }
    .m <- phylolm(.f, .d, .t, model = "OUfixedRoot", boot = .boots)
    .m$call[[2]] <- eval(.f)
    return(.m)
}

# Using chir_cerat always improves the model:
y_vars |>
    map_dfr(\(.y) {
        LL1 <- one_phylolm(.y, x = "chir_cerat", .boots = 0) |>
            logLik() |> getElement(1)
        LL2 <- one_phylolm(.y, x = "chir", .boots = 0) |>
            logLik() |> getElement(1)
        tibble(y = .y, LL_cerat = LL1, LL_chir = LL2,
               cerat_better = LL1 > LL2)
    })




if (!file.exists("_data/family-regressions.rds")) {
    # Takes ~20 sec
    set.seed(126458689)
    regr_df <- tibble(feature = y_vars) |>
        mutate(pl_model = map(feature, one_phylolm)) |>
        mutate(coef = map_dbl(pl_model, \(x) coef(x)[["chir_ceratTRUE"]]),
               coef_lo = map_dbl(pl_model,
                                 \(x) x$bootconfint95[1,"chir_ceratTRUE"]),
               coef_hi = map_dbl(pl_model,
                                 \(x) x$bootconfint95[2,"chir_ceratTRUE"])) |>
        select(feature, starts_with("coef"), pl_model)
    write_rds(regr_df, "_data/family-regressions.rds", compress = "gz")
} else {
    regr_df <- read_rds("_data/family-regressions.rds")
}




regr_df |>
    select(-pl_model) |>
    mutate(signif = ifelse((coef_lo < 0 & coef_hi < 0) |
                               (coef_lo > 0 & coef_hi > 0), "**", "\u00A0"))

# # A tibble: 10 × 5
#    feature                coef coef_lo coef_hi signif
#    <chr>                 <dbl>   <dbl>   <dbl> <chr>
#  1 log_gsize           -0.700  -0.960   -0.437 **
#  2 log_n_genes          0.0649 -0.0143   0.141
#  3 log_sum_interg_len  -0.636  -0.958   -0.314 **
#  4 mean_log_intron_len -0.512  -0.681   -0.339 **
#  5 log_DNA             -1.34   -2.16    -0.548 **
#  6 log_LINE            -1.57   -2.07    -1.07  **
#  7 log_LTR             -1.64   -2.05    -1.22  **
#  8 log_SINE            -1.57   -2.23    -0.878 **
#  9 log_Unclassified    -1.49   -1.99    -0.997 **
# 10 log_non_TE          -1.15   -1.65    -0.647 **




#' ===========================================================================
#' ===========================================================================
#  Correlations - gsize vs ... ----
#' ===========================================================================
#' ===========================================================================



one_cor_phylo <- function(y) {
    # rm(.boots, .f,  .d, .t, y, .cd, m)
    .boots <- 2000
    .f <- as.formula(paste0("~ log_gsize + ", y))
    .d <- gstat_df
    .t <- dip_tr
    if (grepl("intron", y)) {
        .d <- .d |> filter(spp_abbrev != "Cmarin")
        .t <- .t |>
            drop.tip("Cmarin") |>
            reorder("pruningwise") |>
            rescale_phy()
    }
    .cd <- y == "log_SINE"
    m <- cor_phylo(.f, ~spp_abbrev, .t, data = .d, method = "bobyqa",
                   boot = .boots, constrain_d = .cd)
    m$call[[2]] <- eval(.f)
    m$call$constrain_d <- eval(.cd)
    m$call$boot <- eval(.boots)
    # cat(y, "\n")
    return(m)
}


if (!file.exists("_data/gsize-corrs.rds")) {
    # Takes ~50 sec
    set.seed(328382065)
    corr_df <- tibble(feature = y_vars[-1]) |>
        mutate(cp_model = map(feature, one_cor_phylo)) |>
        mutate(corr = map_dbl(cp_model, \(x) x[["corrs"]][1,2]),
               corr_lo = map_dbl(cp_model, \(x) boot_ci(x)[["corrs"]][2,1]),
               corr_hi = map_dbl(cp_model, \(x) boot_ci(x)[["corrs"]][1,2])) |>
        select(feature, starts_with("corr"), cp_model)
    write_rds(corr_df, "_data/gsize-corrs.rds", compress = "gz")
} else {
    corr_df <- read_rds("_data/gsize-corrs.rds")
}


corr_df |>
    select(-cp_model) |>
    mutate(signif = ifelse((corr_lo < 0 & corr_hi < 0) |
                               (corr_lo > 0 & corr_hi > 0), "**", "\u00A0"))

# # A tibble: 9 × 5
#   feature              corr corr_lo corr_hi signif
#   <chr>               <dbl>   <dbl>   <dbl> <chr>
# 1 log_n_genes         0.218  -0.322   0.627
# 2 log_sum_interg_len  0.952   0.639   1.00  **
# 3 mean_log_intron_len 0.866   0.500   0.999 **
# 4 log_DNA             0.789   0.431   0.951 **
# 5 log_LINE            0.676   0.240   0.908 **
# 6 log_LTR             0.727   0.312   0.922 **
# 7 log_SINE            0.721   0.271   0.969 **
# 8 log_Unclassified    0.783   0.420   0.949 **
# 9 log_non_TE          0.710   0.270   0.972 **
