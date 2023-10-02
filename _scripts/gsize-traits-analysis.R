

source("_scripts/00-preamble.R")

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
    mutate(log_sum_interg_len = log10(sum_interg_len),
           log_gsize = log10(gsize),
           log_n_genes = log10(n_genes),
           chir = family == "Chironomidae",
           chir_cerat = family %in% c("Chironomidae", "Ceratopogonidae")) |>
    mutate(across(DNA:Small_RNA, \(x) log10(x))) |>
    rename_with(\(x) paste0("log_", x), DNA:Small_RNA) |>
    select(spp_abbrev, chir_cerat, log_gsize, log_n_genes, log_sum_interg_len,
           mean_log_intron_len, log_DNA:log_Small_RNA) |>
    as.data.frame() |>
    (\(.df) {
        rownames(.df) <- paste(.df[["spp_abbrev"]])
        return(.df)
    })()



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

dip_tr <- read.tree("~/_data/_phylo/chir_mcmctree.nwk") |>
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


if (!file.exists("_data/family-regressions.rds")) {
    # Takes ~20 sec
    set.seed(126458689)
    regr_df <- tibble(feature = colnames(gstat_df)[-1:-2]) |>
        mutate(pl_model = map(feature, \(x) {
            .boots <- 2000
            .f <- as.formula(paste(x, "~ chir_cerat"))
            .d <- gstat_df
            .t <- dip_tr
            if (grepl("intron", x)) {
                .d <- .d |> filter(spp_abbrev != "Cmarin")
                .t <- .t |>
                    drop.tip("Cmarin") |>
                    reorder("pruningwise") |>
                    rescale_phy()
            }
            if (x %in% c("log_Low_complexity", "log_RC", "log_Simple_repeat",
                         "log_Small_RNA")) {
                .m <- phylolm(.f, .d, .t, model = "BM", boot = .boots)
            } else if (x %in% c("mean_log_intron_len", "mean_log_n_introns")) {
                .m <- phylolm(.f, .d, .t, model = "OUfixedRoot",
                              upper.bound = 200, boot = .boots)
            } else {
                .m <- phylolm(.f, .d, .t, model = "OUfixedRoot", boot = .boots)
            }
            .m$call[[2]] <- eval(.f)
            return(.m)
        })) |>
        mutate(coef = map_dbl(pl_model, \(x) coef(x)[["chir_ceratTRUE"]]),
               coef_lo = map_dbl(pl_model,
                                 \(x) x$bootconfint95[1,"chir_ceratTRUE"]),
               coef_hi = map_dbl(pl_model,
                                 \(x) x$bootconfint95[2,"chir_ceratTRUE"])) |>
        select(feature, starts_with("coef"), pl_model)
    write_rds(regr_df, "_data/family-regressions.rds", compress = "xz")
} else {
    regr_df <- read_rds("_data/family-regressions.rds")
}



regr_df |>
    select(-pl_model) |>
    mutate(signif = ifelse((coef_lo < 0 & coef_hi < 0) |
                               (coef_lo > 0 & coef_hi > 0), "**", "\u00A0"))

# # A tibble: 13 × 5
#    feature                coef coef_lo coef_hi signif
#    <chr>                 <dbl>   <dbl>   <dbl> <chr>
#  1 log_gsize           -0.698  -0.934  -0.466  **
#  2 log_n_genes          0.0668 -0.0104  0.148
#  3 log_sum_interg_len  -0.648  -0.949  -0.367  **
#  4 mean_log_intron_len -0.502  -0.665  -0.339  **
#  5 log_DNA             -1.32   -2.15   -0.536  **
#  6 log_LINE            -1.56   -2.04   -1.05   **
#  7 log_LTR             -1.64   -2.04   -1.19   **
#  8 log_Low_complexity  -0.182  -0.825   0.448
#  9 log_RC              -1.53   -4.03    0.989
# 10 log_SINE            -1.56   -2.24   -0.881  **
# 11 log_Satellite       -1.21   -2.48    0.144
# 12 log_Simple_repeat   -0.684  -1.44    0.0782
# 13 log_Small_RNA       -1.36   -2.62   -0.151  **




#' ===========================================================================
#' ===========================================================================
#  Correlations - gsize vs ... ----
#' ===========================================================================
#' ===========================================================================




if (!file.exists("_data/gsize-corrs.rds")) {
    # Takes ~50 sec
    set.seed(328382065)
    corr_df <- tibble(feature = colnames(gstat_df)[-1:-3]) |>
        mutate(cp_model = map(feature, \(x) {
            .boots <- 2000
            .f <- as.formula(paste0("~ log_gsize + ", x))
            .d <- gstat_df
            .t <- dip_tr
            if (grepl("intron", x)) {
                .d <- .d |> filter(spp_abbrev != "Cmarin")
                .t <- .t |>
                    drop.tip("Cmarin") |>
                    reorder("pruningwise") |>
                    rescale_phy()
            }
            .cd <- x %in% c("log_SINE", "log_Satellite", "log_Small_RNA")
            m <- cor_phylo(.f, ~spp_abbrev, .t, data = .d, method = "bobyqa",
                           boot = .boots, constrain_d = .cd)
            m$call[[2]] <- eval(.f)
            # cat(x, "\n")
            return(m)
        })) |>
        mutate(corr = map_dbl(cp_model, \(x) x[["corrs"]][1,2]),
               corr_lo = map_dbl(cp_model, \(x) boot_ci(x)[["corrs"]][2,1]),
               corr_hi = map_dbl(cp_model, \(x) boot_ci(x)[["corrs"]][1,2])) |>
        select(feature, starts_with("corr"), cp_model)

    write_rds(corr_df, "_data/gsize-corrs.rds", compress = "xz")
} else {
    corr_df <- read_rds("_data/gsize-corrs.rds")
}

corr_df |>
    select(-cp_model) |>
    mutate(signif = ifelse((corr_lo < 0 & corr_hi < 0) |
                               (corr_lo > 0 & corr_hi > 0), "**", "\u00A0"))

# # A tibble: 12 × 5
#    feature              corr corr_lo corr_hi signif
#    <chr>               <dbl>   <dbl>   <dbl> <chr>
#  1 log_n_genes         0.160 -0.373    0.583
#  2 log_sum_interg_len  0.962  0.857    1.00  **
#  3 mean_log_intron_len 0.913  0.574    1.00  **
#  4 log_DNA             0.842  0.526    0.975 **
#  5 log_LINE            0.776  0.412    0.946 **
#  6 log_LTR             0.814  0.482    0.956 **
#  7 log_Low_complexity  0.200 -0.429    0.668
#  8 log_RC              0.786  0.437    0.982 **
#  9 log_SINE            0.769  0.407    0.979 **
# 10 log_Satellite       0.432 -0.108    0.868
# 11 log_Simple_repeat   0.487 -0.0583   0.831
# 12 log_Small_RNA       0.250 -0.299    0.730


