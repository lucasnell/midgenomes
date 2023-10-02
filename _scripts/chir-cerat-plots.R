

source("_scripts/00-preamble.R")

library(patchwork)


#' =================================================================
#' =================================================================
#  Info to make labels pretty ----
#' =================================================================
#' =================================================================

#' Repeat element classes in order with formatted names for plotting:
rep_class_map <- list("SINE" = "SINEs",
                      "LINE" = "LINEs",
                      "LTR" = "LTR elements",
                      "DNA" = "DNA transposons",
                      "RC" = "Rolling-circles",
                      "Small_RNA" = "Small RNA",
                      "Satellite" = "Satellites",
                      "Simple_repeat" = "Simple repeats",
                      "Low_complexity" = "Low complexity",
                      "Unclassified" = "Unclassified")



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
    mutate(log_sum_interg_len = log10(sum_interg_len),
           log_gsize = log10(gsize),
           log_n_genes = log10(n_genes),
           chir = family == "Chironomidae",
           chir_cerat = family %in% c("Chironomidae", "Ceratopogonidae")) |>
    mutate(across(DNA:Small_RNA, \(x) log10(x))) |>
    rename_with(\(x) paste0("log_", x), DNA:Small_RNA) |>
    select(spp_abbrev, chir_cerat, log_gsize, log_n_genes, log_sum_interg_len,
           mean_log_intron_len, log_DNA:log_Small_RNA)

