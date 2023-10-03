
#'
#' Create plots of whether families Chironomidae and Ceratopogonidae differ
#' from all others in our phylogeny.
#'
#' File `_data/family-regressions.rds` (created in
#' `_scripts/gsize-traits-analysis.R`) is required for these plots.
#'

source("_scripts/00-preamble.R")

library(phylolm)



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
    mutate(non_TE = rowSums(across(all_of(nonTE_classes))),
           log_sum_interg_len = log10(sum_interg_len),
           log_gsize = log10(gsize),
           log_n_genes = log10(n_genes),
           chir_cerat = family %in% c("Chironomidae", "Ceratopogonidae")) |>
    select(-all_of(nonTE_classes)) |>
    mutate(across(DNA:non_TE, \(x) log10(x))) |>
    rename_with(\(x) paste0("log_", x), DNA:non_TE) |>
    select(spp_abbrev, chir_cerat, log_gsize, log_n_genes, log_sum_interg_len,
           mean_log_intron_len, log_DNA:log_non_TE)




#' =================================================================
#' =================================================================
#  Read CI from phylolm models ----
#' =================================================================
#' =================================================================

regr_ci_df <- read_rds("_data/family-regressions.rds") |>
    select(feature, pl_model) |>
    pmap_dfr(\(feature, pl_model) {
        b0 <- pl_model$coefficients[["(Intercept)"]]
        b1 <- pl_model$coefficients[["chir_ceratTRUE"]]
        others <- quantile(pl_model$bootstrap[,"(Intercept)"], c(0.025, 0.975))
        chir_cerats <- quantile(pl_model$bootstrap[,"(Intercept)"] +
                                    pl_model$bootstrap[,"chir_ceratTRUE"],
                                c(0.025, 0.975))
        tibble(chir_cerat = factor(c(TRUE, FALSE), levels = c(TRUE, FALSE),
                                   labels = c("Chir.+\nCerat.", "others")),
               name = feature,
               value = c(b0+b1, b0),
               lo = c(chir_cerats[["2.5%"]], others[["2.5%"]]),
               hi = c(chir_cerats[["97.5%"]], others[["97.5%"]]))
        }) |>
    mutate(name = pretty$convert(name, to_fct = TRUE, units = TRUE))




#' =================================================================
#' =================================================================
#  Plot ----
#' =================================================================
#' =================================================================

chir_cerat_p <- gstat_df |>
    pivot_longer(- c(spp_abbrev, chir_cerat)) |>
    mutate(name = pretty$convert(name, to_fct = TRUE, units = TRUE),
           chir_cerat = factor(chir_cerat, levels = c(TRUE, FALSE),
                               labels = c("Chir.+\nCerat.", "others"))) |>
    ggplot(aes(chir_cerat, value)) +
    geom_jitter(width = 0.2, height = 0, color = "dodgerblue", na.rm = TRUE) +
    geom_pointrange(data = regr_ci_df, aes(ymin = lo, ymax = hi), shape = 1) +
    ylab(expression({"log"}[10]~"[ Feature ]")) +
    facet_wrap(~ name, scales = "free_y", ncol = 4) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 7, color = "black"),
          strip.text = element_text(size = 7))



save_plot("chir-cerat-traits", chir_cerat_p, 6.5, 4)
