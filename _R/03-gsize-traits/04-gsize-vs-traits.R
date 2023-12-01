
#'
#' Plot genome size versus other traits that may relate to it
#'
#' File `_data/family-regressions.rds` (created in
#' `_R/03-gsize-traits/02-gsize-traits-analysis.R`)
#' is required for these plots.
#'



source("_R/00-preamble.R")

library(phyr)
library(ggtree)
library(treeio)
library(patchwork)
library(grid)



#' =================================================================
#' =================================================================
#  Read genome stats ----
#' =================================================================
#' =================================================================



feature_df <- read_csv("_data/genome-stats.csv", col_types = cols()) |>
    # ordering factors so their levels are in order they show up
    mutate(family = factor(family, levels = unique(family)),
           species = factor(species, levels = species),
           spp_abbrev = factor(spp_abbrev, levels = spp_abbrev)) |>
    mutate(non_TE = rowSums(across(all_of(nonTE_classes))),
           log_sum_interg_len = log10(sum_interg_len),
           log_gsize = log10(gsize),
           log_n_genes = log10(n_genes)) |>
    select(-all_of(nonTE_classes)) |>
    mutate(across(DNA:non_TE, \(x) log10(x))) |>
    rename_with(\(x) paste0("log_", x), DNA:non_TE) |>
    select(species, spp_abbrev, log_gsize, log_n_genes, log_sum_interg_len,
           mean_log_intron_len, log_DNA:log_non_TE)





#' ===========================================================================
#' ===========================================================================
#  Read correlation models ----
#' ===========================================================================
#' ===========================================================================

corr_list <- read_rds("_data/gsize-corrs.rds") |>
    split(~ feature) |>
    map(\(x) {
        if ((x$corr_lo > 0 && x$corr_hi > 0) ||
            (x$corr_lo < 0 && x$corr_hi < 0)) {
            sig_str <- "*"
        } else sig_str <- ""
        sprintf("r = %.2f%s", x$corr, sig_str)
    })



#' ===========================================================================
#' ===========================================================================
#  Plots - gsize vs ----
#' ===========================================================================
#' ===========================================================================


one_gsize_corr_p <- function(.y_var, .y_subtr, .y_lab, .y_breaks, .y_lims) {
    if (inherits(.y_breaks, "waiver")) {
        y_scale <- scale_y_continuous(.y_lab, limits = .y_lims)
    } else {
        y_scale <- scale_y_continuous(.y_lab, breaks = log10(.y_breaks),
                                      labels = .y_breaks, limits = .y_lims)
    }
    if (grepl("n_genes", .y_var)) {
        .hj <- 1;
        .x_fun <- max
    } else {
        .hj <- 0;
        .x_fun <- min
    }
    if (is.null(.y_lims)) {
        .y_max <- max(feature_df[[.y_var]] - .y_subtr, na.rm = TRUE)
    } else .y_max <- max(.y_lims)
    p <- feature_df |>
        ggplot(aes(log_gsize, .data[[.y_var]] - .y_subtr)) +
        geom_point(aes(color = species), size = 2, na.rm = TRUE) +
        geom_text(data = tibble(log_gsize = .x_fun(feature_df$log_gsize),
                                y = .y_max,
                                lab = corr_list[[.y_var]]),
                  aes(y = y, label = lab), size = 6 / 2.83465,
                  hjust = .hj, vjust = 1) +
        y_scale +
        scale_x_continuous("Genome size (Mb)",
                           breaks = log10(100e6 * 2^(0:3)),
                           labels = 100 * 2^(0:3)) +
        scale_fill_manual(NULL, values = full_spp_pal,
                          aesthetics = c("color", "fill"), guide = "none") +
        theme(axis.title.y = element_text(size = 9),
              axis.text = element_text(size = 8))
    if (! .y_var %in% c("log_DNA", "log_non_TE", "log_Unclassified")) {
        p <- p + theme(axis.text.x = element_blank())
    }
    if ( .y_var != "log_non_TE" ) {
        p <- p + theme(axis.title.x = element_blank())
    } else p <- p + theme(axis.title.x = element_text(size = 11,
                                                      margin = margin(0,0,0,t=9)))
    return(p)
}


gsize_corr_ps <- tribble(~.y_var, ~.y_subtr, ~.y_lab, ~.y_breaks, ~.y_lims,
        #' Below, '\u00D7' is unicode for multiply sign, and had to use lowercase
        #' '\u' instead of '\U' to have no space between sign and "1000"
        "log_n_genes", 3, "Protein-coding\ngenes (\u00D71000)", c(12, 15, 19), NULL,
        "log_sum_interg_len", 6, "Total intergenic\nDNA (Mb)", 50 * 3^(0:2), NULL,
        "mean_log_intron_len", 0, "Mean intron\nlength (bp)", 200 * 2^(0:2), NULL,
        "log_SINE", 6, "SINE (Mb)", 5 * 10^(-2:0), NULL,
        "log_LINE", 6, "LINE (Mb)", 10^(0:2), log10(c(0.3, 400)),
        "log_LTR", 6, "LTR (Mb)", 10^(0:2), log10(c(0.3, 400)),
        "log_DNA", 6, "DNA (Mb)", 10^(0:2), log10(c(0.3, 400)),
        "log_non_TE", 6, "non-TE\nrepeats (Mb)", 5^(0:2), log10(c(0.17, 53)),
        "log_Unclassified", 6, "Unclassified\nrepeats (Mb)", 5^(0:2), log10(c(0.17, 53))) |>
    pmap(one_gsize_corr_p)

gsize_corr_p <- do.call(wrap_plots, c(list(nrow = 3), gsize_corr_ps))




#' ========================================================================
#' ========================================================================
# save plot ----
#' ========================================================================
#' ========================================================================


# save_plot("gsize-corrs", gsize_corr_p, 5.5, 4)





# #' Not sure if this is necessary. It's a larger dataset just on number of
# #' proteins vs genome size for a bunch of Dipteran species.
# #' I collected these two variables for all dipterans with assemblies on
# #' InsectBase (on 29 Sep 2023), except for those in genera with many assemblies,
# #' at which point I randomly chose 5 species to represent that genus.
# #' I don't have the phylogenetic info to analyze this properly, but it shows
# #' that compared to a wide range of dipterans, chironomids have an especially
# #' steep slope between genome size and number of proteins.
# #'
# gs_df <- "~/Stanford_Drive/UW/midgenomes/midgenomes.xlsx" |>
#     readxl::read_excel(sheet = "gsize-proteins") |>
#     filter(accession != "x") |>
#     mutate(chir_cerat = family %in% c("Chironomidae", "Ceratopogonidae"))
#
# gs_df |>
#     arrange(chir_cerat) |>
#     ggplot(aes(log10(gsize), log10(n_genes / 1e3))) +
#     geom_point(aes(color = chir_cerat)) +
#     scale_y_continuous("Protein-coding\ngenes (\u00D71000)",
#                        breaks = log10(12 * 1.5^(0:2)),
#                        labels = 12 * 1.5^(0:2)) +
#     scale_x_continuous("Genome size (Mb)",
#                        breaks = log10(100e6 * 2^(0:3)),
#                        labels = 100 * 2^(0:3)) +
#     scale_color_manual(values = c("gray70", "red"), guide = "none")
