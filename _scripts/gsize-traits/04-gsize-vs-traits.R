

source("_scripts/00-preamble.R")

library(ggtree)
library(treeio)
library(patchwork)



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
    p <- feature_df |>
        ggplot(aes(log_gsize, .data[[.y_var]] - .y_subtr)) +
        geom_point(aes(color = species), size = 2, na.rm = TRUE) +
        y_scale +
        scale_x_continuous("Genome size (Mb)",
                           breaks = log10(100e6 * 2^(0:3)),
                           labels = 100 * 2^(0:3)) +
        scale_fill_manual(NULL, values = full_spp_pal,
                          aesthetics = c("color", "fill"), guide = "none")
    if (.y_var != "log_non_TE") {
        p <- p + theme(axis.title.x = element_blank())
    } else p <- p +
        theme(axis.title.x = element_text(size = 14, vjust = 0,
                                          margin = margin(0,0,0,t=1,"lines")))
    if (! .y_var %in% c("log_DNA", "log_non_TE", "log_Unclassified")) {
        p <- p + theme(axis.text.x = element_blank())
    }
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
do.call(wrap_plots, c(list(nrow = 3), gsize_corr_ps))


for (x in c("log_SINE", "log_LINE", "log_LTR", "log_DNA", "log_non_TE", "log_Unclassified")){
    cat(sprintf("%s = %.3f -- %.3f\n", x,
                round(10^min(feature_df[[x]]) / 1e6, 3),
                round(10^max(feature_df[[x]]) / 1e6, 3)))
}
# log_SINE = 0.011 -- 9.848
# log_LINE = 0.351 -- 246.300
# log_LTR = 0.479 -- 294.027
# log_DNA = 1.083 -- 399.441
# log_non_TE = 0.365 -- 52.163
# log_Unclassified = 0.174 -- 45.254





gsize_traits_panels_p <- feature_df |>
    pivot_longer(- c(species, spp_abbrev, log_gsize)) |>
    mutate(name = pretty$convert(name, to_fct = TRUE, units = TRUE)) |>
    ggplot(aes(log_gsize, value)) +
    geom_point(aes(color = species), size = 2, na.rm = TRUE) +
    ylab(expression({"log"}[10]~"[ Feature ]")) +
    scale_x_continuous("Genome size (Mb)",
                       breaks = log10(100e6 * 2^(0:3)),
                       labels = 100 * 2^(0:3)) +
    scale_fill_manual(NULL, values = full_spp_pal,
                      guide = "none",
                      # guide = guide_legend(label.theme = element_text(face = "italic")),
                      aesthetics = c("color", "fill")) +
    facet_wrap(~ name, scales = "free_y", ncol = 3) +
    theme(strip.text = element_text(size = 7), strip.clip = "off")



#' ========================================================================
#' ========================================================================
# inset tree ----
#' ========================================================================
#' ========================================================================

tree_inset <- read.tree("_data/phylo/time-tree.nwk") |>
    (\(x) {
        x$tip.label <- expand_spp(x$tip.label)
        return(x)
    })() |>
    ggtree() |>
    as.treedata() |>
    mutate(species = factor(label, levels = levels(feature_df$species))) |>
    ggtree() +
    geom_rootedge(0.1) +
    # geom_tippoint(aes(color = species), size = 3) +
    geom_tiplab(aes(color = species), size = 7 / 2.83465, fontface = "italic") +
    scale_color_manual(NULL, values = full_spp_pal, guide = "none") +
    theme_inset() +
    coord_cartesian(xlim = c(0, 6))
# tree_inset


#' ========================================================================
#' ========================================================================
# combine and save ----
#' ========================================================================
#' ========================================================================


gsize_traits_p <- tree_inset + gsize_traits_panels_p +
    plot_layout(nrow = 1, widths = c(1, 1.6))


save_plot("gsize-traits", gsize_traits_p, 6.5, 4)





#' Not sure if this is necessary. It's a larger dataset just on number of
#' proteins vs genome size for a bunch of Dipteran species.
#' I collected these two variables for all dipterans with assemblies on
#' InsectBase (on 29 Sep 2023), except for those in genera with many assemblies,
#' at which point I randomly chose 5 species to represent that genus.
#' I don't have the phylogenetic info to analyze this properly, but it shows
#' that compared to a wide range of dipterans, chironomids have an especially
#' steep slope between genome size and number of proteins.
#'
# gs_df <- "~/Stanford_Drive/UW/midgenomes/midgenomes.xlsx" |>
#     readxl::read_excel(sheet = "gsize-proteins") |>
#     filter(accession != "x") |>
#     mutate(chir_cerat = family %in% c("Chironomidae", "Ceratopogonidae"))
#
# gs_df |>
#     ggplot(aes(log10(gsize), log10(n_genes))) +
#     geom_point(aes(color = chir_cerat)) +
#     scale_y_continuous("Protein-coding genes (K)") +
#     scale_x_continuous("Genome size (Mb)",
#                        breaks = log10(100e6 * 2^(0:3)),
#                        labels = 100 * 2^(0:3)) +
#     scale_color_manual(values = c("gray70", "red"))
