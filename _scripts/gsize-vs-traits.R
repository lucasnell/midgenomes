

source("_scripts/00-preamble.R")

library(ggtree)
library(treeio)
library(patchwork)





#' ===========================================================================
#' ===========================================================================
#  Plots - gsize vs ----
#' ===========================================================================
#' ===========================================================================


feature_df |>
    ggplot(aes(log_gsize, log_intron_len)) +
    geom_point(aes(color = species)) +
    scale_y_continuous("Mean intron length (bp)",
                       breaks = log10(120 * 2^(0:3)),
                       labels = 120 * 2^(0:3)) +
    scale_x_continuous("Genome size (Mb)",
                       breaks = log10(100e6 * 2^(0:3)),
                       labels = 100 * 2^(0:3)) +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none")

feature_df |>
    ggplot(aes(log_gsize, n_introns)) +
    geom_point(aes(color = species)) +
    scale_y_continuous("Mean introns per gene") +
    scale_x_continuous("Genome size (Mb)",
                       breaks = log10(100e6 * 2^(0:3)),
                       labels = 100 * 2^(0:3)) +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none")


feature_df |>
    ggplot(aes(log_gsize, n_prots / 1e3)) +
    geom_point(aes(color = species)) +
    scale_y_continuous("Protein-coding genes (K)") +
    scale_x_continuous("Genome size (Mb)",
                       breaks = log10(100e6 * 2^(0:3)),
                       labels = 100 * 2^(0:3)) +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none")

feature_df |>
    # ggplot(aes(log_gsize, log_sum_interg_len - 6)) +
    ggplot(aes(log_gsize, sum_interg_len / gsize)) +
    geom_point(aes(color = species)) +
    scale_y_continuous("Total intergenic DNA (Mb)") +
    scale_x_continuous("Genome size (Mb)",
                       breaks = log10(100e6 * 2^(0:3)),
                       labels = 100 * 2^(0:3)) +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none")




rep_df |>
    filter(class %in% TE_classes) |>
    group_by(species) |>
    summarize(across(elements:prop, sum)) |>
    mutate(gsize = map_dbl(species, \(x) feature_df$gsize[feature_df$species == x]),
           log_gsize = map_dbl(species, \(x) feature_df$log_gsize[feature_df$species == x])) |>
    ggplot(aes(log_gsize, (prop))) +
    geom_point(aes(color = species)) +
    # scale_y_continuous("") +
    scale_x_continuous("Genome size (Mb)",
                       breaks = log10(100e6 * 2^(0:3)),
                       labels = 100 * 2^(0:3)) +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none")






gs_df <- "~/Stanford_Drive/UW/midgenomes/midgenomes.xlsx" |>
    readxl::read_excel(sheet = "gsize-proteins") |>
    filter(accession != "x") |>
    mutate(chir_cerat = family %in% c("Chironomidae", "Ceratopogonidae"))


gs_df |>
    ggplot(aes(log10(gsize), log10(n_genes))) +
    geom_point(aes(color = chir_cerat)) +
    scale_y_continuous("Protein-coding genes (K)") +
    scale_x_continuous("Genome size (Mb)",
                       breaks = log10(100e6 * 2^(0:3)),
                       labels = 100 * 2^(0:3)) +
    scale_color_manual(values = c("gray70", "red"))
