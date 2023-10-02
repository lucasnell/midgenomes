

library(tidyverse)
library(ggridges)
library(ape)
library(phylolm)
library(ggtree)
library(treeio)

theme_set(theme_minimal())


#' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#'
#' CAN SIMPLIFY BELOW BY USING `_data/genome-stats.csv`
#'
#' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

spp_df  <- tibble(species = c("Aaegyp", "Asteph", "Bantar", "Cquinq", "Cripar",
                              "Csonor", "Ctenta", "Mdomes", "Pakamu", "Ppemba",
                              "Pstein", "Pvande", "Tgraci"),
                  family = c("cu", "cu", "ch", "cu", "ch",
                             "ce", "ch", "mu", "ch", "ch",
                             "ch", "ch", "ch") |>
                      factor(levels = c("ch", "ce", "cu", "mu"),
                             labels = c("Chironomidae", "Ceratopogonidae",
                                        "Culicidae", "Muscidae"))) |>
    mutate(species = factor(species,
                            levels = c("Ctenta", "Cripar", "Pvande", "Ppemba",
                                       "Tgraci", "Bantar", "Cmarin", "Pakamu",
                                       "Pstein", "Csonor", "Cquinq", "Aaegyp",
                                       "Asteph", "Mdomes")))

spp_pal <- turbo(100)[c(70+4*0:7, 60, 15+4*2:0, 30)]


intron_df <- map_dfr(spp_df$species, \(.spp) {
    sprintf("_data/introns/%s.csv.xz", .spp) |>
        read_csv(col_types = cols(), progress = FALSE) |>
        mutate(species = .spp, family = spp_df$family[spp_df$species == .spp],
               intron_len = as.integer(end - start + 1)) |>
        select(family, species, gene, id, intron_len)
}) |>
    mutate(log_intron_len = log10(intron_len))




# ===========================================================================*
# ===========================================================================*
#               PLOTS
# ===========================================================================*
# ===========================================================================*



intron_df |>
    ggplot(aes(x = log_intron_len, y = species, fill = species)) +
    geom_density_ridges2(size = 0.25, vline_size = 1,
                         quantile_lines = TRUE, quantile_fun = mean) +
    scale_x_continuous("Intron length (bp)", breaks = 0:3*2,
                       labels = parse(text = paste0("10^", 0:3*2))) +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill")) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = "black"),
          legend.position = "none")




#'
#' Because it looks like some of these are bimodal,
#' I also tried fitting Gaussian mixture models using
#' `mclust::densityMclust(x$loglen, verbose = FALSE, plot = FALSE)`,
#' but this did not yield useful results.
#' Because the data are highly skewed due to a threshold of length >= 1,
#' the mixture models tended to overestimate the number of modes.
#' This resulted in all species having a similar number of modes.
#'

n_intron_df <- intron_df |>
    group_by(family, species, gene) |>
    summarize(n_introns = n(), .groups = "drop")



#' Showing entire distributions with high values (> 15) shown with points.
n_intron_df |>
    ggplot(aes(x = n_introns, y = species, fill = species)) +
    geom_density_ridges2(size = 0.25, vline_size = 1, quantile_lines = TRUE,
                         quantile_fun = mean) +
    geom_jitter(data = n_intron_df |> filter(n_introns > 15),
                aes(color = species), size = 0.5, alpha = 0.5,
                width = 0, height = 0.2) +
    scale_x_continuous("Introns per gene") +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = "black"))

#' Showing distributions without tails:
n_intron_df |>
    ggplot(aes(x = n_introns, y = species, fill = species)) +
    geom_density_ridges2(size = 0.25, vline_size = 1, quantile_lines = TRUE,
                         quantile_fun = mean) +
    scale_x_continuous("Introns per gene") +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = "black")) +
    coord_cartesian(xlim = c(NA, 15))



sum_intron_df <- intron_df |>
    group_by(family, species, gene) |>
    summarize(sum_introns = sum(intron_len), .groups = "drop")

#' Showing entire distributions with high values (> 15) shown with points.
sum_intron_df |>
    ggplot(aes(x = log10(sum_introns), y = species, fill = species)) +
    geom_density_ridges2(size = 0.25, vline_size = 1, quantile_lines = TRUE,
                         quantile_fun = mean) +
    scale_x_continuous("Total intron length per gene (bp)",
                       breaks = 0:3 * 2,
                       labels = parse(text = paste0("10^", 0:3 * 2))) +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = "black"))



