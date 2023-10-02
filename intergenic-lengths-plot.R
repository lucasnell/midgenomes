
library(tidyverse)
library(ggridges)
library(ape)
library(phylolm)
library(ggtree)
library(treeio)
library(viridisLite)

theme_set(theme_minimal())

#' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#'
#' CAN SIMPLIFY BELOW BY USING `_data/genome-stats.csv`
#'
#' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



spp_df  <- tibble(species = c("Aaegyp", "Asteph", "Bantar", "Cquinq", "Cripar",
                              "Csonor", "Ctenta", "Mdomes", "Pakamu", "Ppemba",
                              "Pstein", "Pvande", "Tgraci", "Cmarin"),
                  family = c("cu", "cu", "ch", "cu", "ch",
                             "ce", "ch", "mu", "ch", "ch",
                             "ch", "ch", "ch", "ch") |>
                      factor(levels = c("ch", "ce", "cu", "mu"),
                             labels = c("Chironomidae", "Ceratopogonidae",
                                        "Culicidae", "Muscidae"))) |>
    mutate(species = factor(species,
                            levels = c("Ctenta", "Cripar", "Pvande", "Ppemba",
                                       "Tgraci", "Bantar", "Cmarin", "Pakamu",
                                       "Pstein", "Csonor", "Cquinq", "Aaegyp",
                                       "Asteph", "Mdomes")))

spp_pal <- turbo(100)[c(70+3*0:8, 60, 15+4*2:0, 30)]
names(spp_pal) <- levels(spp_df$species)

interg_df <- map_dfr(spp_df$species, \(.spp) {
    sprintf("_data/intergenic/%s.csv.xz", .spp) |>
        read_csv(col_types = cols(), progress = FALSE) |>
        mutate(species = .spp, family = spp_df$family[spp_df$species == .spp]) |>
        filter(interg_len > 0) |>
        select(family, species, seqid, interg_len)
}) |>
    mutate(log_interg_len = log10(interg_len))






# ===========================================================================*
# ===========================================================================*
#               PLOTS
# ===========================================================================*
# ===========================================================================*



interg_df |>
    ggplot(aes(x = log_interg_len, y = species, fill = species)) +
    geom_density_ridges2(size = 0.25, vline_size = 1,
                         quantile_lines = TRUE, quantile_fun = mean) +
    scale_x_continuous("Intergenic length (bp)", breaks = 0:3*2,
                       labels = parse(text = paste0("10^", 0:3*2))) +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill")) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = "black"),
          legend.position = "none")


sum_interg_df <- interg_df |>
    group_by(family, species) |>
    summarize(sum_interg_len = sum(interg_len),
              mean_interg_len = mean(interg_len),
              mean_log_interg_len = mean(log_interg_len),
              .groups = "drop")

sum_interg_df |>
    ggplot(aes(x = log10(sum_interg_len / 1e6), y = species, color = species)) +
    geom_point(size = 3) +
    scale_x_continuous("Total intergenic length (Mb)",
                       breaks = log10(50 * 2^(0:3)),
                       labels = 50 * 2^(0:3))+
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = "black"))

