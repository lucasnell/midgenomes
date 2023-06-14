

library(tidyverse)
library(ggridges)
library(ape)
library(phylolm)
library(ggtree)
library(treeio)

theme_set(theme_minimal())
spp_pal <- viridisLite::turbo(100, begin = 0, end = 1)[rev(c(45, 15+4*0:2, 60, 70+4*0:7))]


#' data frame of genome stats:
gstat_df <- read_csv("_data/genome-stats.csv", col_types = cols()) |>
    mutate(family = factor(family,
                           labels = c("Chironomidae", "Ceratopogonidae",
                                      "Culicidae", "Muscidae")),
           species = factor(species,
                            levels = c("Ctenta", "Cripar", "Pvande", "Ppemba",
                                       "Tgraci", "Bantar", "Cmarin", "Pakamu",
                                       "Pstein", "Csonor", "Cquinq", "Aaegyp",
                                       "Asteph", "Mdomes"))) |>
    filter(species != "Cmarin") |>
    mutate(chir = family == "Chironomidae",
           chir_cerat = family %in% c("Chironomidae", "Ceratopogonidae")) |>
    as.data.frame()

rownames(gstat_df) <- paste(gstat_df$species)


dip_tr <- read.tree("~/_data/_phylo/chir_mcmctree.nwk") |>
    drop.tip("Cmarin")

dip_tr |> plot(no.margin = TRUE)

phylolm(loglen ~ log10(gsize) + chir_cerat, gstat_df, dip_tr) |>
    summary()

gstat_df |>
    ggplot(aes(log10(gsize), loglen)) +
    geom_point(aes(color = species)) +
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"))

gstat_df |>
    # filter(species != "Pstein") |>
    ggplot(aes(log10(gsize), log_n_introns)) +
    geom_point(aes(color = species)) +
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    # scale_fill_manual(NULL, values = spp_pal[-6], aesthetics = c("color", "fill"))
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"))


# plot(n_introns ~ gsize, gstat_df)



intron_lens <- gstat_df[dip_tr$tip.label, "length"]
n_introns <- gstat_df[dip_tr$tip.label, "n_introns"]
names(intron_lens) <- names(n_introns) <- dip_tr$tip.label
l_recon <- reconstruct(intron_lens, dip_tr, method = "GLS_OUS", CI = TRUE)
n_recon <- reconstruct(n_introns, dip_tr, method = "GLS_OUS", CI = TRUE)


# dip_tr$node.label <-

p1 <- ggtree(dip_tr) %<+% df_tip_data +
    geom_rootedge(0.04) +
    geom_tiplab(size = 10 / 2.83465, fontface = "italic") +
    geom_nodelab(size = 8 / 2.83465, nudge_x = -0.025, nudge_y = 0.25) +
    theme_tree2()

p2 <- p1 |>
    revts() +
    scale_x_continuous("Million years ago", limits = c(-3.53367, 1.6), expand = c(0,0),
                       breaks = -3:0, labels = 3:0 * 100)



