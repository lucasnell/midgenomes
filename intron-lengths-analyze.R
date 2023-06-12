

library(tidyverse)
library(ggridges)
library(ape)
library(phylolm)
library(ggtree)
library(treeio)

theme_set(theme_minimal())


spp_df  <- tibble(species = c("Aaegyp", "Asteph", "Bantar", "Cquinq", "Cripar",
                              "Csonor", "Ctenta", "Mdomes", "Pakamu", "Ppemba",
                              "Pstein", "Pvande", "Tgraci"),
                  family = c("cu", "cu", "ch", "cu", "ch",
                             "ce", "ch", "mu", "ch", "ch",
                             "ch", "ch", "ch") |>
                      factor(levels = c("ch", "ce", "cu", "mu"),
                             labels = c("Chironomidae", "Ceratopogonidae",
                                        "Culicidae", "Muscidae"))) |>
    arrange(family, species) |>
    mutate(species = factor(species, levels = species))


intron_df <- map_dfr(spp_df$species, \(.spp) {
    sprintf("_data/introns-%s.csv.xz", .spp) |>
        read_csv(col_types = cols(), progress = FALSE) |>
        mutate(species = .spp, family = spp_df$family[spp_df$species == .spp],
               length = as.integer(end - start + 1)) |>
        select(family, species, gene, id, length)
}) |>
    mutate(loglen = log10(length))



# ===========================================================================*
# ===========================================================================*
#               PLOTS
# ===========================================================================*
# ===========================================================================*



intron_df |>
    ggplot(aes(species, loglen, color = family, fill = family)) +
    geom_violin() +
    # geom_jitter(alpha = 0.1) +
    stat_summary(fun = "median", colour = "red", size = 2, geom = "point") +
    # stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
    scale_y_continuous("Intron length (bp)", breaks = 0:3*2,
                       labels = parse(text = paste0("10^", 0:3*2))) +
    scale_color_viridis_d(aesthetics = c("color", "fill")) +
    coord_flip()





intron_df |>
    ggplot(aes(x = loglen, y = species, fill = family)) +
    geom_density_ridges(size = 0.25) +
    stat_summary(fun = "mean", colour = "black", size = 2, geom = "point") +
    scale_x_continuous("Intron length (bp)", breaks = 0:3*2,
                       labels = parse(text = paste0("10^", 0:3*2))) +
    scale_color_brewer(NULL, aesthetics = c("color", "fill"), palette = "Dark2") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = "black"))




#'
#' Because it looks like some of these are bimodal,
#' I also tried fitting Gaussian mixture models using
#' `mclust::densityMclust(x$loglen, verbose = FALSE, plot = FALSE)`,
#' but this did not yield useful results.
#' Because the data highly skewed due to a threshold of length >= 1,
#' the mixture models tended to overestimate the number of modes.
#' This resulted in all species having a similar number of modes.
#'

n_intron_df <- intron_df |>
    group_by(family, species, gene) |>
    summarize(n_introns = n(), .groups = "drop")

spp_pal <- viridisLite::turbo(100, begin = 0, end = 1)[rev(c(45, 15+4*0:2, 60, 70+4*0:7))]


#' Showing entire distributions with high values (> 15) shown with points.
n_intron_df |>
    ggplot(aes(x = n_introns, y = species, fill = species)) +
    geom_density_ridges2(size = 0.5, color = "white",
                         quantile_lines = TRUE,
                         quantile_fun = function(x,...) mean(x)) +
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
    geom_density_ridges2(size = 0.5, color = "white",
                         quantile_lines = TRUE,
                         quantile_fun = function(x,...) mean(x)) +
    scale_x_continuous("Introns per gene") +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = "black")) +
    coord_cartesian(xlim = c(NA, 15))



sum_intron_df <- intron_df |>
    group_by(family, species, gene) |>
    summarize(sum_introns = sum(length), .groups = "drop")

#' Showing entire distributions with high values (> 15) shown with points.
sum_intron_df |>
    ggplot(aes(x = log10(sum_introns), y = species, fill = species)) +
    geom_density_ridges2(size = 0.5, color = "white",
                         quantile_lines = TRUE,
                         quantile_fun = function(x,...) mean(x)) +
    # geom_jitter(data = n_intron_df |> filter(n_introns > 15),
    #             aes(color = species), size = 0.5, alpha = 0.5,
    #             width = 0, height = 0.2) +
    scale_x_continuous("Total intron length per gene") +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(color = "black"))



# ===========================================================================*
# ===========================================================================*
#               ANALYSES
# ===========================================================================*
# ===========================================================================*





#' data frame of genome sizes:
gsize_df <- "species	gsize
Tgraci	91827299
Cripar	191837449
Ctenta	213462749
Pvande	118968711
Ppemba	122915671
Bantar	89583723
Pakamu	85836937
Pstein	143569098
Csonor	155941109
Asteph	243475713
Aaegyp	1278732104
Cquinq	573230032
Mdomes	750403944" |>
    read_tsv(col_type = cols()) |>
    mutate(species = factor(species, levels = levels(intron_df$species)))

# Average intron length per species:
intron_len_summ_df <- intron_df |>
    group_by(family, species) |>
    summarize(length = mean(length),
              loglen = mean(loglen),
              .groups = "drop") |>
    mutate(chir = ifelse(family == "Chironomidae", "yes", "no") |>
               factor(levels = c("no", "yes")),
           chir_cerat = ifelse(family %in% c("Chironomidae", "Ceratopogonidae"),
                               "yes", "no") |>
               factor(levels = c("no", "yes"))) |>
    left_join(gsize_df, by = "species") |>
    arrange(species) |>
    as.data.frame()

# Average number of introns per gene per species:
n_intron_summ_df <- intron_df |>
    group_by(family, species, gene) |>
    summarize(n_introns = n(),
              log_n_introns = log10(n()),
              .groups = "drop") |>
    group_by(family, species) |>
    summarize(n_introns = mean(n_introns),
              log_n_introns = mean(log_n_introns),
              .groups = "drop") |>
    mutate(chir = ifelse(family == "Chironomidae", "yes", "no") |>
               factor(levels = c("no", "yes")),
           chir_cerat = ifelse(family %in% c("Chironomidae", "Ceratopogonidae"),
                               "yes", "no") |>
               factor(levels = c("no", "yes"))) |>
    left_join(gsize_df, by = "species") |>
    arrange(species) |>
    as.data.frame()

rownames(intron_len_summ_df) <- paste(intron_len_summ_df$species)
rownames(n_intron_summ_df) <- paste(n_intron_summ_df$species)


dip_tr <- read.tree("~/_data/_phylo/chir_mcmctree.nwk") |>
    drop.tip("Cmarin")

dip_tr |> plot(no.margin = TRUE)

phylolm(loglen ~ log10(gsize) * chir_cerat, intron_len_summ_df, dip_tr) |>
    summary()

intron_len_summ_df |>
    ggplot(aes(log10(gsize), loglen)) +
    geom_point(aes(color = species)) +
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"))

n_intron_summ_df |>
    filter(species != "Pstein") |>
    ggplot(aes(log10(gsize), log_n_introns)) +
    geom_point(aes(color = species)) +
    stat_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    scale_fill_manual(NULL, values = spp_pal[-6], aesthetics = c("color", "fill"))


# plot(n_introns ~ gsize, n_intron_summ_df)



intron_lens <- intron_len_summ_df[dip_tr$tip.label, "length"]
n_introns <- n_intron_summ_df[dip_tr$tip.label, "n_introns"]
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



