
library(tidyverse)
library(viridisLite)
library(ggridges)
library(ape)
library(phylolm)
library(phyr)
library(ggtree)
library(treeio)
library(patchwork)
library(future)
library(future.apply)

theme_set(theme_classic())
spp_pal <- turbo(100)[c(70+3*0:8, 60, 15+4*2:0, 30)] |> as.list()
names(spp_pal) <- c("Ctenta", "Cripar", "Pvande", "Ppemba",
                    "Tgraci", "Bantar", "Cmarin", "Pakamu",
                    "Pstein", "Csonor", "Cquinq", "Aaegyp",
                    "Asteph", "Mdomes")

# For phylolm bootstrapping:
plan(multisession)

#' To convert back to full species names:
spp_name_map <- list("Aaegyp" = "A. aegypti",
                     "Asteph" = "A. stephensi",
                     "Bantar" = "B. antarctica",
                     "Cripar" = "C. riparius",
                     "Ctenta" = "C. tentans",
                     "Cmarin" = "C. marinus",
                     "Cquinq" = "C. quinquefasciatus",
                     "Csonor" = "C. sonorensis",
                     "Mdomes" = "M. domestica",
                     "Pstein" = "P. steinenii",
                     "Ppemba" = "P. pembai",
                     "Pvande" = "P. vanderplanki",
                     "Pakamu" = "P. akamusi",
                     "Tgraci" = "T. gracilentus")
# So it's the same order as spp_pal:
spp_name_map <- spp_name_map[names(spp_pal)]

full_spp_pal <- spp_pal
names(full_spp_pal) <- map_chr(names(spp_pal), \(z) spp_name_map[[z]])



#' =================================================================
#' =================================================================
#  Read genome stats ----
#' =================================================================
#' =================================================================

gstat_df <- read_csv("_data/genome-stats.csv", col_types = cols()) |>
    mutate(family = factor(family,
                           levels = c("Chironomidae", "Ceratopogonidae",
                                      "Culicidae", "Muscidae")),
           species = factor(species,
                            levels = c("Ctenta", "Cripar", "Pvande", "Ppemba",
                                       "Tgraci", "Bantar", "Cmarin", "Pakamu",
                                       "Pstein", "Csonor", "Cquinq", "Aaegyp",
                                       "Asteph", "Mdomes"))) |>
    filter(species != "Cmarin") |>
    mutate(log_sum_interg_len = log10(sum_interg_len),
           full_spp = map_chr(paste(species), \(z) spp_name_map[[z]]) |>
               factor(levels = unlist(spp_name_map)),
           chir = family == "Chironomidae",
           chir_cerat = family %in% c("Chironomidae", "Ceratopogonidae")) |>
    as.data.frame() |>
    (\(x) { rownames(x) <- paste(x$species); return(x) })()



#' =================================================================
#' =================================================================
#  Read repeat-element data ----
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
                      "Low_complexity" = "Low complexity")

# non-TE element classes:
nonTE_classes <- c("RC", "Small_RNA", "Satellite", "Simple_repeat", "Low_complexity")
# TE element classes:
TE_classes <- names(rep_class_map)[!names(rep_class_map) %in% nonTE_classes]

rep_pal <- as.list(plasma(9, begin = 0.2, end = 0.9)[c(5,1,6,2,7,3,8,4,9)])
names(rep_pal) <- paste(rep_class_map)

rep_df <- read_csv("_data/repeats-summary.csv", col_types = cols()) |>
    filter(class != "Unclassified") |>
    # filter(!class %in% nonTE_classes) |>
    mutate(gsize = map_dbl(species, \(s) gstat_df$gsize[gstat_df$species == s]),
           prop = length / gsize,
           plot_class = map_chr(class, \(x) rep_class_map[[x]]) |>
               factor(levels = unlist(rep_class_map)),
           class = factor(class, levels = names(rep_class_map)),
           species = factor(species, levels = levels(gstat_df$species)),
           full_spp = map_chr(paste(species), \(z) spp_name_map[[z]]) |>
               factor(levels = unlist(spp_name_map)),
           family = map_vec(species, \(s) gstat_df$family[gstat_df$species == s]),
           chir = family == "Chironomidae",
           chir_cerat = family %in% c("Chironomidae", "Ceratopogonidae")) |>
    select(-gsize)



rep_elem_df <- rep_df |>
    select(-length, -prop, -plot_class) |>
    pivot_wider(names_from = class, values_from = elements) |>
    mutate(total = rowSums(across(all_of(names(rep_class_map)))),
           total_TE = rowSums(across(all_of(TE_classes))),
           gsize = map_dbl(species, \(x) gstat_df$gsize[gstat_df$species == x]),
           log_gsize = map_dbl(species, \(x) gstat_df$log_gsize[gstat_df$species == x])) |>
    as.data.frame()
rep_len_df <- rep_df |>
    select(-elements, -prop, -plot_class) |>
    pivot_wider(names_from = class, values_from = length) |>
    mutate(total = rowSums(across(all_of(names(rep_class_map)))),
           total_TE = rowSums(across(all_of(TE_classes))),
           gsize = map_dbl(species, \(x) gstat_df$gsize[gstat_df$species == x]),
           log_gsize = map_dbl(species, \(x) gstat_df$log_gsize[gstat_df$species == x])) |>
    as.data.frame()
rep_prop_df <- rep_df |>
    select(-elements, -length, -plot_class) |>
    pivot_wider(names_from = class, values_from = prop) |>
    mutate(total = rowSums(across(all_of(names(rep_class_map)))),
           total_TE = rowSums(across(all_of(TE_classes))),
           gsize = map_dbl(species, \(x) gstat_df$gsize[gstat_df$species == x]),
           log_gsize = map_dbl(species, \(x) gstat_df$log_gsize[gstat_df$species == x])) |>
    as.data.frame()
rownames(rep_elem_df) <- paste(rep_elem_df$species)
rownames(rep_len_df) <- paste(rep_len_df$species)
rownames(rep_prop_df) <- paste(rep_prop_df$species)



#' ===========================================================================
#' ===========================================================================
#  Phylogenetic analyses ----
#' ===========================================================================
#' ===========================================================================

#'
#' This section tests whether traits differ between families
#' Chironomidae and Ceratopogonidae, and all others.
#'


dip_tr <- read.tree("~/_data/_phylo/chir_mcmctree.nwk") |>
    drop.tip("Cmarin") |>
    reorder("pruningwise")
#' Ives and Garland (2010; https://doi.org/10.1093/sysbio/syp074)
#' recommend scaling phylogeny to max depth of 1 for easier interpretation
#' (because alpha scales with phylogeny depth):
dip_tr$edge.length <- dip_tr$edge.length / max(node.depth.edgelength(dip_tr))
max(node.depth.edgelength(dip_tr))

dip_tr |> plot(no.margin = TRUE)


#' genome size ~ family
set.seed(531684999)
gs_mod <- phylolm(log_gsize ~ chir_cerat, gstat_df, dip_tr,
                  model = "OUfixedRoot", boot = 2000)
gs_mod |> summary()
gs_mod$bootconfint95
#       (Intercept) chir_ceratTRUE     sigma2     alpha
# 2.5%     8.587211     -0.9153641 0.05574064  1.017443
# 97.5%    8.960914     -0.4351406 3.73589239 50.000021


#' number of protein-coding genes ~ family
set.seed(741269662)
np_mod <- phylolm(n_prots ~ chir_cerat, gstat_df, dip_tr,
                  model = "OUfixedRoot", boot = 2000)
np_mod |> summary()
np_mod$bootconfint95
#       (Intercept) chir_ceratTRUE    sigma2      alpha
# 2.5%     12075.25      -700.9569   9117071  0.9514159
# 97.5%    16483.06      4593.7824 655126569 50.0000210


#' intron length ~ family
set.seed(1979093900)
il_mod <- phylolm(log_intron_len ~ chir_cerat, gstat_df, dip_tr,
                  model = "OUfixedRoot", boot = 2000)
il_mod |> summary()
il_mod$bootconfint95
#       (Intercept) chir_ceratTRUE     sigma2      alpha
# 2.5%     2.523898     -0.6707964 0.02535029  0.7385094
# 97.5%    2.786116     -0.3464602 2.06566365 50.0000210

#' n introns ~ family
set.seed(502219923)
in_mod <- phylolm(n_introns ~ chir_cerat, gstat_df, dip_tr,
                  model = "OUfixedRoot", upper.bound = 500, boot = 2000)
in_mod |> summary()
in_mod$bootconfint95
#       (Intercept) chir_ceratTRUE      sigma2      alpha
# 2.5%     4.403924     -0.8643251   0.5264592   1.280122
# 97.5%    5.307306      0.2404377 208.1284428 500.000000

#' **no significant effect of family, but when we drop Pstein there is one:
set.seed(81582400)
in_mod2 <- phylolm(n_introns ~ chir_cerat, filter(gstat_df, species != "Pstein"),
                  drop.tip(dip_tr, "Pstein"),
                  model = "OUfixedRoot", upper.bound = 500, boot = 2000)
in_mod2 |> summary()
in_mod2$bootconfint95
#       (Intercept) chir_ceratTRUE     sigma2      alpha
# 2.5%     4.600993     -0.8086019  0.1309948   1.125725
# 97.5%    5.129142     -0.1493447 71.0557109 500.000000


#' intergenic length ~ family
set.seed(502219923)
ig_mod <- phylolm(log_sum_interg_len ~ chir_cerat, gstat_df, dip_tr,
                  model = "OUfixedRoot", boot = 2000)
ig_mod |> summary()
ig_mod$bootconfint95
#       (Intercept) chir_ceratTRUE     sigma2      alpha
# 2.5%     8.212508     -0.8988728 0.06907438  0.5533094
# 97.5%    8.639594     -0.3441207 4.97920148 50.0000210



#' repeat classes ~ family
#' note: the Brownian motion (BM) model fit better for these regressions
set.seed(195060244)
rep_mods <- map(c(levels(rep_df$class), "total", "total_TE"), \(cl) {
    cl_mod <- phylolm(as.formula(paste(cl, "~ chir_cerat")),
                      rep_prop_df, dip_tr, model = "BM", boot = 2000)
    cl_mod[["call"]][[2]] <- as.formula(paste(cl, "~ chir_cerat"))
    return(cl_mod)
})
names(rep_mods) <- c(levels(rep_df$class), "total", "total_TE")

# map(rep_mods, summary)
map(rep_mods, \(z) z$bootconfint95)








#' ===========================================================================
#' ===========================================================================
#  Plots - tip traits ----
#' ===========================================================================
#' ===========================================================================


#' I'm re-reading tree here because the previous version had its branch
#' lengths re-scaled so are no longer in units of mya.
#' I also want full species names.
p_dip_tr <- read.tree("~/_data/_phylo/chir_mcmctree.nwk") |>
    drop.tip("Cmarin") |>
    (\(p) {
        p$tip.label <- map_chr(p$tip.label, \(z) spp_name_map[[z]])
        return(p)
    })()


tree_p <- (p_dip_tr |>
               ggtree() %<+%
               data.frame(taxa = paste(gstat_df$full_spp), spp = gstat_df$full_spp) +
               geom_rootedge(0.04) +
               geom_tiplab(aes(color = spp), size = 8 / 2.83465, fontface = "italic") +
               scale_color_manual(values = full_spp_pal, guide = "none")) |>
    revts() +
    scale_x_continuous("Million years ago", limits = c(-3.1, 1.6),
                       expand = c(0,0), breaks = -3:0, labels = 3:0 * 100) +
    theme_tree2()
tree_p


# Simpler version:
simp_tree_p <- (p_dip_tr |>
               ggtree() +
               geom_rootedge(0.04)) |>
    revts() +
    scale_x_continuous("mya", limits = c(-3.1, 0),
                       expand = c(0,0), breaks = -3:0, labels = 3:0 * 100) +
    theme_tree2()
simp_tree_p

# chironomidae node = 19
# cucidae node = 16
# simp_tree_p +
#     geom_cladelabel(node = 19, label = "Chironomidae", offset = 0.8, align = TRUE)
# simp_tree_p +
#     geom_hilight(node = 19, fill = "gray70")

# # alternatively:
# groupClade(ggtree(p_dip_tr), .node = 19) +
#     aes(color=group) +
#     theme(legend.position="right")

one_phy_trait_p <- function(.y_var, .y_mult, .y_lab, .y_breaks, .y_lims) {
    gstat_df |>
        mutate(species = factor(species, levels = rev(levels(species)))) |>
        select(species, all_of(.y_var)) |>
        ggplot(aes(species, .data[[.y_var]] * .y_mult, color = species)) +
        geom_segment(aes(yend = 0, xend = species), linewidth = 0.5, linetype = 3) +
        geom_point(size = 3) +
        scale_y_continuous(ifelse(grepl("\\*|'|~", .y_lab),
                                  parse(text = .y_lab), .y_lab),
                           breaks = .y_breaks) +
        scale_color_manual(values = rev(spp_pal), guide = "none") +
        coord_flip(ylim = .y_lims) +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.title.x = element_text(size = 10))
}



phy_trait_ps <- tribble(~.y_var, ~.y_mult, ~.y_lab, ~.y_breaks, ~.y_lims,
        "gsize", 1e-6, "Genome size (Mb)", 0:3*400, NULL,
        "intron_len", 1e-3, "Mean intron length (kb)", 0:3 * 4, NULL,
        "n_introns", 1, "Mean introns per gene", 4:6, c(3.95, NA),
        "n_prots", 1e-3, "'Protein-coding genes (' %*% '1000)'", 12 + 0:2 * 4,
            c(10.8, NA),
        "sum_interg_len", 1e-6, "Total intergenic DNA (Mb)", 0:3*200, NULL) |>
    pmap(one_phy_trait_p)

names(phy_trait_ps) <- c("gsize", "intron_len", "n_introns", "n_prots",
                         "sum_interg_len")


# If you want them all together with a single tree:
do.call(wrap_plots, c(list(nrow = 1, tree_p), phy_trait_ps))

# If you instead want each trait plot to have a tree, do this for each:
tree_p +
    phy_trait_ps[["gsize"]] +
    plot_layout(widths = c(0.5, 1), nrow = 1)


rep_p <- rep_df |>
    mutate(species = factor(species, levels = rev(levels(species)))) |>
    ggplot(aes(species, prop * 100, color = plot_class, fill = plot_class)) +
    geom_col(position = "stack", width = 0.7, color = NA) +
    scale_y_continuous("Percent of genome") +
    scale_fill_manual(values = rep_pal, aesthetics = c("color", "fill")) +
    coord_flip() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 9),
          axis.title.x = element_text(size = 10),
          legend.title = element_blank())

# Alternative using simpler tree with Chironomidae highlighted
(simp_tree_p +
        geom_hilight(node = 19, fill = "gray70")) +
    rep_p +
    plot_layout(widths = c(0.5, 1), nrow = 1)






#' ===========================================================================
#' ===========================================================================
#  Plots - gsize vs ----
#' ===========================================================================
#' ===========================================================================


set.seed(1812535342)
gs_nl_cor <- cor_phylo(~ log_gsize + log_intron_len, ~species, dip_tr,
                       data = gstat_df, boot = 2000)
# Pearson correlation coefficient estimate:
gs_nl_cor[["corrs"]][1,2]
# Bootstrapped 95% CI:
with(list(x = boot_ci(gs_nl_cor)[["corrs"]]), {c(x[2,1], x[1,2])})

set.seed(256063708)
gs_ni_cor <- cor_phylo(~ log_gsize + n_introns, ~species, dip_tr,
                       data = gstat_df, boot = 2000)
# Pearson correlation coefficient estimate:
gs_ni_cor[["corrs"]][1,2]
# Bootstrapped 95% CI:
with(list(x = boot_ci(gs_ni_cor)[["corrs"]]), {c(x[2,1], x[1,2])})
# note: removing Pstein doesn't change this!

set.seed(1180981201)
gs_np_cor <- cor_phylo(~ log_gsize + n_prots, ~species, dip_tr,
                       data = gstat_df, boot = 2000)
# Pearson correlation coefficient estimate:
gs_np_cor[["corrs"]][1,2]
# Bootstrapped 95% CI:
with(list(x = boot_ci(gs_np_cor)[["corrs"]]), {c(x[2,1], x[1,2])})

set.seed(538650984)
gs_ig_cor <- cor_phylo(~ log_gsize + log_sum_interg_len, ~species, dip_tr,
                       data = gstat_df, constrain_d = TRUE, method = "bobyqa",
                       boot = 2000)
# Pearson correlation coefficient estimate:
gs_ig_cor[["corrs"]][1,2]
# Bootstrapped 95% CI:
with(list(x = boot_ci(gs_ig_cor)[["corrs"]]), {c(x[2,1], x[1,2])})


gstat_df |>
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

gstat_df |>
    ggplot(aes(log_gsize, n_introns)) +
    geom_point(aes(color = species)) +
    scale_y_continuous("Mean introns per gene") +
    scale_x_continuous("Genome size (Mb)",
                       breaks = log10(100e6 * 2^(0:3)),
                       labels = 100 * 2^(0:3)) +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none")


gstat_df |>
    ggplot(aes(log_gsize, n_prots / 1e3)) +
    geom_point(aes(color = species)) +
    scale_y_continuous("Protein-coding genes (K)") +
    scale_x_continuous("Genome size (Mb)",
                       breaks = log10(100e6 * 2^(0:3)),
                       labels = 100 * 2^(0:3)) +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none")

gstat_df |>
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
    mutate(gsize = map_dbl(species, \(x) gstat_df$gsize[gstat_df$species == x]),
           log_gsize = map_dbl(species, \(x) gstat_df$log_gsize[gstat_df$species == x])) |>
    ggplot(aes(log_gsize, (prop))) +
    geom_point(aes(color = species)) +
    # scale_y_continuous("") +
    scale_x_continuous("Genome size (Mb)",
                       breaks = log10(100e6 * 2^(0:3)),
                       labels = 100 * 2^(0:3)) +
    scale_fill_manual(NULL, values = spp_pal, aesthetics = c("color", "fill"),
                      guide = "none")
