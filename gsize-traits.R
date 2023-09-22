
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
spp_name_map <- list("Aaegyp" = "Aedes aegypti",
                     "Asteph" = "Anopheles stephensi",
                     "Bantar" = "Belgica antarctica",
                     "Cripar" = "Chironomus riparius",
                     "Ctenta" = "Chironomus tentans",
                     "Cmarin" = "Clunio marinus",
                     "Cquinq" = "Culex quinquefasciatus",
                     "Csonor" = "Culicoides sonorensis",
                     "Mdomes" = "Musca domestica",
                     "Pstein" = "Parochlus steinenii",
                     "Ppemba" = "Polypedilum pembai",
                     "Pvande" = "Polypedilum vanderplanki",
                     "Pakamu" = "Propsilocerus akamusi",
                     "Tgraci" = "Tanytarsus gracilentus")
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
           prop_interg_len = sum_interg_len / gsize,
           prop_intron_len = (intron_len * n_introns * n_prots) / gsize,
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
    filter(species != "Cmarin",
           class != "Unclassified") |>
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
           chir_cerat = family %in% c("Chironomidae", "Ceratopogonidae"),
           chir_cerat_fct = factor(chir_cerat, levels = c(TRUE, FALSE),
                                   labels = c("Chir. + Cerat.", "Others"))) |>
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



#' =================================================================
#' =================================================================
#  Read phylogeny ----
#' =================================================================
#' =================================================================


dip_tr <- read.tree("~/_data/_phylo/chir_mcmctree.nwk") |>
    drop.tip("Cmarin") |>
    reorder("pruningwise")
#' Ives and Garland (2010; https://doi.org/10.1093/sysbio/syp074)
#' recommend scaling phylogeny to max depth of 1 for easier interpretation
#' (because alpha scales with phylogeny depth):
dip_tr$edge.length <- dip_tr$edge.length / max(node.depth.edgelength(dip_tr))
max(node.depth.edgelength(dip_tr))

# dip_tr |> plot(no.margin = TRUE)




#' ===========================================================================
#' ===========================================================================
#  Regressions - family versus... ----
#' ===========================================================================
#' ===========================================================================

#'
#' This section tests whether traits differ between families
#' Chironomidae and Ceratopogonidae, and all others.
#'


if (!file.exists("_data/family-regressions.rds")) {
    # Takes ~36 sec
    set.seed(1269567689)
    regr_df <- tibble(feature = c("log_gsize", "n_prots", "log_intron_len",
                                  "n_introns", "n_introns_no_Pstein",
                                  "log_sum_interg_len",
                                  "prop_intron_len", "prop_interg_len")) |>
        mutate(pl_model = map(feature, \(x) {
            .boots <- 2000
            f <- as.formula(paste(x, "~ chir_cerat"))
            if (x %in% c("prop_intron_len", "prop_interg_len")) {
                m <- phylolm(f, gstat_df, dip_tr, model = "OUfixedRoot",
                             upper.bound = 200, boot = .boots)
            } else if (x == "n_introns") {
                m <- phylolm(f, gstat_df, dip_tr, model = "OUfixedRoot",
                             upper.bound = 500, boot = .boots)
            } else if (x == "n_introns_no_Pstein") {
                f <- as.formula("n_introns ~ chir_cerat")
                m <- phylolm(f, filter(gstat_df, species != "Pstein"),
                             drop.tip(dip_tr, "Pstein"),
                             model = "OUfixedRoot", upper.bound = 500,
                             boot = .boots)
            } else {
                m <- phylolm(f, gstat_df, dip_tr, model = "OUfixedRoot",
                             boot = .boots)
            }
            m$call[[2]] <- eval(f)
            return(m)
        })) |>
        mutate(coef = map_dbl(pl_model, \(x) coef(x)[["chir_ceratTRUE"]]),
               coef_lo = map_dbl(pl_model,
                                 \(x) x$bootconfint95[1,"chir_ceratTRUE"]),
               coef_hi = map_dbl(pl_model,
                                 \(x) x$bootconfint95[2,"chir_ceratTRUE"])) |>
        select(feature, starts_with("coef"), pl_model)
    write_rds(regr_df, "_data/family-regressions.rds", compress = "xz")
} else {
    regr_df <- read_rds("_data/family-regressions.rds")
}

regr_df
# # A tibble: 8 × 5
#   feature                  coef   coef_lo  coef_hi pl_model
#   <chr>                   <dbl>     <dbl>    <dbl> <list>
# 1 log_gsize             -0.681    -0.921    -0.459 <phylolm>
# 2 n_prots             2043.     -502.     4646.    <phylolm>
# 3 log_intron_len        -0.502    -0.665    -0.341 <phylolm>
# 4 n_introns             -0.305    -0.896     0.264 <phylolm>
# 5 n_introns_no_Pstein   -0.475    -0.796    -0.158 <phylolm>
# 6 log_sum_interg_len    -0.622    -0.901    -0.341 <phylolm>
# 7 prop_intron_len       -0.312    -0.468    -0.159 <phylolm>
# 8 prop_interg_len        0.0599   -0.0135    0.134 <phylolm>


#' repeat classes ~ family
#' note: the Brownian motion (BM) model fit better for these regressions
set.seed(195060244)
rep_mods <- map(c(levels(rep_df$class), "total", "total_TE"), \(cl) {
    f <- as.formula(paste(cl, "~ chir_cerat"))
    cl_mod <- phylolm(f, rep_len_df, dip_tr, model = "BM", boot = 2000)
    cl_mod$call[[2]] <- eval(f)
    return(cl_mod)
})
names(rep_mods) <- c(levels(rep_df$class), "total", "total_TE")

# map(rep_mods, summary)
map(rep_mods, \(z) z$bootconfint95)


rep_regr_df <- tibble(class = c(levels(rep_df$class), "total", "total_TE")) |>
    mutate(pl_model1 = map(class, \(x) {
        f <- as.formula(sprintf("log10(%s) ~ chir_cerat", x))
        if (x == "LINE") {
            cl_mod <- phylolm(f, rep_len_df, dip_tr, model = "OUfixedRoot",
                              upper.bound = 200, boot = 0)
        } else {
            cl_mod <- phylolm(f, rep_len_df, dip_tr, model = "OUfixedRoot",
                              boot = 0)
        }
        cl_mod$call[[2]] <- eval(f)
        return(cl_mod)
    }),
    pl_model2 = map(class, \(x) {
        f <- as.formula(sprintf("log10(%s) ~ chir_cerat", x))
        cl_mod <- phylolm(f, rep_len_df, dip_tr, model = "BM", boot = 0)
        cl_mod$call[[2]] <- eval(f)
        return(cl_mod)
    })) |>
    # mutate(coef = map_dbl(pl_model, \(x) coef(x)[["chir_ceratTRUE"]]),
    #        coef_lo = map_dbl(pl_model,
    #                          \(x) x$bootconfint95[1,"chir_ceratTRUE"]),
    #        coef_hi = map_dbl(pl_model,
    #                          \(x) x$bootconfint95[2,"chir_ceratTRUE"])) |>
    # select(class, starts_with("coef"), pl_model)
    mutate(aic1 = map_dbl(pl_model1, \(x) AIC(x)),
           aic2 = map_dbl(pl_model2, \(x) AIC(x)),
           daic = aic1 - aic2)
rep_regr_df


rep_len_pglmm_df <- rep_len_df |>
    as_tibble() |>
    select(-total, -total_TE) |>
    pivot_longer(DNA:Small_RNA, names_to = "class", values_to = "len") |>
    mutate(log_len = log10(len)) |>
    group_by(class) |>
    mutate(z_len = (len - mean(len)),
           z_log_len = (log_len - mean(log_len))) |>
    ungroup()

rep_len_pglmm_df |>
    ggplot(aes(log_gsize, z_log_len)) +
    geom_point() +
    facet_wrap(~ class)

rep_len_pglmm_df |>
    ggplot(aes(chir_cerat_fct, z_log_len)) +
    geom_hline(yintercept = 0, color = "gray70", linewidth = 0.75) +
    geom_jitter(aes(color = full_spp), height = 0, width = 0.1, size = 2) +
    facet_wrap(~ class) +
    scale_color_manual(values = full_spp_pal) +
    theme(axis.text.x = element_text(size = 9, color = "black"),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(face = "italic"))










#' ===========================================================================
#' ===========================================================================
#  Correlations - gsize vs ... ----
#' ===========================================================================
#' ===========================================================================




m <- cor_phylo(~ log_gsize + log_intron_len + n_introns + n_prots +
                   log_sum_interg_len + prop_intron_len + prop_interg_len,
               ~species, dip_tr, data = gstat_df,
               # constrain_d = TRUE,
               method = "bobyqa",
               max_iter = 10e3,
               boot = 0)
m



if (!file.exists("_data/gsize-corrs.rds")) {
    # Takes ~1 min
    set.seed(328380265)
    corr_df <- tibble(feature = c("log_intron_len", "n_introns", "n_prots",
                                  "log_sum_interg_len", "prop_intron_len",
                                  "prop_interg_len")) |>
        mutate(cp_model = map(feature, \(x) {
            f <- as.formula(paste0("~ log_gsize + ", x))
            if (x == "log_sum_interg_len" || x == "prop_intron_len") {
                m <- cor_phylo(f, ~species, dip_tr, data = gstat_df,
                               constrain_d = TRUE, method = "bobyqa",
                               boot = 2000)
            } else {
                m <- cor_phylo(f, ~species, dip_tr, data = gstat_df,
                               boot = 2000)
            }
            m$call[[2]] <- eval(f)
            # cat(x, "\n")
            return(m)
        })) |>
        mutate(corr = map_dbl(cp_model, \(x) x[["corrs"]][1,2]),
               corr_lo = map_dbl(cp_model, \(x) boot_ci(x)[["corrs"]][2,1]),
               corr_hi = map_dbl(cp_model, \(x) boot_ci(x)[["corrs"]][1,2])) |>
        select(feature, starts_with("corr"), cp_model)

    write_rds(corr_df, "_data/gsize-corrs.rds", compress = "xz")
} else {
    corr_df <- read_rds("_data/gsize-corrs.rds")
}

corr_df





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
