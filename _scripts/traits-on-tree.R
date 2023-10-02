

source("_scripts/00-preamble.R")

library(ggtree)
library(treeio)
library(patchwork)




#' =================================================================
#' =================================================================
#  Info to make labels pretty ----
#' =================================================================
#' =================================================================

#' Repeat element classes in order with formatted names for plotting,
#' where all non-TE classes are grouped together:
rep_class_map <- list("SINE" = "SINEs",
                      "LINE" = "LINEs",
                      "LTR" = "LTR elements",
                      "DNA" = "DNA transposons",
                      "non_TE" = "non-TE elements",
                      "Unclassified" = "Unclassified")

# non-TE (but known) element classes:
nonTE_classes <- c("RC", "Small_RNA", "Satellite", "Simple_repeat",
                   "Low_complexity")
# TE element classes:
TE_classes <- c("SINE", "LINE", "LTR", "DNA")




#' =================================================================
#' =================================================================
#  Read genome stats ----
#' =================================================================
#' =================================================================


#' We're going to make separate plots for repeats and other features,
#' so they'll be read into separate data frames.


feature_df <- read_csv("_data/genome-stats.csv", col_types = cols()) |>
    # ordering factors so their levels are in order they show up
    mutate(family = factor(family, levels = unique(family)),
           species = factor(species, levels = species),
           spp_abbrev = factor(spp_abbrev, levels = spp_abbrev)) |>
    mutate(log_sum_interg_len = log10(sum_interg_len),
           log_gsize = log10(gsize),
           log_n_genes = log10(n_genes)) |>
    select(species, spp_abbrev, contains("gsize"), contains("n_genes"),
           contains("interg"), contains("intron"))



# Same as above, but all non-TE classes grouped together
repeat_df <- read_csv("_data/genome-stats.csv", col_types = cols()) |>
    # ordering factors so their levels are in order they show up
    mutate(family = factor(family, levels = unique(family)),
           species = factor(species, levels = species),
           spp_abbrev = factor(spp_abbrev, levels = spp_abbrev)) |>
    mutate(non_TE = rowSums(across(all_of(nonTE_classes)))) |>
    select(species, spp_abbrev, all_of(names(rep_class_map))) |>
    pivot_longer(all_of(names(rep_class_map)),
                 names_to = "class", values_to = "len") |>
    mutate(class = factor(class, levels = names(rep_class_map),
                          labels = unlist(rep_class_map)))






#' =================================================================
#' =================================================================
#  Read phylogeny ----
#' =================================================================
#' =================================================================


p_dip_tr <- read.tree("_data/phylo/time-tree.nwk") |>
    (\(p) {
        p$tip.label <- expand_spp(p$tip.label)
        return(p)
    })()




#' ===========================================================================
#' ===========================================================================
#  Inset tree and set aesthetics ----
#' ===========================================================================
#' ===========================================================================


spp_pal <- full_spp_pal <- turbo(100)[c(70+3*0:8, 60, 15+4*2:0, 30)] |> as.list()
names(spp_pal) <- levels(feature_df$spp_abbrev)
names(full_spp_pal) <- levels(feature_df$species)

#' Start with color palette from https://www.nature.com/articles/nmeth.1618
#' then remove those not needed:
rep_pal <- c(c("#000000", "#2271B2", "#3DB7E9", "#F748A5", "#359B73", "#d55e00",
               "#e69f00", "#f0e442")[c(2:6)], "gray70") |>
    as.list()
names(rep_pal) <- paste(rep_class_map)

# An alternative starting point if you want:
# c("#000000", "#AA0DB4", "#FF54ED", "#00B19F", "#EB057A", "#F8071D", "#FF8D1A", "#9EFF37")



shared_theme <- theme(panel.background = element_rect(fill = "transparent", color = NA),
                      plot.background = element_rect(fill = "transparent", color = NA),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.line.y = element_blank(),
                      axis.text.x = element_text(size = 8),
                      axis.title.x = element_text(size = 9))



# color lines by whether it's Chironomidae:
tree_p <- ggtree(p_dip_tr, linewidth = 1, lineend = "square") |>
    groupClade(.node = 20) +
    aes(color = group) +
    theme(legend.position = "none") +
    scale_color_manual(values = c("gray70", "black"))
# tree_p


# # Simpler version:
# simp_tree_p <- (p_dip_tr |>
#                ggtree() +
#                geom_rootedge(0.04)) |>
#     revts() +
#     theme_inset()
# # simp_tree_p






#' ===========================================================================
#' ===========================================================================
#  Feature plot ----
#' ===========================================================================
#' ===========================================================================


one_phy_feature_p <- function(.y_var, .y_mult, .y_lab, .y_breaks, .y_labels, .y_lims) {

    feature_df |>
        mutate(species = factor(species, levels = rev(levels(species)))) |>
        select(species, all_of(.y_var)) |>
        filter(!is.na(.data[[.y_var]])) |>
        # ggplot(aes(species, .data[[.y_var]] * .y_mult,
        #            color = species, fill = species)) +
        ggplot(aes(species, .data[[.y_var]] * .y_mult)) +
        geom_hline(yintercept = 0, color = "gray70", linewidth = 0.75) +
        # geom_segment(aes(yend = 0, xend = species), linewidth = 0.5,
        #              linetype = "22") +
        # geom_point(size = 3) +
        geom_col(width = 0.6, fill = "#00B19F", color = NA) +
        scale_y_continuous(.y_lab, breaks = .y_breaks) +
        scale_x_discrete(drop = FALSE) +
        scale_color_manual(values = rev(full_spp_pal), guide = "none",
                           aesthetics = c("color", "fill")) +
        coord_flip(ylim = .y_lims) +
        shared_theme
}


phy_feature_ps <- tribble(~.y_var, ~.y_mult, ~.y_lab, ~.y_breaks, ~.y_lims,
     "gsize", 1e-6, "Genome size\n(Mb)", 0:2*500, NULL,
     #' Below, '\u00D7' is unicode for multiply sign, and had to use lowercase
     #' '\u' instead of '\U' to have no space between sign and "1000"
     "n_genes", 1e-3, "Protein-coding\ngenes (\u00D71000)", 0:2*10, NULL, # c(10.8, NA),
     "sum_interg_len", 1e-6, "Total intergenic\nDNA (Mb)", 0:2*300, NULL,
     "mean_intron_len", 1e-3, "Mean intron\nlength (kb)", 0:2*6, NULL) |>
    pmap(one_phy_feature_p)

names(phy_feature_ps) <- c("gsize", "n_genes", "sum_interg_len",
                           "mean_intron_len")


# All together with a single tree:
# do.call(wrap_plots, c(list(nrow = 1, tree_p), phy_feature_ps))


# # If you instead want each trait plot to have a tree, do this for each:
# tree_p +
#     phy_feature_ps[["gsize"]] +
#     plot_layout(widths = c(0.5, 1), nrow = 1)



#' ===========================================================================
#' ===========================================================================
#  Repeats plot ----
#' ===========================================================================
#' ===========================================================================




repeat_p <- repeat_df |>
    mutate(species = factor(species, levels = rev(levels(species)))) |>
    ggplot(aes(species, len / 1e9, color = class, fill = class)) +
    geom_col(position = "stack", width = 0.6, color = NA) +
    scale_y_continuous("Genome content (Gb)", breaks = 0:2*0.5) +
    scale_fill_manual(values = rep_pal) +
    coord_flip() +
    shared_theme +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 6),
          legend.key.size = unit(1, "lines"),
          legend.background = element_blank(),
          legend.position = c(0.95, 0.95),
          legend.justification = c(1, 1))





#' ===========================================================================
#' ===========================================================================
#  Combine and save  ----
#' ===========================================================================
#' ===========================================================================


traits_p <- do.call(wrap_plots, c(list(nrow = 1, widths = c(rep(1, 5), 2), tree_p),
                      phy_feature_ps, list(repeat_p)))
# traits_p

# save_plot("phylo-traits", traits_p, 6.5, 4, .png = FALSE)


