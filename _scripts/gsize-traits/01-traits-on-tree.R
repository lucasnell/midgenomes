

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




#' =================================================================
#' =================================================================
#  Read CI from phylolm models ----
#' =================================================================
#' =================================================================

regr_ci_list <- read_rds("_data/family-regressions.rds") |>
    select(feature, pl_model) |>
    pmap(\(feature, pl_model) {
        b0 <- pl_model$coefficients[["(Intercept)"]]
        b1 <- pl_model$coefficients[["chir_ceratTRUE"]]
        others <- quantile(pl_model$bootstrap[,"(Intercept)"], c(0.025, 0.975))
        chir_cerats <- quantile(pl_model$bootstrap[,"(Intercept)"] +
                                    pl_model$bootstrap[,"chir_ceratTRUE"],
                                c(0.025, 0.975))
        tibble(chir_cerat = c(TRUE, FALSE),
               name = feature,
               value = c(b0+b1, b0),
               lo = c(chir_cerats[["2.5%"]], others[["2.5%"]]),
               hi = c(chir_cerats[["97.5%"]], others[["97.5%"]]))
        bind_rows(tibble(x = c(0.75, 4.25),
                         chir_cerat = FALSE,
                         mid = b0,
                         lo = others[["2.5%"]],
                         hi = others[["97.5%"]]),
                  tibble(x = c(4.75, 14.25),
                         chir_cerat = TRUE,
                         mid = b0+b1,
                         lo = chir_cerats[["2.5%"]],
                         hi = chir_cerats[["97.5%"]]))
    })
names(regr_ci_list) <- read_rds("_data/family-regressions.rds")[["feature"]]





#' ===========================================================================
#' ===========================================================================
#  Inset tree and set aesthetics ----
#' ===========================================================================
#' ===========================================================================



#' Start with color palette from https://www.nature.com/articles/nmeth.1618
#' then remove those not needed:
rep_pal <- c(c("#000000", "#2271B2", "#3DB7E9", "#F748A5", "#359B73", "#d55e00",
               "#e69f00", "#f0e442")[c(2:4, 6:7)], "gray70") |>
    as.list()
names(rep_pal) <- repeat_classes

# An alternative starting point if you want:
c("#000000", "#AA0DB4", "#FF54ED", "#00B19F", "#EB057A", "#F8071D", "#FF8D1A", "#9EFF37")



shared_theme <- theme(panel.background = element_rect(fill = "transparent", color = NA),
                      plot.background = element_rect(fill = "transparent", color = NA),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      # axis.line.y = element_blank(),
                      axis.text.x = element_text(size = 8),
                      axis.title.x = element_text(size = 9))


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
    geom_tippoint(aes(color = species), size = 2) +
    # geom_tiplab(aes(color = species), size = 7 / 2.83465, fontface = "italic") +
    scale_color_manual(NULL, values = full_spp_pal, guide = "none") +
    coord_cartesian(xlim = c(0, 4)) +
    theme_inset()
# tree_inset




#' ===========================================================================
#' ===========================================================================
#  Feature plot ----
#' ===========================================================================
#' ===========================================================================


one_phy_feature_p <- function(.y_var, .y_subtr, .y_lab, .y_breaks, .y_lims) {

    .ymin <- min(feature_df[[.y_var]], na.rm = TRUE) - .y_subtr
    if (is.null(.y_lims)) {
        .y_lims <- c(.ymin, NA)
    } else {
        stopifnot(is.numeric(.y_lims) && length(.y_lims) == 2 &&
                      .y_lims[1] < .y_lims[2])
        if (.y_lims[1] < .ymin) .ymin <- .y_lims[1]
    }

    regr_df <- regr_ci_list[[.y_var]]

    .ribbons <- lapply(c(TRUE, FALSE), \(z) {
        geom_ribbon(data = filter(regr_df, chir_cerat == z),
                    aes(x = x, ymin = lo - .y_subtr, ymax = hi - .y_subtr),
                    color = NA, fill = "gray80")
    })

    feature_df |>
        mutate(species = factor(species, levels = rev(levels(species)))) |>
        select(species, all_of(.y_var)) |>
        filter(!is.na(.data[[.y_var]])) |>
        ggplot() +
        .ribbons +
        geom_line(data = regr_df, aes(x, mid - .y_subtr, group = chir_cerat),
                  linetype = 1, color = "gray50") +
        geom_point(aes(species, .data[[.y_var]] - .y_subtr, color = species),
                   size = 2) +
        scale_y_continuous(.y_lab, breaks = log10(.y_breaks),
                           labels = .y_breaks) +
        scale_x_discrete(drop = FALSE) +
        scale_color_manual(values = full_spp_pal, guide = "none",
                           aesthetics = c("color", "fill")) +
        coord_flip(ylim = .y_lims) +
        shared_theme
}


#' phy_feature_ps <- tribble(~.y_var, ~.y_mult, ~.y_lab, ~.y_breaks, ~.y_lims,
#'      "gsize", 1e-6, "Genome size\n(Mb)", 0:2*500, NULL,
#'      #' Below, '\u00D7' is unicode for multiply sign, and had to use lowercase
#'      #' '\u' instead of '\U' to have no space between sign and "1000"
#'      "n_genes", 1e-3, "Protein-coding\ngenes (\u00D71000)", 0:2*10, NULL, # c(10.8, NA),
#'      "sum_interg_len", 1e-6, "Total intergenic\nDNA (Mb)", 0:2*300, NULL,
#'      "mean_intron_len", 1e-3, "Mean intron\nlength (kb)", 0:2*6, NULL) |>
#'     pmap(one_phy_feature_p)
#' names(phy_feature_ps) <- c("gsize", "n_genes", "sum_interg_len",
#'                            "mean_intron_len")

phy_feature_ps <- tribble(~.y_var, ~.y_subtr, ~.y_lab, ~.y_breaks, ~.y_lims,
                           "log_gsize", 6, "Genome size\n(Mb)", 100 * 2^(0:3), NULL,
                         #' Below, '\u00D7' is unicode for multiply sign, and had to use lowercase
                         #' '\u' instead of '\U' to have no space between sign and "1000"
                         "log_n_genes", 3, "Protein-coding\ngenes (\u00D71000)", c(12, 15, 19), NULL,
                         "log_sum_interg_len", 6, "Total intergenic\nDNA (Mb)", 50 * 3^(0:2), NULL,
                         "mean_log_intron_len", 0, "Mean intron\nlength (bp)", 200 * 2^(0:2), NULL,
                         "log_SINE", 6, "SINE (Mb)", 5 * 10^(-2:0), NULL,
                         "log_LINE", 6, "LINE (Mb)", 10^(0:2), log10(c(0.3, 400)),
                         "log_LTR", 6, "LTR (Mb)", 10^(0:2), log10(c(0.3, 400)),
                         "log_DNA", 6, "DNA (Mb)", 10^(0:2), log10(c(0.3, 400)),
                         "log_non_TE", 6, "non-TE\n(Mb)", 5^(0:2), log10(c(0.17, 53)),
                         "log_Unclassified", 6, "Unclassified\n(Mb)", 5^(0:2), log10(c(0.17, 53))
                         ) |>
    pmap(one_phy_feature_p)


#' ===========================================================================
#' ===========================================================================
#  Combine and save  ----
#' ===========================================================================
#' ===========================================================================


traits_p <- do.call(wrap_plots, c(list(nrow = 1), phy_feature_ps[1:4]))
repeats_p <- do.call(wrap_plots, c(list(nrow = 1), phy_feature_ps[5:10]))

# traits_p / repeats_p

save_plot("cc-tree-inset", tree_inset, 1, 1.255, .png = FALSE)
save_plot("cc-traits", traits_p, 5.5, 2, .png = FALSE)
save_plot("cc-repeats", repeats_p, 5.5, 2, .png = FALSE)


