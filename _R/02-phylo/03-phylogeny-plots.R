
#'
#' Plot maximum likelihood and time-calibrated phylogenetic trees
#'

source("_R/00-preamble.R")

library(ggtree)
library(treeio)
library(deeptime)




#' ===========================================================================
#' ===========================================================================
# ML Tree with support ----
#' ===========================================================================
#' ===========================================================================


#' Tree from RAxML-NG with bootstrap (n = 1000) branch support
#' using Felsenstein bootstrap
boot_tr <- read.tree(paste0(dirs$raxml_boot, "/chir_supp/",
                            "chir_raxml_supp.raxml.supportFBP"))
ml_tr <- read.tree("_data/phylo/chir_ml.tree")
ml_tr$node.label <- boot_tr$node.label
ml_tr$tip.label <- expand_spp(ml_tr$tip.label)

ml_tr_p <- ggtree(ml_tr) +
    geom_rootedge(0.01) +
    geom_tiplab(size = 9 / 2.83465, fontface = "italic") +
    geom_nodelab(size = 8 / 2.83465, nudge_x = -0.025, nudge_y = 0.25) +
    theme_tree2() +
    scale_x_continuous("Mean substitutions per site",
                       breaks = 0:4 * 0.2) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(0,0,0,r=1.5, unit = "in"))
# ml_tr_p


# save_plot("ML-tree", ml_tr_p, 6.5, 4)




#' ===========================================================================
#' ===========================================================================
# Check MCMCTree consistency ----
#' ===========================================================================
#' ===========================================================================


tt1 <- read.mcmctree(paste0(dirs$mcmctree, "/mcmc_1/FigTree.tre"))@data
tt2 <- read.mcmctree(paste0(dirs$mcmctree, "/mcmc_2/FigTree.tre"))@data
tt3 <- read.mcmctree(paste0(dirs$mcmctree, "/mcmc_3/FigTree.tre"))@data
tt4 <- read.mcmctree(paste0(dirs$mcmctree, "/mcmc_4/FigTree.tre"))@data

# All rows correspond to the same nodes:
all(tt1$node == tt2$node) && all(tt1$node == tt3$node) && all(tt1$node == tt4$node)

# All times are consistent across runs:
sprintf("%.3f", cor(tt1$reltime, tt2$reltime))
sprintf("%.3f", cor(tt1$reltime, tt3$reltime))
sprintf("%.3f", cor(tt1$reltime, tt4$reltime))
sprintf("%.3f", cor(tt2$reltime, tt3$reltime))
sprintf("%.3f", cor(tt2$reltime, tt4$reltime))
sprintf("%.3f", cor(tt3$reltime, tt4$reltime))


# For each run, calculate minimum effective sample size for any estimate
# (takes ~7 sec with >=4 threads):

mclapply(1:4, \(i) {
    ess <- paste0(dirs$mcmctree, sprintf("/mcmc_%i/chir_mcmctree_%i.txt", i, i)) |>
        read_tsv(col_types = paste(c("i", rep("d", 15)), collapse = "")) |>
        select(-Gen) |>
        summarize(across(everything(), coda::effectiveSize)) |>
        as.numeric()
    return(tibble(run = i, min_ess = min(ess)))
}) |>
    do.call(what = bind_rows)
# # A tibble: 4 × 2
#     run min_ess
#   <int>   <dbl>
# 1     1    642.
# 2     2    652.
# 3     3    559.
# 4     4    567.




#' ===========================================================================
#' ===========================================================================
# Plot time-calibrated tree ----
#' ===========================================================================
#' ===========================================================================

time_tr <- read.mcmctree(paste0(dirs$mcmctree, "/mcmc_1/FigTree.tre")) |>
    (\(tr) {
        # use the pre-processed (ultrametric) tree for phylogeny:
        tr@phylo <- read.tree("_data/phylo/time-tree.nwk")
        # convert units from 100 My to 1 My:
        tr@phylo$edge.length <- tr@phylo$edge.length * 100
        tr@data[["CI"]] <- tr@data[["0.95HPD"]] |>
            map(\(x) 100 * as.numeric(x))
        tr@data[["reltime"]] <- tr@data[["reltime"]] * 100
        # expand species names:
        tr@phylo$tip.label <- expand_spp(tr@phylo$tip.label)
        return(tr)
    })()




time_tr_p0 <- time_tr |>
    mutate(species = factor(label, levels = names(full_spp_pal))) |>
    ggtree() +
    geom_rootedge(0.04) +
    geom_tiplab(size = 9 / 2.83465, fontface = "bold.italic", aes(color = species)) +
    geom_range("CI", color = "gray50", alpha = 0.5, size = 3, center = "reltime") +
    scale_color_manual(NULL, values = full_spp_pal, guide = "none") +
    theme_tree2()

time_tr_p <- time_tr_p0 |>
    revts() +
    scale_x_continuous("Million years ago", breaks = -3:0 * 100,
                       labels = 3:0 * 100) +
    coord_geo(dat = "period",
              xlim = c(-350, 0), ylim = c(0, Ntip(time_tr)+1), neg = TRUE,
              abbrv = FALSE, clip = "off",
              fill = plasma(8, begin = 0.2), color = NA,
              height = unit(1.5, "line"), size = 7/2.8,
              center_end_labels = TRUE,
              skip = c("Quaternary", "Neogene")) +
    theme(plot.margin = margin(0,0,0,r=1.7, unit = "in"))

time_tr_p

# save_plot("time-tree", time_tr_p, 6, 4, .png = FALSE)

