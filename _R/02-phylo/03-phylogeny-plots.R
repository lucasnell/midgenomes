
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
boot_tr <- read.tree(paste0(dirs$raxml_supp, "/chir_supp.raxml.supportFBP"))
ml_tr <- read.tree("_data/phylo/chir_ml.nwk")
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
sprintf("%.9f", cor(tt1$reltime, tt2$reltime))
sprintf("%.9f", cor(tt1$reltime, tt3$reltime))
sprintf("%.9f", cor(tt1$reltime, tt4$reltime))
sprintf("%.9f", cor(tt2$reltime, tt3$reltime))
sprintf("%.9f", cor(tt2$reltime, tt4$reltime))
sprintf("%.9f", cor(tt3$reltime, tt4$reltime))


# For each run, calculate minimum effective sample size for any estimate:

map_dfr(1:4, \(i) {
    ess <- paste0(dirs$mcmctree, sprintf("/mcmc_%i/chir_mcmctree_mcmc.txt", i)) |>
        read_tsv(col_types = paste(c("i", rep("d", 15)), collapse = "")) |>
        select(-Gen) |>
        summarize(across(everything(), coda::effectiveSize)) |>
        as.numeric()
    return(tibble(run = i, min_ess = min(ess)))
})





#' ===========================================================================
#' ===========================================================================
# Plot time-calibrated tree ----
#' ===========================================================================
#' ===========================================================================

#' Force tree to be ultrametric. Taken from
#' https://github.com/PuttickMacroevolution/MCMCtreeR/blob/2330e7a9916c3929513ee217d3854be993965f6b/R/readMCMCTree.R#L53-L70
#'
force_ultrametric <- function(treedata_obj) {
    phy <- treedata_obj@phylo
    outer <- phy$edge[,2]
    inner <- phy$edge[,1]
    totalPath <- c()
    for(i in which(outer<=Ntip(phy))) {
        start <- i
        end <- inner[start]
        edgeTimes <- phy$edge.length[start]
        while(end != inner[1]) {
            start <- which(outer == end)
            end <- inner[start]
            edgeTimes <- c(edgeTimes, phy$edge.length[start])
        }
        totalPath <- c(totalPath, sum(edgeTimes))
    }
    addLength <- max(totalPath) - totalPath
    phy$edge.length[which(outer <= Ntip(phy))] <- phy$edge.length[
        which(outer <= Ntip(phy))] + addLength
    treedata_obj@phylo <- phy
    return(treedata_obj)
}

time_tr <- read.mcmctree(paste0(dirs$mcmctree, "/mcmc_1/FigTree.tre")) |>
    force_ultrametric() |>
    (\(tr) {
        # convert units from 100 My to 1 My:
        tr@phylo$edge.length <- tr@phylo$edge.length * 100
        tr@data[["CI"]] <- tr@data[["0.95HPD"]] |>
            map(\(x) 100 * as.numeric(x))
        tr@data[["reltime"]] <- tr@data[["reltime"]] * 100
        # expand species names:
        tr@phylo$tip.label <- expand_spp(tr@phylo$tip.label)
        return(tr)
    })()

if (!file.exists("_data/phylo/time-tree.nwk")) {
    read.mcmctree(paste0(dirs$mcmctree, "/mcmc_1/FigTree.tre")) |>
        force_ultrametric() |>
        getElement("phylo") |>
        write.tree("_data/phylo/time-tree.nwk")
}




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
    scale_x_continuous("Million years ago", breaks = -3:0 * 100) +
    coord_geo(dat = "period",
              xlim = c(-325, 0), ylim = c(-1, Ntip(time_tr)+1), neg = TRUE,
              abbrv = FALSE, clip = "off",
              fill = plasma(8, begin = 0.2), color = NA,
              skip = c("Quaternary", "Neogene", "Carboniferous")) +
    theme(plot.margin = margin(0,0,0,r=1.5, unit = "in"))

time_tr_p

# save_plot("time-tree", time_tr_p, 6, 4, .png = FALSE)

