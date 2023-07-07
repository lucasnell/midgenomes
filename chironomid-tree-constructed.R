
library(tidyverse)
library(ape)
library(ggtree)
library(treeio)
library(viridis)

if (file.exists(".Rprofile")) source(".Rprofile")

save_plot <- function(n, p, w, h, ...) {
    cairo_pdf(sprintf("~/Stanford_Drive/UW/midgenomes/%s.pdf", n),
              width = w, height = h, bg = NA, ...)
    plot(p)
    dev.off()
}
#' To convert back to full names:
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


#' Tree from RAxML-NG with bootstrap (n = 1000) branch support
#' using Felsenstein bootstrap
boot_tr <- read.tree("~/_data/_phylo/chir_raxml_supp/chir_supp.raxml.supportFBP")
ml_tr <- read.tree("~/_data/_phylo/chir_ml.tree")
ml_tr$node.label <- boot_tr$node.label
ml_tr$tip.label <- map_chr(ml_tr$tip.label, ~ spp_name_map[[.x]])
#' Use this to set x axis bounds and labels
node.depth.edgelength(ml_tr) |> max()

ml_tr_p <- ggtree(ml_tr) +
    geom_rootedge(0.01) +
    geom_tiplab(size = 10 / 2.83465, fontface = "italic") +
    geom_nodelab(size = 8 / 2.83465, nudge_x = -0.025, nudge_y = 0.25) +
    theme_tree2() +
    scale_x_continuous("Mean substitutions per site", limits = c(-0.01, 1.05),
                       expand = c(0,0), breaks = 0:4 * 0.2) +
    NULL
ml_tr_p


# save_plot("tree_ML", ml_tr_p, 6, 3.5)


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

#' I plotted all 4 reps, and they produce the same tree (14 March 2023)
time_tr <- read.mcmctree("~/_data/_phylo/chir_mcmctree/mcmc_1/FigTree.tre") |>
    force_ultrametric()

# write.tree(time_tr@phylo, "~/_data/_phylo/chir_mcmctree.nwk")




# map(1:4, function(i) {
#     mcmc <- read_tsv(paste0("~/_data/_phylo/chir_mcmctree/mcmc_", i,
#                             "/chir_mcmctree_mcmc.txt"))
#     mcmc |>
#         select(starts_with("t_n"))
# })


time_tr@data[["CI"]] <- time_tr@data[["0.95HPD"]] |> map(as.numeric)
time_tr@phylo$tip.label <- map_chr(time_tr@phylo$tip.label, ~ spp_name_map[[.x]])



# Use this to find nodes you want to label:
ggtree(time_tr) + geom_text(aes(label=node), hjust=-.3)

# then to label the nodes, use...
# geom_cladelabel(node = X, label = "CLADE NAME", offset = 0.8, align = TRUE)
# OR
# geom_hilight(node = X, fill = "COLOR")


time_tr_p0 <- ggtree(time_tr) +
    geom_rootedge(0.04) +
    geom_tiplab(size = 10 / 2.83465, fontface = "italic") +
    geom_range("CI", color = "dodgerblue", size = 3, alpha = 0.5,
               center = "reltime") +
    theme_tree2()

time_tr_p <- time_tr_p0 |>
    revts() +
    scale_x_continuous("Million years ago", limits = c(-3.53367, 1.6), expand = c(0,0),
                       breaks = -3:0, labels = 3:0 * 100)

time_tr_p

# save_plot("tree_time", time_tr_p, 6, 3.5)





tree1 <- time_tr@phylo |>
    (function(x) {
        x_genera <- x$tip.label |>
            str_split(" ") |>
            map_chr(~ .x[[1]])
        x <- keep.tip(x, x$tip.label[!duplicated(x_genera)])
        x$tip.label <- x$tip.label |>
            str_split(" ") |>
            map_chr(~ .x[[1]])
        return(x)
    })()

tree2 <- read.tree("Cranston_2012.nwk") |>
    (function(x) {
        # genera in `time_tr@phylo`:
        my_genera <- time_tr@phylo$tip.label |>
            str_split(" ") |>
            map_chr(~ .x[[1]]) |>
            unique()
        x_genera <- x$tip.label |>
            str_split("_") |>
            map_chr(~ .x[[1]])
        lgl1 <- str_detect(x$tip.label, str_c("^", my_genera, collapse = "|"))
        x <- keep.tip(x, x$tip.label[lgl1 & !duplicated(x_genera)])
        x$tip.label <- x$tip.label |>
            str_split("_") |>
            map_chr(~ .x[[1]])
        return(x)
    })()

# plot(tree2); nodelabels(); axisPhylo()


mrca_time <- function(.tree, .sp1, .sp2) {
    .node <- getMRCA(.tree, c(.sp1, .sp2))
    max(node.depth.edgelength(.tree)) - node.depth.edgelength(.tree)[.node]
}


mrca_time(tree2, "Culicoides", "Tanytarsus")
mrca_time(tree2, "Parochlus", "Tanytarsus")
mrca_time(tree2, "Propsilocerus", "Tanytarsus")
mrca_time(tree2, "Clunio", "Tanytarsus")
mrca_time(tree2, "Polypedilum", "Tanytarsus")
mrca_time(tree2, "Polypedilum", "Chironomus")



for (spp in list(c("Culicoides", "Tanytarsus"),
c("Parochlus", "Tanytarsus"),
c("Propsilocerus", "Tanytarsus"),
c("Clunio", "Tanytarsus"),
c("Polypedilum", "Tanytarsus"),
c("Polypedilum", "Chironomus"))) {
    t1 <- mrca_time(tree1, spp[[1]], spp[[2]]) * 100
    t2 <- mrca_time(tree2, spp[[1]], spp[[2]])
    msg1 <- sprintf("%s -- %s = ",
            spp[[1]], spp[[2]]) |>
        str_pad(30, "right")
    msg2 <- sprintf("%7.2f |%7.2f |%5.2f (us | them | rel. diff.)\n",
                    t1, t2, abs(t1 - t2) / mean(c(t1, t2)))
    cat(paste0(msg1, msg2))
}






# =============================================================================
# =============================================================================
# =============================================================================

# BELOW IS EXPLORING DIFFERENT TREES WITH AND WITHOUT OUTGROUPS

# =============================================================================
# =============================================================================
# =============================================================================







bt_lines <- c("(((((Chironomus_riparius:0.028788,Chironomus_tepperi:0.025600):0.008002,Chironomus_tentans:0.019421):0.108075,(Polypedilum_pembai:0.031548,Polypedilum_vanderplanki:0.029075):0.107787):0.033845,Tanytarsus_gracilentus:0.150658):0.085242,(Clunio_marinus:0.114118,Belgica_antarctica:0.120344):0.043675,(Propsilocerus_akamusi:0.124211,((Culicoides_sonorensis:0.309562,Anopheles_stephensi:0.286965):0.130299,Parochlus_steinenii:0.188874):0.150195):0.026016);",
"((((Chironomus_tentans:0.019063,(Chironomus_riparius:0.028818,Chironomus_tepperi:0.025412):0.007965):0.107509,(Polypedilum_pembai:0.031416,Polypedilum_vanderplanki:0.029336):0.106407):0.034020,Tanytarsus_gracilentus:0.150539):0.085533,(Clunio_marinus:0.113747,Belgica_antarctica:0.119879):0.044267,((Parochlus_steinenii:0.189637,(Culicoides_sonorensis:0.308262,Anopheles_stephensi:0.285912):0.130360):0.150050,Propsilocerus_akamusi:0.124203):0.026159);",
"((Propsilocerus_akamusi:0.124500,((Anopheles_stephensi:0.286971,Culicoides_sonorensis:0.312491):0.130820,Parochlus_steinenii:0.190219):0.149749):0.026774,(((Chironomus_tentans:0.019329,(Chironomus_riparius:0.028877,Chironomus_tepperi:0.025466):0.007742):0.106546,(Polypedilum_pembai:0.031411,Polypedilum_vanderplanki:0.029394):0.107079):0.034006,Tanytarsus_gracilentus:0.151206):0.085562,(Belgica_antarctica:0.120046,Clunio_marinus:0.114474):0.043482);",
"((((Parochlus_steinenii:0.190057,(Culicoides_sonorensis:0.310311,Anopheles_stephensi:0.287772):0.131149):0.151462,Propsilocerus_akamusi:0.124546):0.026706,(Clunio_marinus:0.114560,Belgica_antarctica:0.120289):0.043150):0.085973,(((Chironomus_riparius:0.029034,Chironomus_tepperi:0.024908):0.007858,Chironomus_tentans:0.019469):0.107791,(Polypedilum_vanderplanki:0.029858,Polypedilum_pembai:0.031382):0.107674):0.034189,Tanytarsus_gracilentus:0.150574);",
"((((Parochlus_steinenii:0.189321,(Anopheles_stephensi:0.285715,Culicoides_sonorensis:0.310995):0.130366):0.149742,Propsilocerus_akamusi:0.123982):0.026463,(((Chironomus_tentans:0.018972,(Chironomus_tepperi:0.025268,Chironomus_riparius:0.028808):0.007860):0.106977,(Polypedilum_pembai:0.031248,Polypedilum_vanderplanki:0.029236):0.108174):0.033470,Tanytarsus_gracilentus:0.150447):0.086420):0.043931,Clunio_marinus:0.114298,Belgica_antarctica:0.119618);",
"(((Chironomus_tentans:0.019418,(Chironomus_tepperi:0.025086,Chironomus_riparius:0.029073):0.007958):0.107041,(Polypedilum_vanderplanki:0.029433,Polypedilum_pembai:0.031444):0.108209):0.034213,((((Culicoides_sonorensis:0.310853,Anopheles_stephensi:0.284779):0.131698,Parochlus_steinenii:0.189661):0.151380,Propsilocerus_akamusi:0.124252):0.026188,(Belgica_antarctica:0.119511,Clunio_marinus:0.113914):0.043126):0.085895,Tanytarsus_gracilentus:0.151090);")

write_lines(bt_lines, "~/_data/_phylo/boot_noog.bootstraps")


map_lgl(bt_lines, ~ read.tree(text = .x) |> is.ultrametric())


bt_lines0 <- read_lines("/Users/lucasnell/_data/_phylo/chir_raxml_boot_0/chir_phy.raxml.bootstraps")

map_lgl(bt_lines0, ~ read.tree(text = .x) |> is.ultrametric())
read.tree(text = bt_lines0[1])







tree <- read.tree(text = "((((((((Chironomus_tentans:0.019199,(Chironomus_riparius:0.028791,Chironomus_tepperi:0.025378):0.007889):0.107427,(Polypedilum_vanderplanki:0.029482,Polypedilum_pembai:0.031424):0.107480):0.034059,Tanytarsus_gracilentus:0.150572):0.086032,(Clunio_marinus:0.114155,Belgica_antarctica:0.120113):0.043496):0.026316,Propsilocerus_akamusi:0.124744):0.150301,Parochlus_steinenii:0.189277):0.130557,Culicoides_sonorensis:0.310669):0.143039,Anopheles_stephensi:0.143039);
")

plot(tree); # axisPhylo()

# ggtree(tree) +
#     geom_tiplab(size = 3)


tree2 <- read.tree("Cranston_2012.nwk") |>
    (function(x) {
        # genera in `tree`:
        genera <- tree$tip.label |>
            str_split("_") |>
            map_chr(~ .x[[1]]) |>
            unique()
        lgl <- str_detect(x$tip.label, str_c("^", genera, collapse = "|"))
        keep.tip(x, x$tip.label[lgl])
    })()

plot(tree2); nodelabels(); axisPhylo()


mrca_time <- function(.tree, .sp1, .sp2) {
    .node <- getMRCA(.tree, c(.sp1, .sp2))
    max(node.depth.edgelength(.tree)) - node.depth.edgelength(.tree)[.node]
}


mrca_time(tree2, "Culicoides_leechi", "Tanytarsus_sp._1_SRM-2010")
mrca_time(tree2, "Parochlus_spinosus", "Tanytarsus_sp._1_SRM-2010")
mrca_time(tree2, "Polypedilum_vanderplanki", "Tanytarsus_sp._1_SRM-2010")


# ggtree(tree2) +
#     geom_tiplab(size = 3)


