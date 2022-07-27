

#' Make two trees, each based on the ML tree, but with fossil evidence and
#' no branch lengths.
#' They'll be used in MCMCtree and CODEML.

library(ape)
library(tidyverse)

ml_tr <- read.tree("~/_data/phylo/chir_ml.tree")

#'
#' Nodes to calibrate
#' (I'm starting with B bc A is the root calibration that I'll do in the
#'  MCMCtree control file)
#'
nodeB <- getMRCA(ml_tr, tip = c("Culicoides_sonorensis", "Tanytarsus_gracilentus"))
nodeC <- getMRCA(ml_tr, tip = c("Parochlus_steinenii", "Tanytarsus_gracilentus"))
nodeD <- getMRCA(ml_tr, tip = c("Chironomus_riparius", "Tanytarsus_gracilentus"))
nodeE <- getMRCA(ml_tr, tip = c("Clunio_marinus", "Belgica_antarctica"))


#'
#' First create tree with no branch lengths but some fixed time calibrations
#' for use in CODEML.
#'
codeml_tr <- ml_tr
codeml_tr$edge.length <- NULL
codeml_tr$node.label <- rep("", codeml_tr$Nnode)
codeml_tr$node.label[nodeC - Ntip(codeml_tr)] <- "'@2.327'"   # from Cranston et al. (2012)
codeml_tr$node.label[nodeD - Ntip(codeml_tr)] <- "'@1.37573'" # from Cranston et al. (2012)
codeml_tr$node.label[nodeE - Ntip(codeml_tr)] <- "'@0.76'"      # from timetree.org
write.tree(codeml_tr, "~/_data/phylo/codeml_in_tree.nwk")


#'
#' Next create a similar tree but with more information on calibrations, and
#' these calibrations are to describe distributions.
#'

calib_tr <- ml_tr
calib_tr$edge.length <- NULL
calib_tr$node.label <- rep("", calib_tr$Nnode)

#' The extra complexity below is for two reasons:
#'
#' 1. `ape` removes parentheses and commas when writing to a newick file, so
#'    I'm assigning each character a "safe" alternative, then replacing them
#'    in the output string from `ape` before writing to file.
#' 2. I need to add the number of species and trees at the top of the input
#'    file to mcmctree.
#'

calib_tr$node.label[nodeB - Ntip(calib_tr)] <- "'L__2.420-0.10-0.2-0.1--'"  # B
calib_tr$node.label[nodeC - Ntip(calib_tr)] <- "'L__2.013-0.16-0.5--'"      # C
calib_tr$node.label[nodeD - Ntip(calib_tr)] <- "'L__0.935-0.47-0.5--'"      # D
calib_tr$node.label[nodeE - Ntip(calib_tr)] <- "'L__0.339-1.24-1.0--'"      # E

write_file(paste(Ntip(calib_tr), 1), "~/_data/phylo/mcmctree_in_tree.nwk")
write_file("\n", "~/_data/phylo/mcmctree_in_tree.nwk", append = TRUE)
write.tree(calib_tr) %>%
    str_replace_all("__", "(") %>%
    str_replace_all("--", ")") %>%
    str_replace_all("-", ", ") %>%
    write_file("~/_data/phylo/mcmctree_in_tree.nwk", append = TRUE)
write_file("\n", "~/_data/phylo/mcmctree_in_tree.nwk", append = TRUE)





# ((((((((Chironomus_tentans: 0.225825, (Chironomus_riparius: 0.191673, Chironomus_tepperi: 0.191673) [&95%={0.167724, 0.231869}]: 0.034152) [&95%={0.197619, 0.273171}]: 0.737943, (Polypedilum_vanderplanki: 0.234349, Polypedilum_pembai: 0.234349) [&95%={0.205083, 0.283425}]: 0.729419) [&95%={0.843625, 1.16573}]: 0.169057, Tanytarsus_gracilentus: 1.132824) [&95%={0.991739, 1.37068}]: 0.361535, (Clunio_marinus: 1.010433, Belgica_antarctica: 1.010433) [&95%={0.884667, 1.2226}]: 0.483926) [&95%={1.30827, 1.80844}]: 0.058194, Propsilocerus_akamusi: 1.552554) [&95%={1.35938, 1.87877}]: 0.697730, Parochlus_steinenii: 2.250284) [&95%={1.97052, 2.72394}]: 0.670012, Culicoides_sonorensis: 2.920296) [&95%={2.55699, 3.53357}]: 0.000076, Anopheles_stephensi: 2.920372) [&95%={2.55702, 3.53366}];


read_lines("~/_data/phylo/chir_mcmctree/mcmc_1/FigTree.tre") %>%
    keep(~ grepl("^\tUTREE", .x))


mcmc <- read_tsv("~/_data/phylo/chir_mcmctree/mcmc_1/chir_mcmctree_mcmc.txt")


mcmc %>%
    select(starts_with("t_n")) %>%
    summarize_all(list(mean = mean, median = median,
                       low = ~ quantile(.x, 0.025),
                       high = ~ quantile(.x, 0.975))) %>%
    pivot_longer(everything()) %>%
    mutate(variable = str_split(name, "_") %>%
               map_chr(~ paste(.x[[1]], .x[[2]], sep = "_")),
           measure = str_split(name, "_") %>%
               map_chr(~ .x[[3]])) %>%
    select(variable, measure, value) %>%
    pivot_wider(names_from = measure, values_from = value) %>%
    mutate(summary = sprintf("%.4f %.4f (%.4f, %.4f)", mean, median, low, high)) %>%
    select(variable, summary)


mcmc %>%
    filter(lnL == max(lnL)) %>%
    .[["t_n13"]]


mcmc %>%
    select(Gen, starts_with("t_n")) %>%
    pivot_longer(starts_with("t_n")) %>%
    ggplot(aes(Gen, value, color = name)) +
    geom_line() +
    scale_color_viridis_d(guide = "none") +
    facet_wrap(~ name, ncol = 3, scales = "free_y")




#' Below looks for differences in credible intervals among MCMCtree runs

mcmc_ests <- map(1:4, function(i) {
    fn <- paste0("~/_data/phylo/chir_mcmctree/mcmc_", i, "/chir_mcmctree.out")
    lns <- read_lines(fn)
    lns <- lns[(which(grepl("^Posterior mean", lns))+2):length(lns)]
    cis <- str_split(lns, "\\s+") %>%
        map_dbl(~ as.numeric(.x[[2]]))
    return(cis)
})

r2 <- numeric(0)
for (i in 1:3) {
    for (j in (i+1):4) {
        x <- mcmc_ests[[i]]
        y <- mcmc_ests[[j]]
        r2_ij <- cor(x, y)
        r2 <- c(r2, r2_ij)
    }
}; rm(i, j, x, y, r2_ij)



#' I'm calculating the difference between the max and min values
#' for low and high interval estimates relative to the mean value:
map(1:nrow(mcmc_cis[[1]]),
       function(i) {
           # z <- function(x) diff(range(x))
           z <- function(x) diff(range(x)) / abs(mean(x))
           a <- map_dbl(1:4, ~ mcmc_cis[[.x]][i,1]) %>%
               z()
           b <- map_dbl(1:4, ~ mcmc_cis[[.x]][i,2]) %>%
               z()
           return(c(a, b))
       }) %>%
    do.call(what = rbind)

#'
#' This is what was output on 10 July 2022:
#'
#' ```
#'               [,1]         [,2]
#'  [1,] 3.016680e-03 4.227490e-03
#'  [2,] 3.055917e-03 4.199207e-03
#'  [3,] 3.254099e-03 4.676166e-03
#'  [4,] 3.021760e-03 4.536963e-03
#'  [5,] 2.986446e-03 4.214320e-03
#'  [6,] 3.131471e-03 3.949823e-03
#'  [7,] 2.849426e-03 3.868555e-03
#'  [8,] 3.041440e-03 3.669052e-03
#'  [9,] 2.986412e-03 3.890630e-03
#' [10,] 3.419638e-03 3.535130e-03
#' [11,] 3.510659e-03 4.838743e-03
#' [12,] 4.039384e-03 2.192982e-03
#' [13,] 2.322206e-02 7.848523e-03
#' [14,] 1.471054e-06 6.347454e-07
#' ```
#'
#' This looks pretty consistent.
#'


rbind(c(3.016680e-03, 4.227490e-03),
      c(3.055917e-03, 4.199207e-03),
      c(3.254099e-03, 4.676166e-03),
      c(3.021760e-03, 4.536963e-03),
      c(2.986446e-03, 4.214320e-03),
      c(3.131471e-03, 3.949823e-03),
      c(2.849426e-03, 3.868555e-03),
      c(3.041440e-03, 3.669052e-03),
      c(2.986412e-03, 3.890630e-03),
      c(3.419638e-03, 3.535130e-03),
      c(3.510659e-03, 4.838743e-03),
      c(4.039384e-03, 2.192982e-03),
      c(2.322206e-02, 7.848523e-03),
      c(1.471054e-06, 6.347454e-07)) %>%
    mean()

