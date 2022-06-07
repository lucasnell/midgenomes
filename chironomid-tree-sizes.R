
#'
#' Phylogeny of family Chironomidae by species alongside those with estimated
#' genome sizes.
#'


library(ape)
library(tidyverse)
library(ggtree)


# library(grid)
# library(phytools)
# library(gtable)
# library(ggstance) # for horizontal versions of geoms

if (file.exists(".Rprofile")) source(".Rprofile")





# -------------*
# flow cytometry ----
# -------------*
# Mb per C-value (pg)
mb_per_pg <- 978

cyto_mb <- tibble(# From Cornette et al. (2015)
    spp = c(
        "Polypedilum vanderplanki", "Polypedilum tamanigrum", "Polypedilum nubifer",
        "Paratanytarsus grimmii", "Sergentia kizakiensis", "Tanytarsus takahashii",
        "Paraborniella tonnoiri", "Stictochironomus akizukii",
        "Chironomus nippodorsalis", "Chironomus sulfurosus", "Chironomus salinarius",
        "Chironomus plumosus", "Chironomus fusciceps", "Chironomus acerbiphilus",
        "Glyptotendipes tokunagai", "Dicrotendipes pelochloris"),
    pg = c(0.10, 0.11, 0.10, 0.10, 0.10, 0.09, 0.20, 0.13, 0.12, 0.15, 0.16,
           0.15, 0.16, 0.17, 0.16, 0.12)) %>%
    mutate(mb = pg * mb_per_pg) %>%
    select(-pg) %>%
    # From doi: 10.1007/s00427-009-0281-0:
    add_row(spp = "Chironomus riparius", mb = (196.2 + 194.3) / 2) %>%
    mutate(method = "cyto")


# -------------*
# assemblies ----
# -------------*

genome_mb <- list(list("Tanytarsus gracilentus", 184 / (1 - 0.098)), # My newest assembly,
                                        # adjusted for incompleteness measured by BUSCO
                                        # NOTE: this BUSCO info is for the older assembly
                                        # version, so this is just a rough idea.
        # From doi: 10.1038/s41598-019-41549-8 (original genome version
        # doi: 10.1093/gigascience/giw009)
        list("Parochlus steinenii", 145),
        # From doi: 10.1038/ncomms5784
        list("Polypedilum vanderplanki", 104),
        list("Polypedilum nubifer", 107),
        # From doi: 10.1038/ncomms5611
        list("Belgica antarctica", 99),
        # From doi: 10.1534/g3.119.400710
        list("Chironomus riparius", 178)) %>%
    map_dfr(~ tibble(spp = .x[[1]], mb = .x[[2]])) %>%
    mutate(method = "genome")



mb_df <- bind_rows(cyto_mb, genome_mb) %>%
    spread(method, mb) %>%
    mutate(species = spp,
           spp = str_split(spp, " "),
           genus = map_chr(spp, ~ .x[[1]])) %>%
    select(-spp)


mb_df %>%
    group_by(genus) %>%
    summarize(N = n())


# -------------*
# Chironomid genus tree:
# -------------*

tree <- read.tree("chironomidae_genus.nwk")
tree <- drop.tip(tree, tree$tip.label[tree$tip.label != "Buchonomyia" &
                                          !tree$tip.label %in% mb_df$genus])


for (.df in split(mb_df, mb_df$genus)) {

    .gen <- paste(.df$genus[[1]])
    .spp <- paste(.df$species)
    tree$tip.label[tree$tip.label == .gen] <- .spp[[1]]
    tree$tip.label <- gsub("_", " ", tree$tip.label)

    if (length(.spp) > 1) {
        for (i in 2:length(.spp)) {
            j <- which(tree$tip.label == .spp[[1]])[1]
            tree <- bind.tip(tree, .spp[i], edge.length = 0.5,
                             where = j, position = 0)
            tree$tip.label <- gsub("_", " ", tree$tip.label)
        }
    }

}; rm(.df, .gen, .spp)


plot(tree, cex = 0.5, no.margin = TRUE)






g_df <- mb_df %>%
    gather("method", "mb", cyto:genome) %>%
    filter(!is.na(mb)) %>%
    select(-genus) %>%
    rename(taxon = species) %>%
    add_row(taxon = "Buchonomyia", mb = NA)

p <- ggtree(tree) +
    geom_tiplab(size = 3)

# p2 <- facet_plot(p + xlim_tree(350), "Genome size (Mb)", g_df, geom_segment,
#                  aes(x = 0, xend = mb, y = y, yend = y, color = method),
#                  size = 3, na.rm = TRUE) +
p2 <- facet_plot(p + xlim_tree(350), "Genome size (Mb)", g_df, geom = geom_colh,
                 mapping = aes(x = mb, fill = method),
                 position = "dodge",
                 width = 0.5,
                 size = 3, na.rm = TRUE) +
    scale_fill_manual(values = c("dodgerblue", "firebrick"))

fp <- p2 +
    theme_tree2() +
    theme(strip.background = element_blank()) +
    NULL


gt = ggplot_gtable(ggplot_build(fp))
# gtable_show_layout(gt) # will show you the layout - very handy function
# gt # see plot layout in table format
gt$layout$l[grep('panel-2', gt$layout$name)] # you want to find the column specific to panel-2
gt$widths[7] = 0.4 * gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
grid.draw(gt) # plot with grid draw


cairo_pdf("figs/tree_genomes.pdf", width = 6, height = 4)
grid.draw(gt)
dev.off()
