
library(tidyverse)
library(ape)
library(ggtree)
library(viridis)

if (file.exists(".Rprofile")) source(".Rprofile")

save_plot <- function(n, p, w, h, ...) {
    cairo_pdf(sprintf("~/Box Sync/midgenomics/genome-report/%s.pdf", n),
              width = w, height = h, bg = NA, ...)
    plot(p)
    dev.off()
}

#' Previously sequenced species plus Tanytarsus, in the order genera
#' appear in the tree from top to bottom:
prev_spp <- c("Chironomus riparius", "Chironomus tentans",
              "Polypedilum pembai", "Polypedilum vanderplanki",
              "Tanytarsus gracilentus",
              "Belgica antarctica", "Clunio marinus",
              "Propsilocerus akamusi")
#' Unique genera
prev_gen <- str_split(prev_spp, " ") |>
    map_chr(~ .x[[1]]) |>
    unique()

#' Phylogeny of family Chironomidae to genus level from timetree.org
#' downloaded on 6 June 2022
chir_tr <- read.tree("chironomidae_genus.nwk")


d <- tibble(label = chir_tr$tip.label,
            trait = factor(label, levels = prev_gen))


# Note: restart R session if tip labels don't show up. Not sure what's going on.


cols <- turbo(101)[c(15, 80, 5, 25, 90, 35)]
# Make Tanytarsus stand out a bit:
cols[3] <- "gold"

chir_tr_p <- full_join(chir_tr, d, by = "label") |>
    ggtree() +
    geom_tippoint(aes(color = trait), size = 1) +
    geom_tiplab(aes(color = trait), size = 3, offset = 1, fontface = "bold") +
    theme_tree2() +
    scale_color_manual(values = cols, na.value = NA, guide = "none") +
    scale_x_continuous("Time (MYA)", limits = c(0, 250))

save_plot("chir_tree", chir_tr_p, 6, 6)

