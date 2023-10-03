
#'
#' Shared code for all or most scripts in this directory.
#'


#'
#' Update these on your system (do NOT include trailing '/'):
#'
dirs <- new.env()
dirs$raxml_supp <- "~/_data/_phylo/chir_raxml_supp"
dirs$mcmctree <- "~/_data/_phylo/chir_mcmctree"
dirs$cafe <- "~/_data/chir_cafe"
dirs$orthofinder <- "~/_data/chir_orthofinder/orthofinder-output"
dirs$orthofinder_extr <- "~/_data/_orthofinder-extraction"
dirs$go <- "~/_data/_go-terms"
dirs$assembly <- "~/_data/_assemblies"
dirs$proteins <- "~/_data/_proteins"
dirs$features <- "~/_data/_features"
dirs$repeats <- "~/_data/_repeats"
dirs$hyphy_busted <- "~/_data/chir_hyphy_busted"
dirs$hyphy_relax <- "~/_data/chir_hyphy_relax"




library(tidyverse)
library(viridisLite)
library(ape)
library(parallel)

options(mc.cores = max(1L, detectCores()-2L))


if (file.exists(".Rprofile")) source(".Rprofile")


theme_set(theme_classic() +
              theme(strip.background = element_blank()))


save_plot <- function(n, p, w, h, .pdf = TRUE, .png = TRUE, ...) {
    stopifnot(is.character(n) && length(n) == 1)
    stopifnot(is.ggplot(p))
    stopifnot(is.numeric(w) && length(w) == 1)
    stopifnot(is.numeric(h) && length(h) == 1)
    stopifnot(is.logical(.pdf) && length(.pdf) == 1)
    stopifnot(is.logical(.png) && length(.png) == 1)
    old_warn <- getOption("warn")
    options(warn = -1)
    if (.pdf) {
        cairo_pdf(sprintf("_figures/%s.pdf", n), width = w, height = h,
                  bg = NA, ...)
        plot(p)
        tmp <- dev.off()
    }
    if (.png) {
        ggsave(sprintf("_figures/%s.png", n), p, width = w, height = h)
    }
    options(warn = old_warn)
    invisible(NULL)
}

#' Functions to expand and abbreviate species names.
#' Because they both involve reading `_data/species-names-families.csv`,
#' they're best used on vectors.
expand_spp <- function(x) {
    map_list <- read_csv("_data/species-names-families.csv", col_types = cols(),
                         progress = FALSE) |>
        (\(y) {z <- as.list(y$species); names(z) <- y$spp_abbrev; return(z)})()
    map_chr(x, \(.x) map_list[[.x]])
}
abbrev_spp <- function(x) {
    map_list <- read_csv("_data/species-names-families.csv", col_types = cols(),
                         progress = FALSE) |>
        (\(y) {z <- as.list(y$spp_abbrev); names(z) <- y$species; return(z)})()
    map_chr(x, \(.x) map_list[[.x]])
}





