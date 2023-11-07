
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

spp_pal <- full_spp_pal <- turbo(100)[c(70+3*0:8, 60, 15+4*2:0, 30)] |> as.list()
names(spp_pal) <- read_csv("_data/species-names-families.csv", col_types = cols(),
                           progress = FALSE)[["spp_abbrev"]]
names(full_spp_pal) <- read_csv("_data/species-names-families.csv", col_types = cols(),
                                progress = FALSE)[["species"]]


#' Make names of genome features and repeat element classes pretty for plotting:
pretty <- new.env()
pretty$map <- list("gsize" = "Genome size",
                   "n_genes" = "Protein-coding genes",
                   "sum_interg_len" = "Intergenic DNA",
                   "mean_intron_len" = "Intron length",
                   "SINE" = "SINE",
                   "LINE" = "LINE",
                   "LTR" = "LTR",
                   "DNA" = "DNA",
                   "non_TE" = "non-TE",
                   "Unclassified" = "Unclassified")
pretty$convert <- function(yucks, to_fct = FALSE, units = FALSE) {
    yucks <- gsub("log_", "", yucks)
    if (any(! yucks %in% names(pretty$map))) {
        stop("These inputs are not in the pretty map: ",
             paste(yucks[! yucks %in% names(pretty$map)], collapse = ", "))
    }
    pretties <- map_chr(yucks, \(x) pretty$map[[x]])
    if (units) {
        lgl <- pretties != "Protein-coding genes"
        pretties[lgl] <- paste(pretties[lgl], "(bp)")
    }
    if (to_fct) {
        idx <- map_int(gsub(" \\(bp\\)", "", unique(pretties)),
                       \(x) which(paste(pretty$map) == x))
        pretties <- factor(pretties, levels = unique(pretties)[rank(idx)])
    }
    return(pretties)
}

# non-TE (but known) element classes:
nonTE_classes <- c("RC", "Small_RNA", "Satellite", "Simple_repeat",
                   "Low_complexity")
# Repeat classes to keep (in order of plotting):
repeat_classes <- c("SINE", "LINE", "LTR", "DNA", "non_TE", "Unclassified")


save_plot <- function(n, p, w, h, .pdf = TRUE, .png = TRUE, ...) {
    stopifnot(is.character(n) && length(n) == 1)
    stopifnot(is.ggplot(p) || is.function(p))
    stopifnot(is.numeric(w) && length(w) == 1)
    stopifnot(is.numeric(h) && length(h) == 1)
    stopifnot(is.logical(.pdf) && length(.pdf) == 1)
    stopifnot(is.logical(.png) && length(.png) == 1)
    old_warn <- getOption("warn")
    options(warn = -1)
    if (.pdf) {
        cairo_pdf(sprintf("_figures/%s.pdf", n), width = w, height = h,
                  bg = NA, ...)
        if (is.function(p)) {
            p()
        } else {
            plot(p)
        }
        tmp <- dev.off()
    }
    if (.png) {
        .fn <- sprintf("_figures/%s.png", n)
        if (is.function(p)) {
            png(filename = .fn, width = w, height = h, units = "in",
                bg = NA, res = 300)
            p()
            tmp <- dev.off()
        } else {
            ggsave(.fn, p, width = w, height = h)
        }
    }
    options(warn = old_warn)
    invisible(NULL)
}

#' Functions to expand and abbreviate species names.
#' Because they both involve reading `_data/species-names-families.csv`,
#' they're best used on vectors.
expand_spp <- function(x, to_fct = FALSE) {
    map_list <- read_csv("_data/species-names-families.csv", col_types = cols(),
                         progress = FALSE) |>
        (\(y) {z <- as.list(y$species); names(z) <- y$spp_abbrev; return(z)})()
    z <- map_chr(x, \(.x) map_list[[.x]])
    if (to_fct) {
        lvls <- paste(map_list)[paste(map_list) %in% z]
        z <- factor(z, levels = lvls)
    }
    return(z)
}
abbrev_spp <- function(x, to_fct = FALSE) {
    map_list <- read_csv("_data/species-names-families.csv", col_types = cols(),
                         progress = FALSE) |>
        (\(y) {z <- as.list(y$spp_abbrev); names(z) <- y$species; return(z)})()
    z <- map_chr(x, \(.x) map_list[[.x]])
    if (to_fct) {
        lvls <- paste(map_list)[paste(map_list) %in% z]
        z <- factor(z, levels = lvls)
    }
    return(z)
}





