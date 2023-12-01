
#'
#' Shared code for all or most scripts in this directory.
#'


#'
#' Update these on your system (do NOT include trailing '/'):
#'
dirs <- new.env()
#' You should only need to change this one:
dirs$parent <- "~/_data/__to-deposit"
#' ^^^^^^^^
#'

dirs$raxml_supp <- paste0(dirs$parent, "/phylo/chir_raxml_supp")
dirs$mcmctree <- paste0(dirs$parent, "/phylo/chir_mcmctree")
dirs$cafe <- paste0(dirs$parent, "/chir_cafe")
dirs$orthofinder <- paste0(dirs$parent, "/chir_orthofinder/orthofinder-output")
dirs$orthofinder_extr <- paste0(dirs$parent, "/orthofinder-extraction")
dirs$hyphy_busted <- paste0(dirs$parent, "/chir_hyphy_busted")
dirs$hyphy_relax <- paste0(dirs$parent, "/chir_hyphy_relax")

dirs$assembly <- paste0(dirs$parent, "/assemblies")
dirs$go <- paste0(dirs$parent, "/go-terms")
dirs$proteins <- paste0(dirs$parent, "/proteins")
dirs$features <- paste0(dirs$parent, "/features")
dirs$repeats <- paste0(dirs$parent, "/repeats")



suppressPackageStartupMessages({
    library(viridisLite)
    library(ape)
    library(parallel)
    library(tidyverse)
})

#'
#' Alphabetic list of other packages needed:
#'
#' - AnnotationDbi
#' - clusterProfiler
#' - future
#' - future.apply
#' - ggtext
#' - ggtree
#' - GO.db
#' - grid
#' - gt
#' - jsonlite
#' - knitr
#' - org.Dm.eg.db
#' - patchwork
#' - phylolm
#' - phyr
#' - rrvgo
#' - treeio
#' - treemap
#'


options(mc.cores = max(1L, detectCores()-2L))


if (file.exists(".Rprofile")) source(".Rprofile")


theme_set(theme_classic() +
              theme(strip.background = element_blank()))


#' Palette for chironomid species:
spp_pal <- full_spp_pal <- turbo(100)[c(70+3*0:8, 60, 15+4*2:0, 30)] |> as.list()
names(spp_pal) <- read_csv("_data/species-names-families.csv", col_types = cols(),
                           progress = FALSE)[["spp_abbrev"]]
names(full_spp_pal) <- read_csv("_data/species-names-families.csv", col_types = cols(),
                                progress = FALSE)[["species"]]


#' non-TE (but known) element classes:
nonTE_classes <- c("RC", "Small_RNA", "Satellite", "Simple_repeat",
                   "Low_complexity")



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


#' Safe call to `mcmapply` that reverts to `mapply` on Windows:
safe_mcmapply <- function(FUN, ..., MoreArgs = NULL) {
    if (.Platform$OS.type == "unix") {
        out <- mcmapply(FUN, ..., MoreArgs = MoreArgs,
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    } else {
        out <- mapply(FUN, ..., MoreArgs = MoreArgs,
                      SIMPLIFY = FALSE, USE.NAMES = FALSE)
    }
    return(out)
}
#' Safe call to `mclapply` that reverts to `lapply` on Windows:
safe_mclapply <- function(X, FUN, ...) {
    if (.Platform$OS.type == "unix") {
        out <- mclapply(X, FUN, ...)
    } else {
        out <- lapply(X, FUN, ...)
    }
    return(out)
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





