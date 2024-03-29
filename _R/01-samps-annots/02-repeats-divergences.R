

#'
#' Make repeat divergence landscape plots.
#'

source("_R/00-preamble.R")

library(ggtree)
library(patchwork)


if (!file.exists("_data/gsizes-list.rds")) {
    # Takes ~1 min
    gsizes <- read_csv("_data/species-names-families.csv", col_types = "ccc") |>
        getElement("spp_abbrev") |>
        set_names() |>
        safe_mclapply(\(s) {
            .fn <- paste0(dirs$assembly, "/", s, "_assembly.fasta.gz")
            sl <- read_lines(.fn, progress = FALSE)
            hl <- which(grepl("^>", sl))
            ll <- c(hl[-1] - 1L, length(sl))
            stopifnot(length(hl) == length(ll))
            lens <- map_int(1:length(hl), \(i) sum(nchar(sl[(hl[i]+1L):(ll[i])])))
            return(sum(lens))
        })
    write_rds(gsizes, "_data/gsizes-list.rds")
} else {
    gsizes <- read_rds("_data/gsizes-list.rds")
}




simplify_class <- function(.class, combine_nonTE = TRUE) {
    .class <- str_split(.class, "\\/") |>
        map_chr(\(x) x[1])
    # Remove '?' because the RepeatMasker's *.tbl summary includes
    # these with the class minus the '?' modifier
    .class <- str_remove_all(.class, "\\?")
    # rRNA, tRNA, and snRNA are lumped as "Small RNA"
    .class <- ifelse(.class %in% c("rRNA", "tRNA", "snRNA"), "Small_RNA", .class)
    # Unclassified is Unknown or ARTEFACT:
    .class <- ifelse(.class %in% c("Unknown", "ARTEFACT"), "Unclassified", .class)
    if (combine_nonTE) {
        # Combine all non-TE (but known) classes into "non_TE"
        .class <- ifelse(.class %in% nonTE_classes, "non_TE", .class)
    }
    return(.class)
}


#' Make names of genome features and repeat element classes pretty for plotting:
pretty_repeats <- function(yucks, to_fct = FALSE, units = FALSE) {
    pretty_map <- list("gsize" = "Genome size",
                       "n_genes" = "Protein-coding genes",
                       "sum_interg_len" = "Intergenic DNA",
                       "mean_intron_len" = "Intron length",
                       "SINE" = "SINE",
                       "LINE" = "LINE",
                       "LTR" = "LTR",
                       "DNA" = "DNA",
                       "non_TE" = "non-TE",
                       "Unclassified" = "Unclassified")
    yucks <- gsub("log_", "", yucks)
    if (any(! yucks %in% names(pretty_map))) {
        stop("These inputs are not in the pretty map: ",
             paste(yucks[! yucks %in% names(pretty_map)], collapse = ", "))
    }
    pretties <- map_chr(yucks, \(x) pretty_map[[x]])
    if (units) {
        lgl <- pretties != "Protein-coding genes"
        pretties[lgl] <- paste(pretties[lgl], "(bp)")
    }
    if (to_fct) {
        idx <- map_int(gsub(" \\(bp\\)", "", unique(pretties)),
                       \(x) which(paste(pretty_map) == x))
        pretties <- factor(pretties, levels = unique(pretties)[rank(idx)])
    }
    return(pretties)
}



one_spp_ds <- function(.spp) {
    # .spp = "Tgraci"

    ds_lines <- paste0(dirs$repeats, "/repeats_diverg/divsum/",
                       .spp, ".divsum") |>
        read_lines()
    ds_lines <- ds_lines[ds_lines != ""]

    tbl_inds <- c("^Weighted average Kimura divergence for each repeat family",
                  "^Coverage for each repeat class and divergence") |>
        map_int(\(s) which(grepl(s, ds_lines)) + 1L)

    cn <- str_split(ds_lines[tbl_inds[2]], "\\s+")[[1]]
    cn[cn == ""] <- paste0("ZZZZ_", 1:sum(cn == ""))

    ds_df <- ds_lines[(tbl_inds[2]+1L):length(ds_lines)] |>
        read_table(col_types = cols(), col_names = cn) |>
        select(-starts_with("ZZZ")) |>
        pivot_longer(-Div, names_to = "class", values_to = "bp") |>
        mutate(class = simplify_class(class)) |>
        group_by(Div, class) |>
        summarize(bp = sum(bp), .groups = "drop") |>
        mutate(perc = bp / gsizes[[.spp]] * 100) |>
        filter(class != "Unclassified") |>
        mutate(plot_class = pretty_repeats(class, to_fct = TRUE),
               spp_abbrev = .spp)

    return(ds_df)
}


ds_df <- names(gsizes) |>
    map_dfr(one_spp_ds) |>
    mutate(spp_abbrev = factor(spp_abbrev, levels = names(gsizes))) |>
    arrange(spp_abbrev) |>
    mutate(species = expand_spp(spp_abbrev, to_fct = TRUE))


gsize_p <- gsizes |>
    (\(x) {
        tibble(spp_abbrev = factor(names(x), levels = names(x)),
               gsize = log10(as.numeric(x)))
    })() |>
    arrange(spp_abbrev) |>
    mutate(species = expand_spp(spp_abbrev, to_fct = TRUE)) |>
    mutate(species = factor(species, levels = rev(levels(species)))) |>
    ggplot() +
    geom_segment(aes(x = species, y = gsize - 6, xend = species, yend = 0),
                 linewidth = 1, linetype = "22", color = "gray60") +
    geom_point(aes(species, gsize - 6), size = 4) +
    scale_y_continuous("Genome size\n(Mb)", breaks = log10(100 * 2^(0:3)),
                       labels = 100 * 2^(0:3)) +
    scale_x_discrete(drop = FALSE) +
    coord_flip(ylim = c(1.92, NA)) +
    theme(panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          # axis.text.x = element_text(size = 8),
          # axis.title.x = element_text(size = 9),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())



divergence_panels_p <- ds_df |>
    ggplot(aes(Div, perc, fill = plot_class)) +
    geom_bar(stat = "identity", color = NA) +
    scale_fill_viridis_d(NULL, begin = 0.0) +
    xlab("Sequence divergence\n(Kimura distance)") +
    scale_y_continuous("Percent of genome", position = "right", n.breaks = 3) +
    facet_wrap(~ species, ncol = 1, scales = "free_y") +
    theme(strip.background = element_blank(),
          strip.text = element_blank())




# color lines by whether it's Chironomidae:
tree_p <- read.tree("_data/phylo/time-tree.nwk") |>
    (\(x) {
        x$tip.label <- expand_spp(x$tip.label) |>
            str_split(" ") |>
            map_chr(\(x) paste0(str_sub(x[1], end = 1), ". ", x[2]))
        return(x)
    })() |>
    ggtree(linewidth = 1, lineend = "square") |>
    groupClade(.node = 20) +
    aes(color = group) +
    geom_tiplab(size = 8 / 2.83465, fontface = "italic", offset = 0.2) +
    theme(legend.position = "none") +
    scale_color_manual(values = c("gray70", "black")) +
    coord_cartesian(xlim = c(0, 7))



rep_div_p <- tree_p + gsize_p + divergence_panels_p +
    plot_layout(nrow = 1, widths = c(1, 1, 1.5))


# save_plot("repeat-divergences", rep_div_p, w = 7.5, h = 8)
