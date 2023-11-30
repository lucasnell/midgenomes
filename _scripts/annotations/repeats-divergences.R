

#'
#' Make repeat divergence landscape plots.
#'

source("_scripts/00-preamble.R")

library(ggtree)
library(patchwork)


gsizes <- "_data/genome-stats.csv" |>
        read_csv(col_types = cols()) |>
        select(spp_abbrev, gsize) |>
        (\(x) {
            z <- as.list(x[["gsize"]])
            names(z) <- x[["spp_abbrev"]]
            return(z)
        })()


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
        mutate(plot_class = pretty$convert(class, to_fct = TRUE),
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
