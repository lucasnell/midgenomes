library(tidyverse)
library(here)

divsum_dir <- "~/_data/_repeats/repeats_diverg/divsum"

species <- list.files(divsum_dir) |> str_remove("\\.divsum")

theme_set(theme_minimal())

gsizes <- here("_data/genome-stats.csv") |>
    read_csv(col_types = cols()) |>
    select(species, gsize) |>
    (\(x) {
        z <- as.list(x[["gsize"]])
        names(z) <- x[["species"]]
        return(z)
    })()

#' Repeat element classes in order with formatted names for plotting:
class_map <- list("SINE" = "SINEs",
                  "LINE" = "LINEs",
                  "LTR" = "LTR elements",
                  "DNA" = "DNA transposons",
                  "RC" = "Rolling-circles",
                  "Small_RNA" = "Small RNA",
                  "Satellite" = "Satellites",
                  "Simple_repeat" = "Simple repeats",
                  "Low_complexity" = "Low complexity")

simplify_class <- function(.class) {
    .class <- str_split(.class, "\\/") |>
        map_chr(\(x) x[1])
    # Remove '?' because the RepeatMasker's *.tbl summary includes
    # these with the class minus the '?' modifier
    .class <- str_remove_all(.class, "\\?")
    # rRNA, tRNA, and snRNA are lumped as "Small RNA"
    .class <- ifelse(.class %in% c("rRNA", "tRNA", "snRNA"), "Small_RNA", .class)
    # Unclassified is Unknown or ARTEFACT:
    .class <- ifelse(.class %in% c("Unknown", "ARTEFACT"), "Unclassified", .class)
    return(.class)
}


one_spp_ds <- function(.spp) {
    # .spp = "Tgraci"

    ds_lines <- paste0(divsum_dir, "/", .spp, ".divsum") |>
        read_lines()
    ds_lines <- ds_lines[ds_lines != ""]

    tbl_inds <- c("^Weighted average Kimura divergence for each repeat family",
                  "^Coverage for each repeat class and divergence") |>
        map_int(\(s) which(grepl(s, ds_lines)) + 1L)

    ds_lines[tbl_inds[1]:(tbl_inds[2]-2L)] |>
        I() |>
        read_tsv(col_types = cols()) |>
        filter(!grepl("--", `Kimura%`))

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
        mutate(plot_class = map_chr(class, \(x) class_map[[x]]) |>
                   factor(levels = unlist(class_map)),
               species = .spp)

    return(ds_df)
}


ds_df <- species |>
    map_dfr(one_spp_ds) |>
    mutate(species = factor(species,
                            levels = c("Ctenta", "Cripar", "Pvande", "Ppemba",
                                       "Tgraci", "Bantar", "Cmarin", "Pakamu",
                                       "Pstein", "Csonor", "Cquinq", "Aaegyp",
                                       "Asteph", "Mdomes")))


ds_df |>
    filter(Div <= 50) |>
    # Remove non-TE elements:
    filter(!class %in% c("RC", "Small_RNA", "Satellite", "Simple_repeat")) |>
    ggplot(aes(Div, perc, fill = plot_class)) +
    geom_bar(stat = "identity") +
    # scale_fill_brewer(NULL, palette = "Dark2") +
    scale_fill_viridis_d(NULL, begin = 0.0) +
    xlab("Sequence divergence (Kimura distance)") +
    ylab("Percent of genome") +
    facet_wrap(~ species, scales = "free_y")

