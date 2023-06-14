library(tidyverse)


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

one_spp_repeats <- function(.spp) {
    .fn <- sprintf("~/_data/_repeats/%s_repeats.tsv", .spp)
    repeats <- read_table(.fn,
                          skip = 3, col_types = cols(), progress = FALSE,
                          col_names = c("SW", "p_div", "p_del", "p_ins", "query_seq",
                                        "begin", "end", "left", "unk1", "matching_repeat",
                                        "class", "begin2", "end2", "left2", "id", "star"))
    repeat_sums <- repeats |>
        mutate(class = simplify_class(class)) |>
        mutate(species = .spp) |>
        group_by(species, class) |>
        summarize(elements = length(unique(id)),
                  length = sum((end - begin + 1) * (1 - (p_del+p_ins) / 100)),
                  .groups = "drop") |>
        arrange(class)
    cat(.spp, "finished\n")
    return(repeat_sums)
}




# (The warnings about Mdomes are not an issue)
repeats_df <- c("Aaegyp", "Asteph", "Bantar", "Cquinq", "Cripar",
                "Csonor", "Ctenta", "Mdomes", "Pakamu", "Ppemba",
                "Pstein", "Pvande", "Tgraci") |>
    map_dfr(one_spp_repeats)


repeats_df |>
    group_by(class) |>
    summarize(n_species = length(unique(species)),
              species = paste0(sort(unique(species)), collapse = " "))


write_csv(repeats_df, "_data/repeats-summary.csv")

