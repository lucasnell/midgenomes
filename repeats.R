library(tidyverse)

one_spp_repeats <- function(.spp) {
    .fn <- sprintf("~/_data/_repeats/%s_repeats/masker/%s_assembly.fasta.out",
                   .spp, .spp)
    if (!file.exists(.fn)) .fn <- gsub("_assembly.fasta", "_contigs.fasta", .fn)
    repeats <- read_table(.fn,
                          skip = 3, col_types = cols(), progress = FALSE,
                          col_names = c("SW", "p_div", "p_del", "p_ins", "query_seq",
                                        "begin", "end", "left", "unk1", "matching_repeat",
                                        "class", "begin2", "end2", "left2", "id", "star"))
    repeat_sums <- repeats |>
        mutate(class = str_split(class, "\\/") |> map_chr(\(x) x[1]),
               # rRNA and tRNA are lumped as "Small RNA"
               class = ifelse(class %in% c("rRNA", "tRNA"), "Small_RNA", class)) |>
        mutate(species = .spp) |>
        group_by(species, class) |>
        summarize(elements = length(unique(id)),
                  length = sum((end - begin + 1) * (1 - (p_del+p_ins) / 100)),
                  .groups = "drop") |>
        arrange(class)
    cat(.spp, "finished\n")
    return(repeat_sums)
}


repeats_df <- c("Aaegyp", "Asteph", "Bantar", "Cquinq", "Cripar",
                "Csonor", "Ctenta", "Mdomes", "Pakamu", "Ppemba",
                "Pstein", "Pvande", "Tgraci") |>
    map_dfr(one_spp_repeats)


write_csv(repeats_df, "_data/repeats-summary.csv")




