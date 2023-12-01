
#'
#' Summarize output from mantis
#'


source("_R/00-preamble.R")


#' -------------------------------------------------
#' Count number of proteins predicted by each database.
#' -------------------------------------------------


#' Function to properly label databases:
lab_dbs <- function(x) {
    case_when(
        grepl("^Pfam", x) ~ "Pfam",
        grepl("^kofam", x) ~ "KOfam",
        grepl("^NOG", x) ~ "eggNOG",
        grepl("^NCBI", x) ~ "NPFM",
        grepl("^tcdb", x) ~ "TCDB",
        TRUE ~ NA_character_)
}


map_dfr(c("Cripar", "Csonor", "Pstein", "Tgraci"),
        \(spp) {
            anno_file <- paste0(dirs$go, "/zzz-mantis-full/", spp,
                                "_mantis/output_annotation.tsv")
            anno_df <- read_tsv(anno_file, col_types = cols())
            # Total number of proteins annotated:
            total_prots <- anno_df[["Query"]] |> unique() |> length()
            anno_df |>
                mutate(db = lab_dbs(Ref_file)) |>
                group_by(db) |>
                summarize(n_proteins= length(unique(Query))) |>
                add_row(db = "total", n_proteins = total_prots) |>
                mutate(species = spp)
            }) |>
    pivot_wider(names_from = species, values_from = n_proteins)

# # A tibble: 6 × 5
#   db     Cripar Csonor Pstein Tgraci
#   <chr>   <int>  <int>  <int>  <int>
# 1 KOfam   13011  14071  12592  12161
# 2 NPFM     4039   4276   3851   3772
# 3 Pfam    12951  14011  12516  12090
# 4 TCDB     1753   1980   1756   1543
# 5 eggNOG  10736  12692  11024  10182
# 6 total   14304  15829  14003  13424



#' -------------------------------------------------
#' Count number of proteins with various tags.
#' -------------------------------------------------

#' (Doing this by lines bc the number of columns isn't the same among rows)
#'
#' Tags I'm interested in:
tag_interest <- c("kegg_ko", "go", "cog", "enzyme_ec", "kegg_pathway")


map_dfr(c("Cripar", "Csonor", "Pstein", "Tgraci"),
        \(spp) {
            anno_file <- paste0(dirs$go, "/zzz-mantis-full/", spp,
                                "_mantis/consensus_annotation.tsv")
            anno_tags <- read_lines(anno_file) |>
                str_split("\t") |>
                base::`[`(-1) |>
                map(function(x) {
                    x[-1:-6] |>
                        str_split(":") |>
                        map_chr(~ .x[[1]]) |>
                        unique()
                })
            # Total number of proteins tagged with any of the tags of interest
            total_prots <- map_lgl(anno_tags, \(x) any(tag_interest %in% x)) |>
                sum()
            tibble(term = do.call(c, anno_tags)) |>
                filter(term %in% tag_interest) |>
                mutate(term = factor(term, levels = tag_interest)) |>
                group_by(term) |>
                summarize(n_proteins = n()) |>
                add_row(term = "total", n_proteins = total_prots) |>
                mutate(species = spp)
        }) |>
    pivot_wider(names_from = species, values_from = n_proteins)

# # A tibble: 5 × 5
#   term         Cripar Csonor Pstein Tgraci
#   <fct>         <int>  <int>  <int>  <int>
# 1 kegg_ko       10576  10760   9675   9753
# 2 go             7414   7907   7095   6834
# 3 cog            3937   4411   4009   3785
# 4 enzyme_ec      3698   3870   3517   3566
# 5 kegg_pathway   3080   3427   3074   2927
# 6 total         12249  12971  11645  11369



