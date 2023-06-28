
#' Summarize output from mantis
#' Updated on 2023-06-23


library(tidyverse)



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

tany_anno_df <- read_tsv("~/_data/_go-terms/zzz-mantis-full/Tgraci_mantis/output_annotation.tsv")

# Total number of proteins annotated:
tany_anno_df[["Query"]] |> unique() |> length()
# [1] 13424


tany_anno_df |>
    mutate(db = lab_dbs(Ref_file)) |>
    group_by(db) |>
    summarize(n_proteins= length(unique(Query)))
# 1 KOfam       12161
# 2 NPFM         3772
# 3 Pfam        12090
# 4 TCDB         1543
# 5 eggNOG      10182

pstein_anno_df <- read_tsv("~/_data/_go-terms/zzz-mantis-full/Pstein_mantis/output_annotation.tsv")

# Total number of proteins annotated:
pstein_anno_df[["Query"]] |> unique() |> length()
# [1] 14003

pstein_anno_df |>
    mutate(db = lab_dbs(Ref_file)) |>
    group_by(db) |>
    summarize(n_proteins= length(unique(Query)))
# 1 KOfam       12592
# 2 NPFM         3851
# 3 Pfam        12516
# 4 TCDB         1756
# 5 eggNOG      11024


#' -------------------------------------------------
#' Count number of proteins with various tags.
#' -------------------------------------------------

#' (Doing this by lines bc the number of columns isn't the same among rows)
#'
#' Tags I'm interested in:
tag_interest <- c("go", "kegg_ko", "cog", "enzyme_ec", "kegg_pathway")


tany_anno_tags <- read_lines(paste0("~/_data/_go-terms/zzz-mantis-full/Tgraci_mantis/",
                                     "consensus_annotation.tsv")) |>
    str_split("\t") |>
    base::`[`(-1) |>
    map(function(x) {
        x[-1:-6] |>
            str_split(":") |>
            map_chr(~ .x[[1]]) |>
            unique()
    })

tibble(term = do.call(c, tany_anno_tags)) |>
    filter(term %in% tag_interest) |>
    group_by(term) |>
    summarize(n_proteins = n()) |>
    arrange(desc(n_proteins))
# 1 kegg_ko            9753
# 2 go                 6834
# 3 cog                3785
# 4 enzyme_ec          3566
# 5 kegg_pathway       2927


pstein_anno_tags <- read_lines(paste0("~/_data/_go-terms/zzz-mantis-full/Pstein_mantis/",
                                    "consensus_annotation.tsv")) |>
    str_split("\t") |>
    base::`[`(-1) |>
    map(function(x) {
        x[-1:-6] |>
            str_split(":") |>
            map_chr(~ .x[[1]]) |>
            unique()
    })

tibble(term = do.call(c, pstein_anno_tags)) |>
    filter(term %in% tag_interest) |>
    group_by(term) |>
    summarize(n_proteins = n()) |>
    arrange(desc(n_proteins))
# 1 kegg_ko            9675
# 2 go                 7095
# 3 cog                4009
# 4 enzyme_ec          3517
# 5 kegg_pathway       3074




