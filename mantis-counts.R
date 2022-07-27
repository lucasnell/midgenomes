
#' Summarize output from mantis

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

tany_anno_df <- read_tsv("~/_data/annotation/tany_mantis_new/output_annotation.tsv")

# Total number of proteins annotated:
tany_anno_df[["Query"]] %>% unique() %>% length()

tany_anno_df %>%
    mutate(db = lab_dbs(Ref_file)) %>%
    group_by(db) %>%
    summarize(n_proteins= length(unique(Query)))
# 1 eggNOG      10242
# 2 KOfam       12225
# 3 NPFM         3764
# 4 Pfam        12163
# 5 TCDB         1551


pstein_anno_df <- read_tsv("~/_data/annotation/Pstein_mantis/output_annotation.tsv")

# Total number of proteins annotated:
pstein_anno_df[["Query"]] %>% unique() %>% length()

pstein_anno_df %>%
    mutate(db = lab_dbs(Ref_file)) %>%
    group_by(db) %>%
    summarize(n_proteins= length(unique(Query)))
# 1 eggNOG      10731
# 2 KOfam       12517
# 3 NPFM         3851
# 4 Pfam        12443
# 5 TCDB         1710



#' -------------------------------------------------
#' Count number of proteins with various tags.
#' -------------------------------------------------

#' (Doing this by lines bc the number of columns isn't the same among rows)
#'
#' Tags I'm interested in:
tag_interest <- c("go", "kegg_ko", "cog", "enzyme_ec", "kegg_pathway")


tany_anno_tags <- read_lines(paste0("~/_data/annotation/tany_mantis_new/",
                                     "consensus_annotation.tsv")) %>%
    str_split("\t") %>%
    .[-1] %>%
    map(function(x) {
        x[-1:-6] %>%
            str_split(":") %>%
            map_chr(~ .x[[1]]) %>%
            unique()
    })

tibble(term = do.call(c, tany_anno_tags)) %>%
    filter(term %in% tag_interest) %>%
    group_by(term) %>%
    summarize(n_proteins = n()) %>%
    arrange(desc(n_proteins))
# 1 kegg_ko            9804
# 2 go                 6855
# 3 cog                3800
# 4 enzyme_ec          3555
# 5 kegg_pathway       2946


pstein_anno_tags <- read_lines(paste0("~/_data/annotation/Pstein_mantis/",
                                    "consensus_annotation.tsv")) %>%
    str_split("\t") %>%
    .[-1] %>%
    map(function(x) {
        x[-1:-6] %>%
            str_split(":") %>%
            map_chr(~ .x[[1]]) %>%
            unique()
    })

tibble(term = do.call(c, pstein_anno_tags)) %>%
    filter(term %in% tag_interest) %>%
    group_by(term) %>%
    summarize(n_proteins = n()) %>%
    arrange(desc(n_proteins))
# 1 kegg_ko            9660
# 2 go                 7011
# 3 cog                3933
# 4 enzyme_ec          3505
# 5 kegg_pathway       2962




