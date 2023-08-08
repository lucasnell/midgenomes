# *****************************************************************************
# *****************************************************************************

#' This script is now defunct. See `./chtc/03-phylo/cafe.sh`

# *****************************************************************************
# *****************************************************************************



#'
#' Create tab-delimited family "counts" file for CAFE5 from OrthoFinder
#' output.
#'

library(readr)
count_df <- read_tsv("~/_data/chir_orthofinder/orthofinder-output/Phylogenetic_Hierarchical_Orthogroups/orthofinder-hogs-n0.tsv.gz", col_types = cols())
count_df <- count_df[,-which(colnames(count_df) %in%
                                 c("OG", "Gene Tree Parent Clade"))]
#' This species only had ~58% of genes match to orthogroups, so
#' I'm removing it from the analyses:
count_df <- count_df[,-which(colnames(count_df) == "Cmarin")]

spp_names <- colnames(count_df)[-1]

for (sp in spp_names) {
    z <- rep(0L, nrow(count_df))
    na_lgl <- !is.na(count_df[[sp]])
    z[na_lgl] <- sapply(strsplit(count_df[[sp]][na_lgl], ", "), length)
    count_df[[sp]] <- z
}
count_df[["Family ID"]] <- count_df[["HOG"]]
count_df[["Desc"]] <- "(null)"
count_df <- count_df[,c("Desc", "Family ID", spp_names)]
write_tsv(count_df, "${COUNTS_FILE}")


read_tsv("~/_data/chir_orthofinder/orthofinder-output/Phylogenetic_Hierarchical_Orthogroups/N0.tsv", col_types = cols()) |>
    filter(!is.na(Mdomes)) |>
    nrow()
# [1] 9485

read_tsv("~/_data/chir_orthofinder/orthofinder-output/Phylogenetic_Hierarchical_Orthogroups/N0.tsv", col_types = cols()) |>
    mutate(nouts = pmap_int(list(Aaegyp, Asteph, Cquinq),
                            \(x, y, z) {
                                sum(c(is.na(x), is.na(y), is.na(z)))
                            })) |>
    filter(nouts >= 2) |>
    nrow()
# [1] 9675

read_tsv("~/_data/chir_orthofinder/orthofinder-output/Phylogenetic_Hierarchical_Orthogroups/N0.tsv", col_types = cols()) |>
    filter(!is.na(Csonor)) |>
    nrow()
# [1] 9940

read_tsv("~/_data/chir_orthofinder/orthofinder-output/Phylogenetic_Hierarchical_Orthogroups/N1.tsv", col_types = cols()) |>
    mutate(nouts = pmap_int(list(Aaegyp, Asteph, Cquinq),
                        \(x, y, z) {
                            sum(c(is.na(x), is.na(y), is.na(z)))
                        })) |>
    filter(nouts >= 2) |>
    nrow()
# [1] 9435

read_tsv("~/_data/chir_orthofinder/orthofinder-output/Phylogenetic_Hierarchical_Orthogroups/N3.tsv", col_types = cols()) |>
    filter(!is.na(Csonor)) |>
    nrow()
# [1] 10214

read_tsv("~/_data/chir_orthofinder/orthofinder-output/Phylogenetic_Hierarchical_Orthogroups/N5.tsv", col_types = cols()) |>
    filter(!is.na(Pstein)) |>
    nrow()
# [1] 10403


