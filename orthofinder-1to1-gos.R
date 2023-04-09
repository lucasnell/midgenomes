
library(tidyverse)
library(parallel)

options(mc.cores = max(1L, detectCores()-2L))


#' OrthoFinder directory. It's assumed this directory is unchanged from the
#' output from OrthoFinder.
ofd <- function(...) {
    dots <- list(...)
    base_dir <- "~/_data/chir_orthofinder/orthofinder-output"  ## << change as necessary
    base_dir <- str_remove(base_dir, "/$")
    if (length(dots) == 0) return(base_dir)
    do.call(paste0, c(list(base_dir, "/"), dots))
}
#' GO terms directory
god <- function(...) {
    dots <- list(...)
    base_dir <- "~/_data/_go-terms"  ## << change as necessary
    base_dir <- str_remove(base_dir, "/$")
    if (length(dots) == 0) return(base_dir)
    do.call(paste0, c(list(base_dir, "/"), dots))
}
#' This removes transcript id that's an extension separated by
#' either the last '-' or '.'
#' This converts the transcript id to a gene id.
tran2gene <- function(tran) {
    loc <- tail(str_locate_all(tran, "\\.|-")[[1]][,"end"], 1)
    if (length(loc) == 0) return(tran)
    return(str_sub(tran, 1, loc-1))
}


for (.node in c("N0", "N1", "N3", "N5")) {

    node_df <- tibble(node = .node,
                      hog_files = ofd("Single_Copy_HOG_Sequences/", .node) |>
                          list.files("*.faa", full.names = TRUE),
                      hog = hog_files |> basename() |> str_remove(".faa"),
                      hog_lines = mclapply(hog_files, \(hf) {
                          read_lines(hf) |>
                              keep(\(x) startsWith(x, ">"))
                      }),
                      species = hog_lines |>
                          map(\(x) str_remove_all(x, ">|__.*")),
                      gene = hog_lines |>
                          map(\(x) {
                              str_remove_all(x, ".*__") |>
                                  map_chr(tran2gene)
                          })) |>
        select(-hog_files, -hog_lines) |>
        unnest(species:gene) |>
        mutate(go = "")


    # Now fill in the GO terms:
    for (sp in unique(node_df$species)) {

        sp_genes <- node_df |> filter(species == sp) |> getElement("gene")

        go_file <- list.files(god(), sp, full.names = TRUE)
        if (length(go_file) != 1) {
            stop("Number of files found for species ", sp, " = ", length(go_file))
        }

        if (grepl(".gaf$", go_file)) {

            go_df <- read_tsv(go_file, comment = c("!"), col_types = cols(),
                              col_names = c("DB", "DB Object ID",
                                            "DB Object Symbol",
                                            "Qualifier", "GO ID",
                                            "DB:Reference (|DB:Reference)",
                                            "Evidence Code", "With (or) From",
                                            "Aspect", "DB Object Name",
                                            "DB Object Synonym (|Synonym)",
                                            "DB Object Type", "Taxon(|taxon)",
                                            "Date", "Assigned By",
                                            "Annotation Extension",
                                            "Gene Product Form ID")) |>
                select(`DB Object ID`, `GO ID`) |>
                set_names(c("gene", "go")) |>
                filter(gene %in% sp_genes) |>
                group_by(gene) |>
                summarize(go = paste(tolower(trimws(go)), collapse = ";"))

        } else if (grepl(".tab$", go_file)) {

            go_df <- read_tsv(go_file, col_types = cols(),
                              col_names = c("gene", "desc", "uniprot",
                                            "go", "kegg", "pfam")) |>
                # These ones use transcript IDs instead of genes:
                mutate(gene = tran2gene(gene),
                       go = tolower(go) |> str_replace_all(", ", ";")) |>
                select(gene, go) |>
                filter(!is.na(go))

        } else if (grepl(".tsv$", go_file)) {

            go_df <- read_lines(go_file) |>
                base::`[`(-1) |>
                str_split("\t") |>
                mclapply(\(x) {
                    gox <- x[startsWith(x, "go")]
                    tibble(gene = x[[1]],
                           go = ifelse(length(gox) == 0, NA_character_,
                                       tolower(paste(gox, collapse = ";"))))
                }) |>
                do.call(what = bind_rows) |>
                mutate(gene = tran2gene(gene)) |>
                filter(!is.na(go))

        } else {
            stop("GO term file extension for file '",
                 basename(go_file), "' not recognized")
        }

        sp_lgl <- node_df[["species"]] == sp
        for (i in 1:nrow(go_df)) {
            gene_lgl <- node_df[["gene"]] == go_df[["gene"]][[i]]
            node_df[["go"]][sp_lgl & gene_lgl] <- go_df[["go"]][[i]]
        }; rm(i)

    }


    if (!dir.exists(ofd("Single_Copy_HOG_GO"))) {
        dir.create(ofd("Single_Copy_HOG_GO"), recursive = TRUE)
    }

    node_summ_df <- node_df |>
        group_by(node, hog) |>
        summarize(go = paste(go[go != ""], collapse = ";"), .groups = "drop") |>
        mutate(go = str_split(go, ";") |>
                   map_chr(\(x) paste(unique(x), collapse = ";")))

    fn1 <- paste0(.node, "-GO-by-species-genes.tsv")
    fn2 <- paste0(.node, "-GO-by-HOG.tsv")
    write_tsv(node_df, ofd("Single_Copy_HOG_GO/", fn1))
    write_tsv(node_summ_df, ofd("Single_Copy_HOG_GO/", fn2))
    cat("Wrote to\n", fn1, "\nand\n", fn2, "\ninside\n",
        ofd("Single_Copy_HOG_GO"), "\n\n")

}
