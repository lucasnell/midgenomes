
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
trans2genes <- function(trans) {
    locs <- str_locate_all(trans, "\\.|-") |>
        map_int(\(x) {
            if (nrow(x) == 0) return(NA_integer_)
            x[nrow(x),"end"]
        })
    out <- character(length(trans))
    out[is.na(locs)] <- trans[is.na(locs)]
    out[!is.na(locs)] <- str_sub(trans[!is.na(locs)], 1,
                                 locs[!is.na(locs)]-1)
    return(out)
}


for (.node in c("N0", "N1", "N3", "N5")) {

    node_df <- ofd("Phylogenetic_Hierarchical_Orthogroups/", .node, ".tsv") |>
        read_tsv(col_types = cols()) |>
        mutate(across(Aaegyp:Tgraci, \(x) {
            z <- rep(list(NULL), length(x))
            na_lgl <- !is.na(x)
            z[na_lgl] <- str_split(x[na_lgl], ", ")
            return(z)
        })) |>
        pivot_longer(Aaegyp:Tgraci, names_to = "species", values_to = "gene") |>
        select(-`Gene Tree Parent Clade`, -OG) |>
        rename(hog = HOG) |>
        unnest(gene)


    #' Now find all the GO terms:
    #' Took ~6 sec for N0 w/ 6 threads
    go_df <- mclapply(unique(node_df$species), \(sp) {

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
                group_by(gene) |>
                summarize(go = paste(tolower(trimws(go)), collapse = ";"))

        } else if (grepl(".tab$", go_file)) {

            go_df <- read_tsv(go_file, col_types = cols(),
                              col_names = c("gene", "desc", "uniprot",
                                            "go", "kegg", "pfam")) |>
                # These ones use transcript IDs instead of genes:
                mutate(gene = trans2genes(gene),
                       go = tolower(go) |> str_replace_all(", ", ";")) |>
                select(gene, go) |>
                filter(!is.na(go))

        } else if (grepl(".tsv$", go_file)) {

            go_df <- read_lines(go_file) |>
                base::`[`(-1) |>
                str_split("\t") |>
                mclapply(\(x) {
                    gox <- x[startsWith(x, "go")]
                    # Some weird ones were "go:GO" for some reason:
                    gox <- gox[gsub("[0-9]+", "", gox) == "go:"]
                    if (length(gox) > 0) {
                        goxx <- tolower(paste(gox, collapse = ";"))
                    } else goxx <- NA_character_
                    tibble(gene = x[[1]], go = goxx)
                }) |>
                do.call(what = bind_rows) |>
                mutate(gene = trans2genes(gene)) |>
                filter(!is.na(go))

        } else {
            stop("GO term file extension for file '",
                 basename(go_file), "' not recognized")
        }

        # Make sure all GO terms are sensible:
        go_df$go <- go_df$go |>
            str_remove_all("go:") |>
            str_split(";") |>
            map_chr(\(x) {
                # GO terms are 7 digits:
                x <- x[!grepl("\\D", x) & nchar(x) == 7]
                if (length(x) == 0) return(NA_character_)
                x <- paste0("go:", x, collapse = ";")
                return(x)
            })

        go_df <- go_df |>
            filter(!is.na(go)) |>
            mutate(species = sp)

        return(go_df)

    }) |>
        do.call(what = bind_rows)


    node_go_df <- left_join(node_df, go_df, by = c("species", "gene"))


    if (!dir.exists(ofd("All_HOG_GO"))) {
        dir.create(ofd("All_HOG_GO"), recursive = TRUE)
    }

    #' summarize by just HOG, not species:
    node_go_hog_df <- node_go_df |>
        mutate(go = str_split(go, ";"),
               go = map(go, \(x) {
                   if (all(is.na(x))) return(NULL) else return(x)
               })) |>
        group_by(hog) |>
        summarize(go = do.call(c, go) |>
                      unique() |>
                      paste(collapse = ";"),
                  .groups = "drop")

    fn1 <- paste0(.node, "-GO-by-species-genes.tsv")
    fn2 <- paste0(.node, "-GO-by-HOG.tsv")
    write_tsv(node_go_df, ofd("All_HOG_GO/", fn1))
    write_tsv(node_go_hog_df, ofd("All_HOG_GO/", fn2))
    cat("Wrote to\n", fn1, "\nand\n", fn2, "\ninside\n",
        ofd("All_HOG_GO"), "\n\n")

}
