#'
#' Extract GO terms for 1-to-1 OrthoFinder HOGs.
#' Output from here is used for HyPhy.
#'

source("_scripts/00-preamble.R")




#' OrthoFinder directory. It's assumed this directory is unchanged from the
#' output from OrthoFinder.
ofd <- function(...) {
    dots <- list(...)
    if (length(dots) == 0) return(dirs$orthofinder)
    do.call(paste0, c(list(dirs$orthofinder, "/"), dots))
}
#' OrthoFinder Extraction directory. This is where files live that contain
#' information we extracted from the OrthoFinder information.
oed <- function(...) {
    dots <- list(...)
    if (length(dots) == 0) return(dirs$orthofinder_extr)
    do.call(paste0, c(list(dirs$orthofinder_extr, "/"), dots))
}
#' GO terms directory
god <- function(...) {
    dots <- list(...)
    if (length(dots) == 0) return(dirs$go)
    do.call(paste0, c(list(dirs$go, "/"), dots))
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

#'
#' Get GO terms for a dataframe that contains only one species.
#'
get_gos <- function(d) {

    sp <- d[["species"]][[1]]

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
                                        "Gene Product Form ID"),
                          progress = FALSE) |>
            select(`DB Object ID`, `GO ID`) |>
            set_names(c("gene", "go")) |>
            filter(gene %in% d[["gene"]]) |>
            group_by(gene) |>
            summarize(go = paste(tolower(trimws(go)), collapse = ";"))

    } else if (grepl(".tab$", go_file)) {

        go_df <- read_tsv(go_file, col_types = cols(),
                          col_names = c("gene", "desc", "uniprot",
                                        "go", "kegg", "pfam"),
                          progress = FALSE) |>
            # These ones use transcript IDs instead of genes:
            mutate(gene = trans2genes(gene),
                   go = tolower(go) |> str_replace_all(", ", ";")) |>
            select(gene, go) |>
            filter(!is.na(go))

    } else if (grepl(".tsv$", go_file)) {

        go_df <- read_lines(go_file, progress = FALSE) |>
            base::`[`(-1) |>
            str_split("\t") |>
            lapply(\(x) {
                gox <- x[startsWith(x, "go")]
                tibble(gene = x[[1]],
                       go = ifelse(length(gox) == 0, NA_character_,
                                   tolower(paste(gox, collapse = ";"))))
            }) |>
            do.call(what = bind_rows) |>
            mutate(gene = trans2genes(gene)) |>
            filter(!is.na(go))

    } else {
        stop("GO term file extension for file '",
             basename(go_file), "' not recognized")
    }

    #' Make doubly sure that there aren't duplicate entries for each gene
    #' or duplicate GO terms:
    go_df <- go_df |>
        split(~ gene) |>
        map_dfr(\(z) {
            if (nrow(z) == 1) return(z)
            z |>
                summarize(go = paste(go, collapse = ";") |>
                              str_split(";") |>
                              map_chr(\(x) {
                                  paste(unique(x), collapse = ";")
                                  }))
        })

    return(left_join(d, go_df, by = "gene"))

}


all_nodes <- c("N0", "N1", "N3", "N5")
#' Species to not consider for each node, most of which we're removing
#' because they do not fall within a particular node.
#' Species 'Cmarin' only had ~58% of genes match to orthogroups, so
#' I'm removing it from all the analyses.
spp_rm <- list(N0 = "Cmarin",
               N1 = "Mdomes",
               N3 = c("Asteph", "Cquinq", "Aaegyp"),
               N5 = "Csonor")
for (i in 2:length(all_nodes)) spp_rm[[i]] <- c(spp_rm[[i-1]], spp_rm[[i]])


for (.node in all_nodes) {

    hog_gene_df <- read_tsv(ofd("Phylogenetic_Hierarchical_Orthogroups/",
                                .node, ".tsv"),
                            col_types = cols()) |>
        #' Remove columns that are no longer necessary, including species not
        #' on this node or with poor OrthoFinder matching.
        select(-all_of(spp_rm[[.node]]), -OG, -`Gene Tree Parent Clade`) |>
        rename(hog = HOG) |>
        mutate(across(-hog, \(x) {
            z <- rep(list(character(0)), length(x))
            z[!is.na(x)] <- str_split(x[!is.na(x)], ", ")
            return(z)
        })) |>
        #' Filter for just 1-to-1 orthogroups:
        filter(if_all(-hog, \(x) map_lgl(x, \(y) length(y) == 1))) |>
        # Because all species have only 1 gene here, then we can unnest these cols:
        mutate(across(-hog, unlist)) |>
        pivot_longer(-hog, names_to = "species", values_to = "gene") |>
        arrange(species, gene, hog) |>
        select(species, gene, hog) |>
        # Now fill in the GO terms:
        split(~ species) |>
        safe_mclapply(get_gos) |>
        do.call(what = bind_rows)

    if (!dir.exists(oed("Single_Copy_HOG_GO"))) {
        dir.create(oed("Single_Copy_HOG_GO"), recursive = TRUE)
    }

    # Summarizing just by HOG:
    hog_summ_df <- hog_gene_df |>
        group_by(hog) |>
        summarize(go = paste(go[!is.na(go)], collapse = ";"), .groups = "drop") |>
        mutate(go = str_split(go, ";") |>
                   map_chr(\(x) paste(unique(x), collapse = ";")),
               go = ifelse(go == "NA", NA, go))

    fn1 <- paste0(.node, "-GO-by-species-genes.tsv")
    fn2 <- paste0(.node, "-GO-by-HOG.tsv")
    write_tsv(hog_gene_df, oed("Single_Copy_HOG_GO/", fn1))
    write_tsv(hog_summ_df, oed("Single_Copy_HOG_GO/", fn2))
    cat("Wrote to", fn1, "and", fn2, "inside\n  ",
        oed("Single_Copy_HOG_GO"), "\n\n")

}

