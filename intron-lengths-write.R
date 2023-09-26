library(tidyverse)
library(parallel)

#' DEFUNCT -- EVERYTHING HERE IS DONE INSIDE `make-genome-stats.R`


options("mc.cores" = max(1L, detectCores() - 2L))

#' Safe call to `mcmapply` that reverts to `mapply` on Windows:
safe_mcmapply <- function(FUN, ..., MoreArgs = NULL) {
    if (.Platform$OS.type == "unix") {
        out <- mcmapply(FUN, ..., MoreArgs = MoreArgs,
                        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    } else {
        out <- mapply(FUN, ..., MoreArgs = MoreArgs,
                      SIMPLIFY = FALSE, USE.NAMES = FALSE)
    }
    return(out)
}


#' HOG file for N0 node from OrthoFinder:
hog_file_n0 <- paste0("~/_data/chir_orthofinder/orthofinder-output/",
                      "Phylogenetic_Hierarchical_Orthogroups/N0.tsv")

#'
#' Find names of "... genes with orthologues across [all]... species ... and
#' orthogroups with less than four paralogues per species."
#' (Martín-Durán et al. 2021).
#' We're using output from OrthoFinder at the N0 node so that our outgroups
#' are also captured.
#'

gnames_df <- hog_file_n0 |>
    read_tsv(col_types = cols()) |>
    #' This species only had ~58% of genes match to orthogroups, so
    #' I'm removing it from the analyses:
    select(-Cmarin) |>
    mutate(across(Aaegyp:Tgraci, \(x) {
        z <- rep(list(character(0)), length(x))
        z[!is.na(x)] <- str_split(x[!is.na(x)], ", ")
        return(z)
    })) |>
    filter(if_all(Aaegyp:Tgraci,
                  \(x) map_lgl(x, \(y) length(y) >= 1 & length(y) < 4))) |>
    select(HOG, Aaegyp:Tgraci)



#' Get transcript bounds for a given data frame read from a gff file,
#' a string indicating the source database, and a string vector indicating
#' the genes to keep.
#' This is used below in `write_introns` function.
get_trans <- function(.gff_df, .source, .genes) {
    # .gff_df = gff_df

    if (.source == "GenBank") {
        #' I have to first use CDS to map gene names onto transcript IDs
        cds_df <- .gff_df |>
            filter(type == "CDS") |>
            # This one doesn't have necessary info and isn't in my list anyway:
            filter(!grepl("ID=cds-CHIRRI_LOCUS9184", attributes)) |>
            mutate(attributes = str_split(attributes, ";"),
                   gene = attributes |>
                       map_chr(\(x) x[grepl("^protein_id=", x)]) |>
                       str_remove_all("^protein_id="),
                   # Remove redundant ".1" at end of each gene name:
                   gene = sub("\\.[^\\.]+$", "", gene),
                   trans = attributes |>
                       map_chr(\(x) x[grepl("^Parent=", x)]) |>
                       str_remove_all("^Parent=") |>
                       str_replace("^rna-", "rna_") |>
                       sub(pattern = "-[^-]+$", replacement = "") |>
                       str_replace("^rna_", "rna-")) |>
            distinct(gene, trans)
        gt_map <- cds_df[["gene"]] |> as.list()
        names(gt_map) <- cds_df[["trans"]]
        #' Now I can look for transcripts and add gene names from `gt_map`
        trans_df <- .gff_df |>
            filter(type == "mRNA") |>
            # This one doesn't have necessary info and isn't in my list anyway:
            filter(!grepl("ID=rna-CHIRRI_LOCUS9184", attributes)) |>
            mutate(trans = str_split(attributes, ";") |>
                       map_chr(\(x) x[grepl("^ID=", x)]) |>
                       str_remove_all("^ID="),
                   gene = trans |>
                       str_replace("^rna-", "rna_") |>
                       sub(pattern = "-[^-]+$", replacement = "") |>
                       str_replace("^rna_", "rna-") |>
                       map_chr(\(.t) gt_map[[.t]]))
    } else if (.source %in% c("VectorBase", "InsectBase")) {
        trans_df <- .gff_df |>
            filter(type == "mRNA") |>
            mutate(attributes = str_split(attributes, ";"),
                   gene = attributes |>
                       map_chr(\(x) x[grepl("^Parent=", x)]) |>
                       str_remove_all("^Parent="),
                   trans = attributes |>
                       map_chr(\(x) x[grepl("^ID=", x)]) |>
                       str_remove_all("^ID="))
    } else {
        trans_df <- .gff_df |>
            filter(type == "transcript") |>
            mutate(trans = str_split(attributes, ";") |>
                       map_chr(\(x) x[grepl("^ID=", x)]) |>
                       str_remove_all("^ID="),
                   #' Gene names are weird here, so it's easier to remove the
                   #' ".tX" suffix to transcript names:
                   gene = sub("\\.[^\\.]+$", "", trans))
    }
    trans_df <- trans_df |>
        select(seqid, gene, trans, start, end, strand, phase) |>
        filter(gene %in% .genes) |>
        mutate(length = end - start + 1L) |>
        group_by(gene) |>
        #' Get the longest isoform for each gene, and simply take the
        #' first one if there's a tie:
        filter(length == max(length)) |>
        filter(trans == trans[[1]]) |>
        ungroup() |>
        arrange(seqid, start)
    return(trans_df)
}

#' Get exon bounds for a given data frame read from a gff file and a string
#' vector indicating the transcripts you want to keep.
#' This is used below in `write_introns` function.
get_exons <- function(.gff_df, .trans_ids) {
    exon_df <- .gff_df |>
        filter(type == "exon") |>
        mutate(trans = str_split(attributes, ";") |>
                   map_chr(\(x) x[grepl("^Parent=", x)]) |>
                   str_remove_all("^Parent="))
    if (any(grepl(",", exon_df[["trans"]]))) {
        exon_df <- exon_df |>
            mutate(trans = str_split(trans, ",") |>
                       map_chr(\(x) {
                           xx <- x[x %in% .trans_ids]
                           if (length(xx) == 0) xx <- NA_character_
                           return(xx)
                       })) |>
            filter(!is.na(trans))
    } else {
        exon_df <- exon_df |>
            filter(trans %in% .trans_ids)
    }
    exon_df <- exon_df |>
        select(trans, start, end) |>
        arrange(trans, start)
    exon_list <- exon_df |>
        select(trans, start, end) |>
        (\(x) split(x, x$trans))()
    trans_names <- map_chr(1:length(exon_list),
                           \(i) exon_list[[i]][["trans"]][[1]])
    exon_mat_list <- map(exon_list, \(x) as.matrix(x[,c("start", "end")]))
    names(exon_mat_list) <- trans_names
    #' So it's the same order as transcript data frame:
    exon_mat_list <- exon_mat_list[.trans_ids]
    return(exon_mat_list)
}

#' Get intron bounds for transcript `i` based on a list of transcript and exon
#' bounds that are assumed to be in the same order.
#' This is used below in `write_introns` function.
get_one_intron_set <- function(exon_mat, trans_row, .label_format) {

    tstart <- trans_row[["start"]]
    tend <- trans_row[["end"]]

    n_exons <- nrow(exon_mat)
    stopifnot(n_exons > 0)
    #' `n_introns` means all possible introns, so it includes before
    #' and after first and last exons, respectively.
    #' These may not be present and will be filtered out later,
    #' as will any length-0 gaps between exons.
    n_introns <- n_exons + 1
    ibounds <- matrix(0L, n_introns, 2)
    colnames(ibounds) <- c("start", "end")
    if (exon_mat[1,"start"] > tstart) {
        ibounds[1,"start"] <- tstart
        ibounds[1,"end"] <- exon_mat[1,"start"] - 1L
    }
    if (exon_mat[n_exons,"end"] < tend) {
        ibounds[n_introns,"start"] <- exon_mat[n_exons,"end"] + 1L
        ibounds[n_introns,"end"] <- tend
    }
    if (n_exons > 1) {
        for (j in 2:n_exons) {
            ibounds[j,"start"] <- exon_mat[(j-1),"end"] + 1L
            ibounds[j,"end"] <- exon_mat[j,"start"] - 1L
        }
    }
    ibounds <- ibounds[ibounds[,1] > 0,,drop=FALSE]
    if (nrow(ibounds) == 0) return(NULL)
    ibounds <- ibounds[ibounds[,2] >= ibounds[,1],,drop=FALSE]
    if (nrow(ibounds) == 0) return(NULL)

    intron_df_i <- trans_row[,c("seqid", "gene", "trans")] |>
        mutate(start = list(ibounds[,"start"]),
               end = list(ibounds[,"end"])) |>
        unnest(start:end) |>
        mutate(id = sprintf(.label_format, 1:n()))

    return(intron_df_i)
}


#' Get intron sizes and write to csv file based on species abbreviation,
#' source database, and list of genes to use.
write_introns <- function(.spp, .source, .genes) {

    .fn <- paste0("~/_data/_features/", .spp, "_features.gff3.gz")
    gff_df <- read_tsv(.fn, comment = "#",
                       col_names = c("seqid", "source", "type", "start",
                                     "end", "score", "strand", "phase",
                                     "attributes"),
                       col_types = "ccciicccc")

    stopifnot(.source %in% c("VectorBase", "GenBank", "InsectBase", "here"))

    trans_df <- get_trans(gff_df, .source, .genes)
    exon_mats <- get_exons(gff_df, trans_df[["trans"]])

    stopifnot(identical(trans_df[["trans"]], names(exon_mats)))
    stopifnot(all(trans_df$end > trans_df$start))
    stopifnot(all(do.call(rbind, exon_mats)[,"end"] >
                      do.call(rbind, exon_mats)[,"start"]))

    # Used for labelling introns below:
    max_introns <- max(map_int(exon_mats, \(x) nrow(x) + 2L))
    label_format <- sprintf("intron-%%0%ii", ceiling(log10(max_introns)))

    intron_df <- safe_mcmapply(get_one_intron_set,
                               exon_mats, split(trans_df, 1:nrow(trans_df)),
                               MoreArgs = list(.label_format = label_format)) |>
        do.call(what = bind_rows) |>
        select(seqid, gene, trans, id, start, end)

    out_fn <- sprintf("_data/introns-%s.csv.xz", .spp)
    cat("Writing to ", out_fn, "\n")
    write_csv(intron_df, out_fn)

    invisible(NULL)
}



spp_df <- tibble(spp = c("Aaegyp", "Asteph", "Bantar", "Cquinq",
                         "Cripar", "Csonor", "Ctenta", "Mdomes",
                         "Pakamu", "Ppemba", "Pstein", "Pvande",
                         "Tgraci"),
                 source = c("VectorBase", "VectorBase", "InsectBase", "VectorBase",
                            "GenBank", "here", "InsectBase", "InsectBase",
                            "InsectBase", "InsectBase", "here", "InsectBase",
                            "here"),
                 genes = map(spp, \(s) unique(do.call(c, gnames_df[[s]]))))

for (i in 1:nrow(spp_df)) {
    write_introns(spp_df[["spp"]][[i]],
                  spp_df[["source"]][[i]],
                  spp_df[["genes"]][[i]])
}; rm(i)



