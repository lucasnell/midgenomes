library(tidyverse)
library(parallel)


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





#' `"here"`:
#' * for each gene name, find all transcripts (`type == "transcript"`)
#'   that have `geneID=<gene name>` in `attributes`, extract its bounds and
#'   name; the name is under `attributes` column with `ID=<transcript name>`
#' * extract all exons (`type == "exon"`) that have `Parent=<mRNA name>`
#'   in `attributes`
#' * introns are gaps between exons within transcript bounds
#'
#' `"InsectBase"`:
#' * for each gene name, find all with `type == "mRNA"` and
#'   `attribute == "*Parent=<gene name>"`; extract those mRNA bounds and
#'   names, where names are under `attributes` column with `ID=<mRNA name>`
#'   (I checked, and there are 1:1 matches between genes and mRNAs in all
#'   files from InsectBase I downloaded)
#' * extract all exons (`type == "exon"`) that have `Parent=<mRNA name>`
#'   in `attributes`
#' * introns are gaps between exons within mRNA bounds
#'
#' `"VectorBase"`:
#' * for each gene name, find all with `type == "mRNA"` and
#'   `attribute == "*Parent=<gene name>"`
#' * When there is >1 overlapping mRNA, choose the longest isoform
#' * For each mRNA, extract the bounds and name, where names are under
#'   `attributes` column with `ID=<mRNA name>`
#' * extract all exons (`type == "exon"`) that have `Parent=<mRNA name>`
#'   in `attributes` (can have multiple parents, separated by commas
#'   like this: `Parent=<mRNA name 1>, <mRNA name 2>, ...`)
#' * introns are gaps between exons within mRNA bounds
#'
#' `"GenBank"`:
#' * Same as for VectorBase
#'


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
    # .spp = spp_df[["spp"]][[5]]; .source = spp_df[["source"]][[5]]
    # .genes = spp_df[["genes"]][[5]]
    # rm(.spp, .source, .fn, .genes, gff_df, trans_df, exon_mats, max_introns,
    #    label_format, intron_df, out_fn)
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





# ----------------------------------------------------------------------------*
# ----------------------------------------------------------------------------*
# ----------------------------------------------------------------------------*
# ----------------------------------------------------------------------------*
# ----------------------------------------------------------------------------*


tany_gff_df <- read_tsv("~/_data/annotation/tany_tsebra.gff3.gz", skip = 3,
                        col_names = c("seqid", "source", "type", "start",
                                      "end", "score", "strand", "phase",
                                      "attributes"),
                        col_types = "ccciicccc")
pstein_gff_df <- read_tsv("~/_data/annotation/Pstein_tsebra.gff3.gz", skip = 3,
                          col_names = c("seqid", "source", "type", "start",
                                        "end", "score", "strand", "phase",
                                        "attributes"),
                          col_types = "ccciicccc")

#' Just exons, to be added to transcript dataframes later.
tany_exon_df <- tany_gff_df %>%
    filter(type == "exon") %>%
    mutate(attributes = str_remove_all(attributes, "Parent=")) %>%
    rename(protein = attributes) %>%
    select(seqid, start, end, protein)
pstein_exon_df <- pstein_gff_df %>%
    filter(type == "exon") %>%
    mutate(attributes = str_remove_all(attributes, "Parent=")) %>%
    rename(protein = attributes) %>%
    select(seqid, start, end, protein)



#' Extract gene and protein from transcript attributes string.
#' NOTE: only works for transcripts!
get_gene_protein <- function(.df) {
    attr_str <- .df[["attributes"]]
    if (any(grepl("Parent", attr_str))) {
        stop("`get_gene_protein` only works for transcripts!")
    }
    gp_list <- attr_str %>%
        str_split(";") %>%
        map(function(.str) {
            g <- .str[which(grepl("^geneID=", .str))]
            p <- .str[which(grepl("^ID=", .str))]
            c(gsub("geneID=", "", g), gsub("ID=", "", p))
        })
    .df[["gene"]] = map_chr(gp_list, ~ .x[[1]])
    .df[["protein"]] = map_chr(gp_list, ~ .x[[2]])
    return(.df)
}

#' Add exon starts and ends for each transcript. It also orders exons by
#' start position. (This is probably already done in the GFF3 file, but better
#' safe than sorry.)
get_exons <- function(.df, .exon_df) {
    .df$exons = mcmapply(                                   # <<<<<<<<<<<<<<<
        FUN = function(p, s) {
            exd <- .exon_df[.exon_df$protein == p,]
            stopifnot(all(exd$seqid == s))
            exd <- exd[order(exd$start),]
            return(cbind(exd$start, exd$end))
        },
        .df$protein, .df$seqid,
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    return(.df)
}
#' Get intron lengths for transcripts.
#' Dataframe must have `exons`, `start`, and `end` columns.
get_intron_lens <- function(.df) {
    #' For one transcript:
    intron_lens_one_trans <- function(exon_mat, tran_start, tran_end) {
        n_exons <- nrow(exon_mat)
        stopifnot(n_exons > 0)
        #' `n_introns` means all possible introns, so it includes before
        #' and after first and last exons, respectively.
        #' These may not be present and will be filtered out later,
        #' as will any length-0 gaps between exons.
        n_introns <- n_exons + 1
        introns <- numeric(n_introns)
        if (exon_mat[1,1] > tran_start) {
            introns[1] <- exon_mat[1,1] - tran_start
        }
        if (exon_mat[n_exons,2] < tran_end) {
            introns[n_introns] <- tran_end - exon_mat[n_exons,2]
        }
        if (n_exons > 1) {
            for (j in 2:n_exons) {
                introns[j] <- exon_mat[j,1] - exon_mat[(j-1),2] - 1
            }
        }
        introns <- introns[introns > 0]
        return(introns)
    }
    .df$intron_lens = mcmapply(                              # <<<<<<<<<<<<<<<
        FUN = intron_lens_one_trans,
        .df$exons, .df$start, .df$end,
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    return(.df)
}



#' Takes ~2.3 sec on my machine (6 threads)
tany_trans_df <- tany_gff_df %>%
    filter(type == "transcript") %>%
    get_gene_protein() %>%
    select(seqid, start, end, gene, protein) %>%
    get_exons(tany_exon_df) %>%
    get_intron_lens()

tany_trans_df$intron_lens %>% do.call(what = c) %>% mean()



#' Takes ~2.4 sec
pstein_trans_df <- pstein_gff_df %>%
    filter(type == "transcript") %>%
    get_gene_protein() %>%
    select(seqid, start, end, gene, protein) %>%
    get_exons(pstein_exon_df) %>%
    get_intron_lens()

pstein_trans_df$intron_lens %>% do.call(what = c) %>% mean()











#' ============================================================================
#' ============================================================================
#'
#' Numbers of unique proteins and genes
#'
#' ============================================================================
#' ============================================================================

tany_trans_df %>% distinct(protein) %>% nrow()
tany_trans_df %>% distinct(gene) %>% nrow()

pstein_trans_df %>% distinct(protein) %>% nrow()
pstein_trans_df %>% distinct(gene) %>% nrow()

