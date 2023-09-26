library(tidyverse)

#' DEFUNCT -- EVERYTHING HERE IS DONE INSIDE `make-genome-stats.R`


#' Get transcript bounds for a given data frame read from a gff file and
#' a string indicating the source database.
#' This is used below in `write_intergenic` function.
get_trans <- function(.gff_df, .source) {
    # .gff_df = gff_df

    if (.source == "GenBank") {
        #' I have to first use CDS to map gene names onto transcript IDs
        cds_df <- .gff_df |>
            filter(type == "CDS") |>
            # This one doesn't have necessary info and isn't in listed proteins:
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
        group_by(gene) |>
        #' Get the longest isoform for each gene, and simply take the
        #' first one if there's a tie:
        filter((end - start + 1L) == max(end - start + 1L)) |>
        filter(trans == trans[[1]]) |>
        ungroup() |>
        select(seqid, start, end) |>
        arrange(seqid, start)
    return(trans_df)
}


#' Get intron sizes and write to csv file based on species abbreviation,
#' source database, and list of genes to use.
write_intergenic <- function(.spp, .source, .seq_len_df) {

    stopifnot(.source %in% c("VectorBase", "GenBank", "InsectBase", "here"))

    .fn <- paste0("~/_data/_features/", .spp, "_features.gff3.gz")
    gff_df <- read_tsv(.fn, comment = "#",
                       col_names = c("seqid", "source", "type", "start",
                                     "end", "score", "strand", "phase",
                                     "attributes"),
                       col_types = "ccciicccc", progress = FALSE)

    trans_df <- get_trans(gff_df, .source)

    stopifnot(all(trans_df$end > trans_df$start))

    igl_df <- trans_df |>
        (\(x) split(x, x$seqid))() |>
        map_dfr(\(.d) {
            .L <- filter(.seq_len_df, species == .spp, seqid == .d$seqid[[1]]) |>
                getElement("length")
            .o <- tibble(seqid = .d$seqid[[1]],
                         interg_len = c(.d$start - lag(.d$end, 1L, 1L) - 1L,
                                        .L - tail(.d$end, 1)))
            return(.o)
        })

    out_fn <- sprintf("_data/intergenic-%s.csv.xz", .spp)
    cat("Writing to ", out_fn, "\n")
    write_csv(igl_df, out_fn)

    invisible(NULL)
}



spp_df <- tibble(spp = c("Aaegyp", "Asteph", "Bantar", "Cmarin", "Cquinq",
                         "Cripar", "Csonor", "Ctenta", "Mdomes",
                         "Pakamu", "Ppemba", "Pstein", "Pvande",
                         "Tgraci"),
                 source = c("VectorBase", "VectorBase", "InsectBase", "InsectBase", "VectorBase",
                            "GenBank", "here", "InsectBase", "InsectBase",
                            "InsectBase", "InsectBase", "here", "InsectBase",
                            "here"))


seq_len_df <- read_csv("_data/seq_lens.csv", col_types = "cci")


for (i in 1:nrow(spp_df)) {
    write_intergenic(spp_df[["spp"]][[i]],
                  spp_df[["source"]][[i]],
                  seq_len_df)
}; rm(i)



