
#' Make `genome-stats.csv`, `seq_lens.csv`, `intergenic-%s.csv.xz`,
#' `introns-%s.csv.xz`, `repeats-summary.csv` files all inside `_data`.


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
#' Safe call to `mclapply` that reverts to `lapply` on Windows:
safe_mclapply <- function(X, FUN, ...) {
    if (.Platform$OS.type == "unix") {
        out <- mclapply(X, FUN, ...)
    } else {
        out <- lapply(X, FUN, ...)
    }
    return(out)
}



#' --------------------
#' Set these directories:
#' --------------------
assembly_dir <- "~/_data/_assemblies"
proteins_dir <- "~/_data/_proteins"
features_dir <- "~/_data/_features"
orthofinder_dir <- "~/_data/chir_orthofinder/orthofinder-output"
repeats_dir <- "~/_data/_repeats"




#' ===========================================================================
# Read starting info on species and families ----
#' ===========================================================================


gstats_df <- read_csv("_data/species-names-families.csv", col_types = "ccc") |>
    #' This is the order I want to keep these species bc it matches
    #' with the phylogeny:
    mutate(species = factor(species, levels = species),
           spp_abbrev = factor(spp_abbrev, levels = spp_abbrev))

# Manually adding the source for annotations:
gstats_df$annot_source <- ""
gstats_df$annot_source[gstats_df$species == "Tanytarsus gracilentus"] <- "here"
gstats_df$annot_source[gstats_df$species == "Chironomus riparius"] <- "GenBank"
gstats_df$annot_source[gstats_df$species == "Chironomus tentans"] <- "InsectBase"
gstats_df$annot_source[gstats_df$species == "Polypedilum vanderplanki"] <- "InsectBase"
gstats_df$annot_source[gstats_df$species == "Polypedilum pembai"] <- "InsectBase"
gstats_df$annot_source[gstats_df$species == "Belgica antarctica"] <- "InsectBase"
gstats_df$annot_source[gstats_df$species == "Clunio marinus"] <- "InsectBase"
gstats_df$annot_source[gstats_df$species == "Propsilocerus akamusi"] <- "InsectBase"
gstats_df$annot_source[gstats_df$species == "Parochlus steinenii"] <- "here"
gstats_df$annot_source[gstats_df$species == "Culicoides sonorensis"] <- "here"
gstats_df$annot_source[gstats_df$species == "Anopheles stephensi"] <- "VectorBase"
gstats_df$annot_source[gstats_df$species == "Aedes aegypti"] <- "VectorBase"
gstats_df$annot_source[gstats_df$species == "Culex quinquefasciatus"] <- "VectorBase"
gstats_df$annot_source[gstats_df$species == "Musca domestica"] <- "InsectBase"



#' ===========================================================================
# Genome sizes ----
#' ===========================================================================


#'
#' Note: `seq_len_df` is necessary for intergenic sequences below
#'
if (file.exists("_data/seq_lens.csv")) {
    seq_len_df <- read_csv("_data/seq_lens.csv", col_types = "cci")
} else {
    seq_len_df <- gstats_df$spp_abbrev |>
        map_dfr(\(s) {
            .fn <- paste0(assembly_dir, "/", s, "_assembly.fasta.gz")
            sl <- read_lines(.fn, progress = FALSE)
            hl <- which(grepl("^>", sl))
            ll <- c(hl[-1] - 1L, length(sl))
            stopifnot(length(hl) == length(ll))
            lens <- map_int(1:length(hl), \(i) sum(nchar(sl[(hl[i]+1L):(ll[i])])))
            ids <- sl[hl] |> str_remove_all(">") |> str_remove("\ .*")
            tibble(spp_abbrev = s, seqid = ids, length = lens)
        })
    write_csv(seq_len_df, "_data/seq_lens.csv")
}


gsize_df <- seq_len_df |>
    group_by(spp_abbrev) |>
    summarize(gsize = sum(length)) |>
    mutate(spp_abbrev = factor(spp_abbrev, levels = levels(gstats_df$spp_abbrev))) |>
    arrange(spp_abbrev)



#' ===========================================================================
# Number of genes ----
#' ===========================================================================



ngenes_df <- gstats_df |>
    select(spp_abbrev, annot_source) |>
    rename(.spp = spp_abbrev, .source = annot_source) |>
    pmap_dfr(\(.spp, .source) {
        stopifnot(.source %in% c("VectorBase", "GenBank", "InsectBase", "here"))
        faa_lines <- read_lines(paste0(proteins_dir, "/", .spp,
                                       "_proteins.faa.gz"), progress = FALSE)
        faa_lines <- faa_lines[grepl("^>", faa_lines)]
        .n_prots <- length(faa_lines)
        if (.source == "VectorBase") {
            stopifnot(sum(grepl("gene=", faa_lines)) == length(faa_lines))
            .n_genes <- faa_lines |>
                str_split(" \\| ") |>
                map_chr(\(x) x[grepl("^gene=", x)]) |>
                str_remove("^gene=") |>
                unique() |>
                length()
        } else if (.source == "InsectBase" || .source == "GenBank") {
            .n_genes <- faa_lines |>
                str_split(" ") |>
                map_chr(\(x) x[1]) |>
                str_remove("^>") |>
                sub(pattern = "\\.[^\\.]+$", replacement = "") |>
                unique() |>
                length()
        } else {
            .n_genes <- faa_lines |>
                str_remove("^>") |>
                sub(pattern = "\\.[^\\.]+$", replacement = "") |>
                unique() |>
                length()
        }
        tibble(spp_abbrev = .spp, n_prots = .n_prots, n_genes = .n_genes)
    })


# ngenes_df
# Updated on 2023-09-22
# # A tibble: 14 × 3
#    spp_abbrev n_prots n_genes
#    <fct>        <int>   <int>
#  1 Ctenta       20615   20615
#  2 Cripar       16522   16522
#  3 Pvande       17631   17631
#  4 Ppemba       17339   17339
#  5 Tgraci       15561   15499
#  6 Bantar       10853   10853
#  7 Cmarin       21259   21259
#  8 Pakamu       14557   14557
#  9 Pstein       16259   16105
# 10 Csonor       18212   18080
# 11 Cquinq       24531   15094
# 12 Aaegyp       28391   14718
# 13 Asteph       29673   12705
# 14 Mdomes       14215   14215




#' ===========================================================================
# Introns ----
#' ===========================================================================

#' HOG file for N0 node from OrthoFinder:
hog_file_n0 <- paste0(orthofinder_dir, "/Phylogenetic_Hierarchical_Orthogroups/N0.tsv")


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
#' This is used below in `write_introns` and `write_intergenic` functions.
get_trans <- function(.gff_df, .source, .genes = NULL) {
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
        select(seqid, gene, trans, start, end, strand, phase)

    if (!is.null(.genes)) {
        trans_df <- trans_df |>
            filter(gene %in% .genes)
    }

    trans_df <- trans_df |>
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


#' Get intron sizes and optionally write to csv file based on species abbreviation,
#' source database, and list of genes to use.
#' Then return tibble to summarize by species.
get_introns <- function(.spp, .source, .genes, .write = FALSE) {

    # .spp = "Tgraci"; .source = "here"
    # .genes = unique(do.call(c, gnames_df[[.spp]])); .write = FALSE

    .fn <- paste0(features_dir, "/", .spp, "_features.gff3.gz")
    gff_df <- read_tsv(.fn, comment = "#",
                       col_names = c("seqid", "source", "type", "start",
                                     "end", "score", "strand", "phase",
                                     "attributes"),
                       col_types = "ccciicccc", progress = FALSE)

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

    if (.write) {
        out_fn <- sprintf("_data/introns-%s.csv.xz", .spp)
        cat("Writing to ", out_fn, "\n")
        write_csv(intron_df, out_fn)
    }

    intron_df <- intron_df |>
        mutate(spp_abbrev = .spp,
               intron_len = as.integer(end - start + 1),
               log_intron_len = log10(intron_len)) |>
        group_by(spp_abbrev, gene) |>
        summarize(n_introns = n(),
                  log_n_introns = log10(n_introns),
                  intron_len = list(intron_len),
                  log_intron_len = list(log_intron_len),
                  .groups = "drop") |>
        group_by(spp_abbrev) |>
        summarize(mean_intron_len = mean(do.call(c, intron_len)),
                  sum_intron_len = sum(do.call(c, intron_len)),
                  mean_log_intron_len = mean(do.call(c, log_intron_len)),
                  mean_n_introns = mean(n_introns),
                  mean_log_n_introns = mean(log_n_introns))

    return(intron_df)
}


intron_df <- gstats_df |>
    #' This species only had ~58% of genes match to orthogroups, so
    #' I'm removing it from the analyses:
    filter(spp_abbrev != "Cmarin") |>
    select(spp_abbrev, annot_source) |>
    rename(.spp = spp_abbrev, .source = annot_source) |>
    mutate(.genes = map(.spp, \(s) unique(do.call(c, gnames_df[[s]])))) |>
    pmap_dfr(get_introns) |>
    # Have to add this species back in so that this dataframe matches others
    add_row(spp_abbrev = "Cmarin") |>
    mutate(spp_abbrev = factor(spp_abbrev,
                               levels = levels(gstats_df$spp_abbrev))) |>
    arrange(spp_abbrev)


#' ===========================================================================
# Intergenic sequences ----
#' ===========================================================================



#' Get intron sizes and optionally write to csv file based on species abbreviation,
#' source database, and list of genes to use.
#' Then return tibble to summarize by species.
get_intergenic <- function(.spp, .source,
                           .seq_len_df = seq_len_df, .write = FALSE) {

    stopifnot(.source %in% c("VectorBase", "GenBank", "InsectBase", "here"))

    .fn <- paste0(features_dir, "/", .spp, "_features.gff3.gz")
    gff_df <- read_tsv(.fn, comment = "#",
                       col_names = c("seqid", "source", "type", "start",
                                     "end", "score", "strand", "phase",
                                     "attributes"),
                       col_types = "ccciicccc", progress = FALSE)

    trans_df <- get_trans(gff_df, .source) |>
        select(seqid, start, end)

    stopifnot(all(trans_df$end > trans_df$start))

    igl_df <- trans_df |>
        split(~ seqid) |>
        map_dfr(\(.d) {
            .L <- filter(.seq_len_df, spp_abbrev == .spp, seqid == .d$seqid[[1]]) |>
                getElement("length")
            .o <- tibble(seqid = .d$seqid[[1]],
                         interg_len = c(.d$start - lag(.d$end, 1L, 1L) - 1L,
                                        .L - tail(.d$end, 1)))
            return(.o)
        })

    if (.write) {
        out_fn <- sprintf("_data/intergenic-%s.csv.xz", .spp)
        cat("Writing to ", out_fn, "\n")
        write_csv(igl_df, out_fn, progress = FALSE)
    }

    igl_df <- igl_df |>
        mutate(spp_abbrev = .spp) |>
        filter(interg_len > 0) |>
        select(spp_abbrev, seqid, interg_len) |>
        group_by(spp_abbrev) |>
        summarize(sum_interg_len = sum(interg_len), .groups = "drop")

    return(igl_df)
}



interg_df <- gstats_df |>
    select(spp_abbrev, annot_source) |>
    rename(.spp = spp_abbrev, .source = annot_source) |>
    pmap_dfr(get_intergenic)

# interg_df





#' ===========================================================================
# Repeat elements ----
#' ===========================================================================



# make sure this section works

simplify_class <- function(.class) {
    .class <- str_split(.class, "\\/") |>
        map_chr(\(x) x[1])
    # Remove '?' because the RepeatMasker's *.tbl summary includes
    # these with the class minus the '?' modifier
    .class <- str_remove_all(.class, "\\?")
    # rRNA, tRNA, and snRNA are lumped as "Small RNA"
    .class <- ifelse(.class %in% c("rRNA", "tRNA", "snRNA"), "Small_RNA", .class)
    # Unclassified is Unknown or ARTEFACT:
    .class <- ifelse(.class %in% c("Unknown", "ARTEFACT"), "Unclassified", .class)
    return(.class)
}

one_spp_repeats <- function(.spp) {
    .col_names <- c("SW", "p_div", "p_del", "p_ins", "query_seq",
                    "begin", "end", "left", "unk1", "matching_repeat",
                    "class", "begin2", "end2", "left2", "id", "star")
    .fn <- sprintf("%s/%s_repeats.tsv", repeats_dir, .spp)
    repeats <- read_table(.fn, skip = 3, col_types = cols(), progress = FALSE,
                          col_names = .col_names)
    repeat_sums <- repeats |>
        mutate(class = simplify_class(class),
               spp_abbrev = .spp) |>
        group_by(spp_abbrev, class) |>
        summarize(elements = length(unique(id)),
                  length = round(sum((end - begin + 1) *
                                         (1 - (p_del+p_ins) / 100))),
                  .groups = "drop") |>
        arrange(class)
    cat(.spp, "finished\n")
    return(repeat_sums)
}


repeats_df <- c("Aaegyp", "Asteph", "Bantar", "Cmarin", "Cquinq", "Cripar",
                "Csonor", "Ctenta", "Mdomes", "Pakamu", "Ppemba",
                "Pstein", "Pvande", "Tgraci") |>
    safe_mclapply(one_spp_repeats) |>
    do.call(what = bind_rows)

write_csv(repeats_df, "_data/repeats-summary.csv")


repeats_len_df <- repeats_df |>
    select(spp_abbrev, class, length) |>
    pivot_wider(names_from = class, values_from = length) |>
    mutate(spp_abbrev = factor(spp_abbrev, levels = levels(gstats_df$spp_abbrev))) |>
    arrange(spp_abbrev)




#' ===========================================================================
# Add columns to `gstats_df` and write to CSV
#' ===========================================================================

if (all(gstats_df$spp_abbrev == gsize_df$spp_abbrev)) {
    gstats_df$gsize <- gsize_df$gsize
} else stop("not all(gstats_df$spp_abbrev == gsize_df$spp_abbrev)")

if (all(gstats_df$spp_abbrev == ngenes_df$spp_abbrev)) {
    gstats_df$n_genes <- ngenes_df$n_genes
} else stop("not all(gstats_df$spp_abbrev == ngenes_df$spp_abbrev)")


if (all(gstats_df$spp_abbrev == interg_df$spp_abbrev)) {
    gstats_df$sum_interg_len <- interg_df$sum_interg_len
} else stop("not all(gstats_df$spp_abbrev == interg_df$spp_abbrev)")


if (all(gstats_df$spp_abbrev == intron_df$spp_abbrev)) {
    for (x in c("mean_intron_len", "sum_intron_len", "mean_log_intron_len",
                "mean_n_introns", "mean_log_n_introns")) {
        gstats_df[[x]] <- intron_df[[x]]
    }
} else stop("not all(gstats_df$spp_abbrev == intron_df$spp_abbrev)")


if (all(gstats_df$spp_abbrev == repeats_len_df$spp_abbrev)) {
    for (x in colnames(repeats_len_df)[-1]) {
        gstats_df[[x]] <- repeats_len_df[[x]]
    }
} else stop("not all(gstats_df$spp_abbrev == repeats_len_df$spp_abbrev)")

write_csv(gstats_df, "_data/genome-stats.csv")

