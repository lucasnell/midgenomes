
library(tidyverse)
library(ape)



#' OrthoFinder directory. It's assumed this directory is unchanged from the
#' output from OrthoFinder, except for separate folders added to the parent
#' directory.
ofd <- function(...) {
    dots <- list(...)
    base_dir <- "~/_data/chir_orthofinder/orthofinder-output"  ## << change as necessary
    base_dir <- str_remove(base_dir, "/$")
    if (length(dots) == 0) return(base_dir)
    do.call(paste0, c(list(base_dir, "/"), dots))
}
#' Protein sequence directory. It's assumed the names within follow the pattern:
#' <species name>.faa
#' and that they have duplicated proteins removed
psd <- function(...) {
    dots <- list(...)
    base_dir <- "~/_data/chir_orthofinder/chir_proteins"  ## << change as necessary
    base_dir <- str_remove(base_dir, "/$")
    if (length(dots) == 0) return(base_dir)
    do.call(paste0, c(list(base_dir, "/"), dots))
}

orig_gnames <- read_tsv(ofd("Phylogenetic_Hierarchical_Orthogroups/N0.tsv"),
                    col_types = cols()) |>
    #' This species only had ~58% of genes match to orthogroups, so
    #' I'm removing it from the analyses:
    select(-Cmarin) |>
    mutate(across(Aaegyp:Tgraci, \(x) {
        z <- rep(list(character(0)), length(x))
        z[!is.na(x)] <- str_split(x[!is.na(x)], ", ")
        return(z)
    }))


orig_gnums <- orig_gnames |>
    mutate(across(Aaegyp:Tgraci, \(x) {
        z <- map_int(x, length)
        return(z)
    }))

#' Filter for just 1-to-1 orthogroups:
gnames11 <- orig_gnames |>
    filter(if_all(Aaegyp:Tgraci, \(x) map_lgl(x, \(y) length(y) == 1))) |>
    # Because all species have only 1 gene here, then we can unnest these cols:
    unnest(Aaegyp:Tgraci)
gnums11 <- orig_gnums |>
    filter(if_all(Aaegyp:Tgraci, \(x) x == 1))



#' ===========================================================================
#' ===========================================================================
#'
#' Can I use the Orthogroup_Sequences folder?
#'
#' ===========================================================================
#' ===========================================================================


#' The Orthogroup_Sequences contains sequences for orthogroups, but I *think*
#' these sequences might be from the "normal" orthogroups, not the
#' "Phylogenetic Hierarchical Orthogroups" that OrthoFinder now recommends
#' using because they're more accurate.

#' When we filter for exactly 1 gene for all species, 3 OGs are repeated.
#' Other OGs are unique.
gnums11 |> nrow()
gnums11 |>
    distinct(OG, .keep_all = TRUE) |>
    nrow()
gnums11 |>
    filter(OG %in% OG[duplicated(OG)]) |>
    select(HOG, OG, starts_with("Gene"))


#' Number of sequences and sequence names for each OG from the
#' Orthogroup_Sequences folder
os_df <- gnums11 |>
    distinct(OG) |>
    mutate(os_names = map(OG,
                       \(og) {
                           fn <- ofd("Orthogroup_Sequences/", og, ".fa")
                           read_lines(fn) |>
                               keep(\(x) grepl("^>", x)) |>
                               str_remove_all("^>")
                       }),
           os_n = map_int(os_names, length))

#' I could go on, but we can see here that many of the Orthogroups from
#' Orthogroup_Sequences that should be 1-to-1 based on the HOG categorization
#' contain a number of sequences (i.e., genes) that are not equal to
#' the number of species:

# Numbers of genes:
os_df |>
    arrange(desc(os_n))
# Number of species:
orig_gnames |>
    select(Aaegyp:Tgraci) |>
    ncol()




#' ===========================================================================
#' ===========================================================================
#'
#' Extracting sequences for 1-to-1 HOGs:
#'
#' ===========================================================================
#' ===========================================================================



#' Write HOG sequences to directory.
#'
#' @param hog_df Data frame from one of the tsv files within
#'        `Phylogenetic_Hierarchical_Orthogroups`.
#' @param hog_seq_dir Directory where sequence faa files should be written.
#' @param proteins List of `AAbin` objects containing all proteins for all
#'        species.
#'
write_hog_seqs <- function(hog_df, hog_seq_dir, proteins) {

    stopifnot(inherits(hog_df, "data.frame"))
    stopifnot(inherits(hog_seq_dir, "character"))
    stopifnot(inherits(proteins, "list"))
    stopifnot(!is.null(names(proteins)))
    stopifnot(all(sapply(proteins, inherits, what = "AAbin")))

    if ("OG" %in% colnames(hog_df)) hog_df[["OG"]] <- NULL
    if ("Gene Tree Parent Clade" %in% colnames(hog_df)) {
        hog_df[["Gene Tree Parent Clade"]] <- NULL
    }

    spp_names <- hog_df |>
        select(-HOG) |>
        colnames()

    stopifnot(all(spp_names %in% names(proteins)))

    if (!dir.exists(hog_seq_dir)) dir.create(hog_seq_dir, recursive = TRUE)

    # Convert raw vector to single string:
    to_str <- \(x) paste(rawToChar(x), collapse = "")

    for (i in 1:nrow(hog_df)) {

        hog <- hog_df$HOG[[i]]
        n_genes <- length(unlist(hog_df[i,spp_names]))
        hog_lines <- character(2 * n_genes)
        j <- 1

        for (spp in spp_names) {

            gene_vec <- hog_df[[spp]][[i]]

            for (gene in gene_vec) {

                hog_lines[j] <- sprintf(">%s__%s", spp, gene)
                j <- j + 1

                hog_lines[j] <- to_str(proteins[[spp]][[gene]])
                j <- j + 1

            }

        }

        write_lines(hog_lines, paste0(hog_seq_dir, "/", hog, ".faa"))

    }

    cat(sprintf("Finished writing to %s\n", hog_seq_dir))

    invisible(NULL)

}



#' All species names (minus Cmarin bc of its poor matching)
all_spp_names <- orig_gnames |> select(Aaegyp:Tgraci) |> colnames() |> sort()

# All proteins from all species
proteins <- map(all_spp_names, \(n) {
    rfa <- read.FASTA(psd(n, ".faa"), "AA")
    stopifnot(all(!duplicated(names(rfa))))
    return(rfa)
})
names(proteins) <- all_spp_names


#' ============================================================================
#' ============================================================================
#'
#' Root node POGs
#'
#' ============================================================================
#' ============================================================================

write_hog_seqs(gnames11, ofd("Single_Copy_HOG_Sequences/N0"), proteins)




#' ============================================================================
#' ============================================================================
#'
#' Non-root POGs
#'
#' ============================================================================
#' ============================================================================

# Species tree with node labels:
tr <- read.tree(ofd("Species_Tree/SpeciesTree_rooted_node_labels.txt"))

#' To convert back to full names:
spp_name_map <- list("Aaegyp" = "Aedes aegypti",
                     "Asteph" = "Anopheles stephensi",
                     "Bantar" = "Belgica antarctica",
                     "Cripar" = "Chironomus riparius",
                     "Ctenta" = "Chironomus tentans",
                     "Cmarin" = "Clunio marinus",
                     "Cquinq" = "Culex quinquefasciatus",
                     "Csonor" = "Culicoides sonorensis",
                     "Mdomes" = "Musca domestica",
                     "Pstein" = "Parochlus steinenii",
                     "Ppemba" = "Polypedilum pembai",
                     "Pvande" = "Polypedilum vanderplanki",
                     "Pakamu" = "Propsilocerus akamusi",
                     "Tgraci" = "Tanytarsus gracilentus")

tr$tip.label <- map_chr(tr$tip.label, \(a) spp_name_map[[a]])

plot(tr, show.node.label = TRUE, no.margin = TRUE, label.offset = 0.1)

# pdf("orthofinder-species-tree.pdf", width = 8, height = 5)
# plot(tr, show.node.label = TRUE, no.margin = TRUE, label.offset = 0.1)
# dev.off()


#'
#' From this, other nodes that may be of interest, in order from most to least
#' inclusive, are
#' - N1: Culicomorpha (no Musca domestica)
#' - N3: just families Chironomidae and Ceratopogonidae (no Musca domestica
#'       or family Culicidae)
#' - N5: family Chironomidae
#'

for (node in c("N1", "N3", "N5")) {

    hogdf <- read_tsv(ofd("Phylogenetic_Hierarchical_Orthogroups/", node, ".tsv"),
             col_types = cols()) |>
        #' Remove species not represented by this node:
        select(\(x) any(!is.na(x))) |>
        #' This species only had ~58% of genes match to orthogroups, so
        #' I'm removing it from the analyses:
        select(-Cmarin) |>
        select(-OG, -starts_with("Gene")) |>
        mutate(across(-HOG, \(x) {
            z <- rep(list(character(0)), length(x))
            z[!is.na(x)] <- str_split(x[!is.na(x)], ", ")
            return(z)
        })) |>
        filter(if_all(-HOG, \(x) map_lgl(x, \(y) length(y) == 1)))

    write_hog_seqs(hogdf, ofd("Single_Copy_HOG_Sequences/", node), proteins)
}





#' ------------------------------------------------------
#' ------------------------------------------------------
#'
#' Number of single-copy HOGs by node:
#'
#' N0 = 1,735
#' N1 = 1,922
#' N3 = 2,460
#' N5 = 3,443
#'

for (node in c("N0", "N1", "N3", "N5")) {
    hogdf <- read_tsv(ofd("Phylogenetic_Hierarchical_Orthogroups/", node, ".tsv"),
                      col_types = cols()) |>
        #' Remove species not represented by this node:
        select(\(x) any(!is.na(x))) |>
        #' This species only had ~58% of genes match to orthogroups, so
        #' I'm removing it from the analyses:
        select(-Cmarin) |>
        select(-OG, -starts_with("Gene")) |>
        mutate(across(-HOG, \(x) {
            z <- rep(list(character(0)), length(x))
            z[!is.na(x)] <- str_split(x[!is.na(x)], ", ")
            return(z)
        })) |>
        filter(if_all(-HOG, \(x) map_lgl(x, \(y) length(y) == 1)))
    cat(sprintf("%s = %s\n", node, prettyNum(nrow(hogdf), big.mark = ",")))
}

#' ------------------------------------------------------
#' ------------------------------------------------------
#'
#' Number of total HOGs by node:
#'
#' N0 = 19,575
#' N1 = 19,534
#' N3 = 17,918
#' N5 = 16,794
#'
#'

for (node in c("N0", "N1", "N3", "N5")) {
    hogdf <- read_tsv(ofd("Phylogenetic_Hierarchical_Orthogroups/", node, ".tsv"),
                      col_types = cols()) |>
        #' This species only had ~58% of genes match to orthogroups, so
        #' I'm removing it from all analyses.
        select(-Cmarin) |>
        #' The rest of this section is to remove HOGs only containing genes
        #' from Cmarin:
        select(-OG, -starts_with("Gene")) |>
        mutate(across(-HOG, \(x) {
            z <- integer(length(x))
            z[!is.na(x)] <- str_count(x[!is.na(x)], ", ") + 1L
            return(z)
        })) |>
        filter(if_any(-HOG, \(x) x > 0))
    cat(sprintf("%s = %s\n", node, prettyNum(nrow(hogdf), big.mark = ",")))
}

