
library(tidyverse)
library(ape)



#' OrthoFinder directory. It's assumed this directory is unchanged from the
#' output from OrthoFinder.
ofd <- function(...) {
    base_dir <- "~/_data/chir_orthofinder/orthofinder-output"  ## << change as necessary
    base_dir <- paste0(str_remove(base_dir, "/$"), "/")
    dots <- list(...)
    do.call(paste0, c(list(base_dir), dots))
}
#' Protein sequence directory. It's assumed the names within follow the pattern:
#' <species name>.faa
#' and that they have duplicated proteins removed
psd <- function(...) {
    base_dir <- "~/_data/chir_orthofinder/chir_proteins"  ## << change as necessary
    base_dir <- paste0(str_remove(base_dir, "/$"), "/")
    dots <- list(...)
    do.call(paste0, c(list(base_dir), dots))
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


#' When we filter for exactly 1 gene for all species, 3 OGs are repeated.
#' Other OGs are unique.
gnums11 |> nrow()
gnums11 |>
    distinct(OG, .keep_all = TRUE) |>
    nrow()
gnums11 |>
    filter(OG %in% OG[duplicated(OG)]) |>
    select(HOG, OG, starts_with("Gene"))



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


# Convert raw vector to single string:
to_str <- \(x) paste(rawToChar(x), collapse = "")

spp_names <- orig_gnames |> select(Aaegyp:Tgraci) |> colnames() |> sort()

n = spp_names[1]

# for (n in spp_names) {

rfa <- read.FASTA(psd(n, ".faa"), "AA")
stopifnot(!any(duplicated(names(rfa))))

lens <- map_int(rfa, length) |> as.integer()
sum(duplicated(lens))
mean(duplicated(lens))

#' Unique lengths that aren't duplicated:
nond_lens <- lens[duplicated(lens)] |> unique()
length(nond_lens)

#' Now go through each duplicated length and see if any proteins themselves
#' are duplicated.
dup_seq <- logical(length(rfa))

for (len in nond_lens) {
    len_inds <- which(lens == len)
    seqs <- rfa[len_inds] |> map_chr(to_str) |> as.character()
    seq_dups <- duplicated(seqs)
    if (any(seq_dups)) dup_seq[len_inds[seq_dups]] <- TRUE
}

# cat(sprintf("%s\n  total   = %i\n  percent = %.1f\n\n", n, sum(dup_seq),
#             100 * mean(dup_seq)))

nond_rfa <- rfa[!dup_seq]

write.FASTA(nond_rfa, psd(n, "_proteins_nond.faa"))

# }






