
#'
#' Create table with summary stats about alignments of proteins matching
#' OrthoDB Diptera proteins.
#'

source("_R/00-preamble.R")

read_fasta <- function(fn) {
    z <- read_lines(fn, progress = FALSE)
    heads <- which(grepl("^>", z))
    nh <- length(heads)
    sl <- c(map(2:nh, \(i) (heads[i-1]+1L):(heads[i]-1L)),
            list((heads[nh]+1L):length(z)))
    zz <- map_chr(sl, \(inds) paste0(z[inds], collapse = ""))
    names(zz) <- z[heads] |>
        str_remove_all(">") |>
        str_split("__") |>
        map_chr(\(x) x[[1]])
    return(zz)
}

# Number of genes in species' annotations matching OrthoDB single-copy genes:
gene_df <- paste0(dirs$odb, "/shared_odb/single-copy-gene-counts.tsv") |>
    read_tsv(col_types = "ci") |>
    mutate(species = expand_spp(species, to_fct = TRUE))

# Percent missing in alignments (before trimming):
miss_df <- list.files(dirs$mafft, "*.faa", full.names = TRUE) |>
    safe_mclapply(\(f) {
        faa <- read_fasta(f)
        tibble(species = expand_spp(names(faa), TRUE),
               gene = str_remove(basename(f), "\\.faa"),
               total = nchar(faa),
               missing = str_count(faa, "-"))
    }) |>
    do.call(what = bind_rows) |>
    group_by(species) |>
    summarize(missing_p = sum(missing) / sum(total) * 100)

full_join(gene_df, miss_df, by = "species") |>
    arrange(species)



