
#'
#' Read an assembly from a FASTA file and look at its contig / scaffold
#' size stats.
#'

library(tidyverse)
library(parallel)

options(mc.cores = max(parallel::detectCores()-2L, 1L))



read_fasta <- function(fn) {

    fasta <- read_lines(fn)

    splits <- tibble(start = which(str_detect(fasta, "^>")),
                     end = c(tail(start, -1) - 1, length(fasta))) %>%
        mutate(start = start + 1)

    seqs <- pmap_chr(splits,
                     function(start, end) str_c(fasta[start:end],
                                                collapse = ""))

    nms <- map_chr(splits$start - 1, ~ fasta[.x]) %>%
        str_remove("^>")

    names(seqs) <- nms
    class(seqs) <- "fasta"

    return(seqs)
}

print.fasta <- function(x, ...) {
    cat(sprintf("fasta file with %i sequences", length(x)))
}

len_str <- function(len, .digits) {
    fmt <- paste0("%.", .digits, "f")
    if (len > 1e9) {
        z <- sprintf(paste(fmt, "Gb"), len / 1e9)
    } else if (len > 1e6) {
        z <- sprintf(paste(fmt, "Mb"), len / 1e6)
    } else if (len > 1e3) {
        z <- sprintf(paste(fmt, "kb"), len / 1e3)
    } else {
        z <- sprintf(paste(fmt, "b"), len)
    }
    return(z)
}


# Report genome size, # contigs, N50
report_stats <- function(seqs, .digits = 2) {

    lens <- sort(nchar(seqs), decreasing = TRUE)
    i <- which(cumsum(lens) >= (sum(lens) / 2))[1]

    pn <- function(x) prettyNum(x, big.mark = ",", scientific = FALSE)

    cat(sprintf("size = %s\n", len_str(sum(lens), .digits)))
    cat(sprintf("%s contigs\n", pn(length(lens))))
    cat(sprintf("N50 = %s\n", len_str(lens[i], .digits)))
    cat(sprintf("min = %s\n", len_str(min(lens), .digits)))
    cat(sprintf("max = %s\n", len_str(max(lens), .digits)))
    n_count <- len_str(sum(str_count(seqs, "N|n")), .digits)
    cat(sprintf("total N = %s\n", n_count))

    invisible(NULL)

}

