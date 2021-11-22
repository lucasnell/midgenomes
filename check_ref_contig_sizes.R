
library(tidyverse)




read_fasta <- function(fn) {

    fasta <- read_lines(fn)

    splits <- tibble(start = which(str_detect(fasta, "^>")),
                     end = c(tail(start, -1) - 1, length(fasta))) %>%
        mutate(start = start + 1)

    seqs <- pmap_chr(splits,
                     function(start, end) str_c(fasta[start:end], collapse = ""))

    return(seqs)
}


# Report genome size, # contigs, N50
report_stats <- function(seqs) {

    lens <- sort(nchar(seqs), decreasing = TRUE)
    i <- which(cumsum(lens) >= (sum(lens) / 2))[1]

    cat(sprintf("size = %.2f Mb\n", sum(lens) / 1e6))
    cat(sprintf("%s contigs\n", prettyNum(length(lens),
                                        big.mark = ",", scientific = FALSE)))
    cat(sprintf("N50 = %.2f kb\n", lens[i] / 1e3))

    invisible(NULL)

}



assembly <- read_fasta("~/_data/assembly_shasta/Assembly.fasta")
polish1 <- read_fasta("~/_data/polished_pepper/polished_hap1.fasta")
polish2 <- read_fasta("~/_data/polished_pepper/polished_hap2.fasta")


report_stats(assembly)
report_stats(polish1)
report_stats(polish2)


#' CANU > REDUNDANS version:
#'
#' size = 112 Mb
#' 1,176 contigs
#' N50 = 169.6 kb
#'
#'
