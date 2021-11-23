
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


# Results from just SHASTA > purge_dups
hap_np <- read_fasta("/Users/lucasnell/_data/haploid_purge_dups_nopolish/seqs/contigs_shasta.hap.fa")
hap_p_np <- read_fasta("/Users/lucasnell/_data/haploid_purge_dups_nopolish/seqs/contigs_shasta.purged.fa")

report_stats(hap_np)
report_stats(hap_p_np)


# Results from SHASTA > PEPPER > purge_dups (haplotype 1)
hap1 <- read_fasta("/Users/lucasnell/_data/haploid_purge_dups_pepper1/seqs/polished_hap1.hap.fa")
hap_p1 <- read_fasta("/Users/lucasnell/_data/haploid_purge_dups_pepper1/seqs/polished_hap1.purged.fa")

report_stats(hap1)
report_stats(hap_p1)


# SHASTA > PEPPER > purge_dups (haplotype 1) > purge_dups (purged contigs + haplotype 2)
hap <- read_fasta("/Users/lucasnell/_data/haploid_purge_dups/seqs/pepper_haps.hap.fa")
hap_p <- read_fasta("/Users/lucasnell/_data/haploid_purge_dups/seqs/pepper_haps.purged.fa")

report_stats(hap)
report_stats(hap_p)




all_haps <- read_fasta("/Users/lucasnell/_data/haploid_purge_dups/pepper_haps.fasta")
report_stats(all_haps)

