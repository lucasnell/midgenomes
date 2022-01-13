
library(tidyverse)




read_fasta <- function(fn) {

    fasta <- read_lines(fn)

    splits <- tibble(start = which(str_detect(fasta, "^>")),
                     end = c(tail(start, -1) - 1, length(fasta))) %>%
        mutate(start = start + 1)

    seqs <- pmap_chr(splits,
                     function(start, end) str_c(fasta[start:end],
                                                collapse = ""))

    nms <- map_chr(splits$start - 1, ~ fasta[.x])

    names(seqs) <- nms

    return(seqs)
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
    cat(sprintf("total N = %s\n", len_str(sum(str_count(seqs, "N")), .digits)))

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


# Results from SHASTA > PEPPER > purge_dups
hap <- read_fasta("~/_data/haploid_purge_dups/seqs/polished_hap1.purged.fa")
report_stats(hap)

# Results from SHASTA > PEPPER > purge_dups > LongStitch
scaffs_ls <- read_fasta(paste0("~/_data/longstitch-then-rna/scaffold_longstitch/",
                            "scaffold_longstitch.fasta"))
report_stats(scaffs_ls)

# Results from SHASTA > PEPPER > purge_dups > LongStitch > BESST_RNA
scaffs_ls_besst <- read_fasta(paste0("~/_data/longstitch-then-rna/",
                                   "scaffold_besst/pass1/Scaffolds-pass1.fa"))
report_stats(scaffs_ls_besst)


# Results from SHASTA > PEPPER > purge_dups > BESST_RNA
scaffs_besst <- read_fasta(paste0("~/_data/scaffold_besst.fasta.gz"))
report_stats(scaffs_besst)

# Results from SHASTA > PEPPER > purge_dups > P_RNA_Scaffolder
scaffs_p_rna <- read_fasta("~/_data/scaffolds_p_rna/P_RNA_scaffold.fasta")
report_stats(scaffs_p_rna)


# Results from SHASTA > PEPPER > purge_dups > BESST_RNA > LongStitch
scaffs_besst_ls <- read_fasta(paste0("~/_data/scaffold_longstitch_besst/",
                                     "scaffold_longstitch_besst.fasta"))
report_stats(scaffs_besst_ls)


# Results from SHASTA > PEPPER > purge_dups > P_RNA_Scaffolder > LongStitch
scaffs_p_rna_ls <- read_fasta("~/_data/scaffolds_longstitch_p_rna.fasta.gz")
report_stats(scaffs_p_rna_ls)


