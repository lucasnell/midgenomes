
# Looking at p-vals of SNAPE output and masked files

library(tidyverse)

# z <- read_lines("~/_data/Vik-19_S41_snape_masked_P-MINOR.txt.gz") %>%
#     as.numeric()
# sum(is.na(z))
# hist(z)
# zz <- pmin(z, 1 - z)
# hist(zz)
# range(zz)

# scaffold-1	1335	C	1	0	70	70	C	0.03412	2.032e-06	0.004035
x <- read_tsv("~/_data/Vik-19_S41_snape.txt.gz",
              col_names = c("seq", "pos", "ref",
                            "n_ref", "n_alt",
                            "q_ref", "q_alt",
                            "nt",
                            "p1", "p2", "p3"),
              col_types = "ccccccccccc")
set.seed(1)
x <- x[sort(sample.int(nrow(x), nrow(x) %/% 2)),]


xx <- x %>%
    filter(n_ref != "*") %>%
    select(n_ref, n_alt, nt, p1, p2, p3) %>%
    mutate(across(c(p1, p2, p3), as.numeric)) %>%
    mutate(across(c(n_ref, n_alt), as.integer)) %>%
    filter(n_alt > 0 | n_ref > 0)



xx %>%
    filter(p1 > 0.9, p2 > 0.1) %>%
    filter(n_ref > 0)


xx %>%
    # filter(p1 > 0.9) %>%
    filter(p1 > 0.9, p2 < 0.1) %>%
    # filter(!(p1 <= 0.9 & p1 >= (1 - 0.9))) %>%
    # mutate(is_zero = ifelse(n_ref == 0 | n_alt == 0, 1, 0) %>%
    #            factor(levels = 0:1, labels = c("no", "yes"))) %>%
    # ggplot(aes(n_ref, n_alt)) +
    # geom_point(aes(color = is_zero), alpha = 0.5) +
    # scale_color_manual(values = c("black", "red")) +
    # theme_minimal()
    .[["p3"]] %>%
    hist()


# info <= float(options.maxsnape) and info >= float(1 - float(options.maxsnape))


library(tidyverse)
library(readxl)

# Excel sheet provided by biotech center
all_names <- read_excel("~/Box Sync/midges/fastq-stats.xlsx", "trimmed") %>%
    filter(to_use == 1) %>%
    .[["Sample Name"]] %>%
    str_remove_all("trimmed_|_L002_R1_001|_L002_R2_001") %>%
    unique()

# Spatial samples:
spat_filt <- ! Reduce(`|`, lapply(c("^SN-", "^KS-", "^H-"), grepl, x = all_names))
spat_samples <- all_names[spat_filt]
paste0(spat_samples, "*", collapse = " ")
paste0(spat_samples, collapse = " ")
# Temporal samples:
temp_samples <- all_names[!spat_filt]
paste0(temp_samples, "*", collapse = " ")
paste0(temp_samples, collapse = " ")


# To get sample names and # adults for use in npstat
read_csv("~/Box Sync/midges/full-DNA-info.csv",
         col_types = "cfcddidcccddidiDccldd") %>%
    filter(to_use == 1) %>%
    distinct(biotech_id, n_adults) %>%
    mutate(win = map(1:n(), ~ c(1L, 10L))) %>%
    unnest(win) %>%
    arrange(desc(win), biotech_id) %>%
    format_delim(", ") %>%
    str_replace_all(",", ", ") %>%
    cat()


# And for quickmerge:




crossing(ass1 = c("contigs_flye_ont-polish_purgedups.fasta",
                  "contigs_necat_ont-polish_purgedups.fasta",
                  "contigs_next_ont-polish.fasta",
                  "contigs_smart_ont-polish_purgedups.fasta"),
         ass2 = ass1,
         # 2 for saving all output, 1 for saving just FASTA, 0 for no saving
         save_out = 1L) %>%
    filter(ass1 != ass2) %>%
    arrange(ass1, ass2) %>%
    mutate(nm1 = str_remove_all(ass1, "^contigs_|_ont-polish|_purgedups|.fasta$"),
           nm2 = str_remove_all(ass2, "^contigs_|_ont-polish|_purgedups|.fasta$"),
           out_name = sprintf("contigs_merged_%s_%s", nm1, nm2)) %>%
    select(ass1, ass2, out_name, save_out) %>%
    format_delim(", ") %>%
    str_replace_all(",", ", ") %>%
    cat()

crossing(ass1 = c("contigs_merged_next_flye.fasta",
                  "contigs_merged_next_smart.fasta",
                  "contigs_merged_smart_necat.fasta",
                  "contigs_merged_smart_next.fasta"),
         ass2 = ass1,
         # 2 for saving all output, 1 for saving just FASTA, 0 for no saving
         save_out = 1L) %>%
    filter(ass1 != ass2 & (grepl("merged_next_flye", ass1) |
                               grepl("merged_next_flye", ass2))) %>%
    arrange(ass1, ass2) %>%
    mutate(nm1 = str_remove_all(ass1, "^contigs_merged_|.fasta$") %>%
               str_replace("_", "-"),
           nm2 = str_remove_all(ass2, "^contigs_merged_|.fasta$") %>%
               str_replace("_", "-"),
           out_name = sprintf("contigs_merged_%s_%s", nm1, nm2)) %>%
    select(ass1, ass2, out_name, save_out) %>%
    format_delim(", ") %>%
    str_replace_all(",", ", ") %>%
    cat()




# second round where we tried both combos of smart + smart--necat assemblies
crossing(ass1 = c("merged_smart_necat.fasta",
                  "contigs_smart_ont-polish_nextpolish_purgedups.fasta"),
         ass2 = ass1,
         # 1 for saving output, 0 for not
         save_out = 0L) %>%
    filter(ass1 != ass2) %>%
    mutate(nm1 = str_remove_all(ass1, "^merged_|^contigs_|_ont-polish_nextpolish_purgedups|.fasta$"),
           nm2 = str_remove_all(ass2, "^merged_|^contigs_|_ont-polish_nextpolish_purgedups|.fasta$"),
           out_name = sprintf("merged_%s__%s", nm1, nm2)) %>%
    select(-nm1, -nm2) %>%
    select(ass1, ass2, out_name, save_out) %>%
    format_delim(", ") %>%
    str_replace_all(",", ", ") %>%
    cat()





# Round 2: NextDenovo and best from previous version:
crossing(ass1 = c("merged_smart__smart_necat.fasta",
                  "contigs_next_ont-polish_nextpolish.fasta"),
         ass2 = ass1) %>%
    filter(ass1 != ass2) %>%
    mutate(nm1 = str_remove_all(ass1, "^contigs_|^merged_|_ont-polish_nextpolish|.fasta$"),
           nm2 = str_remove_all(ass2, "^contigs_|^merged_|_ont-polish_nextpolish|.fasta$"),
           nm = sprintf("merged_%s___%s.fasta", nm1, nm2)) %>%
    select(-nm1, -nm2) %>%
    mutate(save_out = 1L) %>%
    format_delim(", ") %>%
    str_replace_all(",", ", ") %>%
    cat()




# LongStitch:

crossing(ass = c("contigs_merged_next_smart"),
         k_ntLink = c(24, 32, 40),
         w = c(100, 250, 500),
         ark= 0:1,
         save_out = 0L) %>%
    format_delim(", ") %>%
    str_replace_all(",", ", ") %>%
    cat()

