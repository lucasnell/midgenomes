
#'
#' ```bash
#' # Beforehand:
#' # Download tar files from `fastp-dna` step into the `~/_data/trimmed_reports`
#' # directory.
#' # Install `multiqc` onto working conda environment (v 1.12).
#'
#'
#' # Go to this directory of reports and un-tar all files
#' cd ~/_data/trimmed_reports
#' for f in *.tar.gz; do tar -xzf $f; done
#' rm *.tar.gz
#'
#' conda activate working-env
#'
#' # Adjust names for fastp reports so multiqc can parse sample names:
#' for f in trimmed_*; do
#'     cd $f
#'     mv fastp.html ${f}.html
#'     mv fastp.json ${f}.json
#'     cd ..
#' done
#'
#' # Run MultiQC
#' multiqc .
#'
#' # If no issues above, open the report
#' open multiqc_report.html
#' ```
#'


library(tidyverse)
library(readxl)


# Heatmap of pass / fail / warning output from multiqc report:
heatmap <- read_csv("~/Box Sync/midges/fastqc-status-check-heatmap.csv") %>%
    rename(sample = `Series 1 (y)`,
           stat = `Section Name`,
           flag = `Series 1 (value)`) %>%
    mutate(flag = case_when(flag == 1 ~ "pass",
                            flag == 0.5 ~ "warn",
                            flag == 0.25 ~ "fail",
                            TRUE ~ NA_character_)) %>%
    mutate(stat = case_when(stat == "Sequence Duplication Levels" ~ "duplicate_flag",
                            stat == "Overrepresented Sequences" ~ "overrep_flag",
                            stat == "Adapter Content" ~ "adapter_flag",
                            TRUE ~ stat)) %>%
    filter(stat %in% c("duplicate_flag", "overrep_flag", "adapter_flag")) %>%
    pivot_wider(names_from = stat, values_from = flag)

# Main table of fastq info:
stat_df <- read_excel("~/Box Sync/midges/fastq-stats.xlsx", sheet = "trimmed") %>%
    rename(sample = `Sample Name`)

# They *should* be the same order, but better to check:
identical(stat_df$sample, heatmap$sample)


# write_csv(left_join(stat_df, heatmap, by = "sample"),
#           "~/Box Sync/midges/fastq-stats-with-flags.csv")



#' ```bash
#' export OUT_CSV=mpileup_summs.csv
#' for f in *_mpileup.txt.gz; do
#'     echo -n ${f} "started..."
#'     echo -n ${f/_mpileup.txt.gz/}"," >> ${OUT_CSV}
#'     gunzip -c ${f} | cut -f 4 | awk '{s+=$1} END {printf "%.0f\n", s}' \
#'         >> ${OUT_CSV}
#'     echo " and finished"
#' done
#' ```


# Total sequencing output across entire genome based on mpileup files
final_out_df <- read_csv("~/_data/final_coverage/mpileup_summs.csv",
                         col_names = c("biotech_id", "total_seq"),
                         col_types = "cd")

# Add biotech_id col to match
stat_df <- stat_df %>%
    mutate(biotech_id = str_remove_all(sample, "trimmed_|_L002_R2_001|_L002_R1_001"))

identical(sort(unique(stat_df$biotech_id)), final_out_df$biotech_id)

# Dividing by two bc the `stat_df` is split by FASTQ file (2 per sample),
# while `final_out_df` is split by BAM file (1 per sample)
stat_final_out_df <- left_join(stat_df,
                               mutate(final_out_df, total_seq = total_seq / 2),
                               by = "biotech_id")
# write_csv(stat_final_out_df,
#           "~/Box Sync/midges/fastq-stats-with-final-out.csv")


# Now combine these in the excel file.
# To see which are to be analyzed:

stat_df <- read_excel("~/Box Sync/midges/fastq-stats.xlsx", sheet = "trimmed") %>%
    rename(sample = `Sample Name`) %>%
    mutate(biotech_id = str_remove_all(sample, "trimmed_|_L002_R2_001|_L002_R1_001")) %>%
    filter(to_use == 1)

# Make sure there aren't any loners:
stat_df %>%
    group_by(biotech_id) %>%
    summarize(N = n(), .groups = "drop") %>%
    .[["N"]] %>% unique()

an_names <- stat_df %>%
    .[["biotech_id"]] %>%
    unique()

temp_filt <- Reduce(`|`, lapply(c("^SN-", "^KS-", "^H-"), grepl, x = an_names))
temp_samples <- an_names[temp_filt]
paste0(temp_samples, "*", collapse = " ")


