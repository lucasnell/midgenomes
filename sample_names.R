

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


# To get sample names and # adults for use in SNAPE-pooled
read_csv("~/Box Sync/midges/full-DNA-info.csv") %>%
    distinct(biotech_id, n_adults) %>%
    format_delim(", ") %>%
    str_replace_all(",", ", ") %>%
    cat()


