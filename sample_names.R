

library(tidyverse)
library(readxl)

# Excel sheet provided by biotech center
all_names <- read_excel("~/_data/biotech-files/DNA/coverage.xlsx") %>%
    .[["Sample Name"]] %>%
    str_remove_all("_L002_R1_001|_L002_R2_001") %>%
    unique()

# All spatial samples:
spat_filt <- ! Reduce(`|`, lapply(c("^SN-", "^KS-", "^H-"), grepl, x = all_names))
spat_samples <- all_names[spat_filt]
paste0(spat_samples, "*", collapse = " ")
# Spatial samples started to be processed:
proc_rn <- c("Blik-19_S6", "Lys-19_S14", "MyBR-19_S17", "MyKS-19-B_S18", "MySN-19_S19")
paste0(proc_rn, "*", collapse = " ")
# Spatial samples not yet processed:
spat_nop_filt <- spat_filt & ! all_names %in% proc_rn
paste0(all_names[spat_nop_filt], "*", collapse = " ")



