
library(tidyverse)
library(readxl)

# Excel sheet provided by biotech center
all_names <- read_excel("~/Box Sync/midges/fastq-stats.xlsx", "trimmed") %>%
    .[["Sample Name"]] %>%
    str_remove_all("trimmed_|_L002_R1_001|_L002_R2_001") %>%
    unique()

# Easy lookup from simpler names to full ones (e.g., "MyBR-19" to "MyBR-19_S17")
simple_names <- all_names %>%
    as.list() %>%
    set_names(str_split(all_names, "_S") %>% map_chr(~ .x[[1]]))


# emergence trap locations:
trap_df <- read_csv("~/Box Sync/midges/location_data/trap_locations.csv",
         col_types = cols()) %>%
    # Simplify site names:
    mutate(trap = case_when(trap == "Syðri Neslönd" ~ "SN",
                            trap == "Kálfaströnd" ~ "KS",
                            trap == "Haganes" ~ "H",
                            trap == "Vindbelgur" ~ "V",
                            TRUE ~ NA_character_))


dd_ts <- read_excel("~/Box Sync/midges/tany_identifying/tany_from_archive.xlsx",
           sheet = "data-entry") %>%
    mutate(date = as.Date(sprintf("%i-%i-%i", col_year, col_month, col_day)),
           contam_pilot = ifelse(is.na(contam_pilot), FALSE, contam_pilot),
           lake = "Mývatn") %>%
    filter(!is.na(biotech_id))

#' These samples were prepped together:
#'
#' KS-2-90 and KS-2-90-B
#' SN-2-91 and SN-2-91-B
#' KS-2-83 and KS-2-83-B

for (x in c("KS-2-90", "SN-2-91", "KS-2-83")) {
    nm <- dd_ts$num_males[dd_ts$biotech_id == paste0(x, "-B")]
    nm2 <- dd_ts$num_males[dd_ts$biotech_id == x]
    dd_ts$num_males[dd_ts$biotech_id == paste0(x, "-B")] <- nm + nm2
}

dd_ts <- dd_ts %>%
    filter(biotech_id %in% names(simple_names)) %>%
    mutate(biotech_id = map_chr(biotech_id, ~ simple_names[[.x]])) %>%
    select(biotech_id, num_males, date, lake, site, contam_pilot) %>%
    rename(n_adults = num_males) %>%
    mutate(lat = map_dbl(site, ~ trap_df$lat[trap_df$trap == .x]),
           lon = map_dbl(site, ~ trap_df$lon[trap_df$trap == .x]))

# Should be FALSE
any(map_lgl(dd_ts, ~ any(is.na(.))))



dd_spat <- read_excel("~/Box Sync/midges/tany_identifying/tany_other_lakes.xlsx",
           sheet = "data-entry") %>%
    mutate(date = as.Date(sprintf("2019-%i-%i", col_month, col_day)),
           contam_pilot = ifelse(is.na(contam_pilot), FALSE, contam_pilot)) %>%
    rename(month = col_month, day = col_day) %>%
    filter(n_tany > 0) %>%
    filter(!is.na(biotech_id)) %>%
    filter(sequenced) %>%
    filter(biotech_id %in% names(simple_names)) %>%
    select(biotech_id, date, lake, sample, month, day, n_tany, contam_pilot) %>%
    # Had to add this bc in the "all-samples" sheet, before June 20th,
    # when only a single sample was taken from a site on a particular day,
    # I called it sample # 1.
    # After June 20th, I just left it blank.
    # In the "data-entry" sheet, I just always left these blank.
    mutate(sample = ifelse(is.na(sample) & month == 6 & day < 20,
                           1, sample)) %>%
    left_join(read_excel("~/Box Sync/midges/tany_identifying/tany_other_lakes.xlsx",
                         sheet = "all-samples", na = c("", "NA")) %>%
                  select(lake, sample, month, day, lat, lon) %>%
                  mutate(month = map_int(month, ~ which(month.name == .x))),
              by = c("lake", "sample", "month", "day")) %>%
    mutate(biotech_id = map_chr(biotech_id, ~ simple_names[[.x]])) %>%
    select(biotech_id, n_tany, date, lake, contam_pilot, lat, lon) %>%
    mutate(site = map_chr(str_split(lake, " - "), ~ tail(.x, 1)),
           lake = map_chr(str_split(lake, " - "), ~ .x[[1]]),
           site = case_when(lake != "Mývatn" ~ NA_character_,
                            site == "Mývatn" ~ "BR",
                            TRUE ~ site)) %>%
    # "Blonduos" is the lake's nearest town; the lake is "Grafarvatn"
    mutate(lake = ifelse(lake == "Blonduos", "Grafarvatn", lake)) %>%
    rename(n_adults = n_tany)




dd_other_info <- bind_rows(dd_ts, dd_spat) %>%
    mutate(n_adults = as.integer(n_adults))


dd_seq_info <- read_excel("~/Box Sync/midges/fastq-stats.xlsx",
                          sheet = "trimmed") %>%
    rename(sample_name = `Sample Name`, perc_dups = `% Dups`,
           perc_gc = `% GC`, mill_seqs = `M Seqs`,
           read_length = `Read Length`) %>%
    mutate(biotech_id = str_remove_all(sample_name,
                                       "trimmed_|_L002_R1_001|_L002_R2_001"),
           read_length = str_remove_all(read_length, " bp") %>%
               as.integer(),
           read = str_split(sample_name, "_R") %>%
               map_chr(~ .x[[2]]) %>%
               str_remove_all("_001") %>% as.integer()) %>%
    select(biotech_id, read, everything(), -coverage)


full_dd <- left_join(dd_seq_info, dd_other_info, by = "biotech_id")

write_csv(full_dd, "~/Box Sync/midges/full-DNA-info.csv")


# dd_other_info %>%
#     filter(biotech_id %in% proc_rn) %>%
#     select(biotech_id, n_adults)

