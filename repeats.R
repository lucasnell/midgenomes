library(tidyverse)

repeats <- read_table("~/_data/annotation/tany_repeats/masker/tany_contigs.fasta.out",
                      skip = 3, col_types = cols(),
                      col_names = c("score", "div", "del", "ins", "sequence",
                                    "q_begin", "q_end", "q_left",
                                    "sign", "repeat_name", "repeat_class",
                                    "r_begin", "r_end", "r_left",
                                    "id", "aster"))

repeats %>%
    filter(repeat_class == "Unknown") %>%
    select(sequence, starts_with("q_"), "repeat_name") %>%
    .[["repeat_name"]] %>% table() %>% sort()


repeats %>%
    filter(repeat_name == "rnd-5_family-1225") %>%
    distinct(sequence)



