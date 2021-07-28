
library(tidyverse)


# Looking at N50 for reduced assembly
fasta <- read_lines("~/_data/redundans/contigs.reduced.fa")

splits <- tibble(start = which(str_detect(fasta, "^>")),
                 end = c(tail(start, -1) - 1, length(fasta))) %>%
    mutate(start = start + 1)

seqs <- pmap_chr(splits,
                 function(start, end) str_c(fasta[start:end], collapse = ""))

lens <- sort(nchar(seqs), decreasing = TRUE)
sum(lens)

i <- which(cumsum(lens) >= (sum(lens) / 2))[1]

# N50:
lens[i]
