# midgenomes

<!-- [![DOI](https://zenodo.org/badge/632674820.svg)](https://zenodo.org/badge/latestdoi/632674820) -->

Code and small datasets related to the shared features underlying compact 
genomes and extreme habitat use in chironomid midges.

This repository was created and is maintained by
[Lucas A. Nell](https://github.com/lucasnell).

# Organization

Folder contents:

```
.
├── LICENSE.md
├── README.md
├── midgenomes.Rproj
├── _data
├── _R
└── chtc
```

The following files/folders should be present:

- `LICENSE.md`: file containing the CC-BY license for this repository
- `README.md`: this file
- `midgenomes.Rproj`: file saving this RStudio Project's preferences
- `_data`: Small datasets associated with this paper.
- `_R`: R scripts mostly used to create figures and analyze data output 
  from the programs used inside the `chtc` folder's scripts, but some create
  necessary output used in `chtc`.

- `chtc`: 


The folders `_R`, `_data`, and `chtc` have separate
`README.md` files that have more information on each.



# Replicating the R environment

My session info:

```
 setting  value
 version  R version 4.3.1 (2023-06-16)
 os       macOS Sonoma 14.1.1
 system   aarch64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/Los_Angeles
 date     2023-11-30
 rstudio  2023.06.0+421 Mountain Hydrangea (desktop)
 pandoc   3.1.9 @ /opt/homebrew/bin/ (via rmarkdown)
```

Alphabetized list of all loaded packages and their version numbers:

```
AnnotationDbi@1.62.2
ape@5.7.1
clusterProfiler@4.8.2
future@1.33.0
future.apply@1.11.0
ggtext@0.1.2
ggtree@3.8.2
GO.db@3.17.0
grid@4.3.1
gt@0.10.0
jsonlite@1.8.7
knitr@1.45
org.Dm.eg.db@3.17.0
parallel@4.3.1
patchwork@1.1.3
phylolm@2.6.2
phyr@1.1.2
rrvgo@1.12.2
tidyverse@2.0.0
treeio@1.24.3
treemap@2.4.4
viridisLite@0.4.2
```

To install all these packages in the versions I used:

```r
pkgs <- c("AnnotationDbi@1.62.2", "ape@5.7.1", "clusterProfiler@4.8.2",
          "future@1.33.0", "future.apply@1.11.0", "ggtext@0.1.2",
          "ggtree@3.8.2", "GO.db@3.17.0", "grid@4.3.1", "gt@0.10.0",
          "jsonlite@1.8.7", "knitr@1.45", "org.Dm.eg.db@3.17.0",
          "parallel@4.3.1", "patchwork@1.1.3", "phylolm@2.6.2", "phyr@1.1.2",
          "rrvgo@1.12.2", "tidyverse@2.0.0", "treeio@1.24.3", "treemap@2.4.4",
          "viridisLite@0.4.2")
install.packages(pkgs)
```

Note that these only install the proper versions of the packages I manually 
installed, so dependencies might vary from what I used.


