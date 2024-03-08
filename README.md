# midgenomes

[![DOI](https://zenodo.org/badge/499583048.svg)](https://zenodo.org/badge/latestdoi/499583048)

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
└── _bash
```

The following files/folders should be present:

- `LICENSE.md`: file containing the CC-BY license for this repository
- `README.md`: this file
- `midgenomes.Rproj`: file saving this RStudio Project's preferences
- `_data`: Small datasets associated with this paper.
- `_R`: R scripts mostly used to create figures and analyze data output 
  from the programs used inside the `_bash` folder's scripts, but some create
  necessary output used in `_bash`.
- `_bash`: Bash scripts used to run much of the bioinformatic
  software for this project. Most were run on UW--Madison's
  [Center for High Throughput Computing (CHTC)](https://chtc.cs.wisc.edu/),
  using the [HTCondor Software Suite](http://htcondor.org/).


The folders `_bash`, `_data`, and `_R` have separate
`README.md` files that have more information on each.


# Required dataset

The smaller datasets necessary to run the code herein is located inside
the `_data` folder.
See `_R/README.md` for how to download and use the larger datasets.



# Replicating the R environment

My session info:

```
 setting  value
 version  R version 4.3.2 (2023-10-31)
 os       macOS Sonoma 14.3.1
 system   aarch64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/Los_Angeles
 date     2024-03-05
 rstudio  2023.12.1+402 Ocean Storm (desktop)
 pandoc   3.1.9 @ /opt/homebrew/bin/pandoc
```

Alphabetized list of all loaded packages and their version numbers:

```
AnnotationDbi@1.62.2
ape@5.7.1
clusterProfiler@4.8.2
coda@0.19.4
deeptime@1.0.1
future@1.33.0
future.apply@1.11.0
ggtext@0.1.2
ggtree@3.8.2
GO.db@3.17.0
grid@4.3.2
gt@0.10.0
jsonlite@1.8.7
knitr@1.45
org.Dm.eg.db@3.17.0
parallel@4.3.2
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
          "coda@0.19.4", "deeptime@1.0.1",
          "future@1.33.0", "future.apply@1.11.0", "ggtext@0.1.2", 
          "ggtree@3.8.2", "GO.db@3.17.0", "grid@4.3.2", "gt@0.10.0", 
          "jsonlite@1.8.7", "knitr@1.45", "org.Dm.eg.db@3.17.0", 
          "parallel@4.3.2", "patchwork@1.1.3", "phylolm@2.6.2", 
          "phyr@1.1.2", "rrvgo@1.12.2", "tidyverse@2.0.0", "treeio@1.24.3", 
          "treemap@2.4.4", "viridisLite@0.4.2")
install.packages(pkgs)
```

Note that these only install the proper versions of the packages I manually 
installed, so dependencies might vary from what I used.



# Replicating the bash environment

Most scripts in `_bash` should be run using the Docker image found at
<https://hub.docker.com/r/lucasnell/midgenomes>, so to replicate this
environment you can pull the latest version:

```bash
docker pull lucasnell/midgenomes:v1.0.12
```


The main exception to this is `_bash/01-contigs/backmap.sh`,
which uses the image found at <https://hub.docker.com/r/lucasnell/tany_backmap>,
which you can pull using the following:

```bash
docker pull lucasnell/tany_backmap:v0.0.1
```
