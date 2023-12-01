
# `_R`

This folder contains R scripts mostly used to create figures and analyze
data output from the programs used inside the `chtc` folder's scripts, but 
some create necessary output used in `chtc`.
Folders are sorted by the order they're presented in the text, and
scripts within are sorted by when they're used.


Folder contents:

```
.
├── 00-preamble.R
├── 01-samps-annots
│   ├── 01-sample-maps.R
│   ├── 02-repeats-divergences.R
│   ├── 03-cds-busco.R
│   └── 04-mantis-counts.R
├── 02-phylo
│   ├── 01-ddBD.R
│   └── 02-phylogeny-plots.R
├── 03-gsize-traits
│   ├── 01-make-genome-stats.R
│   ├── 02-gsize-traits-analysis.R
│   ├── 03-traits-on-tree.R
│   └── 04-gsize-vs-traits.R
├── 04-cafe
│   ├── 01-cafe-orthofinder-gos.R
│   └── 02-cafe-gos.R
├── 05-hyphy
│   ├── 01-hyphy-orthofinder-gos.R
│   ├── 02-hyphy-genes.R
│   └── 03-hyphy-summarize.R
└── README.md
```


## Content descriptions

- `00-preamble.R`: shared code used in all R scripts.


### `01-samps-annots` folder:

- `01-sample-maps.R`: Maps of *Tanytarsus gracilentus* sample locations.
- `02-repeats-divergences.R`: Make repeat divergence landscape plots.
- `03-cds-busco.R`: Make plot of BUSCO scores for all species' transcriptomes.
- `04-mantis-counts.R`: Summarize output from mantis



### `02-phylo` folder:

- `01-ddBD.R`: Use ddBD method to inform priors for the speciation birth--death 
  process in MCMCTree
- `02-phylogeny-plots.R`: Plot maximum likelihood and time-calibrated 
  phylogenetic trees


### `03-gsize-traits` folder:

- `01-make-genome-stats.R`: Make `midgenomes/_data/genome-stats.csv` file 
  based on assemblies, annotations, repeats, and OrthoFinder output.
- `02-gsize-traits-analysis.R`: Analyze the following:
  (1) Do genome traits differ for families Chironomidae and Ceratopogonidae
  versus all other dipterans in the phylogeny?
  (2) Do genome traits significantly correlate with genome size when accounting
  for phylogenetic autocorrelation?
- `03-traits-on-tree.R`: Create plots of traits along phylogeny and confidence 
  intervals for traits for either (a) families Chironomidae and Ceratopogonidae
  or (b) all others.
- `04-gsize-vs-traits.R`: Plot genome size versus other traits that may relate to it



### `04-cafe` folder:

- `01-cafe-orthofinder-gos.R`: Extract GO terms for all OrthoFinder HOGs. 
  Output is used for CAFE.
- `02-cafe-gos.R`: Filters for HOGs that significantly expanded in Chironomidae,
  does overrepresentation test on GO terms in expanded HOGs, and
  produces treemap plot of overrepresented GO terms.


### `05-hyphy` folder:

- `01-hyphy-orthofinder-gos.R`: Extract GO terms for 1-to-1 OrthoFinder HOGs. 
  Output is used for HyPhy.
- `02-hyphy-genes.R`: Produces files necessary to use in HyPhy only 1-to-1 HOGs 
  that are linked to GO terms related to stressful environments.
- `03-hyphy-summarize.R`: Summarize significant results for HyPhy tests BUSTED 
  and RELAX, including describing HOGs that were significant for both tests.




