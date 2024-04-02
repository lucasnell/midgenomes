
# `_data`

This folder contains small datasets associated with this paper.


Folder contents:

```
.
├── README.md
├── family-regressions.rds**
├── genome-stats.csv
├── gsize-corrs.rds**
├── gsizes-list.rds**
├── species-names-families.csv
├── hyphy
│   ├── focal-hog-go.csv
│   └── hyphy-hog-genes.csv
├── phylo
│   ├── chir_mega_relTimes.nwk
│   ├── chir_ml.tree
│   ├── time-tree-noCmarin.nwk
│   └── time-tree.nwk
└── spatial
    ├── Myvatn_WSGUTM28.geojson
    ├── Tgraci-sample-locations.csv
    └── gadm41_ISL_0.geojson
```

`**` indicates files that won't be here when you first download this repo,
but will be generated when running the scripts inside `midgenomes/_R`.


## Content descriptions

#### `family-regressions.rds`:
RDS file for a list of `phylolm` objects for the 
analysis of genome features in Chironimidae and Ceratopogonidae versus 
other dipterans.
Created in `_R/03-gsize-traits/02-gsize-traits-analysis.R`.

#### `genome-stats.csv`: 
Table of genome features for all species in the phylogeny.
Created in `_R/03-gsize-traits/01-make-genome-stats.R`.<br>
Column descriptions:

- `family`: family species belongs to
- `species`: full species name
- `spp_abbrev`: abbreviated species name
- `annot_source`: source of genome annotation
- `gsize`: size of genome in bp
- `n_genes`: number of protein-coding genes
- `sum_interg_len`: total intergenic sequence in the genome (bp)
- `mean_intron_len`: mean intron length across all genes that match to HOGs 
  with representative genes from all species and with < 4 paralogs per species
- `sum_intron_len`: total intron length for genes that meet conditions listed
  for `mean_intron_len`
- `mean_log_intron_len`: mean of log10-transformed intron lengths for genes
  that meet conditions listed for `mean_intron_len`
- `mean_n_introns`: mean number of introns per gene for genes
  that meet conditions listed for `mean_intron_len`
- `mean_log_n_introns`: mean of log10-transformed number of introns per gene
  for genes that meet conditions listed for `mean_intron_len`
- `DNA`: total length of genome comprised of repeat-element type `"DNA"`
- `LINE`: total length of genome comprised of repeat-element type `"LINE"`
- `LTR`: total length of genome comprised of repeat-element type `"LTR"`
- `Low_complexity`: total length of genome comprised of repeat-element type `"Low_complexity"`
- `RC`: total length of genome comprised of repeat-element type `"RC"`
- `SINE`: total length of genome comprised of repeat-element type `"SINE"`
- `Satellite`: total length of genome comprised of repeat-element type `"Satellite"`
- `Simple_repeat`: total length of genome comprised of repeat-element type `"Simple_repeat"`
- `Small_RNA`: total length of genome comprised of repeat-element type `"Small_RNA"`
- `Unclassified`: total length of genome comprised of repeat-element type `"Unclassified"`


#### `gsize-corrs.rds`:
RDS file for a list of `cor_phylo` objects for the analysis
of whether genome traits correlate with genome size when accounting
for phylogenetic autocorrelation.
Created in `_R/03-gsize-traits/02-gsize-traits-analysis.R`.


#### `gsizes-list.rds`:
RDS file for a named list of genome sizes with abbreviated 
species names as the names of the list elements.
Created in `_R/01-samps-annots/02-repeats-divergences.R`.

#### `species-names-families.csv`: 
Table of species info for all species in the phylogeny.<br>
Column descriptions:

- `family`: family species belongs to
- `species`: full species name
- `spp_abbrev`: abbreviated species name



### `hyphy` folder:

#### `focal-hog-go.csv`:
Table of focal GO terms related to stressful environments, descriptions, and
HOGs that match.
Created in `_R/05-hyphy/02-hyphy-genes.R`<br>
Column descriptions:

- `go`: parent GO term that relates to stressful environments
- `term`: description of GO term in `go`
- `offspring`: all offspring GO terms of the GO terms in `go`; 
  multiples separated by `";"`
- `hogs`: all OrthoFinder 1-to-1 HOGs that match to the GO term in `go` or 
  to any of its offspring in `offspring`; multiples separated by `";"`

#### `hyphy-hog-genes.csv`:
All genes across the phylogeny that are part of 1-to-1 HOGs that match to our
list of GO terms related to stressful environments.
Created in `_R/05-hyphy/02-hyphy-genes.R`<br>
Column descriptions:

- `species`: abbreviated species name
- `gene`: gene name
- `hog`: HOG name



### `phylo` folder:

#### `chir_mega_relTimes.nwk`:
Tree output from RelTime-ML in MEGA-CC in newick format.

#### `chir_ml.tree`:
Maximum likelihood tree from RAxML-NG in newick format.

#### `time-tree.nwk`:
Time tree from MCMCTree made to be ultrametric and in newick format
Created in `_bash/03-phylo/mcmctree.sh`.

#### `time-tree-noCmarin.nwk`:
Ultrametric time-tree without species *Clunio marinus* because of its
poor matching in OrthoFinder.


### `spatial` folder:

#### `Myvatn_WSGUTM28.geojson`:
Outline of Lake Myvatn, Iceland.

#### `Tgraci-sample-locations.csv`:
Locations of where *Tanytarsus gracilentus* were sampled.

#### `gadm41_ISL_0.geojson`:
Outline of Iceland.
