

- `orthofinder-all-gos.R`: Extract GO terms for all OrthoFinder HOGs. 
  Output is used for CAFE.
  
- `orthofinder-1to1-gos.R`: Extract GO terms for 1-to-1 OrthoFinder HOGs. 
  Output is used for HyPhy.

- `phylogeny-plots.R`: Plot maximum likelihood and time-calibrated 
  phylogenetic trees

- `ddBD.R`: Use ddBD method to inform priors for the speciation birth--death 
  process in MCMCTree

- `01-sample-maps.R`: Maps of *Tanytarsus gracilentus* sample locations.

- `cafe-gos.R`: Filters for HOGs that significantly expanded in Chironomidae,
  does overrepresentation test on GO terms in expanded HOGs, and
  produces treemap plot of overrepresented GO terms.

- `hyphy-genes.R`: Produces files necessary to use in HyPhy only 1-to-1 HOGs 
  that are linked to GO terms related to stressful environments.

- `hyphy-summarize.R`: Summarize significant results for HyPhy tests BUSTED 
  and RELAX, including describing HOGs that were significant for both tests.

- `mantis-counts.R`: Summarize output from mantis

- `cds-busco.R`: Make plot of BUSCO scores for all species' transcriptomes.

- `repeats-divergences.R`: Make repeat divergence landscape plots.

- `01-gsize-traits-analysis.R`: Analyze the following:
  (1) Do genome traits differ for families Chironomidae and Ceratopogonidae
  versus all other dipterans in the phylogeny?
  (2) Do genome traits significantly correlate with genome size when accounting
  for phylogenetic autocorrelation?

- `02-traits-on-tree.R`: Create plots of traits along phylogeny and confidence 
  intervals for traits for either (a) families Chironomidae and Ceratopogonidae
  or (b) all others.

- `03-gsize-vs-traits.R`: Plot genome size versus other traits that may relate to it

- `02-make-genome-stats.R`: Make `midgenomes/_data/genome-stats.csv` file 
  based on assemblies, annotations, repeats, and OrthoFinder output.

- `00-preamble.R`: shared code used in all R scripts.

