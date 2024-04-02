# `03-phylo`

Construct time-calibrated phylogeny for nine chironomids and five dipteran 
outgroups.

Folder contents:

```
.
├── README.md
├── codeml.sh
├── mafft.sh
├── mcmctree-modelfinder.sh
├── mcmctree.sh
├── mcmctree.sub
├── mega.sh
├── modeltest-combine.sh
├── modeltest-sep.sh
├── modeltest-sep.sub
├── orthodb.sh
├── orthodb.sub
├── phylo.sub
├── prequal.sh
├── raxml-boot.sh
├── raxml-boot.sub
├── raxml-mt-combine.sh
├── raxml-mt-combine.sub
├── raxml-mt-run.sh
└── raxml-mt-run.sub
```


The bash scripts in this folder are as follows:

- `codeml.sh`: Use codeml in PAML to estimate amino acid substitution rate from
  concatenated alignments.
- `mafft.sh`: Sequence alignment using MAFFT, then concatenation of alignments.
- `mcmctree-modelfinder.sh`: Use ModelFinder to merge partitions and do model
  selection for use in MCMCTree (not RAxML-NG).
- `mcmctree.sh`: Ultrametric timetree based on ML tree and fossil records using
  MCMCTree in PAML.
- `mega.sh`: Create time tree for input to ddBD using RelTime-ML in MEGA.
- `modeltest-combine.sh`: Combine output from the 10 runs of `modeltest-sep.sh`.
- `modeltest-sep.sh`: Run ModelTest-NG to do model selection on dataset that
  was partitioned by gene and not merged. This task was split into 10 jobs,
  and each run of this script did 1/10 of the work.
- `orthodb.sh`: Use BUSCO to find single-copy genes from OrthoDB (diptera_odb10) 
  within each assembly, then find the genes shared amongst all species.
- `prequal.sh`: Pre-alignment quality filter with prequal.
- `raxml-boot.sh`: Run up to 1,000 bootstrap replicates for ML tree using 
  RAxML-NG, then evaluate ML tree branch support.
- `raxml-mt-combine.sh`: Combine output from the 10 runs of `raxml-mt-run.sh`.
- `raxml-mt-run.sh`: Create ML tree using RAxML-NG. This task was split into
  10 jobs, and each run of this script did 1/10 of the work.
