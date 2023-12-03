# `03-phylo`

Construct time-calibrated phylogeny for nine chironomids and five dipteran 
outgroups.

Folder contents:

```
.
├── README.md
├── codeml.sh
├── mafft.sh
├── mcmctree.sh
├── mega.sh
├── orthodb.sh
├── orthodb.sub
├── phylo.sub
├── prequal.sh
├── raxml-boot.sh
├── raxml-boot.sub
├── raxml-supp.sh
├── raxml.sh
└── raxml.sub
```


The bash scripts in this folder are as follows:

- `codeml.sh`: Use codeml in PAML to estimate amino acid substitution rate from
  concatenated alignments.
- `mafft.sh`: Sequence alignment using MAFFT, then concatenation of alignments.
- `mcmctree.sh`: Ultrametric timetree based on ML tree and fossil records using
  MCMCTree in PAML.
- `mega.sh`: Create time tree for input to ddBD using RelTime-ML in MEGA.
- `orthodb.sh`: Use BUSCO to find single-copy genes from OrthoDB (diptera_odb10) 
  within each assembly, then find the genes shared amongst all species.
- `prequal.sh`: Pre-alignment quality filter with prequal.
- `raxml-boot.sh`: Run 100 bootstrap replicates for ML tree using RAxML-NG.
- `raxml-supp.sh`: Use bootstrapped trees to evaluate ML tree branch support 
  using RAxML-NG.
- `raxml.sh`: Create ML tree using RAxML-NG.

