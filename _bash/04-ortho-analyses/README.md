# `04-ortho-analyses`


This folder contains scripts for identifying orthogroups using OrthoFinder,
and for analyses resulting from OrthoFinder output: 
gene family evolution (using CAFE) and positive selection (using HyPhy).


Folder contents:

```
.
├── README.md
├── cafe.sh
├── hyphy-align.sh
├── hyphy-busted.sh
├── hyphy-busted.sub
├── hyphy-relax.sh
├── hyphy-relax.sub
└── orthofinder.sh
```


The bash scripts in this folder are as follows:

- `cafe.sh`: Gene family evolution analysis using CAFE.
- `hyphy-align.sh`: Codon-aware alignments used in HyPhy selection tests.
- `hyphy-busted.sh`: Selection tests using HyPhy BUSTED method.
- `hyphy-relax.sh`: Selection tests using HyPhy RELAX method.
- `orthofinder.sh`: Find orthogroups using OrthoFinder.
