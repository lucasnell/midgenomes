# `01-contigs`

Create genome assembly contigs for *Tanytarsus gracilentus* and 
*Parochlus steinii*.

Folder contents:

```
.
├── README.md
├── backmap.sh
├── flye.sh
├── flye.sub
├── necat.sh
├── necat.sub
├── nextdenovo-Pstein.sh
├── nextdenovo-Pstein.sub
├── nextdenovo.sh
├── nextdenovo.sub
├── nextpolish.sh
├── nextpolish.sub
├── ont-polish.sh
├── ont-polish.sub
├── purgedups.sh
├── quickmerge.sh
├── quickmerge.sub
├── simplify-final.sh
├── smartdenovo.sh
└── smartdenovo.sub
```


The bash scripts in this folder are as follows:

- `backmap.sh`: Use backmap.sh to estimate size of *T. gracilentus* genome by 
  back-mapping ONT reads to assembly.
- `flye.sh`: Assemble *T. gracilentus* genome using Flye.
- `necat.sh`: Assemble *T. gracilentus* genome using NECAT.
- `nextdenovo-Pstein.sh`: Assemble *P. steinii* genome using NextDenovo.
- `nextdenovo.sh`: Assemble *T. gracilentus* genome using NextDenovo.
- `nextpolish.sh`: Polish assembly using Illumina reads with NextPolish.
- `ont-polish.sh`: Polish assembly using ONT reads with Racon and medaka.
- `purgedups.sh`: Use purge_dups to remove duplicate regions from the
   NECAT and SMARTdenovo *T. gracilentus* assemblies
- `quickmerge.sh`: Combine assemblies using quickmerge.
- `simplify-final.sh`: Simplify the final reference assembly, specifically the 
  contig names and the file name.
- `smartdenovo.sh`: Assemble *T. gracilentus* genome using SMARTdenovo.
