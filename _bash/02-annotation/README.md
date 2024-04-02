# `02-annotation`

Create genome annotations for *Tanytarsus gracilentus*, *Parochlus steinii*,
*Culicoides sonorensis*, and (only functional annotations) *Chironomus riparius*.

Folder contents:

```
.
├── README.md
├── braker-prot.sh
├── braker-prot.sub
├── braker-rna.sh
├── braker-rna.sub
├── cds-busco.sh
├── mantis-download.sh
├── mantis-download.sub
├── mantis.sh
├── mantis.sub
├── repeats-Aaegyp-mask.sh
├── repeats-Aaegyp-mask.sub
├── repeats-Aaegyp-model.sh
├── repeats-Aaegyp-model.sub
├── repeats-divergences.sh
├── repeats.sh
├── repeats.sub
├── sendsketch.sh
├── tsebra.sh
└── tsebra.sub
```


The bash scripts in this folder are as follows:

- `braker-prot.sh`: Use BRAKER2 for ".. prediction of protein coding gene 
  structures" using proteins from a database (pipeline C)
- `braker-rna.sh`: Use BRAKER2 for ".. prediction of protein coding gene 
  structures" using RNAseq alignments (pipeline B)
- `cds-busco.sh`: Run BUSCO on CDS from our assemblies and others.
- `mantis-download.sh`: Initial database downloads for functional annotations using mantis.
- `mantis.sh`: Functional annotations using mantis.
- `repeats-Aaegyp-mask.sh`: Use RepeatMasker to count repeat elements and to 
  softmask the *Aedes aegypti* assembly. This is a separate script because
  using the combined script (`repeats.sh`) on *A. aegypti* exceeded CHTC's
  time limit.
- `repeats-Aaegyp-model.sh`: Use RepeatModeler to create a repeat library for 
  the *Aedes aegypti* assembly. This is a separate script because
  using the combined script (`repeats.sh`) on *A. aegypti* exceeded CHTC's
  time limit.
- `repeats-divergences.sh`: Calculate repeat divergences.
- `repeats.sh`: Use RepeatModeler to create a repeat library for an assembly,
  use TEclass to classify many of the unknowns from RepeatModeler, and
  use RepeatMasker to count repeat elements and to softmask the assembly.
- `sendsketch.sh`: Use BBmap's sendsketch.sh to look for contamination in the
  *Tanytarsus gracilentus* and *Parochlus steinii* assemblies.
- `tsebra.sh`: Combine 2 BRAKER annotations (RNAseq + OrthoDB proteins) 
  using TSEBRA.
