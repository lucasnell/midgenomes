# `00-prep-reads`

Prepare (and/or download) raw reads.

Folder contents:

```
.
├── README.md
├── download-Cson-RNA.sh
├── fastp.sh
├── fastp.sub
├── guppy.sh
├── guppy.sub
└── porechop.sh
```


The bash scripts in this folder are as follows:

- `download-Cson-RNA.sh`: download RNAseq files for *Culicoides sonorensis*
  from SRA
- `fastp.sh`: Use fastp to trim paired-end, DNAseq or RNAseq Illumina reads.
- `guppy.sh`: Basecalling for Nanopore reads (used only for
  *Tanytarsus gracilentus*)
- `porechop.sh`: Code to attempt to use PoreChop to remove adapters from 
  *Parochlus steinii* ONT reads, but it turns out they're already removed.
