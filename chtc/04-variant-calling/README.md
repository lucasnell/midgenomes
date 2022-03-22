
# bwa.sh

- Merge overlapping paired reads using `bbmerge.sh`
- Separately align merged and unmerged reads using `bwa mem`
- Combine alignments using `picard MergeSamFiles`


# mpileup.sh

- Mark and remove duplicates using `picard MarkDuplicates`
- Re-align sequences in the proximity of indels with `IndelRealigner` and
  `RealignerTargetCreator` in `GATK`
- Produce pileup text files from BAMs using `samtools mpileup`


# snape.sh

- Estimate allele frequencies using `snape-pooled`
- Produce gSYNC files from the snape output
- Filter gSYNC files based on coverage, indels, and repetitive elements


From Kapun et al. (2021):
"Sites were filtered from gSYNC files if they had: 
1) minimum read depth <10; 
2) maximum read depth >the 95% coverage percentile of a given chromosomal 
   arm and sample; 
3) located within repetitive elements as defined by RepeatMasker;
4) within 5-bp distance up-and downstream of indel polymorphisms identified 
   by the GATK IndelRealigner."








# OLD DESCRIPTION:


- trim and filter reads using
    [`fastp`](https://bioconda.github.io/recipes/fastp/README.html)
    - quality trimming not necessary:
        https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-S8-S8
- map reads to genome using 
    [`bwa mem`](https://bioconda.github.io/recipes/bwa/README.html)
- duplicates marked using `samtools markdup`
- filter reads (Q >= 20, no secondary alignment) using `samtools view`
- add group information using `AddOrReplaceReadGroups` in `picard`
- re-align sequences in the proximity of indels with `IndelRealigner` and
    `RealignerTargetCreator` in
    [`GATK`](https://bioconda.github.io/recipes/gatk/README.html)
- use `samtools mpileup` to generate pileup file
- estimate allele frequencies using SNAPE-pool
    - [This paper](https://doi.org/10.1111/1755-0998.13343) 
      indicates it works really well
    - [This paper](https://doi.org/10.1093/molbev/msab259)
      created a bunch of tools to deal with output
    - See...
        - https://github.com/EmanueleRaineri/snape-pooled
        - https://github.com/DEST-bio/DEST_freeze1/blob/main/mappingPipeline/scripts/Mpileup2Snape.sh
        - https://github.com/DEST-bio/DEST_freeze1/blob/main/mappingPipeline/scripts/SNAPE2SYNC.py
        - https://github.com/DEST-bio/DEST_freeze1/blob/main/mappingPipeline/scripts/MaskSYNC_snape.py
        - https://github.com/DEST-bio/DEST_freeze1/blob/main/mappingPipeline/scripts/Snape_to_VCF.py







- damage profiling using
    [`DamageProfiler`](https://bioconda.github.io/recipes/damageprofiler/README.html)
- for population genomics analyses:
    [`poolfstat`](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13557)


## Variant calling from pileup:

[MAPGD](https://github.com/LynchLab/MAPGD)
- performs well via https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13343
- no indels
- I'm going with this one

[snape](https://github.com/EmanueleRaineri/snape-pooled)
- performs well via https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13343
- what should priors be??
- no indels

[PoPoolation2](https://sourceforge.net/p/popoolation2/wiki/Manual/)
- documentation doesn't look good
- any papers on its performance??
- it only hard-filters from mpileup, not a great method

[CRISP](https://github.com/vibansal/crisp/)
- used a lot, but is it actually good?
- indels!



see this for more info: `https://evodify.com/gatk-in-non-model-organism/`
