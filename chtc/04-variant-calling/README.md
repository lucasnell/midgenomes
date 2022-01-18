

- trim and filter reads using
    [`fastp`](https://bioconda.github.io/recipes/fastp/README.html)
    - quality trimming not necessary:
        https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-S8-S8
- map reads to genome using 
    [`bwa mem`](https://bioconda.github.io/recipes/bwa/README.html)
- filter reads (Q >= 20, no secondary alignment) using `samtools`
- duplicates marked using `MarkDuplicates` from
    [`picard`](https://bioconda.github.io/recipes/picard/README.html)
- add group information using `AddOrReplaceReadGroups` in `picard`
- re-align sequences in the proximity of indels with `IndelRealigner` and
    `RealignerTargetCreator` in
    [`GATK`](https://bioconda.github.io/recipes/gatk/README.html)
- use `samtools mpileup` to generate pileup file
- estimate allele frequencies using `mapgd pool` 
    from [MAPGD](https://github.com/LynchLab/MAPGD)
- damage profiling using
    [`DamageProfiler`](https://bioconda.github.io/recipes/damageprofiler/README.html)
- for population genomics analyses:
    [`poolfstat`](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13557)


## Variant calling from pileup:

[MAPGD](https://github.com/LynchLab/MAPGD)
- performs well via https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13343
- no indels
- I'm ging with this one

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
