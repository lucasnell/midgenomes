
# `_docker`

This folder contains the files used to create the Docker container in which 
most scripts inside `chtc` were run.


Folder contents:

```
.
├── Dockerfile
├── README.md
├── RepBaseRepeatMaskerEdition-20181026.tar.gz**
├── env-annotate.yml
├── env-assembly.yml
├── env-main.yml
├── env-phylo.yml
├── env-repeat.yml
├── gm_key_64.gz**
├── gmes_linux_64_4.tar.gz**
├── helpers.sh
├── longest-isoforms.py
├── pretty-csv.py
└── summ-scaffs.py
```

`**` indicates files that I cannot share here because of the licenses for these
programs


File descriptions:

- `Dockerfile`: File used to construct the Docker container.
- `RepBaseRepeatMaskerEdition-20181026.tar.gz`: Repbase-derived RepeatMasker 
  libraries found at <https://www.girinst.org/server/RepBase/>
- `env-annotate.yml`: YAML file used to construct the `conda` environment for
  genome annotations.
- `env-assembly.yml`: YAML file used to construct the `conda` environment for
  assembling genomes.
- `env-main.yml`: YAML file used to construct the primary `conda` environment
  with some generally useful packages.
- `env-phylo.yml`: YAML file used to construct the `conda` environment for
  constructing phylogenies, plus running CAFE.
- `env-repeat.yml`: YAML file used to construct the `conda` environment for
  constructing repeat libraries and for masking and summarizing repeats
  in genomes.
- `gm_key_64.gz`: License file for GeneMark-ES/ET/EP+ required for BRAKER2.
  Download from here: <http://exon.gatech.edu/GeneMark/license_download.cgi>
- `gmes_linux_64_4.tar.gz`: Compiled software package for GeneMark-ES/ET/EP+
  required to use BRAKER2.
  Download from here: <http://exon.gatech.edu/GeneMark/license_download.cgi>
- `helpers.sh`: A bash script containing helper functions.
- `longest-isoforms.py`: Python script that filters for the longest isoform
  for each gene. See the description inside this file for the assumptions of
  this script.
- `pretty-csv.py`: Combine output from `summ-scaffs.py` and BUSCO into a 
  single CSV file to make it easier to input to a summary table.
- `summ-scaffs.py`: Summarize contigs and scaffolds.


