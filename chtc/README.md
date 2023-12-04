
# `chtc`

This folder contains bash scripts used to run much of the bioinformatic
software for this project.


Folders are sorted by the order they're presented in the text.
All folders except `_docker_backmap` contain `README.md` files with 
descriptions of the files therein.


Folder contents:

```
.
├── README.md
├── 00-prep-reads
├── 01-contigs
├── 02-annotation
├── 03-phylo
├── 04-ortho-analyses
├── _docker
└── _docker_backmap
```


Folder descriptions:

- `00-prep-reads`: Prepare (and/or download) raw reads.
- `01-contigs`: Create genome assembly contigs for *Tanytarsus gracilentus* and 
  *Parochlus steinii*.
- `02-annotation`: Create genome annotations for *Tanytarsus gracilentus*, 
  *Parochlus steinii*, *Culicoides sonorensis*, and (only functional 
  annotations) *Chironomus riparius*.
- `03-phylo`: Construct time-calibrated phylogeny for nine chironomids and 
  five dipteran outgroups.
- `04-ortho-analyses`: This folder contains scripts for identifying orthogroups
  using OrthoFinder, and for analyses resulting from OrthoFinder output: 
  gene family evolution (using CAFE) and positive selection (using HyPhy).
- `_docker`: Files used to create the Docker container in which most scripts
  inside `chtc` were run.
- `_docker_backmap`: Dockerfile for using backmap on *Tanytarsus gracilentus*
  reference assembly. This is in a separate docker container bc the qualimap 
  conda install hangs for some reason.

