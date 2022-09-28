#!/bin/bash


#'
#' IMPORTANT: uses the docker container from BlobToolKit.
#'


export ASSEMBLY_NAME=tany
export BLOB_DIR=${ASSEMBLY_NAME}_blob
export ONT_READS=basecalls_guppy-5.0.11



mkdir -p ${BLOB_DIR}
cd ${BLOB_DIR}
mkdir busco_lineages
mkdir databases
mkdir tmp


cp /staging/lnell/assemblies/tany_contigs.fasta.gz ${ASSEMBLY_NAME}.fasta.gz \
    && gunzip ${ASSEMBLY_NAME}.fasta.gz

cp /staging/lnell/ont/${ONT_READS}.fastq.gz ./


cat << EOF > ${ASSEMBLY_NAME}.yaml
assembly:
  record_type: contig
  scaffold-count: 45
  span: 91827299
  prefix: ${ASSEMBLY_NAME}

taxon:
  taxid: 288803
  name: Tanytarsus gracilentus

busco:
  lineages:
    - diptera_odb10
    - arthropoda_odb10
    - eukaryota_odb10
  lineage_dir: $(pwd)/busco_lineages

reads:
  single:
    -
      - ${ONT_READS}
      - OXFORD_NANOPORE
      - 23082704153
  paired: []

settings:
  blobtools2_path: /blobtoolkit/blobtools2
  taxonomy: $(pwd)/databases/ncbi_2022_06
  tmp: $(pwd)/tmp
  blast_chunk: 100000
  blast_max_chunks: 10
  blast_overlap: 500
  chunk: 1000000

similarity:
  defaults:
    evalue: 1e-25
    max_target_seqs: 10
    root: 1
    mask_ids: []
  databases:
    -
      local: $(pwd)/databases/ncbi_2022_06
      name: nt
      source: ncbi
      tool: blast
      type: nucl
    -
      local: $(pwd)/databases/uniprot_2022_06
      max_target_seqs: 1
      name: reference_proteomes
      source: uniprot
      tool: diamond
      type: prot
  taxrule: bestsumorder

keep_intermediates: true
EOF




