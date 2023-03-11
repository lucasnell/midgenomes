

library(tidyverse)


#' "RepeatMasker type/subtype to Dfam Classification Lineage"
#' From `https://github.com/Dfam-consortium/RepeatModeler/blob/afa4e4fddea85078e3a3ea20e638b7b2b034883b/RepeatClassifier`
#' (accessed on 27 Jan 2023)

rd_table <- '"other" => "Other",
"artefact" => "Artifact",
"segmental" => "Segmental_Duplication",
"low_complexity" => "Low_Complexity",
"other/dna_virus" => "Accidental;Normally_Non-integrating_Virus",
"unknown" => "Interspersed_Repeat;Unknown",
"simple_repeat" => "Tandem_Repeat;Simple",
"satellite" => "Tandem_Repeat;Satellite",
"unknown/centromeric" => "Interspersed_Repeat;Unknown;Centromeric",
"rna" => "Interspersed_Repeat;Pseudogene;RNA",
"satellite/centromeric" => "Tandem_Repeat;Satellite;Centromeric",
"satellite/macro" => "Tandem_Repeat;Satellite;Macro",
"satellite/y-chromosome" => "Tandem_Repeat;Satellite;Y-chromosomal",
"satellite/acromeric" => "Tandem_Repeat;Satellite;Acromeric",
"satellite/w-chromosome" => "Tandem_Repeat;Satellite;W-chromosomal",
"satellite/subtelomeric" => "Tandem_Repeat;Satellite;Subtelomeric",
"dna" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase",
"rc" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Helicase",
"line" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE",
"scrna" => "Interspersed_Repeat;Pseudogene;RNA;scRNA",
"rrna" => "Interspersed_Repeat;Pseudogene;RNA;rRNA",
"trna" => "Interspersed_Repeat;Pseudogene;RNA;tRNA",
"snrna" => "Interspersed_Repeat;Pseudogene;RNA;snRNA",
"dna/crypton" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton",
"dna/p" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;P_Element",
"dna/kolobok" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Kolobok",
"dna/hat" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT",
"dna/ginger" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Ginger",
"dna/zisupton" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Zisupton",
"dna/zator" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Zator",
"dna/pif" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PIF-Harbinger",
"dna/merlin" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Merlin",
"dna/mule" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Mutator-like",
"dna/dada" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Dada",
"dna/novosib" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Novosib",
"dna/tcmar" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner",
"dna/is3eu" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;IS3EU",
"dna/maverick" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;DNA_Polymerase;Maverick",
"rc/helitron-2" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Helicase;Helitron-2",
"rc/helitron" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Helicase;Helitron-1",
"retroposon" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter",
"sine" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE",
"line/penelope" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements",
"ltr" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element",
"unknown/tate" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Tyrosine_Recombinase_Elements;Viper-group;TATE",
"dna/crypton-s" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-S",
"dna/crypton-r" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-R",
"dna/crypton-v" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-V",
"dna/crypton-f" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-F",
"dna/crypton-c" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-C",
"dna/crypton-i" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-I",
"dna/crypton-x" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-X",
"dna/crypton-a" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-A",
"dna/crypton-h" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Tyrosine_Recombinase;Crypton;Crypton-H",
"dna/p-fungi" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;P_Element;Fungi-specific_Branch",
"dna/kolobok-hydra" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Kolobok;Hydra-specific_Branch",
"dna/kolobok-h" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Kolobok;Kolobok-H",
"dna/kolobok-e" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Kolobok;Kolobok-E",
"dna/kolobok-t2" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Kolobok;T2",
"dna/hat-restless" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Restless",
"dna/hat-blackjack" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Blackjack",
"dna/hat-hobo" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hobo",
"dna/hat-tag1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Tag1",
"dna/hat-hatw" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hATw",
"dna/hat-tip100" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Tip100",
"dna/hat-hat6" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hAT6",
"dna/hat-hat19" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hAT19",
"dna/hat-hat1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hAT1",
"dna/hat-hatx" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hATx",
"dna/hat-pegasus" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Pegasus",
"dna/hat-hatm" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hATm",
"dna/hat-hat5" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;hAT5",
"dna/hat-ac" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Activator",
"dna/hat-charlie" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;hAT;Charlie",
"dna/cmc-transib" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;CACTA;Transib",
"dna/cmc" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;CACTA;CMC",
"dna/sola-1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Sola;Sola-1",
"dna/sola-3" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Sola;Sola-3",
"dna/sola-2" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Sola;Sola-2",
"dna/pif-isl2eu" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PIF-Harbinger;ISL2EU",
"dna/pif-harbinger" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PIF-Harbinger;Harbinger",
"dna/pif-spy" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PIF-Harbinger;Spy",
"dna/pif-harbs" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PIF-Harbinger;HarbS",
"dna/academ-h" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Academ;Academ-H",
"dna/academ-2" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Academ;Academ-2",
"dna/academ-1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Academ;Academ-1",
"dna/piggybac-a" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PiggyBac_Group;PiggyBac-A",
"dna/piggybac" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PiggyBac_Group;PiggyBac",
"dna/piggybac-x" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;PiggyBac_Group;PiggyBac-X",
"dna/mule-f" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Mutator-like;F",
"dna/mule-mudr" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Mutator-like;MuDR",
"dna/mule-nof" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Mutator-like;NOF",
"dna/casposons" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;DNA_Polymerase;Casposon",
"dna/tcmar-tc2" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc2-group;Tc2",
"dna/tcmar-tc4" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc4",
"dna/tcmar-mogwai" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Mogwai",
"dna/tcmar-stowaway" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Stowaway",
"dna/tcmar-mariner" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Mariner",
"dna/tcmar-cweed" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Cweed",
"dna/tcmar-sagan" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Sagan",
"dna/tcmar-gizmo" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Gizmo",
"dna/tcmar-ant1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Ant1",
"dna/tcmar-tc1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc1",
"dna/tcmar-isrm11" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;ISRm11",
"dna/tcmar-m44" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;m44",
"retroposon/r4-derived" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;R4-derived",
"retroposon/l1-dep" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;L1-dependent",
"retroposon/l1-derived" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;L1-derived",
"retroposon/l2-derived" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;L2-derived",
"retroposon/rte-derived" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;RTE-derived",
"retroposon/i-derived" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;I-derived",
"sine/7sl" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;7SL-RNA_Promoter",
"sine/5s" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter",
"ltr/trim" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;TRIM",
"ltr/pao" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Bel-Pao",
"ltr/copia" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Ty1-Copia",
"line/cre-odin" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-I;Odin",
"line/cre" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-I;CRE",
"line/cre-ambal" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-I;Ambal",
"line/genie" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Genie",
"dna/cmc-enspm" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;CACTA;CMC;EnSpm",
"dna/cmc-mirage" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;CACTA;CMC;Mirage",
"dna/tcmar-tigger" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc2-group;Tigger",
"dna/tcmar-fot1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc2-group;Fot1",
"dna/tcmar-pogo" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;Tc2-group;Pogo",
"retroposon/sva" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;Lacking_Small_RNA_pol_III_Promoter;L1-dependent;SVA",
"ltr/cassandra" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae",
"ltr/gypsy" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Gypsy",
"ltr/dirs" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Tyrosine_Recombinase_Elements;DIRS",
"ltr/viper" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Tyrosine_Recombinase_Elements;Viper-group;Viper",
"ltr/ngaro" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Tyrosine_Recombinase_Elements;Ngaro",
"line/cre-2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-I;CRE;CRE-2",
"line/cre-1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-I;CRE;CRE-1",
"line/proto1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;Proto-1",
"line/proto2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;Proto-2",
"line/rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;RTE-like",
"dna/cmc-chapaev" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;CACTA;CMC;Chapaev_group;Chapaev",
"dna/cmc-chapaev-3" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;CACTA;CMC;Chapaev_group;Chapaev-3",
"sine/trna-5s" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_and_5S_RNA;No_or_Unknown_Core;Unknown_LINE-dependent",
"sine/ceph" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;Unknown_Promoter;Ceph-core;RTE-end",
"sine/core-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;Unknown_Promoter;MIR-core;RTE-end",
"sine/core" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;Unknown_Promoter;MIR-core;Unknown_LINE-dependent",
"sine/trna-v-cr1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;V-core;CR1-end",
"sine/trna-v" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;V-core;Unknown_LINE-dependent",
"sine/trna-sauria-l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Sauria-core;L2-end",
"sine/trna-sauria-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Sauria-core;RTE-end",
"sine/trna-sauria" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Sauria-core;Unknown_LINE-dependent",
"sine/trna-ceph-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Ceph-core;RTE-end",
"sine/trna-ceph" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Ceph-core;Unknown_LINE-dependent",
"sine/trna-l1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;L1-dependent",
"sine/trna-i" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;I-end",
"sine/trna-l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;L2-end",
"sine/trna-rex" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;Rex-end",
"sine/trna-cr1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;CR1-end",
"sine/trna-jockey" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;Jockey-end",
"sine/r1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;R1-end",
"sine/trna-tad1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;Tad1_End",
"sine/trna-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;RTE-end",
"sine/trna-r2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;R2-end",
"sine/trna" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;Unknown_LINE-dependent",
"sine/rte-bovb" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No_or_Unknown_Core;BovB-end",
"sine/trna-deu-i" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;I-end",
"sine/trna-deu-l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;L2-end",
"sine/trna-deu-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;RTE-end",
"sine/trna-deu-cr1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;CR1-end",
"sine/trna-deu" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Deu-core;Unknown_LINE-dependent",
"sine/id" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No-core;L1-dependent",
"sine/trna-v-core-l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;V_and_MIR-core;L2-end",
"sine/trna-meta" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;Meta-core;Unknown_LINE-dependent",
"sine/mir" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;MIR-core;L2-end",
"sine/trna-mermaid" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;MIR-core;Mermaid",
"sine/trna-core-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;MIR-core;RTE-end",
"sine/trna-core" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;MIR-core;Unknown_LINE-dependent",
"sine/5s-sauria-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;Sauria-core;RTE-end",
"sine/5s-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;No_or_Unknown_Core;RTE-end",
"sine/5s-deu-l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;Deu-core;L2-end",
"sine/5s-deu" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;Deu-core;Unknown_LINE-dependent",
"sine/5s-core-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;5S-RNA_Promoter;MIR-core;RTE-end",
"sine/u" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;U-RNA_Promoter;No_or_Unknown_Core;Unknown_LINE-dependent",
"sine/b4" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_and_7SL_RNA;No-core;L1-dependent",
"sine/trna-7sl" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_and_7SL_RNA;No_or_Unknown_Core;Unknown_LINE-dependent",
"ltr/erv-foamy" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Spumaretrovirinae",
"ltr/erv" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae",
"ltr/caulimovirus" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Pararetroviridae;Caulimoviridae",
"line/r2-hero" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;R2-like;Hero",
"line/r2-nesl" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;R2-like;NeSL",
"line/r2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;R2-like;R2",
"line/l1-dre" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;L1-like;DRE",
"line/l1-zorro" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;L1-like;Zorro",
"line/dualen" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;R4-like;Dualen",
"line/dong-r4" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;R4-like;Dong-R4",
"line/deceiver" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;R4-like;Deceiver",
"line/rte-orte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;RTE-like;ORTE",
"line/tad1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;Tad1",
"sine/alu" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;7SL-RNA_Promoter;No-core;L1-dependent;Alu",
"sine/b2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;No-core;L1-dependent;B2",
"ltr/erv1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV1",
"line/l1-tx1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;L1-like;L1-group;Tx1",
"line/l1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-1;L1-like;L1-group;L1",
"line/rte-x" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;RTE-like;RTE-group;RTE-X",
"line/rte-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;RTE-like;RTE-group;RTE",
"line/rte-bovb" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;RTE-like;RTE-group;BovB",
"line/rex-babar" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;CR1-group;Rex-Babar",
"line/cr1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;CR1-group;CR1",
"ltr/erv-lenti" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;Lenti",
"ltr/erv4" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;ERV4",
"ltr/ervk" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;ERV2",
"ltr/ervl" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;ERV3",
"line/l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;CR1-group;L2-group;L2",
"line/cr1-zenon" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;CR1-group;CR1;Zenon",
"line/r1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;R1-group;R1-subgroup;R1",
"line/r1-loa" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;R1-group;R1-subgroup;LOA",
"line/i" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;R1-group;I-group;I",
"line/i-jockey" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE;Group-II;Group-2;R1-like;R1-group;I-group;Jockey",
"ltr/ervl-malr" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV2-group;ERV3;MaLR",
"dna/mule-ricksha" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Mutator-like;Ricksha",
"ltr/dirs-q" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Tyrosine_Recombinase_Elements;DIRS;Q",
"dna/maverick-mavirus" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;DNA_Polymerase;Maverick;Maverick-Mavirus",
"dna/tcmar-is885" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Tc1-Mariner;IS885",
"sine/trna-v-l2" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;V-core;L2-end",
"sine/trna-v-rte" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;tRNA_Promoter;V-core;RTE-end",
"dna/ginger-1" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Ginger;Ginger-1",
"dna/ginger-2" => "Interspersed_Repeat;Transposable_Element;Class_II_DNA_Transposition;Transposase;Ginger;Ginger-2",
"retroposon/sno" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;snoRNA",
"sine/erv1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Long_Terminal_Repeat_Element;Gypsy-ERV;Retroviridae;Orthoretrovirinae;ERV1;SINE-like",
"sine/u-l1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;SINE;U-RNA_Promoter;No_or_Unknown_Core;L1-dependent",
"retroposon/sno-l1" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;LINE-dependent_Retroposon;snoRNA;L1-derived",
"ple/chlamys" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Chlamys",
"ple/hydra" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Hydra",
"ple/naiad" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Naiad",
"ple/athena" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Athena",
"ple/poseidon" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Poseidon",
"ple/nematis" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Nematis",
"ple/neptune" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Neptune",
"ple/coprina" => "Interspersed_Repeat;Transposable_Element;Class_I_Retrotransposition;Retrotransposon;Penelope-like_Elements;Coprina",'


#' RepeatMaker classes:
rm_classes <- rd_table |>
    str_split("\n") |>
    getElement(1) |>
    str_remove_all('"|,') |>
    str_split(" => ") |>
    map_chr(~ .x[[1]])

#' Dfam classes:
df_classes <- rd_table |>
    str_split("\n") |>
    getElement(1) |>
    str_remove_all('"|,') |>
    str_split(" => ") |>
    map_chr(~ .x[2]) |>
    tolower()

teclass_out <- "/Users/lucasnell/_data/TEclass/Tanytarsus_gracilentus_repeats-lib.txt" |>
    read_lines() |>
    keep(~ grepl("^>", .x))
teclass_out2 <- "/Users/lucasnell/_data/annotation/Tgrac_repeats/tec_out/Tgrac_repeats_lib_unknown.fasta.lib" |>
    read_lines() |>
    keep(~ grepl("^>", .x))


# rm_guesses <-
teclass_out |>
    str_remove(".*#") |>
    str_remove(" \\[.*\\]") |>
    str_remove(" \\(.*\\)") |>
    str_remove("\\|TEclass result") |>
    str_remove("\\|forward|\\|reverse complemented") |>
    str_remove("\\|ORFs.*") |>
    str_split(": ") |>
    do.call(what = rbind) |>
    as.data.frame() |>
    set_names(c("rm", "tec")) |>
    as_tibble() |>
    filter(rm == "Unknown") |>
    getElement("tec") |> table() |>
    # getElement("rm") |> table() |>
    print(n = 99)
teclass_out2 |>
    str_remove(".*#") |>
    str_remove(" \\[.*\\]") |>
    str_remove(" \\(.*\\)") |>
    str_remove("\\|TEclass result") |>
    str_remove("\\|forward|\\|reverse complemented") |>
    str_remove("\\|ORFs.*") |>
    str_split(": ") |>
    do.call(what = rbind) |>
    as.data.frame() |>
    set_names(c("rm", "tec")) |>
    as_tibble() |>
    filter(rm == "Unknown") |>
    getElement("tec") |> table() |>
    # getElement("rm") |> table() |>
    print(n = 99)



z <- read_table("~/_data/_annotations/zzz-pre-2023/Tgrac_repeats/masker/Tgrac_contigs.fasta.out",
                col_names = c("SW", "p_div", "p_del", "p_ins", "query_seq",
                              "begin", "end", "left", "unk1", "matching_repeat",
                              "class", "begin2", "end2", "left2", "id", "star"),
                skip = 3)

# From *.tbl file:
# DNA transposons     392     131109
# LINEs:              97      58363
# Low complexity:     8980    476511
# LTR elements:       613     265060
# Rolling-circles     432     342046
# Small RNA:          18      106799
# Satellites:         697     293417
# Simple repeats:     21373   921658
# Unclassified:       17875   4764645

tbl_nums <- tribble(~ class, ~ elements, ~length,
                    "DNA",               392,     131109,
                    "LINE",              97,      58363,
                    "Low_complexity",    8980,    476511,
                    "LTR",               613,     265060,
                    "RC",                432,     342046,
                    "rRNA",              18,      106799,
                    "Satellite",         697,     293417,
                    "Simple_repeat",     21373,   921658,
                    "Unknown",           17875,   4764645)

# Try to replicate this as closely as possible (`elements` I know is accurate):
z_nums <- z |> mutate(begin2 = as.numeric(ifelse(unk1 == "+", begin2, left2)),
            class = str_remove(class, "\\/.*")) |>
    group_by(class) |>
    summarize(# elements = length(unique(id)),
              length = sum((end - begin + 1) * (1 - (p_del+p_ins) / 100))) |>
    arrange(class)

z_nums |>
    mutate(length = ((length - tbl_nums$length) / tbl_nums$length)) |>
    (function(x) {
        print(x)
        cat("\n")
        x |> select(-class) |> summarize(abs_mean = mean(abs(length)),
                                         mean = mean(length))
    })()


x <- "~/Desktop/AaegL5_TE_repeats.gff3" |>
    read_tsv(col_names = c("seq", "source", "type", "start", "end", "score",
                           "strand", "phase", "attr"), comment = "#")

x$attr |> head() |> str_remove("^Target=") |> str_split(" ") |> map_chr(~ .x[[1]])





#' using `(end - begin + 1)`:
#'  - `abs_mean` =  0.0588
#'  - `mean` =      0.0588
#' using `(end2 - begin2 + 1)`:
#'  - `abs_mean` =  0.0876
#'  - `mean` =      0.0764
#'
#' using `(end - begin + 1) * (1 - p_div / 100)`:
#'  - `abs_mean` =  0.0844
#'  - `mean` =     -0.0436
#' using `(end2 - begin2 + 1) * (1 - p_div / 100)`:
#'  - `abs_mean` =  0.106
#'  - `mean` =     -0.0278
#'
#' using `(end - begin + 1) * (1 - p_del / 100)`:
#'  - `abs_mean` =  0.0505
#'  - `mean` =      0.0271
#' using `(end2 - begin2 + 1) * (1 - p_del / 100)`:
#'  - `abs_mean` =  0.0734
#'  - `mean` =      0.0425
#'
#' using `(end - begin + 1) * (1 - (p_del+p_ins) / 100)`:
#'  - `abs_mean` =  0.0516
#'  - `mean` =      0.00751
#' using `(end2 - begin2 + 1) * (1 - (p_del+p_ins) / 100)`:
#'  - `abs_mean` =  0.0758
#'  - `mean` =      0.0233
#'









all_rbfa_classes <- "~/_data/RepBase/RepBase28.01.fasta/invrep.ref" |>
    read_lines() |>
    keep(~ grepl("^>", .x)) |>
    str_split("\t") |>
    map_chr(~.x[[2]]) |>
    tolower()

rbfa_classes <- unique(all_rbfa_classes)

match_class <- function(from, to, split = TRUE) {
    if (split) {
        to <- str_split(to, "\\/") |>
            map_chr(~ paste0("\\b", .x, "\\b", collapse = "|"))
    }
    str_detect(from, to)
}


mean(map_lgl(rbfa_classes, ~ any(match_class(df_classes, .x))))
mean(map_lgl(rbfa_classes, ~ any(match_class(rm_classes, .x))))
mean(map_lgl(rbfa_classes, ~ any(match_class(df_classes, .x) | match_class(rm_classes, .x))))

sum(map_lgl(rbfa_classes, ~ !any(match_class(df_classes, .x) | match_class(rm_classes, .x))))
# which(map_lgl(rbfa_classes, ~ !any(match_class(df_classes, .x) | match_class(rm_classes, .x))))

# For a custom library, the header format must follow: >Repeat_name#CLASS/Subclass with CLASS in "DNA, LINE, LTR, SINE, MITE, Helitron, Simple Repeat, Satellite"

#' From above, the unique CLASS elements are:
rm_classes |> str_split("\\/") |> map_chr(~ .x[[1]]) |> unique()
# [1] "other"          "artefact"       "segmental"      "low_complexity"
# [5] "unknown"        "simple_repeat"  "satellite"      "rna"
# [9] "dna"            "rc"             "line"           "scrna"
# [13] "rrna"           "trna"           "snrna"          "retroposon"
# [17] "sine"           "ltr"            "ple"

#' Unique Subclass elements:
#' (commented bc length == 214)
# rm_classes |> str_split("\\/") |>
#     map_chr(~ if (length(.x) < 2) { NA } else { .x[[2]] }) |>
#     unique()

#' Classes that do NOT have matches:
rbfa_nomatch <- all_rbfa_classes |>
    (function(.arf) {
        lgl1 <- map_lgl(.arf, ~ any(match_class(df_classes, .x)))
        lgl2 <- map_lgl(.arf, ~ any(match_class(rm_classes, .x)))
        .arf[! (lgl1 | lgl2) ]
    })() |>
    table() |> as.data.frame() |> as_tibble() |>
    set_names(c("class", "freq")) |> arrange(desc(freq))
rbfa_nomatch |> print(n = 20)

#' Some manually created maps:
manual_matches <-  list("ltr retrotransposon" = "ltr",
                        "dna transposon" = "dna",
                        "non-ltr retrotransposon" = "",
                        "sola1" = "dna/sola-1",
                        "sola2" = "dna/sola-2",
                        "sola3" = "dna/sola-3",
                        "rtex" = "line"

)


#' Classes that do have matches:
rbfa_match <- all_rbfa_classes |>
    (function(.arf) {
        lgl1 <- map_lgl(.arf, ~ any(match_class(df_classes, .x)))
        lgl2 <- map_lgl(.arf, ~ any(match_class(rm_classes, .x)))
        .arf[(lgl1 | lgl2) ]
    })() |>
    table() |> as.data.frame() |> as_tibble() |>
    set_names(c("class", "freq")) |> arrange(desc(freq))
rbfa_match |> print(n = 200)


#' Keyword (KW) and ID lines for RepBase embl file
rb_em_kwids <- "~/_data/RepBase/RepBase28.01.embl/invrep.ref" |>
    read_lines() |>
    keep(~ grepl("^ID|^KW", .x))


id_inds <- which(grepl("^ID", rb_em_kwids))
all_rb_em_classes <- map_chr(1:length(id_inds), function(i) {
    if (i == length(id_inds)) {
        .inds <- (id_inds[i]+1):length(rb_em_kwids)
    } else .inds <- (id_inds[i]+1):(id_inds[i+1]-1)
    .str <- rb_em_kwids[.inds] |>
        str_remove_all("KW\\s+|\\.") |>
        str_c(collapse = " ") |>
        tolower() |>
        str_split("; ") |>
        getElement(1) |>
        head(-1)


        identity()
    return(.str)
})

rb_em_classes <- unique(all_rb_em_classes)


nomatch_inds <- which(str_detect(all_rb_em_classes, "transposable element") &
    ! map_lgl(all_rbfa_classes, ~ any(match_class(df_classes, .x) | match_class(rm_classes, .x))))

tibble(em = all_rb_em_classes[nomatch_inds],
       fa = all_rbfa_classes[nomatch_inds]) |>
    # filter(fa == "dna transposon") |>
    group_by(fa) |>
    summarize(n = n()) |>
    filter(n > 1) |>
    arrange(desc(n)) |>
    print(n = 200)




all_rb_em_classes[nomatch_inds] |> head()
all_rbfa_classes[nomatch_inds] |> head()

# em_matcher <- function(em_string) {
# rm(em_string, splt_em_str, .x, rm_matches, df_matches)
em_string <- all_rb_em_classes[7856]
splt_em_str <- str_split(em_string, "; ")[[1]]

.x <- splt_em_str[2]
rm_matches <- rm_classes[str_detect(rm_classes, .x)]
if(length(rm_matches) == 0) {
    df_matches <- df_classes[str_detect(df_classes, .x)]
}





# }

