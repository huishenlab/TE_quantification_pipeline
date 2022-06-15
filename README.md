# Overview

This a pipeline for the quantification of Transposable Elements in the single-cell STORM-seq samples. This pipeline can be used with other technologies like SMART-seq, SMART-seq2 etc. 

The pipeline takes a config_file.tsv and a read_names.tsv as inputs and can produce raw count matrices, TE count matrices and and TE enrichment as outputs. Details are dicussed in the latter section. 

# Prerequisites

## Command Line tools

1) trim_galore
2) cutadapt
3) fastqc
4) HISAT2
5) samtools
6) featureCounts

## R libraries

1) data.table
2) dplyr
3) scales
4) stringr
5) GenomicRanges

The current version of the pipeline is designed to be used with CHM13 v2.0 assembly and hence the supporting files (intergenic_intronic_tes.txt and pc_te_chrM.saf) are generated using the CHM13 assembly's catLiftOffGenesV1.gff3 (https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/genes/catLiftOffGenesV1.gff3.gz), chm13_chrom_sizes.txt (https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/GCA_009914755.4.chrom.sizes.txt and https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/GCA_009914755.4.chromAlias.txt), GCA_009914755.4_CHM13_T2T_v2.0.sorted_chrNames.fa.out (https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/GCA_009914755.4.repeatMasker.out.gz) files. The intergenic_intronic_tes.txt and pc_te_chrM.saf files can be generated for the assembly of choice using the generate_supporting_files.R with the equivalent of the above mentioned required files. 


** More details coming soon **
