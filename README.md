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

# Steps in the pipeline

1) Raw reads are first trimmed using trim-galore with the recommned options (--illumina --trim-n -length 36 --clip_R2 3 --fastqc --cores) for STORM-seq data. If the pipeline is used with other technologies appropriate options can be specified in the config.tsv file provided in the supporting_files directory.
2) The trimmed reads are then aligned to the reference genome (default supporting files are provided for CHM13 v2.0 assembly) using HISAT2 with the "--dta" option which can be updated in the config.tsv file as per the requirement. Same is applicable for the HISAT2 indexed genome (It must be updated as per your own HISAT2 indexed genome location). Also, if the reads are paired-end, duplicate marking is performed using samtools and only the duplicate marked coordinate sorted aligned .bam files are retained.
3) The duplicate marked bam files are used to count reads in each of the featues (genes and TE's) in the pc_te_chrM.saf file using featureCounts and options "--fracOverlap 0.1 -M -s 0 --fraction". 
4) feature counts files for all the cells is then used to build combined raw count, counts per million, raw count for only intergenic and intronic TEs and counts per millions for only intergenic and intronic TEs matrices. 
5) counts per millions for only intergenic and intronic TEs matrix is used for the log enrichment calculation of TE's


** More details coming soon **
