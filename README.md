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
6) Enrichment score is calculated as per the folloing formula.
7) 
enrichment score=((Number of TE subfamilies>1cpm)/(Number of TEs>1cpm))/((Number of TE subfamilies)/(Number of TEs))

# How to use the pipeline

1) Directory structure - For a particular set of samples please create a target_directory with any name of your choice. This directory should include "read_names.tsv" and "config.tsv" files as well as "raw_fastq_files", "trimmed_fastq_files", "aligned_files", "featureCounts". If the above mentioned folders are expected to contain files in the downstream steps of the pipeline, they will be created on their own in the target_dir.
2) Fill the read file (sample in 'supporting_files' directory - read_names.tsv) which contains information about the cells and the steps you would like to implement. You can use the following R commands and make the required changes to generate the read_names.tsv file. 

```
files <- list.files("directory/to/raw/fastq/files", "extension")
cells <- unique(files %>% gsub("suffix.*", "", .))
paired <- rep("TRUE", length(cells))
raw_reads_1 <- paste0(cells, "_R1.fastq.gz")
raw_reads_2 <- paste0(cells, "_R2.fastq.gz")
trim <- rep("TRUE", length(cells))
trimmed_reads_1 <- rep("NULL", length(cells))
trimmed_reads_2 <- rep("NULL", length(cells))
align <-  rep("NULL", length(cells))
fc <- rep("NULL", length(cells))q
fc_files <- rep("NULL", length(cells))
sc <- rep("TRUE", length(cells))

read_names <- data.frame(cells, paired, raw_reads_1, raw_reads_2,
                           trim, trimmed_reads_1, trimmed_reads_2, align,
                           aligned_files, fc, fc_files, sc)
 
fwrite(read_names,
       file="target/directory/read_names.tsv",
       quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
       
```

#### Brief explanation of read_names.tsv file

* The read_names.tsv contains two types of columns - logical and variable. "paired" (whether paired-end), "trim", "align", "fc" (feature counts) and "sc" (whether the cell/sample is single-cell) are logical columns. Please use "TRUE"/"FALSE" for yes/no. "paired" and "sc" columns are the sepcifications for the data type and "trim", "align", "fc" are the specifications if you want to carry out those steps in the pipeline. "raw_reads_1", "raw_reads_2", "trimmed_reads_1", "trimmed_reads_2" and "fc_files" are the variable columns. These store the names for the resepctive steps' files. Note: in case of single-end data, set "raw_reads_2" to be "NULL"
* Depending upon the steps you want the pipeline to run and skip, you will have to fill the read_names.tsv file. For example, if you don't want to carry out trimming step but want all the other steps downstream of it, you will have to fill the "trimmed_reads_1" and "trimmed_reads_2" columns. Even if the read files are raw and you still don't want to carry out trimming step, use the raw read names in "trimmed_reads_1" and "trimmed_reads_2" columns. In case of single-end data, set trimmed_reads_2 to be "NULL".
* If you want to run all the steps in the pipeline, set raw_read_1 and raw_read_2 columns with R1 and R2 files and rest of the variable columns with "NULL". As the pipeline proceeds sample-wise (row-wise) and finishes steps in the pipeline, it replaces variable names with the respective output file names and sets the logical columns, that it has finished, to be "FALSE". So, when you run the pipeline again, it won't run the steps for the samples that it has already finished
* Note: if you don't want run any of the steps prior or post to your step of interest, replace all the respective logical variables to be "FALSE". The variable columns won't matter in that case but for sake of clarity, you can just use "NULL".


3) Fill the config.tsv file as per the following.
*   feature_annotation_file - .saf annotation file (generated using generate_supporting_files.R) including genes, Mt and TE annotation.
*   intronic_intergenic_tes - file (generated using generate_supporting_files.R) including TE id's in feature_annotation_file and their corresponding loci, family name, class etc.
*   trim_galore - direction to trim_galore (can use "which trim_galore" in command line to get the location, copy and paste)
*   cutadapt - direction to cutadapt
*   trim_opts - trimming options you want to use with trim_galore. (just write them as you would write with trim_galore)
*   hisat2 - direction to hisat2
*   samtools - direction to samtools
*   hisat2_genome - direction to hisat2 prepared reference genome
*   hisat2_opts - options you want to use with hisat2 (just write them as you would write with hisat2)
*   featureCounts - direction to featureCounts
*   threads - number of threads

4) You can use the following chunk of code in R and make tweaks as per your data and needs

```
source("directory/to/TE_quantification_pipeline.R")

target_dir <- "main/target/directory"
config_file <- "target_dir/config_file.tsv"
read_names_file <- "target_dir/read_names.tsv"
sample <- "STORM-seq" ## or whatever same you want to use

create_count_matrix <- "TRUE"      ##TRUE if you want to create a raw count matrix (will be stored in the target_dir as count_matrix.tsv), else FALSE
scuttle_filter <- "TRUE"           ##TRUE if it is single cell data and you want to perform a quick filtering using scuttle, else FALSE 
log_enrich <- "TRUE"               ##TRUE if you want to calculate log enrichment of TE's (will be stored in the target_dir as log_enrichment.tsv and filtered_log_enrichment.tsv - if the scuttle filtering was on), else FALSE 

quantify_TEs(config_file, read_names_file, sample, target_dir, create_count_matrix, scuttle_filter, log_enrich)

```

