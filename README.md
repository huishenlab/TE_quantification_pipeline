# Overview

This a pipeline for the quantification of Transposable Elements in the single-cell STORM-seq samples. This pipeline can be used with other technologies like SMART-seq, SMART-seq2 etc. 

The pipeline takes a config_file.tsv and a read_names.tsv as inputs and can produce raw count matrices, TE count matrices and and TE enrichment as outputs. Details are dicussed in the latter section. 

1) trim_galore
2) cutadapt
3) fastqc
4) HISAT2
5) samtools
6) featureCounts
