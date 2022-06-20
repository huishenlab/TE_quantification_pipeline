library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(stringr)

options(scipen=999)

################################################ function to trim reads #####################################################

trim_reads <- function(trim_galore, cutadapt, trim_opts, read_names_file, target_dir){
  
  read_names <- fread(read_names_file, header = TRUE) %>% as.data.frame()
  cells <- read_names$cells
  trim <- read_names$trim
  paired <- read_names$paired
  raw_reads_1 <- read_names$raw_reads_1
  raw_reads_2 <- read_names$raw_reads_2
  mk_trimmed_dir <- paste0("mkdir -p ", target_dir, "/trimmed_fastq_files","\n")
  
  for (i in 1:length(cells)){
    cell <- cells[i]
    
    if (trim[i] == "TRUE"){
      if (paired[i] == "TRUE"){
        trim_reads <- paste0(trim_galore, " ", trim_opts,
                             " --path_to_cutadapt ", cutadapt , " --basename ", cell," -o ", target_dir, "/trimmed_fastq_files",
                             " --paired ",  target_dir, "/raw_fastq_files", "/", raw_reads_1[i], " ",  target_dir, "/raw_fastq_files", "/", raw_reads_2[i], "\n")
        add_trimmed_name <- paste0("sed -Ei \"", as.character(i + 1), " s/[^\t]*/\"", cell, "_val_1.fq.gz\"/6\" ", read_names_file, "\n")
        add_trimmed_name <- paste0(add_trimmed_name, "\nsed -Ei \"", as.character(i + 1), " s/[^\t]*/\"", cell, "_val_2.fq.gz\"/7\" ", read_names_file, "\n")
        add_trim_req <- paste0("sed -Ei \"", as.character(i + 1), " s/[^\t]*/\"FALSE\"/5\" ", read_names_file, "\n")
        }
      else {
        trim_reads <- paste0(trim_galore, " ", trim_opts,
                             " --path_to_cutadapt ", cutadapt , " --basename ", cell," -o ", target_dir, "/trimmed_fastq_files", " ",
                             target_dir, "/raw_fastq_files", "/", raw_reads_1[i], "\n")
        add_trimmed_name <- paste0("sed -Ei \"", as.character(i + 1), " s/[^\t]*/\"", cell, "_trimmed.fq.gz\"/6\" ", read_names_file, "\n")
        add_trim_req <- paste0("sed -Ei \"", as.character(i + 1), " s/[^\t]*/\"FALSE\"/5\" ", read_names_file, "\n")
      }
      system(paste0("set -e\n", mk_trimmed_dir, trim_reads, add_trimmed_name, add_trim_req))
    }
  }
}

################################################ function to align fastq's ###################################################

align_fastqs <- function(hisat2, samtools, hisat2_genome, hisat2_opts, threads, read_names_file, target_dir){

  read_names <- fread(read_names_file, header = TRUE) %>% as.data.frame()
  cells <- read_names$cells
  paired <- read_names$paired
  trimmed_fastq_dir <- paste0(target_dir, "/trimmed_fastq_files")
  trimmed_reads_1 <- read_names$trimmed_reads_1
  trimmed_reads_2 <- read_names$trimmed_reads_2
  align <- read_names$align
  aligned_dir <- paste0(target_dir, "/aligned_files")
  
  mk_aligned_dir <- paste0("mkdir -p ", aligned_dir, "\n")
  
  for (i in (1:length(cells))){
    cell <- cells[i]
    
    if (align[i] == "TRUE"){
      
      cd_to_aligned_dir <- paste0("cd ", aligned_dir, "\n")
      
      if (paired[i] == "TRUE"){
        hisat2_alignment <- paste0("( ", hisat2, " -x ", hisat2_genome, " ", hisat2_opts, " -1 ", 
                                   trimmed_fastq_dir, "/" , trimmed_reads_1[i]," -2 ", trimmed_fastq_dir, "/" , trimmed_reads_2[i], " -p ", as.character(threads), " | ", 
                                   samtools, " sort -@ ", as.character(threads), " -O BAM -o ", cell, ".bam ) > ", cell,"_alignment.log 2>&1\n",
                                   samtools, " index ", cell, ".bam\n", samtools, " sort -@ ", as.character(threads), " -n -O BAM -o ", cell, "_NS.bam ", cell, ".bam\n",
                                   samtools, " fixmate -@ ", as.character(threads), " -O BAM -m ", cell, "_NS.bam ", cell, "_NS_fixmate.bam\n",
                                   samtools, " sort -@ ", as.character(threads), " -O BAM -o ", cell, "_NS_fixmate_CS.bam ", cell, "_NS_fixmate.bam\n",
                                   samtools, " markdup -@ ", as.character(threads), " -O BAM ", cell, "_NS_fixmate_CS.bam ", cell, "_NS_fixmate_CS_dup_marked.bam\n",
                                   "rm ",cell, ".bam ", cell, ".bam.bai ", cell, "_NS.bam ", cell, "_NS_fixmate.bam ", cell, "_NS_fixmate_CS.bam\n")
        add_aligned_name <- paste0("sed -Ei \"", as.character(i + 1), " s/[^\t]*/\"", cell, "_NS_fixmate_CS_dup_marked.bam\"/9\" ", read_names_file, "\n")
        add_alig_req <- paste0("sed -Ei \"", as.character(i + 1), " s/[^\t]*/\"FALSE\"/8\" ", read_names_file, "\n")
      }
      else {
        hisat2_alignment <- paste0("( ", hisat2, " -x ", hisat2_genome, " ", hisat2_opts, " -U ", 
                                   trimmed_fastq_dir, "/" , trimmed_reads_1[i], " -p ", as.character(threads), " | ", 
                                   samtools, " sort -@ ", as.character(threads), " -O BAM -o ", cell, ".bam ) > ", cell,"_alignment.log 2>&1\n",
                                   samtools, " index ", cell, ".bam\n")
        add_aligned_name <- paste0("sed -Ei \"", as.character(i + 1), " s/[^\t]*/\"", cell, "_NS_fixmate_CS_dup_marked.bam\"/9\" ", read_names_file, "\n")
        add_alig_req <- paste0("sed -Ei \"", as.character(i + 1), " s/[^\t]*/\"FALSE\"/8\" ", read_names_file, "\n")
      }
      system(paste0("set -e\n", mk_aligned_dir, cd_to_aligned_dir, hisat2_alignment,
                    add_aligned_name, add_alig_req))
    }
  }
  
}

################################################ function to count features ###################################################

count_features <- function(feature_counts, feature_annotation_file, threads, read_names_file, target_dir){
  
  read_names <- fread(read_names_file, header = TRUE) %>% as.data.frame()
  cells <- read_names$cells
  paired <- read_names$paired
  aligned_dir <- paste0(target_dir, "/aligned_files")
  aligned_files <- read_names$aligned_files
  fc <- read_names$fc
  fc_output_dir <- paste0(target_dir, "/featureCounts")
  
  mk_fc_dir <- paste0("mkdir -p ", fc_output_dir, "\n")
  
  for (i in 1:length(cells)){
    
    cell <- cells[i]
    
    if (fc[i] == "TRUE"){
      if (paired[i] == "TRUE"){
        paired_option = "-p"
      }
      else{
        paired_option = ""
      }
      feature_counting <- paste0("( ", feature_counts, " -F SAF -O -B ", paired_option, " --fracOverlap 0.1 -M -s 0 --fraction -T ", 
                                 as.character(threads), " -a ", feature_annotation_file, " -o ", fc_output_dir, "/", cell, "_fc ", 
                                 aligned_dir, "/", aligned_files[i], " ) > ", fc_output_dir, "/", cell, "_fc.log 2>&1\n")
      
      add_fc_name <- paste0("sed -Ei \"", as.character(i + 1), " s/[^\t]*/\"", cell, "_fc\"/11\" ", read_names_file, "\n")
      add_fc_req <- paste0("sed -Ei \"", as.character(i + 1), " s/[^\t]*/\"FALSE\"/10\" ", read_names_file, "\n")
      
      system(paste0("set -e\n", mk_fc_dir, feature_counting, add_fc_name, add_fc_req))
    }
  }
}

###################################### function to create count matrix ###########################################

create_a_matrix_file <- function(read_names_file, intronic_intergenic_tes, target_dir){
  
  read_names <- fread(read_names_file, header = TRUE) %>% as.data.frame()
  fc_output_dir <- paste0(target_dir, "/featureCounts")
  fc_files <- read_names$fc_files
  cells <- read_names$cells
  single_cell <- read_names$sc
  
  for (i in 1:length(cells)){
    cell <- cells[i]
    
    if (i == 1){
      feat_counts_file <- fread(paste0(fc_output_dir, "/", fc_files[i]), select = c(1,7), 
                                nThread = 6)
      names(feat_counts_file)[ncol(feat_counts_file)] <- cell
      count_matrix <- feat_counts_file
    }
    else {
      feat_counts_file <- fread(paste0(fc_output_dir, "/", fc_files[i]), select = c(7), 
                                nThread = 6)
      names(feat_counts_file)[ncol(feat_counts_file)] <- cell
      count_matrix[, cell] <- feat_counts_file[, cell]
    }
  }
  
  Geneids <- count_matrix$Geneid
  
  colnames(count_matrix) <- c("Geneid", cells)
  fwrite(count_matrix %>% as.data.frame(),
         file=paste0(target_dir, "/count_matrix.tsv"),
         quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
  
  count_matrix <- count_matrix[,2:ncol(count_matrix)]
  
  if (single_cell == "TRUE"){
    qcstats <- scuttle::perCellQCMetrics(count_matrix) ###### quality control #######
    qcfilter <- scuttle::quickPerCellQC(qcstats)
    count_matrix_filtered <- count_matrix[,!qcfilter$discard]
    
    fwrite(cbind(count_matrix_filtered, Geneids) %>% as.data.frame(),
           file=paste0(target_dir, "/count_matrix_filtered.tsv"),
           quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
  }
  else {
    count_matrix_filtered <- count_matrix
  }
  
  intergenic_intronic_tes_df <- fread(intronic_intergenic_tes)  ####### extracting TE's ########
  colnames(intergenic_intronic_tes_df) <- gsub(colnames(intergenic_intronic_tes_df)[length(colnames(intergenic_intronic_tes_df)) - 1], "Geneid", colnames(intergenic_intronic_tes_df))
  
  count_matrix_filtered$Geneid <- Geneids
  count_matrix_filtered_TE <- inner_join(count_matrix_filtered, intergenic_intronic_tes_df[,c(6,7,8,9)], by = "Geneid")
  
  count_matrix_filtered <- count_matrix_filtered[, 1:(ncol(count_matrix_filtered)-1)]
  count_matrix_filtered_cpm <- sweep(count_matrix_filtered, 2, colSums(count_matrix_filtered)/1000000, `/`)
  count_matrix_filtered_cpm$Geneid <- Geneids
  
  count_matrix_filtered_cpm_TE <- inner_join(count_matrix_filtered_cpm, intergenic_intronic_tes_df[,c(6,7,8,9)], by = "Geneid")
  
  fwrite(count_matrix_filtered_TE %>% as.data.frame(),
         file=paste0(target_dir, "/count_matrix_filtered_TE.tsv"),
         quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
  
  fwrite(count_matrix_filtered_cpm_TE %>% as.data.frame(),
         file=paste0(target_dir, "/count_matrix_filtered_cpm_TE.tsv"),
         quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
  
  return(count_matrix_filtered_cpm_TE)
  
}

########################################### function to calculate enrichment ###################################################

log_enrichment <- function(tes_cpm, sample){
  
  cells <- colnames(tes_cpm)[1:(ncol(tes_cpm) - 4)]
  for (i in 1:(ncol(tes_cpm) - 4)){
    cell <- colnames(tes_cpm)[i]
    te_cpm <- tes_cpm[,c(cell, "Geneid", "repName", "repClass", "repFamily")]
    te_cpm_grt_1_cpm <- te_cpm[te_cpm[,cell] > 1, ]
    
    te_cpm_grouped <- te_cpm %>% group_by(repFamily,repClass)
    te_cpm_grouped_family_counts <- te_cpm_grouped %>% count(repFamily) 
    te_cpm_grouped_family_counts <- as.data.frame(te_cpm_grouped_family_counts)
    
    te_cpm_grt_1_grouped <- te_cpm_grt_1_cpm %>% group_by(repFamily,repClass)
    te_cpm_grt_1_grouped_family_counts <- te_cpm_grt_1_grouped %>% count(repFamily) 
    te_cpm_grt_1_grouped_family_counts <- as.data.frame(te_cpm_grt_1_grouped_family_counts)
    
    te_cpm_lst_1_grouped_family_counts <- te_cpm_grouped_family_counts[!te_cpm_grouped_family_counts$repFamily %in% te_cpm_grt_1_grouped_family_counts$repFamily,]
    te_cpm_lst_1_grouped_family_counts$n <- rep(0,nrow(te_cpm_lst_1_grouped_family_counts))
    
    te_cpm_grt_lst_1_grouped_family_counts <- rbind(te_cpm_grt_1_grouped_family_counts, te_cpm_lst_1_grouped_family_counts)
    
    te_cpm_grt_lst_1_grouped_family_counts$enrich_numerator <- (te_cpm_grt_lst_1_grouped_family_counts$n)/sum(te_cpm_grt_lst_1_grouped_family_counts$n) + 0.01
    te_cpm_grouped_family_counts$enrich_denominator <- (te_cpm_grouped_family_counts$n)/sum(te_cpm_grouped_family_counts$n) + 0.01
    
    te_enrichment <- inner_join(te_cpm_grt_lst_1_grouped_family_counts, te_cpm_grouped_family_counts, by="repFamily")
    te_enrichment$enrichment <- (te_enrichment$enrich_numerator)/(te_enrichment$enrich_denominator)
    te_enrichment$log_enrichment <- log2(te_enrichment$enrichment)
    te_enrichment$cells <- rep(cell, nrow(te_enrichment))
    te_enrichment$samples <- rep(sample, nrow(te_enrichment))
    
    if (i == 1){
      combined_Te_enrichment <- te_enrichment
    }
    else{
      combined_Te_enrichment <- rbindlist(list(combined_Te_enrichment, te_enrichment))
    }
  }
  return(combined_Te_enrichment)
}

################################################## input to output ##############################################################

quantify_TEs <- function(config_file, read_names_file, sample, target_dir, create_count_matrix, log_enrich){
  
  configs <- fread(config_file, header = TRUE) %>% as.data.frame()
  
  feature_annotation_file <- configs[configs$option == "feature_annotation_file", 2]
  intronic_intergenic_tes <- configs[configs$option == "intronic_intergenic_tes", 2]
  trim_galore <- configs[configs$option == "trim_galore", 2]
  cutadapt <- configs[configs$option == "cutadapt", 2]
  trim_opts <- configs[configs$option == "trim_opts", 2]
  hisat2 <- configs[configs$option == "hisat2", 2]
  samtools <- configs[configs$option == "samtools", 2]
  hisat2_genome <- configs[configs$option == "hisat2_genome", 2]
  hisat2_opts <- configs[configs$option == "hisat2_opts", 2]
  featureCounts <- configs[configs$option == "featureCounts", 2]
  threads <- configs[configs$option == "threads", 2]
  
  trim_reads(trim_galore, cutadapt, trim_opts, read_names_file, target_dir)
  align_fastqs(hisat2, samtools, hisat2_genome, hisat2_opts, threads, read_names_file, target_dir)
  count_features(featureCounts, feature_annotation_file, threads, read_names_file, target_dir)
  
  if (create_count_matrix == "TRUE"){
    tes_cpm <- create_a_matrix_file(read_names_file, intronic_intergenic_tes, target_dir)
  }
  
  if (log_enrich == "TRUE"){
    if (create_count_matrix == "TRUE"){
      enrich_out <- log_enrichment(tes_cpm, sample) %>% as.data.frame()
    }
    if (create_count_matrix == "FALSE"){
      tes_cpm <- fread(paste0(target_dir, "/count_matrix_filtered_cpm_TE.tsv"), header = TRUE) %>% as.data.frame()
      enrich_out <- log_enrichment(tes_cpm, sample) %>% as.data.frame()
    }
    fwrite(enrich_out %>% as.data.frame(),
           file=paste0(target_dir, "/log_enrichment.tsv"),
           quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
  }
}
