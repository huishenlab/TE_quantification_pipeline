
################################################ Generate supporting files for TE quantification pipeline ####################################################
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(stringr)

options(scipen=999)

generate_supporting_files <-  function(annotation_gtf, repeat_masker_file, chrom_length_file, output_dir, skips){
  
  ref <- fread(annotation_gtf, skip = skips, nThread = 6) 
  #ref$Chr <- gsub("chr","", ref$V1)
  ref_exon <- ref[ref$V3 == "exon"]
  ref_exon$gene_type <- gsub(".*gene_biotype=", "", ref_exon$V9) %>% gsub(";.*", "", .)
  ref_exon$gene_name <- gsub(".*gene_name=", "", ref_exon$V9) %>% gsub(";.*", "", .)
  
  ref_exon_mRNA <- ref_exon[ref_exon$gene_type == "protein_coding", ]
  ref_exon_mRNA_saf <- ref_exon_mRNA[, c(11, 1, 4, 5, 7)] %>% unique()
  colnames(ref_exon_mRNA_saf) <- c("GeneID", "Chr", "Start","End", "Strand")
  ref_exon_mRNA_saf.gr <- makeGRangesFromDataFrame(ref_exon_mRNA_saf, keep.extra.columns = T)
  ref_exon_mRNA_saf.gr_ns <- ref_exon_mRNA_saf.gr
  strand(ref_exon_mRNA_saf.gr_ns) <- "*"
  
  ref_tx <- ref_exon_mRNA_saf[, list(start = min(Start), end = max(End)), by = c("GeneID", "Strand", "Chr")]
  ref_tx.gr <- makeGRangesFromDataFrame(ref_tx, keep.extra.columns = T)
  ref_tx.gr_ns <- ref_tx.gr
  strand(ref_tx.gr_ns) <- "*"
  
  te <- fread(repeat_masker_file, fill = T)
  #te$Chr <- gsub("chr","", te$V5)
  te$class <- gsub("/.*", "", te$V11)
  te$class <- gsub("\\?", "", te$class)
  
  te$family <- gsub(".*/", "", te$V11)
  te$family <- gsub("\\?", "", te$family)
  
  te_sub <- te[, c(5, 6, 7, 9, 10, 16, 17)]
  colnames(te_sub) <- c("chr", "start","end","strand", "repName", "repClass", "repFamily")
  te_sub$strand[te_sub$strand == "C"] <- "-"
  selected_classes <- c("DNA", "LINE", "LTR", "SINE","RC", "Retroposon")
  
  tes_class_family <- te_sub[te_sub$repClass %in% selected_classes,] 
  tes_class_family$name <- paste0("te_", 1:nrow(tes_class_family))
  
  te_sub.gr <- makeGRangesFromDataFrame(te_sub[te_sub$repClass %in% selected_classes], keep.extra.columns = T)
  te_sub.gr$name <- paste0("te_", 1:length(te_sub.gr))
  
  ### exon te chrM saf
  
  te_sub.gr_to_genes <- findOverlapPairs(te_sub.gr, ref_tx.gr_ns, ignore.strand = T)
  intergenic_tes <- te_sub.gr[!te_sub.gr$name %in% te_sub.gr_to_genes@first$name]
  intergenic_tes$annot <- rep("intergenic", length(intergenic_tes))
  te_sub.gr_to_exons <- findOverlapPairs(te_sub.gr, ref_exon_mRNA_saf.gr_ns, ignore.strand = T)
  intronic_tes <- te_sub.gr[!te_sub.gr$name %in% c(te_sub.gr_to_exons@first$name, intergenic_tes$name)]
  intronic_tes$annot <- rep("intronic", length(intronic_tes))
  intergenic_intronic_tes <- c(intergenic_tes, intronic_tes)
  
  write.table(as(intergenic_intronic_tes, "data.frame"),
              file = paste0(output_dir, "/intergenic_intronic_tes.txt"),
              sep ="\t", quote = FALSE,
              col.names = TRUE,
              row.names = FALSE)
  
  tes_saf <- data.table(GeneID = intergenic_intronic_tes$name, 
                        Chr = as.character(seqnames(intergenic_intronic_tes)),
                        Start = as.numeric(start(intergenic_intronic_tes)),
                        End = as.numeric(end(intergenic_intronic_tes)),
                        Strand = as.character(strand(intergenic_intronic_tes)))
  
  chr_lengths <- fread(chrom_length_file)
  #chr_lengths$Chr <- gsub("chr","", chr_lengths$V1)
  chrM <- chr_lengths[chr_lengths$V1 == "chrM",]
  chrM_saf <- data.table(GeneID = "M", Chr = "chrM", Start = 1, End = chrM$V2, Strand = "+")
  
  
  selected_cols <- paste0("chr", c(1:22, "X", "Y"), sep="")
  # selected_cols <- c("1", "2", "3", "4", "5", "6", "7", "8", 
  #                    "9", "10", "11", "12", "13", "14", "15",
  #                    "16", "17", "18", "19", "20", "21", "22", 
  #                    "X", "Y")
  
  ref_exon_mRNA_saf_sub <- ref_exon_mRNA_saf[ref_exon_mRNA_saf$Chr %in% selected_cols]
  ref_exon_mRNA_saf_sub$type <- "pc"
  tes_saf_sub <- tes_saf[tes_saf$Chr %in% selected_cols]
  tes_saf_sub$type <- "te"
  chrM_saf$type <- "M"
  exon_tes_chrM_saf <- rbind(ref_exon_mRNA_saf_sub, tes_saf_sub, chrM_saf)
  
  write.table(as(exon_tes_chrM_saf, "data.frame"), 
              file = paste0(output_dir, "/pc_te_chrM.saf"), 
              row.names = F, col.names = T, 
              quote = F, sep = "\t")
  
}
