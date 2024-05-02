#!/usr/local/bin/Rscript

#Install packages
if(!require(pacman)) install.packages("pacman")

pacman::p_load(
  readr,
  dplyr,
  data.table)

result = read.table("${score_files}", header = TRUE, sep = "\\t")

scoring_file <- data.frame(SNP = result\$SNP, A1 = result\$effect_allele, BETA = result\$effect_weight)
scoring_file_uniq <- scoring_file[!duplicated(scoring_file\$SNP), ]
score_file <- na.omit(scoring_file_uniq)
clean_score_file <- score_file[score_file[, 1] != "", , drop = FALSE]
write.table(clean_score_file, file = "${blood_trait}_${pgs_id}_modified.txt", sep = "\\t", col.names = T, row.names = F, quote = F)

