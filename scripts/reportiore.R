if (!require("scRepertoire")) {
  install.packages("devtools")
  devtools::install_github("ncborcherding/scRepertoire")
}
library(scRepertoire)
library(Seurat)
library(tidyverse)
tcr_files <- list.files(pattern = "full-length_productive_TCR_table.tsv$")
bcr_files <- list.files(pattern = "full-length_productive_BCR_table.tsv$")
# Load TCR
tcr_combined <- combineTCR(
  tcr_list,
  samples = sapply(tcr_list, function(x) unique(x$sample))
)
tcr_combined <- combineTCR(
  tcr_list,
  samples = sapply(tcr_list, function(x) unique(x$sample))
)
#bcr
bcr_list <- lapply(bcr_list, function(df) {
  df$chain <- ifelse(grepl("^IGH", df$v_gene), "IGH",
                     ifelse(grepl("^IGK", df$v_gene), "IGK",
                            ifelse(grepl("^IGL", df$v_gene), "IGL", NA)))
  return(df)
})
# Fix all columns to ensure no list/factor issues
bcr_list <- lapply(bcr_list, function(df) {
  df$barcode     <- as.character(df$barcode)
  df$v_gene      <- as.character(df$v_gene)
  df$j_gene      <- as.character(df$j_gene)
  df$cdr3        <- as.character(df$cdr3)
  df$cdr3_nt     <- as.character(df$cdr3_nt)
  df$sample      <- as.character(df$sample)
  df$chain       <- as.character(df$chain)
  df$celltype    <- as.character(df$celltype)
  df$patient_id  <- as.character(df$patient_id)
  df$tissue_type <- as.character(df$tissue_type)
  df$disease     <- as.character(df$disease)
  df$batch_id    <- as.character(df$batch_id)
  return(df)
})

bcr_combined <- combineBCR(
  bcr_list,
  samples = sapply(bcr_list, function(x) unique(x$sample))
)
bcr_df <- do.call(rbind, bcr_list)
tcr_df <- do.call(rbind, tcr_combined)
dim(tcr_df)

# View basic structure
str(tcr_combined)
str(bcr_df)

# Check sample distribution
table(tcr_combined$sample)
table(bcr_df$sample)

# Check chain types
table(tcr_combined$chain)
table(bcr_df$chain)

# Unique CDR3 sequences (for clonotype diversity)
length(unique(tcr_combined$cdr3))
length(unique(bcr_df$cdr3))


