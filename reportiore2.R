install.packages("devtools")
devtools::install_github("ncborcherding/scRepertoire")
library(scRepertoire)
library(Seurat)
library(tidyverse)
# List all files
tcr_files <- list.files(pattern = "full-length_productive_TCR_table.tsv$")
bcr_files <- list.files(pattern = "full-length_productive_BCR_table.tsv$")

# Read TCR files into a named list
tcr_list <- lapply(tcr_files, function(f) {
  df <- tryCatch(read.delim(f), error = function(e) NULL)
  if (!is.null(df)) df$sample <- sub(".tsv.gz$", "", f)  # Add sample name
  return(df)
})
names(tcr_list) <- tcr_files
tcr_list <- Filter(Negate(is.null), tcr_list)  # Remove failed reads

# Read BCR files into a named list
bcr_list <- lapply(bcr_files, function(f) {
  df <- tryCatch(read.delim(f), error = function(e) NULL)
  if (!is.null(df)) df$sample <- sub(".tsv.gz$", "", f)  # Add sample name
  return(df)
})
names(bcr_list) <- bcr_files
bcr_list <- Filter(Negate(is.null), bcr_list)
# Define required columns
required_cols <- c("barcode", "v_gene", "j_gene", "cdr3", "cdr3_nt", "sample")

# Clean function
clean_df <- function(df, required_cols) {
  for (col in required_cols) {
    if (!col %in% colnames(df)) df[[col]] <- NA
    df[[col]] <- as.character(df[[col]])
  }
  return(df[, required_cols])
}

# Clean and bind TCR
tcr_df <- do.call(rbind, lapply(tcr_list, clean_df, required_cols))

# Clean and bind BCR
bcr_df <- do.call(rbind, lapply(bcr_list, clean_df, required_cols))
bcr_df$chain <- ifelse(grepl("^IGH", bcr_df$v_gene), "IGH",
                       ifelse(grepl("^IGK", bcr_df$v_gene), "IGK",
                              ifelse(grepl("^IGL", bcr_df$v_gene), "IGL", NA)))
# clean
tcr_df <- tcr_df[!is.na(tcr_df$cdr3) & nchar(tcr_df$cdr3) >= 6, ]
bcr_df <- bcr_df[!is.na(bcr_df$cdr3) & nchar(bcr_df$cdr3) >= 6, ]
library(dplyr)

# TCR clonotypes
tcr_df <- tcr_df %>%
  mutate(clonotype_id = as.numeric(factor(cdr3)))

# BCR clonotypes
bcr_df <- bcr_df %>%
  mutate(clonotype_id = as.numeric(factor(cdr3)))
# TCR
tcr_freq <- tcr_df %>%
  group_by(clonotype_id) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

# BCR
bcr_freq <- bcr_df %>%
  group_by(clonotype_id) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
# vis
library(dplyr)

# Count TCR clonotype frequencies
tcr_freq <- tcr_df %>%
  group_by(clonotype_id) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

# View top 10
head(tcr_freq, 10)
library(ggplot2)

ggplot(tcr_freq[1:10, ], aes(x = reorder(as.factor(clonotype_id), -count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Top 10 TCR Clonotypes",
       x = "Clonotype ID",
       y = "Cell Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# same for bcr 
library(dplyr)
library(vegan)  # For diversity calculation

# Function to compute Shannon index for each sample
compute_shannon <- function(df) {
  df %>%
    group_by(sample, clonotype_id) %>%
    summarise(count = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = clonotype_id, values_from = count, values_fill = 0) %>%
    column_to_rownames("sample") %>%
    diversity(index = "shannon")
}

# For TCR
tcr_shannon <- compute_shannon(tcr_df)
tcr_shannon_df <- data.frame(sample = names(tcr_shannon), shannon = tcr_shannon, receptor = "TCR")

# For BCR
bcr_shannon <- compute_shannon(bcr_df)
bcr_shannon_df <- data.frame(sample = names(bcr_shannon), shannon = bcr_shannon, receptor = "BCR")

# Combine
diversity_df <- rbind(tcr_shannon_df, bcr_shannon_df)
library(ggplot2)

ggplot(diversity_df, aes(x = basename(sample), y = shannon, fill = receptor)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Shannon Diversity Index by Sample",
       x = "Sample",
       y = "Shannon Index (Diversity)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#
library(scRepertoire)

# Rebuild from raw list of data frames (tcr_list)
# Make sure each data frame has the expected structure
# Add chain annotation based on v_gene
tcr_list <- lapply(tcr_list, function(df) {
  df$chain <- ifelse(grepl("^TRA", df$v_gene), "TRA",
                     ifelse(grepl("^TRB", df$v_gene), "TRB", NA))
  return(df)
})
tcr_combined <- combineTCR(
  tcr_list,
  samples = names(tcr_list),
  ID = "sample"
)
clonalHomeostasis(tcr_combined, cloneCall = "aa", chain = "both")
clonalProportion(tcr_combined, cloneCall = "aa", chain = "both")
clonalOverlap(tcr_combined, cloneCall = "aa", method = "morisita")
tcr_combined <- lapply(tcr_combined, function(df) {
  df$barcode <- sub(".*_", "", df$barcode)  # remove prefix up to last underscore
  return(df)
})

seurat <- combineExpression(tcr_combined, seurat, cloneCall = "aa", group.by = "sample")
DimPlot(seurat, group.by = "Cluster_id")  # if clone_id added as metadata
bcr_combined <- lapply(bcr_combined, function(df) {
  df$barcode <- sub(".*_", "", df$barcode)  # Keep only the cell barcode
  return(df)
})
# merge bcr files
library(dplyr)

# Step 1: Clean sample names
sample_names <- sub("GSM\\d+_", "", names(bcr_list))
sample_names <- sub("_full-length_productive_BCR_table.tsv", "", sample_names)

# Step 2: Add clean sample names to each df
for (i in seq_along(bcr_list)) {
  bcr_list[[i]]$sample <- sample_names[i]
}

# Step 3: Combine all BCRs into one big dataframe
bcr_df <- bind_rows(bcr_list)

# Step 4: Add chain info (if missing)
if (!"chain" %in% colnames(bcr_df)) {
  bcr_df$chain <- ifelse(grepl("^IGH", bcr_df$v_gene), "IGH",
                         ifelse(grepl("^IGK", bcr_df$v_gene), "IGK",
                                ifelse(grepl("^IGL", bcr_df$v_gene), "IGL", NA)))
}

# Step 5: Keep only productive sequences with valid CDR3
bcr_df <- bcr_df %>%
  filter(!is.na(cdr3), nchar(cdr3) >= 6)

# Step 6: Make barcode unique (remove prefixes if any)
bcr_df$barcode <- sub(".*_", "", bcr_df$barcode)

# Step 7: Assign clonotype IDs based on amino acid sequence
bcr_df$clonotype_id <- as.numeric(factor(bcr_df$cdr3))

# Step 8: Optional summary — view top clonotypes
top_clones <- bcr_df %>%
  count(clonotype_id, sort = TRUE) %>%
  head(10)
print(top_clones)
# add to seurat 
# Let's assume 'seurat' is your existing Seurat object
# We'll join clonotype info into metadata

# Format for merge
bcr_meta <- bcr_df %>%
  distinct(barcode, sample, clonotype_id, .keep_all = TRUE) %>%
  select(barcode, clonotype_id)

# Make sure barcodes match those in Seurat
bcr_meta$barcode <- gsub("-1$", "", bcr_meta$barcode)  # if needed to match Seurat format

# Set barcode as rownames for easy joining
bcr_meta_unique <- bcr_meta[!duplicated(bcr_meta$barcode), ]
rownames(bcr_meta_unique) <- bcr_meta_unique$barcode

# Then merge with Seurat
seurat <- AddMetaData(seurat, metadata = bcr_meta_unique["clonotype_id"])

# Merge into Seurat
library(dplyr)
library(ggplot2)

bcr_freq <- bcr_df %>%
  count(clonotype_id, sort = TRUE)

ggplot(bcr_freq[1:10, ], aes(x = reorder(as.factor(clonotype_id), -n), y = n)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  labs(title = "Top 10 BCR Clonotypes",
       x = "Clonotype ID", y = "Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DimPlot(seurat, group.by = "clonotype_id", label = FALSE, pt.size = 1) +
  ggtitle("UMAP Colored by BCR Clonotype")
expanded_clones <- bcr_df %>%
  count(clonotype_id) %>%
  filter(n >= 3) %>%
  pull(clonotype_id)

seurat$expanded_clone <- ifelse(seurat$clonotype_id %in% expanded_clones, "Expanded", "Not Expanded")

DimPlot(seurat, group.by = "expanded_clone") +
  ggtitle("Expanded BCR Clonotypes")





