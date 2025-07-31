library(scRepertoire) # Make sure scRepertoire is loaded

# --- Step 1: Get and Filter Files ---
# This will find all .tsv and .tsv.gz files in the current directory
all_files <- list.files(pattern = "\\.tsv$|\\.tsv\\.gz$", full.names = TRUE)

# Filter for only the BCR productive tables
# The pattern correctly identifies your BCR files based on the list.files() output
bcr_files <- grep("_full-length_productive_BCR_table\\.tsv(\\.gz)?$", all_files, value = TRUE)

# Print the list to verify (optional, but good for checking)
message("Found BCR files:")
print(bcr_files)

# --- Step 2: Extract Sample Names for 'samples' and 'ID' arguments ---

# Extract sample names (e.g., C9_R, C12_R, C12_pBMC, C17_I)
# This regex captures the part between "GSM[0-9]+_" and "_full-length_productive_BCR_table"
my_samples <- sub(".*/GSM[0-9]+_([A-Za-z0-9_]+)_full-length_productive_BCR_table\\.tsv(\\.gz)?$", "\\1", bcr_files)
# load the tcr 
# load the data 
tcr_files <- list.files(pattern = "*TCR_table.tsv", full.names = TRUE)
# Load and add chain for TCR 
contig.frames <- lapply(files, function(file) {
  df <- read.delim(file)
  df$chain <- ifelse(grepl("TRBC", df$c_gene), "TRB",
                     ifelse(grepl("TRAC", df$c_gene), "TRA", NA))
  df
})
# Smarter sample name extraction
sample_names_tcr <- sub("^[^_]+_([^_]+_[^_]+).*", "\\1", basename(tcr_files))
sample_names_bcr <- sub("^[^_]+_([^_]+_[^_]+).*", "\\1", basename(bcr_files))

# Combine
combined.TCR <- combineTCR(contig.frames, samples = sample_names_tcr, ID = NULL)
# load the bcr
# Extract secondary IDs (e.g., R, pBMC, I) for the 'ID' argument
# This regex extracts the last segment after the last underscore
my_ids <- sub(".*_([A-Za-z0-9]+)$", "\\1", my_samples)

# Print to verify (optional)
message("\nGenerated Sample Names:")
print(my_samples)
message("\nGenerated IDs:")
print(my_ids)


# --- Step 3: Read and Pre-process Dataframes ---
bcr_frames <- lapply(bcr_files, function(file_path) {
  # Read the TSV (handling .gz compressed files)
  if (grepl("\\.gz$", file_path)) {
    df <- read.delim(gzfile(file_path), stringsAsFactors = FALSE)
  } else {
    df <- read.delim(file_path, stringsAsFactors = FALSE)
  }
  
  # Add 'productive' column if missing
  if (!"productive" %in% colnames(df)) {
    df$productive <- TRUE
  }
  
  # Add 'chain' column if missing
  if (!"chain" %in% colnames(df)) {
    df$chain <- ifelse(grepl("^IGK", df$v_gene), "IGK",
                       ifelse(grepl("^IGH", df$v_gene), "IGH",
                              ifelse(grepl("^IGL", df$v_gene), "IGL", NA)))
  }
  
  # --- CRITICAL: Add the 'reads' column ---
  # *** REPLACE "YOUR_ACTUAL_READS_COLUMN_NAME" with the name you find in your files ***
  # If you found a column named 'umis', use: if ("umis" %in% colnames(df)) { df$reads <- df$umis }
  # If you found a column named 'clone_reads', use: if ("clone_reads" %in% colnames(df)) { df$reads <- df$clone_reads }
  # etc.
  if ("UMI_count" %in% colnames(df)) { # Replace "UMI_count" with your actual column name
    df$reads <- df$UMI_count
  } else if ("duplicate_count" %in% colnames(df)) { # Common alternative
    df$reads <- df$duplicate_count
  } else if ("count" %in% colnames(df)) { # Another possible common name
    df$reads <- df$count
  } else if ("nReads" %in% colnames(df)) { # Yet another common name
    df$reads <- df$nReads
  } else {
    warning(paste0("No explicit count column found in ", basename(file_path), ". Setting 'reads' to 1 for this sample."))
    df$reads <- 1 # Fallback if no actual count column found
  }
  df$reads <- as.numeric(df$reads) # Ensure it's numeric
  
  return(df)
})

message("\nAll BCR dataframes read and pre-processed.")

# --- Step 4: Load Contigs with scRepertoire ---
contig_list <- loadContigs(input = bcr_frames)
message("Contigs loaded successfully into 'contig_list'.")

# --- Step 5: Combine BCR data ---
combined.BCR <- combineBCR(
  input.data = contig_list,
  samples = my_samples,
  ID = NULL, # Using the extracted IDs for secondary grouping
  threshold = 0.85
)


message("\nBCR data loading and combining complete!")
message("Your combined BCR data is in 'combined.BCR'.")


 

# 
#clonal call options aa,nt,strict, gene
clonalQuant(combined.TCR, 
            cloneCall="strict", 
            chain = "both", 
            scale = TRUE)
clonalQuant(combined.BCR, 
            cloneCall="aa", 
            chain = "both", 
            scale = TRUE)
clonalAbundance(combined.TCR, 
                cloneCall = "gene", 
                scale = FALSE)
clonalAbundance(combined.BCR, 
                cloneCall = "gene", 
                scale = TRUE)
clonalLength(combined.TCR, 
             cloneCall="aa", 
             chain = "both") 
clonalLength(combined.TCR, 
             cloneCall="aa", 
             chain = "TRA", 
             scale = TRUE) 
clonalScatter(combined.TCR, 
              cloneCall ="gene", 
              x.axis = "C16_R", 
              y.axis = "U4_R",
              dot.size = "total",
              graph = "proportion")
clonalHomeostasis(combined.TCR, 
                  cloneCall = "gene")
clonalHomeostasis(combined.TCR, 
                  cloneCall = "gene",
                  cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1, Large = 0.3, Hyperexpanded =
                                  1))
clonalProportion(combined.TCR, 
                 cloneCall = "gene") 
clonalProportion(combined.TCR, 
                 cloneCall = "nt",
                 clonalSplit = c(1, 5, 10, 100, 1000, 10000)) 
vizGenes(combined.BCR,
         x.axis = "IGH", 
         y.axis = NULL,
         plot = "barplot",
         scale = TRUE)
vizGenes(combined.TCR, 
         x.axis = "TRBV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)
vizGenes(combined.TCR[c(4,5,6,7,8)], # uc samples
         x.axis = "TRBV",
         y.axis = "TRBJ",
         plot = "heatmap",  
         scale = TRUE)
vizGenes(combined.TCR[c(29,30,31,32,33)], # normal 
         x.axis = "TRBV",
         y.axis = "TRBJ",
         plot = "heatmap",  
         scale = TRUE)
percentAA(combined.TCR, 
          chain = "TRB", 
          aa.length = 20)

positionalEntropy(combined.TCR, 
                  chain = "TRB", 
                  aa.length = 20)
percentGenes(combined.TCR, 
             chain = "TRB", 
             gene = "Vgene")
#
df.genes <- percentGenes(combined.TCR, 
                         chain = "TRB", 
                         gene = "Vgene", 
                         exportTable = TRUE)

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes)))
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) + 
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) + 
  theme_classic() 
percentVJ(combined.TCR, 
          chain = "TRB")
#
df.genes <- percentVJ(combined.TCR, 
                      chain = "TRB", 
                      exportTable = TRUE)

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes))) 
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) + 
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) + 
  theme_classic() 
percentKmer(combined.TCR, 
            cloneCall = "aa",
            chain = "TRB", 
            motif.length = 3, 
            top.motifs = 25)
percentKmer(combined.TCR, 
            cloneCall = "nt",
            chain = "TRB", 
            motif.length = 3, 
            top.motifs = 25)
clonalDiversity(combined.TCR, 
                cloneCall = "gene", 
                n.boots = 20)
clonalDiversity(combined.TCR, 
                metrics = c("shannon", "ACE"),
                cloneCall = "strict", 
                n.boots = 20)
clonalRarefaction(combined.BCR,
                  plot.type = 1,
                  hill.numbers = 1,
                  n.boots = 2)
#0 - species-richness
#1 - Shannon
#2 - Simpson
clonalSizeDistribution(combined.TCR, 
                       cloneCall = "aa", 
                       method= "ward.D2")
clonalOverlap(combined.TCR, 
              cloneCall = "strict", 
              method = "morisita")
# add to seurat 

# Remove sample prefix from TCR barcodes
for (i in seq_along(combined.TCR)) {
  combined.TCR[[i]]$barcode <- gsub("^[^_]+_", "", combined.TCR[[i]]$barcode)
}
seurat <- combineExpression(
  combined.TCR,
  seurat,
  cloneCall = "gene",           # or "aa", "nt", or "strict" depending on your goal
  group.by = "sample",          # or another metadata column in seurat@meta.data
  proportion = TRUE             # adds clonal proportion info to metadata
)
# Clean barcodes
# Remove prefixes like C9_R_ from barcodes
for (i in seq_along(combined.BCR)) {
  combined.BCR[[i]]$barcode <- sub("^[^_]+_", "", combined.BCR[[i]]$barcode)
}
seurat <- combineExpression(
  combined.BCR,
  seurat,
  cloneCall = "gene",
  group.by = "combined_sample_temp",
  proportion = TRUE
)
# won't work
for (i in seq_along(combined.BCR)) {
  df <- combined.BCR[[i]]
  
  # Combine heavy and light chain genes (adjust this if different columns are used)
  df$gene <- paste(df$IGH, df$IGLC, sep = "_")
  
  # Add sample label
  df$combined_sample_temp <- names(combined.BCR)[i]
  
  # Update the list
  combined.BCR[[i]] <- df
}
seurat <- combineExpression(
  combined.BCR,
  seurat,
  cloneCall = "gene",
  group.by = "combined_sample_temp",
  proportion = TRUE
)
library(dplyr)

# Extract clone vectors from each sample
clone_lists <- lapply(combined.TCR, function(x) x$CTgene)

# Make a data.frame with clone and sample info
clone_sample_df <- do.call(rbind, lapply(names(clone_lists), function(sample) {
  data.frame(sample = sample, cloneId = clone_lists[[sample]], stringsAsFactors = FALSE)
}))

library(dplyr)

clone_sharing <- clone_sample_df %>%
  group_by(cloneId) %>%
  summarise(n_samples = n_distinct(sample)) %>%
  ungroup()
clone_sharing <- clone_sharing %>%
  mutate(clone_type = ifelse(n_samples > 1, "Public", "Private"))
clone_sample_df <- clone_sample_df %>%
  left_join(clone_sharing %>% dplyr::select(cloneId, clone_type), by = "cloneId")
table_by_sample <- clone_sample_df %>%
  group_by(sample, clone_type) %>%
  summarise(clone_count = n()) %>%
  tidyr::pivot_wider(names_from = clone_type, values_from = clone_count, values_fill = 0)
library(ggplot2)

# Example plot of clone sharing counts
ggplot(clone_sharing, aes(x = n_samples)) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "Distribution of Clone Sharing Across Samples",
    x = "Number of Samples Clone is Found In",
    y = "Count of Clones"
  ) +
  theme_minimal()
clone_sharing$clone_type <- ifelse(clone_sharing$n_samples > 1, "Public", "Private")

ggplot(clone_sharing, aes(x = clone_type, fill = clone_type)) +
  geom_bar() +
  labs(title = "Count of Public vs Private Clones") +
  theme_minimal()




# Add to Seurat
seurat <- combineExpression(
  combined.BCR,
  seurat,
  cloneCall = "gene",       # or "aa", "nt", or "strict"
  group.by = "sample",      # must match your metadata column
  proportion = TRUE
)



#
# Load required package
library(dplyr)

# Define cloneCall: can be "gene", "aa", "nt", or "strict"
clone_call <- "gene"

# Step 1: Extract clone IDs from each sample in combined.TCR
clone_lists <- lapply(combined.TCR, function(df) df[[paste0("CT", clone_call)]])
names(clone_lists) <- names(combined.TCR)

# Step 2: Create a long-format data frame (sample, cloneId)
clone_sample_df <- do.call(rbind, lapply(names(clone_lists), function(sample) {
  data.frame(
    sample = sample,
    cloneId = clone_lists[[sample]],
    stringsAsFactors = FALSE
  )
}))

# Step 3: Add patient ID (e.g., from sample name like "C9_R" → "C9")
clone_sample_df$patient <- sub("_.*", "", clone_sample_df$sample)

# Step 4: Count how many patients each clone appears in (publicity)
publicity_df <- clone_sample_df %>%
  distinct(cloneId, patient) %>%
  group_by(cloneId) %>%
  summarise(n_patients = n_distinct(patient), .groups = "drop") %>%
  mutate(clone_type = ifelse(n_patients > 1, "Public", "Private"))

# Step 5: Add clone type back to the full data
clone_sample_df <- clone_sample_df %>%
  left_join(publicity_df, by = "cloneId")

# Step 6: View summary
table(clone_sample_df$clone_type)

# Optional: Get public clones in a specific patient, e.g., C9
public_C9 <- clone_sample_df %>%
  filter(patient == "C9", clone_type == "Public") %>%
  distinct(cloneId)

# Print result
print(public_C9)
# Filter combined.TCR to keep only public C9 clones
library(dplyr)

public_clone_ids <- public_C9$cloneId

# Get all C9 rows from combined.TCR that match public clones
public_c9_data <- combined.TCR[["C9_R"]] %>%
  filter(CTgene %in% public_clone_ids)

# View common V-J pairings
table(public_c9_data$TCR1, public_c9_data$TCR2)
library(ggplot2)


write.csv(public_C9, "Public_Clones_C9.csv", row.names = FALSE)

# Build clone_sample_df from combined.TCR
clone_sample_df <- do.call(rbind, lapply(names(combined.TCR), function(sample) {
  df <- combined.TCR[[sample]]
  data.frame(patient = sample, cloneId = df$CTgene, stringsAsFactors = FALSE)
}))
library(dplyr)

# Identify public clones (appear in >1 patient)
publicity_df <- clone_sample_df %>%
  group_by(cloneId) %>%
  summarise(n_patients = n_distinct(patient)) %>%
  mutate(clone_type = ifelse(n_patients > 1, "public", "private")) %>%
  ungroup()
# Filter for patient C9
c9_clone_info <- clone_sample_df %>%
  filter(patient == "C9_R") %>%
  left_join(publicity_df, by = "cloneId")
library(ggplot2)

ggplot(c9_clone_info, aes(x = clone_type)) +
  geom_bar(fill = "steelblue") +
  labs(title = "Public vs Private Clones in C9", x = "Clone Type", y = "Count") +
  theme_minimal()



















