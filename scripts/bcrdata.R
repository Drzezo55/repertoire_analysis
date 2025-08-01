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
  ID = my_ids, # Using the extracted IDs for secondary grouping
  threshold = 0.85
)


message("\nBCR data loading and combining complete!")
message("Your combined BCR data is in 'combined.BCR'.")