# load the data 
files <- list.files(pattern = "*TCR_table.tsv", full.names = TRUE)
# Load and add chain for TCR 
contig.frames <- lapply(files, function(file) {
  df <- read.delim(file)
  df$chain <- ifelse(grepl("TRBC", df$c_gene), "TRB",
                     ifelse(grepl("TRAC", df$c_gene), "TRA", NA))
  df
})
# Smarter sample name extraction
sample_names <- sub("^[^_]+_([^_]+_[^_]+).*", "\\1", basename(files))
# Combine
combined.TCR <- combineTCR(contig.frames, samples = sample_names, ID = "TCR")
