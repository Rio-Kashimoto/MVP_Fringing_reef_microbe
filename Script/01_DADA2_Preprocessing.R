#Author: Rio Kashinmoto
#Last edit: April 2nd, 2026
#Project: MVP Coral Microbiome Analysis (Porites & Pocillopora)
#Description: DADA2 pipeline from filtered reads to initial Phyloseq object.

# ---- 1. Load Required Packages ----
library(dada2)
library(phyloseq)
library(Biostrings)
library(tidyverse)

# ---- 2. Define Paths & Sample Names ----
# Set path to the cutadapt filtered fastq files
filteredReadsPath <- "path"
taxaDBPath <- "path/silva_nr99_v138.1_train_set.fa.gz"

# Get list of filtered reads
forwardReads_filt <- sort(list.files(filteredReadsPath, pattern=".fastq.gz", full.names = TRUE))

# Extract sample names and ensure they are unique
sample.names <- sapply(strsplit(basename(forwardReads_filt), ".fastq.gz"), `[`, 1)
sample.names <- make.unique(sample.names)

cat("Found", length(sample.names), "samples in the filtered directory.\n")

# ---- 3. Initial Read Tracking Function ----
# Function to count reads in FASTQ files
countReads <- function(file) {
  if (file.exists(file)) {
    total_lines <- length(readLines(file))
    return(total_lines / 4) # Each read is 4 lines
  } else {
    return(NA)
  }
}

# Calculate input read counts for tracking
filteredFastq <- sapply(forwardReads_filt, countReads)

# ---- 4. DADA2 Core Pipeline ----
# Learn error rates
set.seed(123) # Set seed for reproducibility
error_forward <- learnErrors(forwardReads_filt, multithread = TRUE)

# Perform DADA on forward reads
dadaFs <- dada(forwardReads_filt, err = error_forward, multithread = TRUE)

# Create sequence table
seqtab <- makeSequenceTable(dadaFs)

# Remove chimeric sequences
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
num_chim_removed <- 1 - (sum(seqtab_nochim) / sum(seqtab))
cat("Proportion of chimeric sequences removed:", round(num_chim_removed, 4), "\n")

# Assign unique names to rows
rownames(seqtab_nochim) <- sample.names

# ---- 5. Compile Tracking Table ----
# Build the tracking table for all initial samples
getN <- function(x) sum(getUniques(x))

track <- data.frame(
  input = filteredFastq,                      
  filtered = filteredFastq,                   
  denoisedF = sapply(dadaFs, getN),           
  nonchim = rowSums(seqtab_nochim)            
)
rownames(track) <- sample.names

# Save the tracking table to check for sample failures
write.csv(track, "Data/DADA2_read_tracking.csv")

# ---- 6. Filter Low-Read Samples ----
# Remove samples with fewer than 900 reads to prevent artifacts in diversity analysis
read_threshold <- 1000
keep_samples <- rowSums(seqtab_nochim) >= read_threshold

seqtab_filtered <- seqtab_nochim[keep_samples, ]

cat("Samples kept after filtering (<", read_threshold, "reads):", nrow(seqtab_filtered), "out of", nrow(seqtab_nochim), "\n")

# ---- 7. Assign Taxonomy ----
# Assign taxonomy ONLY to the high-quality samples to save computation time
taxa <- assignTaxonomy(seqtab_filtered, taxaDBPath, multithread = TRUE)

# ---- 8. Build Metadata ----
# Extract remaining valid sample names
final_sample_names <- rownames(seqtab_filtered)

# Generate metadata based on sample names
metadata <- tibble(Sample_names = final_sample_names) %>%
  mutate(
    Year = case_when(
      str_detect(Sample_names, "2017") ~ "2017",
      str_detect(Sample_names, "2018") ~ "2018",
      str_detect(Sample_names, "2019") ~ "2019",
      TRUE ~ NA_character_
    ),
    Host = case_when(
      str_detect(Sample_names, "POC") ~ "Pocillopora spp",
      str_detect(Sample_names, "POR") ~ "Porites lobota",
      TRUE ~ NA_character_
    ),
    Site = str_extract(Sample_names, "(?<=_)[0-9\\.]+(?=_)")
  ) %>%
  column_to_rownames("Sample_names")

# ---- 9. Construct and Save Phyloseq Object ----
# Build the phyloseq object
ps <- phyloseq(otu_table(seqtab_filtered, taxa_are_rows = FALSE), 
               sample_data(metadata), 
               tax_table(taxa))

# Extract DNA sequences into a Biostrings object and rename taxa to ASV1, ASV2...
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Final sanity check
print(ps)

# Save the finalized, clean phyloseq object
saveRDS(ps, "Data/ps_POC_POR_filtered.rds")
cat("Pipeline complete. Filtered Phyloseq object saved to 'Data/ps_POC_POR_filtered.rds'\n")
