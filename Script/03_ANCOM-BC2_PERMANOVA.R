#Author: Rio Kashinmoto
#Last edit: April 6th, 2026
#Project: MVP Coral Microbiome Analysis (Porites & Pocillopora)
#Description: PERMANOVA and ANCOM-BC2 differential abundance testing for Pocillopora (POC) and Porites lobota (POR) across sites for 2017 and 2019.

### ANCOM-BC2 POC only 
# ---- 1. Setup and Load Data ----
library(ANCOMBC)
library(tidyverse)
library(phyloseq)
library(vegan)

# Load the cleaned phyloseq object
ps <- readRDS("Data/ps_POC_POR_filtered_clean.rds")

# Keep only Pocillopora (POC) samples to avoid mixed-host artifacts
ps_poc <- prune_samples(grepl("POC", sample_names(ps)), ps)
ps_poc <- prune_taxa(taxa_sums(ps_poc) > 0, ps_poc)


# ==============================================================================
# ---- PART A: 2017 ANALYSIS (Baseline: Site 0) ----
# ==============================================================================

# ---- 2. Subset & Format 2017 Data ----
ps_2017 <- subset_samples(ps_poc, Year == "2017")
ps_2017 <- prune_taxa(taxa_sums(ps_2017) > 0, ps_2017)

# Parse Site from sample names (e.g., "POCFRG2017_0_1" -> "0")
sample_data(ps_2017)$Site <- as.factor(str_extract(sample_names(ps_2017), "(?<=_)[0-9\\.]+(?=_)"))

# Set Baseline/Reference site "0" for 2017
# Set Baseline/Reference site "1, 3.5, 5 and 5.5" for 2017 separetely
sample_data(ps_2017)$Site <- relevel(sample_data(ps_2017)$Site, ref = "0")
#sample_data(ps_2017)$Site <- relevel(sample_data(ps_2017)$Site, ref = "1")
#sample_data(ps_2017)$Site <- relevel(sample_data(ps_2017)$Site, ref = "3.5")
#sample_data(ps_2017)$Site <- relevel(sample_data(ps_2017)$Site, ref = "5")
#sample_data(ps_2017)$Site <- relevel(sample_data(ps_2017)$Site, ref = "5.5")

cat("✅ 2017 Data Prepared. Baseline Site:", levels(sample_data(ps_2017)$Site)[1], "\n")

# ---- 3. Run PERMANOVA (2017) ----
set.seed(123)
# Calculate Bray-Curtis distance matrix natively via phyloseq
bray_2017 <- distance(ps_2017, method = "bray")
meta_2017 <- as(sample_data(ps_2017), "data.frame")

# Run adonis2
adon_2017 <- adonis2(bray_2017 ~ Site, data = meta_2017, permutations = 999)
print("PERMANOVA Results for 2017:")
print(adon_2017)

# ---- 4. Run ANCOM-BC2 (2017) ----
out_2017 <- ancombc2(
  data = ps_2017,
  tax_level = "Family",       
  fix_formula = "Site",       
  p_adj_method = "holm",
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE
)

# Extract and save FULL results for POC site 0 vs different site in 2017
# Extract and save FULL results for POC  site "1, 3.5, 5 and 5.5" by changing the cvs file name.
res_2017 <- out_2017$res
write.csv(res_2017, "Data/POC_refsite0_2017_ANCOMBC2_Full.csv", row.names = FALSE)
