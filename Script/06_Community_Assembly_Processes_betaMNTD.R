#Author: Rio Kashinmoto
#Last edit: April 5th, 2026
#Project: MVP Coral Microbiome Analysis (Porites & Pocillopora)
#Description: Calculates betaMNTD and betaNTI to determine community assembly processes (Deterministic vs. Stochastic) for both Pocillopora (POC) and Porites (POR) samples.

# ---- 1. Load Packages ----
set.seed(123)
library(phyloseq)
library(dplyr)
library(reshape2)
library(ggplot2)
library(picante) 

# ---- 2. Data Preparation & Cleaning ----
# Load previously merged phyloseq object (must contain phylogenetic tree)
ps <- readRDS("Data/ps_POC_POR_filtered_clean.rds")

# Remove samples with 0 total reads (Causes NaN errors)
ps <- prune_samples(sample_sums(ps) > 0, ps)

# Remove taxa with 0 total abundance
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Prune Tree to match OTU table exactly (Prevents 'out of bounds' error)
target_taxa <- taxa_names(ps)
target_tree <- phy_tree(ps)
pruned_tree <- keep.tip(target_tree, target_taxa)

# Extract cleaned OTU table (Samples as rows)
otu_tab <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) { otu_tab <- t(otu_tab) }

# ---- 3. Calculate betaMNTD & betaNTI ----
# Generate the distance matrix from the tree once (saves time)
phlyo_dist <- cophenetic(pruned_tree)

print("Calculating Observed betaMNTD...")
betaMNTD_obs <- as.matrix(comdistnt(otu_tab, phlyo_dist, abundance.weighted = TRUE))

# Calculate Null Model (Randomization)
# Using 999 iterations for publication-level rigor
set.seed(123)
null_list <- list()
print("Starting Null Model Randomizations (n=999)...")

for(i in 1:999){
  # Shuffle the labels on the tree tips (null model)
  rand_dist <- phlyo_dist
  rand_names <- sample(colnames(rand_dist))
  rownames(rand_dist) <- rand_names
  colnames(rand_dist) <- rand_names
  
  null_list[[i]] <- as.matrix(comdistnt(otu_tab, rand_dist, abundance.weighted = TRUE))
  
  if(i %% 50 == 0) cat("Iteration:", i, "\n") # Progress tracker
}

# Calculate betaNTI
null_array <- array(unlist(null_list), dim = c(nrow(betaMNTD_obs), ncol(betaMNTD_obs), length(null_list)))
null_mean  <- apply(null_array, c(1,2), mean)
null_sd    <- apply(null_array, c(1,2), sd)
betaNTI    <- (betaMNTD_obs - null_mean) / null_sd


# ---- 4. Prepare Metadata & Mapping ----
meta <- data.frame(sample_data(ps))
meta$SampleID <- rownames(meta)

# Extract Species, Year, and Region
meta$Species <- ifelse(grepl("^POC", meta$SampleID), "POC", "POR")
meta$Year <- as.factor(meta$Year)
meta$Site <- sub("[A-Z]+[0-9]{4}_([0-9.]+)_.*", "\\1", meta$SampleID)

meta$Region <- case_when(
  meta$Site %in% c("0","1","2") ~ "North",
  meta$Site %in% c("3","3.5","4") ~ "East",
  meta$Site %in% c("5","5.5","6") ~ "West",
  TRUE ~ NA_character_
)

# Convert betaNTI matrix to long format
betaNTI_long <- melt(betaNTI, varnames = c("Sample1", "Sample2"), value.name = "bNTI")

# Remove self-comparisons
betaNTI_long <- betaNTI_long[betaNTI_long$Sample1 != betaNTI_long$Sample2, ]

# Map metadata to the pairwise comparisons
betaNTI_long$Species1 <- meta$Species[match(betaNTI_long$Sample1, meta$SampleID)]
betaNTI_long$Species2 <- meta$Species[match(betaNTI_long$Sample2, meta$SampleID)]
betaNTI_long$Year1    <- meta$Year[match(betaNTI_long$Sample1, meta$SampleID)]
betaNTI_long$Year2    <- meta$Year[match(betaNTI_long$Sample2, meta$SampleID)]
betaNTI_long$Region1  <- meta$Region[match(betaNTI_long$Sample1, meta$SampleID)]
betaNTI_long$Region2  <- meta$Region[match(betaNTI_long$Sample2, meta$SampleID)]

# Classify Process based on thresholds
# |betaNTI| > 2 is Deterministic (Selection), |betaNTI| < 2 is Stochastic
betaNTI_long$Process <- cut(betaNTI_long$bNTI,
                            breaks = c(-Inf, -2, 2, Inf),
                            labels = c("Homogeneous selection", "Stochastic", "Variable selection"))


# ---- 5. Visualizations (Pie Charts) ----

# Global Plot Settings
process_colors <- c(
  "Homogeneous selection" = "#1f78b4", # Blue
  "Stochastic"            = "#cccccc", # Grey
  "Variable selection"    = "#e31a1c"  # Red
)

pie_theme <- theme_void() + 
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

# -----------------------------------------------------
# Plot 1: Species Assembly Processes by Year
# -----------------------------------------------------
summary_sp_year <- betaNTI_long %>%
  filter(Species1 == Species2 & Year1 == Year2) %>%
  group_by(Year1, Species1, Process) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Year1, Species1) %>%
  mutate(percent = n / sum(n) * 100,
         label = paste0(round(percent, 1), "%"))

p1 <- ggplot(summary_sp_year, aes(x = "", y = percent, fill = Process)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  facet_grid(Year1 ~ Species1) +
  geom_text(aes(label = ifelse(percent > 5, label, "")), 
            position = position_stack(vjust = 0.5), size = 3, fontface = "bold") +
  scale_fill_manual(values = process_colors) +
  pie_theme +
  labs(title = "Species Assembly Processes by Year")
print(p1)

# -----------------------------------------------------
# Plot 2: Regional Assembly Processes by Year (All Species)
# -----------------------------------------------------
summary_reg_year <- betaNTI_long %>%
  filter(Region1 == Region2 & Year1 == Year2 & !is.na(Region1)) %>%
  group_by(Year1, Region1, Process) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Year1, Region1) %>%
  mutate(percent = n / sum(n) * 100,
         label = paste0(round(percent, 1), "%"))

summary_reg_year$Region1 <- factor(summary_reg_year$Region1, levels = c("North", "East", "West"))

p2 <- ggplot(summary_reg_year, aes(x = "", y = percent, fill = Process)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  facet_grid(Year1 ~ Region1) +
  geom_text(aes(label = ifelse(percent > 5, label, "")), 
            position = position_stack(vjust = 0.5), size = 3, fontface = "bold") +
  scale_fill_manual(values = process_colors) +
  pie_theme +
  labs(title = "Regional Assembly Processes by Year")
print(p2)

# -----------------------------------------------------
# Plot 3: Regional Assembly Processes (POC vs POR)
# -----------------------------------------------------
summary_reg_sp <- betaNTI_long %>%
  filter(Region1 == Region2 & Species1 == Species2 & !is.na(Region1)) %>%
  group_by(Region1, Species1, Process) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Region1, Species1) %>%
  mutate(percent = n / sum(n) * 100,
         label = paste0(round(percent, 1), "%"))

summary_reg_sp$Region1 <- factor(summary_reg_sp$Region1, levels = c("North", "East", "West"))
summary_reg_sp$Species1 <- factor(summary_reg_sp$Species1, levels = c("POC", "POR"))

p3 <- ggplot(summary_reg_sp, aes(x = "", y = percent, fill = Process)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  facet_grid(Species1 ~ Region1) +
  geom_text(aes(label = ifelse(percent > 5, label, "")), 
            position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = process_colors) +
  pie_theme + 
  theme(legend.position = "bottom", plot.title = element_text(margin = margin(b=20))) +
  labs(title = "Regional Assembly Processes: POC vs. POR", fill = "Assembly Process")
print(p3)

# -----------------------------------------------------
# Plot 4: POC & POR Separated by Region and Year
# -----------------------------------------------------
summary_final <- betaNTI_long %>%
  filter(Region1 == Region2 & Year1 == Year2 & Species1 == Species2 & !is.na(Region1)) %>%
  group_by(Species1, Year1, Region1, Process) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Species1, Year1, Region1) %>%
  mutate(percent = n / sum(n) * 100,
         label = paste0(round(percent, 1), "%"))

summary_final$Region1 <- factor(summary_final$Region1, levels = c("North", "East", "West"))

# POC Specific Plot
p_poc <- ggplot(subset(summary_final, Species1 == "POC"), aes(x = "", y = percent, fill = Process)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  facet_grid(Year1 ~ Region1) +
  geom_text(aes(label = ifelse(percent > 5, label, "")), position = position_stack(vjust = 0.5), size = 3, fontface = "bold") +
  scale_fill_manual(values = process_colors) + pie_theme +
  labs(title = "Pocillopora (POC): Regional Assembly by Year", fill = "Process")
print(p_poc)

# POR Specific Plot
p_por <- ggplot(subset(summary_final, Species1 == "POR"), aes(x = "", y = percent, fill = Process)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  facet_grid(Year1 ~ Region1) +
  geom_text(aes(label = ifelse(percent > 5, label, "")), position = position_stack(vjust = 0.5), size = 3, fontface = "bold") +
  scale_fill_manual(values = process_colors) + pie_theme +
  labs(title = "Porites (POR): Regional Assembly by Year", fill = "Process")
print(p_por)

# -----------------------------------------------------
# Plot 5: Overall Summaries (Year, Species, Region)
# -----------------------------------------------------
# By Year
summary_year <- betaNTI_long %>% filter(Year1 == Year2) %>% group_by(Year1, Process) %>% summarise(n = n(), .groups = "drop") %>% group_by(Year1) %>% mutate(percent = n / sum(n) * 100, label = paste0(round(percent, 1), "%"))
p_year_all <- ggplot(summary_year, aes(x = "", y = percent, fill = Process)) + geom_bar(stat = "identity", width = 1, color = "white") + coord_polar("y", start = 0) + facet_wrap(~Year1) + geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4, fontface = "bold") + scale_fill_manual(values = process_colors) + pie_theme + labs(title = "Overall Assembly Processes by Year")
print(p_year_all)

# By Species
summary_species_all <- betaNTI_long %>% filter(Species1 == Species2) %>% group_by(Species1, Process) %>% summarise(n = n(), .groups = "drop") %>% group_by(Species1) %>% mutate(percent = n / sum(n) * 100, label = paste0(round(percent, 1), "%"))
p_species_all <- ggplot(summary_species_all, aes(x = "", y = percent, fill = Process)) + geom_bar(stat = "identity", width = 1, color = "white") + coord_polar("y", start = 0) + facet_wrap(~Species1) + geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 5, fontface = "bold") + scale_fill_manual(values = process_colors) + pie_theme + labs(title = "Overall Species Assembly Comparison")
print(p_species_all)

# By Region
summary_region_all <- betaNTI_long %>% filter(Region1 == Region2 & !is.na(Region1)) %>% group_by(Region1, Process) %>% summarise(n = n(), .groups = "drop") %>% group_by(Region1) %>% mutate(percent = n / sum(n) * 100, label = paste0(round(percent, 1), "%"))
summary_region_all$Region1 <- factor(summary_region_all$Region1, levels = c("North", "East", "West"))
p_region_all <- ggplot(summary_region_all, aes(x = "", y = percent, fill = Process)) + geom_bar(stat = "identity", width = 1, color = "white") + coord_polar("y", start = 0) + facet_wrap(~Region1) + geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4, fontface = "bold") + scale_fill_manual(values = process_colors) + pie_theme + labs(title = "Overall Regional Assembly Comparison")
print(p_region_all)
