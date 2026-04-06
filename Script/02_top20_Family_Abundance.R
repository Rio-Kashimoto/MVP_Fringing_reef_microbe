#Author: Rio Kashinmoto
#Last edit: April 5th, 2026
#Project: MVP Coral Microbiome Analysis (Porites & Pocillopora)
#Description: Generates top 20 mean relative abundance stacked bar plots for Pocillopora (POC) and Porites lobota (POR) samples across different years and sites.

###POC Mean relative abundance of top 20
# ---- 1. Load Packages ----
library(phyloseq)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(forcats)

# ---- 2. Load and Filter Data ----
# Use relative paths for GitHub reproducibility
ps <- readRDS("Data/ps_POC_POR_filtered_clean.rds")

# KEEP ONLY POC SAMPLES
ps_poc <- prune_samples(grepl("POC", sample_names(ps)), ps)

# Remove any taxa that are now empty (0 reads) after removing POR samples
ps_poc <- prune_taxa(taxa_sums(ps_poc) > 0, ps_poc)

cat("Total POC samples:", nsamples(ps_poc), "\n")

# ---- 3. Format Metadata ----
metadata <- data.frame(Sample = sample_names(ps_poc)) %>%
  mutate(
    Year = str_extract(Sample, "20[0-9]{2}"),
    Site = str_extract(Sample, "(?<=_)\\d+\\.?\\d*(?=_)")
  ) %>%
  tibble::column_to_rownames("Sample")

sample_data(ps_poc) <- sample_data(metadata)

# ---- 4. Aggregate & Calculate Relative Abundance ----
# Group by family
ps_family <- tax_glom(ps_poc, "Family", NArm = TRUE)

# Convert to relative abundance (%)
ps_rel <- transform_sample_counts(
  ps_family,
  function(x) x / sum(x) * 100
)

# Convert to dataframe
df <- psmelt(ps_rel)

# ---- 5. Calculate Mean Abundance per Site ----
mean_df <- df %>%
  group_by(Year, Site, Family) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

# Remove sites with 0 abundance (empty groups)
mean_df <- mean_df %>%
  group_by(Year, Site) %>%
  filter(sum(Abundance) > 0) %>%
  ungroup()

# ---- 6. Select Top 20 Families ----
top20 <- mean_df %>%
  group_by(Family) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice(1:20) %>%
  pull(Family)

# ---- 7. Group the Rest as "Other" ----
mean_df <- mean_df %>%
  mutate(
    Family = ifelse(Family %in% top20, as.character(Family), "Other")
  ) %>%
  group_by(Year, Site, Family) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# ---- 8. Normalize to 100% ----
mean_df <- mean_df %>%
  group_by(Year, Site) %>%
  mutate(Abundance = Abundance / sum(Abundance) * 100) %>%
  ungroup()

# ---- 9. Format Site Labels ----
mean_df <- mean_df %>%
  mutate(
    SiteLabel = case_when(
      Site %in% c("0","1","2") ~ paste0(Site, "_N"),
      Site %in% c("3","3.5","4") ~ paste0(Site, "_E"),
      Site %in% c("5","5.5","6") ~ paste0(Site, "_W"),
      TRUE ~ Site
    )
  )

mean_df$SiteLabel <- factor(
  mean_df$SiteLabel,
  levels = c("0_N","1_N","2_N","3_E","3.5_E","4_E","5_W","5.5_W","6_W")
)

# ---- 10. Order Families (Top at Top, 'Other' at Bottom) ----
mean_df <- mean_df %>%
  group_by(Family) %>%
  mutate(TotalAbundance = sum(Abundance)) %>%
  ungroup() %>%
  mutate(
    # Order by abundance first
    Family = fct_reorder(Family, TotalAbundance, .desc = TRUE),
    # Force "Other" to be the very last factor level so it plots at the bottom/top
    Family = fct_relevel(Family, "Other", after = Inf) 
  )

# ---- 11. Assign Custom Colors ----
custom_colors <- c(
  "Endozoicomonadaceae" = "forestgreen",
  "Xenococcaceae"       = "cyan",
  "Rhodobacteraceae"    = "deeppink",
  "Methanosarcinaceae"  = "yellow",
  "Halieaceae"          = "blue",
  "Cyanobiaceae"        = "seashell4",
  "Saprospiraceae"      = "orchid1",
  "Cryomorphaceae"      = "chartreuse",
  "Alteromonadaceae"    = "darkolivegreen1",
  "Cyclobacteriaceae"   = "#522477",
  "Flavobacteriaceae"   = "#4682B4",
  "Vibrionaceae"        = "#E4003A",
  "Hyphomonadaceae"     = "#C3A580",
  "Kiloniellaceae"      = "darkorchid1",
  "Phormidesmiaceae"    = "darkgoldenrod1",
  "Sphingomonadaceae"   = "cornsilk",
  "Nostocaceae"         = "lavenderblush1",
  "Sandaracinaceae"     = "grey80",
  "SAR116 clade"        = "sienna3",
  "Other"               = "black" 
)

# Find any top 20 families not in the custom list
other_families <- setdiff(unique(mean_df$Family), names(custom_colors))

# Generate palette for missing families
other_colors <- brewer.pal(
  n = max(3, min(length(other_families), 12)),
  "Paired"
)

# Expand palette if needed
if (length(other_families) > length(other_colors)) {
  other_colors <- rep(other_colors, length.out = length(other_families))
}

other_colors <- setNames(other_colors, other_families)

# Combine custom and generated colors
family_colors <- c(custom_colors, other_colors)

# ---- 12. Plot Mean Abundance ----
p <- ggplot(mean_df, aes(x = SiteLabel, y = Abundance, fill = Family)) +
  geom_bar(
    stat = "identity",
    width = 1,
    color = "gray30",
    linewidth = 0.25 # Updated from 'size' to prevent ggplot warnings
  ) +
  scale_fill_manual(values = family_colors) +
  facet_grid(~Year, scales = "free_x", space = "free_x") +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    title = "Mean Family Abundance (Pocillopora)",
    x = "Site",
    y = "Relative Abundance (%)"
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13.5, angle = 90, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    panel.spacing.x = unit(0.5, "lines"),
    panel.grid.major.x = element_blank()
  )

print(p)


###POR Mean relative abundance of top 20
# KEEP ONLY POR SAMPLES
ps_por <- prune_samples(grepl("POR", sample_names(ps)), ps)

# EXCLUDE SPECIFIC SAMPLES (Low replication, n < 3)
samples_to_exclude <- c("PORFRG2019_2_64", "PORFRG2019_2_65")
ps_por <- prune_samples(!(sample_names(ps_por) %in% samples_to_exclude), ps_por)

# Remove any taxa that are now empty (0 reads) after filtering
ps_por <- prune_taxa(taxa_sums(ps_por) > 0, ps_por)

cat("Total POR samples after filtering:", nsamples(ps_por), "\n")

# ---- 3. Format Metadata ----
metadata <- data.frame(Sample = sample_names(ps_por)) %>%
  mutate(
    Year = str_extract(Sample, "20[0-9]{2}"),
    Site = str_extract(Sample, "(?<=_)\\d+\\.?\\d*(?=_)")
  ) %>%
  tibble::column_to_rownames("Sample")

sample_data(ps_por) <- sample_data(metadata)

# ---- 4. Aggregate & Calculate Relative Abundance ----
# Group by family
ps_family <- tax_glom(ps_por, "Family", NArm = TRUE)

# Convert to relative abundance (%)
ps_rel <- transform_sample_counts(
  ps_family,
  function(x) x / sum(x) * 100
)

# Convert to dataframe
df <- psmelt(ps_rel)

# ---- 5. Calculate Mean Abundance per Site ----
mean_df <- df %>%
  group_by(Year, Site, Family) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

# Remove sites with 0 abundance (empty groups)
mean_df <- mean_df %>%
  group_by(Year, Site) %>%
  filter(sum(Abundance) > 0) %>%
  ungroup()

# ---- 6. Select Top 20 Families ----
top20 <- mean_df %>%
  group_by(Family) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice(1:20) %>%
  pull(Family)

# ---- 7. Group the Rest as "Other" ----
mean_df <- mean_df %>%
  mutate(
    Family = ifelse(Family %in% top20, as.character(Family), "Other")
  ) %>%
  group_by(Year, Site, Family) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# ---- 8. Normalize to 100% ----
mean_df <- mean_df %>%
  group_by(Year, Site) %>%
  mutate(Abundance = Abundance / sum(Abundance) * 100) %>%
  ungroup()

# ---- 9. Format Site Labels ----
mean_df <- mean_df %>%
  mutate(
    SiteLabel = case_when(
      Site %in% c("0","1","2") ~ paste0(Site, "_N"),
      Site %in% c("3","3.5","4") ~ paste0(Site, "_E"),
      Site %in% c("5","5.5","6") ~ paste0(Site, "_W"),
      TRUE ~ Site
    )
  )

mean_df$SiteLabel <- factor(
  mean_df$SiteLabel,
  levels = c("0_N","1_N","2_N","3_E","3.5_E","4_E","5_W","5.5_W","6_W")
)

# ---- 10. Order Families (Top at Top, 'Other' at Bottom) ----
mean_df <- mean_df %>%
  group_by(Family) %>%
  mutate(TotalAbundance = sum(Abundance)) %>%
  ungroup() %>%
  mutate(
    # Order by abundance first
    Family = fct_reorder(Family, TotalAbundance, .desc = TRUE),
    # Force "Other" to be the very last factor level so it plots at the bottom/top
    Family = fct_relevel(Family, "Other", after = Inf) 
  )

# ---- 11. Assign Custom Colors ----
custom_colors <- c(
  "Endozoicomonadaceae" = "forestgreen",
  "Xenococcaceae"       = "cyan",
  "Rhodobacteraceae"    = "deeppink",
  "Methanosarcinaceae"  = "yellow",
  "Halieaceae"          = "blue",
  "Cyanobiaceae"        = "seashell4",
  "Saprospiraceae"      = "orchid1",
  "Cryomorphaceae"      = "chartreuse",
  "Alteromonadaceae"    = "darkolivegreen1",
  "Cyclobacteriaceae"   = "#522477",
  "Flavobacteriaceae"   = "#4682B4",
  "Vibrionaceae"        = "#E4003A",
  "Hyphomonadaceae"     = "#C3A580",
  "Kiloniellaceae"      = "darkorchid1",
  "Phormidesmiaceae"    = "darkgoldenrod1",
  "Sphingomonadaceae"   = "cornsilk",
  "Nostocaceae"         = "lavenderblush1",
  "Sandaracinaceae"     = "grey80",
  "SAR116 clade"        = "sienna3",
  "Other"               = "black" 
)

# Find any top 20 families not in the custom list
other_families <- setdiff(unique(mean_df$Family), names(custom_colors))

# Generate palette for missing families
other_colors <- brewer.pal(
  n = max(3, min(length(other_families), 12)),
  "Paired"
)

# Expand palette if needed
if (length(other_families) > length(other_colors)) {
  other_colors <- rep(other_colors, length.out = length(other_families))
}

other_colors <- setNames(other_colors, other_families)

# Combine custom and generated colors
family_colors <- c(custom_colors, other_colors)

# ---- 12. Plot Mean Abundance ----
p <- ggplot(mean_df, aes(x = SiteLabel, y = Abundance, fill = Family)) +
  geom_bar(
    stat = "identity",
    width = 1,
    color = "gray30",
    linewidth = 0.25 
  ) +
  scale_fill_manual(values = family_colors) +
  facet_grid(~Year, scales = "free_x", space = "free_x") +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    title = "Mean Family Abundance (Porites lobata)",
    x = "Site",
    y = "Relative Abundance (%)"
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13.5, angle = 90, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    panel.spacing.x = unit(0.5, "lines"),
    panel.grid.major.x = element_blank()
  )

print(p)
