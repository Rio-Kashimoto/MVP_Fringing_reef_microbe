#Author: Rio Kashinmoto
#Last edit: April 6th, 2026
#Project: MVP Coral Microbiome Analysis (Porites & Pocillopora)
#Description:Calculates Alpha Diversity (Shannon, PD) for Pocillopora (POC) and Porites lobota (POR) runs statistical testing, and formats environmental data (Temperature & Nutrients) for plotting.

# ---- 1. Load Libraries ----
set.seed(123)
library(phyloseq)
library(picante)
library(microViz)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(data.table)
library(readr)
library(tidyr)

### alpha diversity for POC
# ---- 2. Load & Subset Data (POC Only) ----
ps <- readRDS("Data/ps_POC_POR_filtered_clean.rds")

# Keep ONLY Pocillopora samples and remove empty taxa
ps_poc <- prune_samples(grepl("POC", sample_names(ps)), ps)
ps_poc <- prune_taxa(taxa_sums(ps_poc) > 0, ps_poc)

# Fix taxonomy names formatting
ps_poc <- tax_fix(
  ps_poc,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
)

# ---- 3. Calculate Alpha Diversity (Genus Level) ----
# Agglomerate to Genus
ps_poc_genus <- tax_glom(ps_poc, taxrank = "Genus", NArm = FALSE)

# Calculate standard alpha diversity (Observed, InvSimpson, Shannon)
alphaDiv <- estimate_richness(ps_poc_genus, measures = c("Observed", "InvSimpson", "Shannon"))
alphaDiv$Samples <- rownames(alphaDiv)

# Build Phylogenetic Tree for Faith's PD
seqs <- as.character(refseq(ps_poc))
names(seqs) <- taxa_names(ps_poc)

alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), anchor = NA)
phang.align <- phangorn::phyDat(as.matrix(alignment), type = "DNA")
dm <- phangorn::dist.ml(phang.align)
treeNJ <- phangorn::NJ(dm)
fit <- phangorn::pml(treeNJ, data = phang.align)
fitGTR <- phangorn::optim.pml(fit, model="GTR", optInv=TRUE, optGamma=TRUE)

tree <- fitGTR$tree
tree$tip.label <- names(seqs)

# Merge tree and calculate Faith's PD
ps_poc_tree <- merge_phyloseq(ps_poc, phy_tree(tree))
ps_poc_genus_tree <- tax_glom(ps_poc_tree, taxrank = "Genus", NArm = FALSE)

otu <- as(otu_table(ps_poc_genus_tree), "matrix")
tree.genus <- phy_tree(ps_poc_genus_tree)

pd.res <- picante::pd(otu, tree.genus, include.root = FALSE)
pd.df <- as.data.frame(pd.res)
pd.df$Samples <- rownames(pd.df)

# Merge metadata with alpha diversity metrics
samdf <- data.frame(sample_data(ps_poc_genus), check.names = FALSE)
samdf$Samples <- rownames(samdf)

alphaDiv <- merge(samdf, alphaDiv, by = "Samples")

# Add Region assignments
alphaDiv <- alphaDiv %>%
  mutate(
    Site = sub("POCFRG[0-9]{4}_([0-9.]+)_.*", "\\1", Samples),
    Region = case_when(
      Site %in% c("0","1","2") ~ "North",
      Site %in% c("3","3.5","4") ~ "East",
      Site %in% c("5","5.5","6") ~ "West",
      TRUE ~ NA_character_
    ),
    Region = factor(Region, levels = c("North","East","West")),
    Year = factor(Year, levels = c("2017","2018","2019"))
  )

write.csv(alphaDiv, "Results/POC_2017_2019_alphaDiv_Complete.csv", row.names = FALSE)

# ---- 4. Filter Specific Samples ----
keep_POC_2017 <- c(paste0("POCFRG2017_0_", 1:5), paste0("POCFRG2017_1_", 6:9), paste0("POCFRG2017_3.5_", 12:14), paste0("POCFRG2017_5_", 16:20), paste0("POCFRG2017_5.5_", 21:25), paste0("POCFRG2017_6_", 26:32))
keep_POC_2018 <- c(paste0("POCFRG2018_2_", 35:39), paste0("POCFRG2018_5_", 41:43), paste0("POCFRG2018_5.5_", 44:48), paste0("POCFRG2018_6_", 49:52))
keep_POC_2019 <- c(paste0("POCFRG2019_1_", 54:58), paste0("POCFRG2019_2_", 59:62), paste0("POCFRG2019_3_", 63:68), paste0("POCFRG2019_3.5_", 69:73), paste0("POCFRG2019_4_", 74:77), paste0("POCFRG2019_5_", 78:88), paste0("POCFRG2019_6_", 89:94))

keep_samples_poc <- c(keep_POC_2017, keep_POC_2018, keep_POC_2019)
alphaDiv_POC <- alphaDiv %>% filter(Samples %in% keep_samples_poc)

# ---- 5. Plot Temporal Trend per Region ----
box_w <- 0.6

plot_year_POC <- ggplot(alphaDiv_POC, aes(x = Year, y = Shannon, color = Year)) +
  geom_boxplot(width = box_w, outlier.shape = NA, fill = NA, linewidth = 0.9) +
  geom_point(position = position_jitter(width = 0.15), size = 1.4, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 2.8, shape = 18) +
  stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 1.2, color = "black", alpha = 0.7) +
  scale_color_manual(values = c("2017" = "black", "2018" = "blue", "2019" = "orange")) +
  facet_wrap(~Region, scales = "free_x") +
  theme_bw() + 
  labs(title = "POC Shannon Diversity: Temporal Trend per Region", x = "Year", y = "Shannon Index") +
  theme(strip.text = element_text(face = "bold", size = 12), axis.text = element_text(face = "bold"), aspect.ratio = 1)

print(plot_year_POC)

# ---- 6. Statistical Testing (Tukey HSD) ----
# Tukey Test: Across Years within Regions
stat_test_region_year <- alphaDiv_POC %>%
  filter(!is.na(Region), !is.na(Shannon)) %>%
  group_by(Region) %>%
  filter(n_distinct(Year) > 1) %>%
  tukey_hsd(Shannon ~ Year) %>%
  ungroup()

write.csv(stat_test_region_year, "Results/POC_Shannon_Tukey_Region_Across_Year.csv", row.names = FALSE)


# ==============================================================================
# ---- 7. LTER Temperature Data Processing ----
# ==============================================================================
dat <- fread("Data/MCR_LTER01_BTM_Fringe_20250521.csv")

dat[, time_local := as.POSIXct(time_local, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]
dat[, `:=`(year = as.integer(format(time_local, "%Y")), month = as.integer(format(time_local, "%m")), date = as.Date(time_local))]

# Daily Top 10% Mean Temperature (Example: 2019)
daily_top10 <- dat[year == 2019, {
    cutoff <- quantile(temperature_c, 0.9, na.rm = TRUE)
    .(top10_mean_temp = mean(temperature_c[temperature_c >= cutoff], na.rm = TRUE))
  }, by = .(date, year, month)]

monthly_avg_top10 <- daily_top10[, .(mean_top10_temp = mean(top10_mean_temp, na.rm = TRUE)), by = .(year, month)][order(month)]
fwrite(monthly_avg_top10, "Results/LTER01_2019_MonthlyMean_Top10_Temp.csv")


# ==============================================================================
# ---- 8. Environmental / Nutrient Plot (WITH REGION COLORS) ----
# ==============================================================================
nutrient_data <- read_csv("Data/Fringing_2017_2019ave.csv", show_col_types = FALSE)

nutrient_filtered <- nutrient_data %>%
  mutate(Year = as.character(Year), Site = as.character(Site)) %>%
  filter(Year %in% c("2017","2019")) %>%
  mutate(
    Region = case_when(
      Site %in% c("0","1","2") ~ "North",
      Site %in% c("3","3.5","4") ~ "East",
      Site %in% c("5","5.5","6") ~ "West",
      TRUE ~ NA_character_
    ),
    Region = factor(Region, levels = c("North","East","West"))
  )

#  Nutrient Plot
plot_nutrient <- ggplot(nutrient_filtered, aes(x = Site, y = Si, color = Region)) +
  geom_point(size = 4, alpha = 0.8) +
  facet_wrap(~Year) +
  scale_color_manual(values = c(
    "North" = "dodgerblue",
    "East"  = "goldenrod",
    "West"  = "deeppink"  
  )) +
  theme_bw() +
  labs(title = "Silicate (Si) Concentrations by Region & Year", x = "Site", y = "Silicate Concentration") +
  theme(strip.text = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 11),
        axis.title = element_text(face = "bold", size = 12))

print(plot_nutrient)
