#Author: Rio Kashinmoto
#Last edit: April 6th, 2026
#Project: MVP Coral Microbiome Analysis (Porites & Pocillopora)
#Description:　Beta diversity (PCA via rCLR/CLR), PERMANOVA testing, 3D PCA plotting, and Environmental Vector fitting (EnvFit) 
for Pocillopora (POC) and Porites lobota (POR) samples.

# ---- 1. Load Packages ----
set.seed(123)
library(phyloseq)
library(microbiome)
library(vegan)
library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(gridExtra)
library(plotly)
library(data.table)
library(readr)
library(ggrepel)
library(scales)

### Beta diversity for POC
# ---- 2. Load & Subset Data (POC Only) ----
ps <- readRDS("Data/ps_POC_POR_filtered_clean.rds")
ps@refseq <- NULL
ps@phy_tree <- NULL

# Keep strictly Pocillopora (POC)
ps <- prune_samples(grepl("POC", sample_names(ps)), ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Aggregate to Genus
ps.genus <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)
ps.genus <- prune_taxa(taxa_sums(ps.genus) > 0, ps.genus)

# Define standard color palettes
year_colors <- c("2017" = "#F8766D", "2018" = "#00BE67", "2019" = "#00B8E7")
region_colors <- c("North" = "blue", "East" = "#20A387FF", "West" = "deeppink")

# ==============================================================================
# ---- PART A: rCLR Transformation, PERMANOVA, and 2D PCA Plots ----
# ==============================================================================

# 1. rCLR Transformation
otu <- as(otu_table(ps.genus), "matrix")
if(!taxa_are_rows(ps.genus)) { otu <- t(otu) }

otu <- otu + 1
log_otu <- log(otu)
otu_rclr <- log_otu - rowMeans(log_otu)

# 2. PCA Calculation
pca_res <- rda(t(otu_rclr))
eig_vals <- pca_res$CA$eig
var_exp <- eig_vals / sum(eig_vals)
pc1_var <- round(var_exp[1] * 100, 1)
pc2_var <- round(var_exp[2] * 100, 1)

scores_df <- scores(pca_res, display="sites") %>% as.data.frame()
scores_df$SampleID <- rownames(scores_df)

# 3. Metadata Formatting
meta_df <- data.frame(sample_data(ps.genus), stringsAsFactors = FALSE)
meta_df$SampleID <- rownames(meta_df)
meta_df$Site <- str_extract(meta_df$SampleID, "(?<=_)[0-9]+\\.?[0-9]*(?=_)")
meta_df$Region <- case_when(
  meta_df$Site %in% c("0","1","2") ~ "North",
  meta_df$Site %in% c("3","3.5","4") ~ "East",
  meta_df$Site %in% c("5","5.5","6") ~ "West",
  TRUE ~ NA_character_
)
meta_df$Year <- factor(meta_df$Year, levels=c("2017","2018","2019"))

plot_df <- left_join(scores_df, meta_df, by="SampleID")

# 4. Distance Matrix & PERMANOVA
dist_mat <- vegdist(t(otu_rclr), method="euclidean")

# PERMANOVA: Region within Year
permanova_region_year <- plot_df %>%
  split(.$Year) %>%
  map_df(function(df){
    ids <- df$SampleID
    dist_sub <- dist_mat[ids, ids]
    meta_sub <- meta_df[ids,]
    res <- adonis2(dist_sub ~ Region, data = meta_sub)
    data.frame(
      Year = unique(meta_sub$Year),
      R2 = round(res$R2[1],3),
      p = signif(res$`Pr(>F)`[1],3)
    )
  })

permanova_region_year$label <- paste0("R² = ", permanova_region_year$R2, "\np = ", permanova_region_year$p)

# PERMANOVA: Overall Year
permanova_year <- adonis2(dist_mat ~ Year, data=meta_df)
label_year <- paste0("R² = ", round(permanova_year$R2[1],3), "\np = ", signif(permanova_year$`Pr(>F)`[1],3))

# 5. Plot: PCA by Region (Faceted by Year)
p_region <- ggplot(plot_df, aes(PC1, PC2, color = Region)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = Region), linewidth = 0.8, type = "norm", na.rm = TRUE) +
  facet_wrap(~Year) +
  geom_text(data = permanova_region_year, aes(x = -Inf, y = -Inf, label = label), inherit.aes = FALSE, hjust = -0.1, vjust = -0.3, size = 3.2) +
  scale_color_manual(values = region_colors) +
  scale_x_continuous(labels = label_number(accuracy = 1)) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  coord_fixed() + theme_bw() +
  labs(title = "PCA of Genus-level rCLR (POCFRG)", x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC2 (", pc2_var, "%)"), color = "Region") +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), strip.text = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"))

# 6. Plot: PCA Colored by Year
p_year <- ggplot(plot_df, aes(PC1, PC2, color = Year)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = Year), linewidth = 1, type = "t") +
  scale_color_manual(values = year_colors) +
  annotate("text", x = -Inf, y = -Inf, label = label_year, hjust = -0.1, vjust = -0.3, size = 3.2) +
  coord_fixed() + theme_bw() +
  labs(title = "PCA of Genus-level rCLR (POCFRG)", x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC2 (", pc2_var, "%)"), color = "Year") +
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"))

print(p_region)
print(p_year)


# ==============================================================================
# ---- PART B: 3D PCA Plotting (CLR Transformed) ----
# ==============================================================================

ps_clr <- microbiome::transform(ps.genus, "clr")
otu_clr <- as.data.frame(otu_table(ps_clr))
if(taxa_are_rows(ps_clr)) { otu_clr <- t(otu_clr) }

pca_res_3d <- prcomp(otu_clr, center = TRUE, scale. = FALSE)
var_expl_3d <- summary(pca_res_3d)$importance[2,1:3]
PC1_var <- round(var_expl_3d[1]*100,1)
PC2_var <- round(var_expl_3d[2]*100,1)
PC3_var <- round(var_expl_3d[3]*100,1)

scores_3d <- as.data.frame(pca_res_3d$x[,1:3])
scores_3d$SampleID <- rownames(scores_3d)
scores_3d <- left_join(scores_3d, meta_df, by="SampleID")

pca_3d <- plot_ly(
  data = scores_3d, x = ~PC1, y = ~PC2, z = ~PC3,
  color = ~Year, colors = year_colors,
  type = "scatter3d", mode = "markers", marker = list(size = 5, opacity = 0.9)
) %>% layout(
  scene = list(
    xaxis = list(title = list(text = paste0("<b>PC1 (", PC1_var, "%)</b>"))),
    yaxis = list(title = list(text = paste0("<b>PC2 (", PC2_var, "%)</b>"))),
    zaxis = list(title = list(text = paste0("<b>PC3 (", PC3_var, "%)</b>")))
  ),
  title = list(text = "<b>3D PCA (CLR-transformed POCFRG)</b>", font = list(size = 22))
)

pca_3d


# ==============================================================================
# ---- PART C: PCA with Environmental Vectors (EnvFit) ----
# ==============================================================================

# 1. Load & Format Environmental Data
dat_temp <- fread("Data/LTER01_06_2017_2019_Jan_Dec_MonthlyMean_Temp_12PM.csv")
temp_avg <- dat_temp %>%
  filter(grepl("Wet_NovApr", Time)) %>%
  mutate(Year = as.character(Year), Site = trimws(as.character(Site))) %>%
  group_by(Year, Site) %>%
  summarise(Temperature = mean(Temp, na.rm = TRUE), .groups = "drop")

dat_nut <- read_csv("Data/Fringing_2017_2019ave.csv", show_col_types = FALSE)
nutrient_avg <- dat_nut %>%
  mutate(Year = as.character(Year), Site = trimws(as.character(Site))) %>%
  group_by(Year, Site) %>%
  summarise(C = mean(C, na.rm = TRUE), H = mean(H, na.rm = TRUE), N = mean(N, na.rm = TRUE), CN_ratio = mean(CN_ratio, na.rm = TRUE), .groups = "drop")

# 2. Merge Environment Data with existing plot_df (from Part A)
env_df <- plot_df %>%
  select(SampleID, Year, Site, Region, PC1, PC2) %>%
  mutate(Year = as.character(Year)) %>%
  left_join(temp_avg, by = c("Year","Site")) %>%
  left_join(nutrient_avg, by = c("Year","Site")) %>%
  filter(!is.na(Temperature), !is.na(C), !is.na(H), !is.na(N), !is.na(CN_ratio), !is.na(Region))

# 3. Calculate EnvFit
env_numeric <- env_df %>% select(Temperature, C, H, N, CN_ratio)
pca_scores_filtered <- env_df %>% select(PC1, PC2) %>% as.matrix()

envfit_res <- envfit(pca_scores_filtered, env_numeric, permutations = 999)
env_vectors <- as.data.frame(scores(envfit_res, display = "vectors"))
env_vectors$Variable <- rownames(env_vectors)
env_vectors$pval <- envfit_res$vectors$pvals

# 4. Plot EnvFit by Region
arrow_scale <- 4

p_env_region <- ggplot(env_df, aes(x = PC1, y = PC2, color = Region)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Region), linewidth = 0.8, type = "t") +
  geom_segment(data = env_vectors, aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale), arrow = arrow(length = unit(0.3, "cm")), color = "black", inherit.aes = FALSE) +
  geom_text_repel(data = env_vectors, aes(x = PC1 * arrow_scale, y = PC2 * arrow_scale, label = paste0(Variable, "\n(p=", round(pval, 3), ")")), size = 4.5, color = "black", segment.color = "grey40", box.padding = 0.5, inherit.aes = FALSE) +
  scale_color_manual(values = region_colors) +
  coord_fixed() + theme_bw() +
  labs(title = "POC PCA with Wet Season Env Vectors (By Region)", x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC2 (", pc2_var, "%)")) +
  theme(axis.title = element_text(size = 16, face = "bold"), axis.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"))

print(p_env_region)

# 5. Plot EnvFit by Year
p_env_year <- ggplot(env_df, aes(x = PC1, y = PC2, color = as.factor(Year))) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Year), linewidth = 0.8, type = "t") +
  geom_segment(data = env_vectors, aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale), arrow = arrow(length = unit(0.3, "cm")), color = "black", inherit.aes = FALSE) +
  geom_text_repel(data = env_vectors, aes(x = PC1 * arrow_scale, y = PC2 * arrow_scale, label = paste0(Variable, "\n(p=", round(pval, 3), ")")), size = 4.5, color = "black", segment.color = "grey40", box.padding = 0.5, inherit.aes = FALSE) +
  scale_color_manual(values = year_colors) +
  coord_fixed() + theme_bw() +
  labs(title = "POC PCA with Wet Season Env Vectors (By Year)", x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC2 (", pc2_var, "%)"), color = "Year") +
  theme(axis.title = element_text(size = 16, face = "bold"), axis.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"))

print(p_env_year)
