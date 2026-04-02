
# BSNB: Species composition, amplification effect, composition effect, key taxa

## PACKAGES, DATA, FOLDERS
library(tidyverse)
library(vegan)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(RColorBrewer)
library(scales)
library(reshape2)
library(pheatmap)
library(ggpubr)
library(tidytext)
library(dplyr)
library(stringr)
library(tidyr)
library(ggrepel)
library(indicspecies)
library(ggridges)
library(grid)
library(tibble)


# data
data <- read.csv("data3.csv") %>% janitor::clean_names()

# Output folders

#all plots included are in fungal_diversity_results1 folder
dir.create("fungal_diversity_results1", showWarnings = FALSE)
dir.create("composition_results", showWarnings = FALSE)
dir.create("key_taxa_results", showWarnings = FALSE)
dir.create("ordination_results", showWarnings = FALSE)
dir.create("correlation_plots_alternative", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. DATA PREPARATION
## ---------------------------------------------------------------------------

# Disease severity categories
data <- data %>%
  mutate(
    # standardize text
    Disease_class = str_to_lower(disease_class),
    
    # map to ordered Severity factor
    Severity = case_when(
      disease_class == "low"      ~ "Low",
      disease_class == "moderate" ~ "Moderate",
      disease_class == "high"     ~ "High",
      TRUE                       ~ NA_character_
    ),
    Severity = factor(Severity, levels = c("Low", "Moderate", "High"))
  )

# Taxonomic cleaning
data <- data %>%
  mutate(
    species_clean = if_else(
      is.na(speciesidentity) | speciesidentity %in% c("NA", ""),
      paste0("Unidentified_", coalesce(genus, "sp")),
      speciesidentity
    ),
    order_clean   = if_else(is.na(order) | order %in% c("NA", ""),
                            "Unknown_Order", order),
    lineage_clean = if_else(is.na(lineage_class) | lineage_class %in% c("NA", ""),
                            "Unknown_Class", lineage_class),
    genus_clean   = if_else(is.na(genus) | genus %in% c("NA", ""),
                            "Unknown_Genus", genus)
  )

#a.	Total isolates, species, genera, lineages

cat(sprintf("  Total samples: %d\n", nrow(data)))
cat(sprintf("  Unique species: %d\n", length(unique(data$species_clean))))
cat(sprintf("  Unique genera: %d\n", length(unique(data$genus_clean))))
cat(sprintf("  Unique orders: %d\n", length(unique(data$order_clean))))
cat(sprintf("  Unique lineages/classes: %d\n", length(unique(data$lineage_clean))))
cat(sprintf("  Unique functionalclass: %d\n", length(unique(data$functionalclass))))
print(unique(data$functionalclass))

# Filter to functional guilds of interest (endophytes + pathogens) – reused later
allowed_classes <- c("Plant Pathogen", "Endophyte")
#top pathogens and endophytes proportions

pe_data <- data %>%
  mutate(
    functionalclass = str_trim(functionalclass),
    species_clean   = str_trim(species_clean)
  ) %>%
  filter(
    !is.na(functionalclass),
    functionalclass %in% allowed_classes,
    !is.na(species_clean),
    species_clean != "",
    species_clean != "NA"
  )

# 2. Total number of pathogen + endophyte isolates
total_pe <- nrow(data)

# 3. Species frequencies and percentages
species_summary <- pe_data %>%
  count(species_clean, sort = TRUE) %>%
  mutate(
    perc_of_total_PE = 100 * n / total_pe
  )

# View the full table 
species_summary

# 4. Get the top k most frequently isolated taxa
top_k <- 10  # change to 5, 8, etc. as needed

top_species_summary <- species_summary %>%
  slice_max(order_by = n, n = top_k)

top_species_summary

# Filtered per functional class
data_filtered <- data %>%
  mutate(functionalclass = str_trim(functionalclass)) %>%
  filter(!is.na(functionalclass),
         functionalclass %in% allowed_classes,
         !is.na(species_clean),
         species_clean != "",
         species_clean != "NA")

# 2. SPECIES COMPOSITION 
## 2.1 Lineage/Class composition by severity  (Fig 4A )

composition_severity <- data %>%
  group_by(Severity, lineage_clean) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Severity) %>%
  mutate(proportion = count / sum(count) * 100)
composition_severity
p_lineage_severity <- ggplot(composition_severity,
                             aes(x = Severity, y = proportion,
                                 fill = lineage_clean)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = colorRampPalette(
    brewer.pal(12, "Set3"))(length(unique(composition_severity$lineage_clean)))
  ) +
  labs(title = "Fungal Lineage/Class Composition by Disease Severity",
       x = "Disease Severity",
       y = "Percentage (%)",
       fill = "Lineage/Class") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        legend.text = element_text(size = 10))

ggsave("fungal_diversity_results1/fig_lineage_composition_severity.png",
       p_lineage_severity, width = 12, height = 7, dpi = 300)

## 2.2 Top species (Endophytes + Plant Pathogens) by severity  (Fig 1B)

top_species <- data_filtered %>%
  count(species_clean, sort = TRUE) %>%
  slice_max(n, n = 5) %>%
  pull(species_clean)
top_species


species_composition <- data_filtered %>%
  mutate(species_grouped = ifelse(species_clean %in% top_species,
                                  species_clean, "Other")) %>%
  group_by(functionalclass, Severity, species_grouped) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(proportion = 100 * count / sum(count)) %>%
  ungroup() %>%
  filter(species_grouped != "Other")
print(species_composition)

p_top_species_severity <- ggplot(
  species_composition,
  aes(x = Severity, y = proportion, fill = species_grouped)
) +
  geom_col(width = 0.7, alpha = 1, colour = "black", linewidth = 0.15) +
  labs(
    #title = "Top Species(Endophytes & Plant Pathogens)",
    x = "Disease Severity",
    y = "Percentage (%)",
    fill = "Species"
    
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        legend.text = element_text(size = 10, face = "italic"))
p_top_species_severity
ggsave("fungal_diversity_results1/fig_top_species_composition.png",
       p_top_species_severity, width = 10, height = 7, dpi = 300)

# Combined all composition summaries (all plots Fig 1 used)
p_composition_taxa <- p_lineage_severity | p_top_species_severity 
ggsave("fungal_diversity_results1/Fig4_Composition_severity.png",
       p_composition_taxa, width = 15, height = 5, dpi = 300)
p_composition_taxa


# 3. AMPLIFICATION EFFECT  (richness & α-diversity vs disease severity and region)
## 3.1 Richness by severity and region (your Fig 2)

total_species  <- n_distinct(data$species_clean,  na.rm = TRUE)
total_lineages <- n_distinct(data$lineage_clean, na.rm = TRUE)

# Richness by severity
richness_severity <- data %>%
  filter(!is.na(Severity)) %>%
  group_by(Severity) %>%
  summarise(
    Species = n_distinct(species_clean),
    Lineage = n_distinct(lineage_clean),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(Species, Lineage),
               names_to = "Taxonomic_Level",
               values_to = "Richness") %>%
  mutate(
    Total = case_when(
      Taxonomic_Level == "Species" ~ total_species,
      Taxonomic_Level == "Lineage" ~ total_lineages
    ),
    Percent = (Richness / Total) * 100
  ) %>%
  select(Severity, Taxonomic_Level, Richness, Percent)

p_richness_severity <- ggplot(
  richness_severity,
  aes(x = Severity, y = Percent, fill = Taxonomic_Level)
) +
  geom_bar(
    stat = "identity",
    position = "dodge",
    width = 0.7,
    colour = "black",
    linewidth = 0.2
  ) +
  scale_fill_viridis_d(option = "cividis", end = 0.9) +
  labs(
    x = "Disease Severity",
    y = "Percent of Total Richness",
    fill = "Taxonomic Level"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  ) +
  guides(fill = guide_legend(nrow = 1))

ggsave("fungal_diversity_results1/fig_richness_severity.png",
       p_richness_severity, width = 12, height = 6, dpi = 300)
richness_severity

# Richness by region
richness_region <- data %>%
  group_by(region) %>%
  summarise(
    Species = n_distinct(species_clean),
    Lineage = n_distinct(lineage_clean),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(Species, Lineage),
               names_to = "Taxonomic_Level",
               values_to = "Richness") %>%
  mutate(
    Total = case_when(
      Taxonomic_Level == "Species" ~ total_species,
      Taxonomic_Level == "Lineage" ~ total_lineages
    ),
    Percent = (Richness / Total) * 100
  ) %>%
  select(region, Taxonomic_Level, Richness, Percent)

p_richness_region <- ggplot(
  richness_region,
  aes(x = region, y = Percent, fill = Taxonomic_Level)
) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Fungal Richness by Region",
       x = "Region",
       y = "Percent of Total Richness",
       fill = "Taxonomic Level") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right")

ggsave("fungal_diversity_results1/fig_richness_region.png",
       p_richness_region, width = 12, height = 6, dpi = 300)
richness_region
p_richness_combined <- p_richness_severity | p_richness_region
ggsave("fungal_diversity_results1/Fig2_richness_severity_region.png",
       p_richness_combined, width = 14, height = 5, dpi = 300)




##  Diversity metrics per plot & correlation with severity ( Fig 6)
# 2.	Association between fungal species richness and disease severity
diversity_per_plot <- data %>%
  group_by(plot, proportion, Severity, region, county) %>%
  summarise(
    Species_Richness  = n_distinct(species_clean),
    Genus_Richness    = n_distinct(genus_clean),
    Order_Richness    = n_distinct(order_clean),
    Lineage_Richness  = n_distinct(lineage_clean),
    Total_Isolates    = n(),
    .groups = "drop"
  )

shannon_per_plot <- data %>%
  group_by(plot,species_clean) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(plot) %>%
  summarise(
    Shannon_Diversity = diversity(count, index = "shannon"),
    Simpson_Diversity = diversity(count, index = "simpson"),
    .groups = "drop"
  )

diversity_per_plot <- diversity_per_plot %>%
  left_join(shannon_per_plot, by = "plot")
diversity_per_plot

write.csv(diversity_per_plot,
          "fungal_diversity_results1/diversity_metrics_per_plot.csv",
          row.names = FALSE)

# Encode severity as ordinal for Spearman
diversity_per_plot <- diversity_per_plot %>%
  mutate(Severity_num = as.numeric(Severity))

cor_species <- cor.test(diversity_per_plot$Severity_num,
                        diversity_per_plot$Species_Richness,
                        method = "spearman")
cor_shannon <- cor.test(diversity_per_plot$Severity_num,
                        diversity_per_plot$Shannon_Diversity,
                        method = "spearman")
cor_shannon
cor_species
# Boxplots
p_box_species <- ggplot(diversity_per_plot,
                        aes(x = Severity, y = Species_Richness, fill = Severity)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  scale_fill_manual(values = c("Low" = "#66C2A5",
                               "Moderate" = "#FC8D62",
                               "High" = "#E78AC3")) +
  labs(title = "Species Richness by Disease Severity",
       subtitle = sprintf("Spearman rho = %.3f, p = %.4f",
                          cor_species$estimate, cor_species$p.value),
       x = "Disease Severity Category",
       y = "Species Richness") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ggsave("correlation_plots_alternative/fig_species_richness_boxplot.png",
       p_box_species, width = 8, height = 6, dpi = 300)

p_box_shannon <- ggplot(diversity_per_plot,
                        aes(x = Severity, y = Shannon_Diversity, fill = Severity)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  scale_fill_manual(values = c("Low" = "#66C2A5",
                               "Moderate" = "#FC8D62",
                               "High" = "#E78AC3")) +
  labs(x = "Disease Severity Category",
       y = "Shannon Diversity Index") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ggsave("correlation_plots_alternative/fig_shannon_diversity_boxplot.png",
       p_box_shannon, width = 8, height = 6, dpi = 300)

# Ridgeline distribution of richness
p_ridge_species <- ggplot(
  diversity_per_plot,
  aes(x = Species_Richness, y = Severity, fill = Severity)
) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 1.5,
    colour = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "Low" = "#66C2A5",
      "Moderate" = "#FC8D62",
      "High" = "#E78AC3",
      "Unknown" = "grey70"
    ),
    na.value = "grey70"
  ) +
  labs(
    x = "Species Richness",
    y = "Disease Severity Category"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_ridge_species
ggsave("correlation_plots_alternative/fig_species_richness_ridgeline.png",
       p_ridge_species, width = 10, height = 6, dpi = 300)

# Optional barplots with SE
summary_stats <- diversity_per_plot %>%
  group_by(Severity) %>%
  summarise(
    Mean_Richness = mean(Species_Richness, na.rm = TRUE),
    SE_Richness   = sd(Species_Richness, na.rm = TRUE) / sqrt(n()),
    Mean_Shannon  = mean(Shannon_Diversity, na.rm = TRUE),
    SE_Shannon    = sd(Shannon_Diversity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_bar_species <- ggplot(summary_stats,
                        aes(x = Severity, y = Mean_Richness, fill = Severity)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_Richness - SE_Richness,
                    ymax = Mean_Richness + SE_Richness),
                width = 0.2, size = 0.8) +
  scale_fill_manual(values = c("Low" = "#66C2A5",
                               "Moderate" = "#FC8D62",
                               "High" = "#E78AC3")) +
  labs(
       x = "Disease Severity Category",
       y = "Mean Species Richness") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

p_bar_shannon <- ggplot(summary_stats,
                        aes(x = Severity, y = Mean_Shannon, fill = Severity)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_Shannon - SE_Shannon,
                    ymax = Mean_Shannon + SE_Shannon),
                width = 0.2, size = 0.8) +
  scale_fill_manual(values = c("Low" = "#66C2A5",
                               "Moderate" = "#FC8D62",
                               "High" = "#E78AC3")) +
  labs(
       x = "Disease Severity Category",
       y = "Mean Shannon Diversity Index") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ggsave("correlation_plots_alternative/fig_species_richness_barplot.png",
       p_bar_species, width = 8, height = 6, dpi = 300)
ggsave("correlation_plots_alternative/fig_shannon_diversity_barplot.png",
       p_bar_shannon, width = 8, height = 6, dpi = 300)

# all plots figure 2 used ( Association) plot
p_correlation_combined <- p_richness_severity / p_ridge_species | p_box_species / p_box_shannon  
ggsave("fungal_diversity_results1/Fig4_diversity_severity_combined.png",
       p_correlation_combined, width = 14, height = 10, dpi = 300)

p_correlation_combined



# # Combined species + lineage panel ( not used)
# p_functional_taxa <- p_lineage_severity / p_richness_severity | p_top_species_severity / p_richness_region
# ggsave("fungal_diversity_results1/Fig1_species_lineage_severity.png",
#        p_functional_taxa, width = 15, height = 10, dpi = 300)
# p_functional_taxa

# 4. COMPOSITIONAL EFFECT (community PCA + PERMANOVA + environmental PCA)
# 4.1 Community matrix  fig 3
species_matrix_all <- data %>%
  group_by(plot, species_clean) %>%
  summarise(abundance = n(), .groups = "drop") %>%
  pivot_wider(names_from = species_clean,
              values_from = abundance,
              values_fill = 0)

plot_metadata_all <- data %>%
  distinct(plot, .keep_all = TRUE) %>%
  select(plot, Severity, region, county, proportion) %>%
  filter(!is.na(Severity), !is.na(region))

species_data_all <- species_matrix_all %>%
  filter(plot %in% plot_metadata_all$plot) %>%
  column_to_rownames("plot") %>%
  as.matrix()

plot_metadata_all <- plot_metadata_all %>%
  filter(plot %in% rownames(species_data_all)) %>%
  arrange(match(plot, rownames(species_data_all)))

# 2. Hellinger transform + PCA
species_hell_all <- decostand(species_data_all, method = "hellinger")
comm_pca_all     <- rda(species_hell_all)

pca_var_all <- eigenvals(comm_pca_all) / sum(eigenvals(comm_pca_all)) * 100

pca_scores_all <- scores(comm_pca_all, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("plot") %>%
  left_join(plot_metadata_all, by = "plot") %>%
  mutate(
    Severity = str_trim(Severity),
    Severity = factor(Severity, levels = c("Low", "Moderate", "High")),
    region   = str_trim(region),
    region   = factor(region, levels = c("North", "South"))  # keep 2 regions
  ) %>%
  # keep only rows with non-missing scores and region
  filter(!is.na(PC1), !is.na(PC2), !is.na(region), !is.na(Severity))

table(pca_scores_all$region)
table(pca_scores_all$Severity)


# 3. PCA plot by region with two ellipses
p_pca_region <- ggplot(pca_scores_all,
                       aes(x = PC1, y = PC2)) +
  geom_point(aes(color = region, shape = Severity),
             size = 4, alpha = 0.7) +
  stat_ellipse(aes(color = region, group = region),
               level = 0.95, linewidth = 1) +
  scale_color_brewer(palette = "Set1", name = "Region") +
  labs(
    title    = "Fungal Community  by Region",
    x = sprintf("PC1 (%.1f%%)", pca_var_all[1]),
    y = sprintf("PC2 (%.1f%%)", pca_var_all[2])
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 16)
  )

p_pca_region #didnt include in figure
ggsave("composition_results/fig_pca_by_region.png",
       p_pca_region, width = 12, height = 8, dpi = 300)


# # PCA by Severity
# p_pca_severity <- ggplot(pca_scores_all, 
#                          aes(x = PC1, y = PC2)) +
#   geom_point(aes(color = Severity, shape = region),
#              size = 4, alpha = 0.7) +
#   stat_ellipse(aes(color = Severity, group = Severity),  # ← ADD THIS
#                level = 0.95, linewidth = 1) +
#   scale_color_brewer(palette = "Set1", name = "Severity") + 
#   labs(
#     title = "Fungal Community by Diseasegroup",
#     x = sprintf("PC1 (%.1f%%)", pca_var_all[1]),
#     y = sprintf("PC2 (%.1f%%)", pca_var_all[2])
#   ) +
#   theme_classic(base_size = 14) +
#   theme(
#     legend.position = "right",
#     plot.title = element_text(face = "bold", size = 16)
#   )
# 
# p_pca_severity


library(dplyr)
library(ggplot2)

## 1. Build a data frame of hull points for each Severity
hulls <- pca_scores_all %>%
  group_by(Severity) %>%
  slice(chull(PC1, PC2))   # chull() gives the convex hull indices

## 2. Plot with polygons (hulls) instead of ellipses
p_pca_severity <- ggplot(pca_scores_all,
                         aes(x = PC1, y = PC2)) +
  # hull outline + light fill
  #geom_polygon(data = hulls,
               #aes(x = PC1, y = PC2, group = Severity, fill = Severity),
               #color = NA, alpha = 0.12) +
  #geom_path(data = hulls,
            #aes(x = PC1, y = PC2, group = Severity, color = Severity),
            #linewidth = 1) +
  # points
  geom_point(aes(color = Severity, shape = region),
             size = 4, alpha = 0.8) +
  scale_color_brewer(palette = "Set1", name = "Severity") +
  #scale_fill_brewer(palette = "Set1", guide = "none") +
  labs(
    #title = "Fungal Community by Diseasegroup",
    x = sprintf("PC1 (%.1f%%)", pca_var_all[1]),
    y = sprintf("PC2 (%.1f%%)", pca_var_all[2])
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16)
  )

p_pca_severity

# fig 3 used
ggsave("composition_results/fig1_pca_by_severity.png", 
       p_pca_severity, width = 12, height = 8, dpi = 300)

## Environmental PCA + community PCA by region (for combined figure 6) 

cat("\n=== ENVIRONMENTAL PCA ===\n")

# Select environmental variables (remove NAs)
env_vars <- data %>%
  select(isolates, plot, subplot, Severity, region, county,
         tree_height, tree_dbh, basel_area, slope, elevation) %>%
  distinct(plot, .keep_all = TRUE) %>%  # One row per plot
  drop_na(tree_height, tree_dbh, basel_area, elevation)

# Standardize environmental variables
env_matrix <- env_vars %>%
  select(tree_height, tree_dbh, basel_area, slope, elevation) %>%
  scale()

# Run PCA
env_pca <- prcomp(env_matrix, center = FALSE, scale. = FALSE)

# Extract scores
pca_scores <- as.data.frame(env_pca$x) %>%
  bind_cols(env_vars %>% select(plot, Severity, region, county))

# Variance explained
var_explained <- summary(env_pca)$importance[2,] * 100

cat(sprintf("  PC1 explains %.1f%% of variance\n", var_explained[1]))
cat(sprintf("  PC2 explains %.1f%% of variance\n", var_explained[2]))

# Plot: Environmental PCA
p_env_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Severity, shape = region)) +
  geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = c("Low" = "#66C2A5", "Moderate" = "#FC8D62", "High" = "#E78AC3")) +
  labs(title = "Environmental PCA",
       x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2])) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 16))
p_env_pca
ggsave("ordination_results/fig1_environmental_pca.png", 
       p_env_pca, width = 10, height = 7, dpi = 300)

# PCA biplot with shorter loadings
arrow_mult <- 3.0

loadings <- as.data.frame(env_pca$rotation[, 1:2]) %>%
  mutate(
    Variable = rownames(.),
    Variable = recode(
      Variable,
      tree_dbh    = "Tree DBH",
      basel_area  = "Basal area",
      elevation   = "Elevation",
      slope       = "Slope",
      tree_height = "Tree height"
    ),
    arrow_x = PC1 * arrow_mult,
    arrow_y = PC2 * arrow_mult,
    label_x = PC1 * (arrow_mult + 0.35),
    label_y = PC2 * (arrow_mult + 0.35)
  )

p_env_biplot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Severity), size = 3, alpha = 0.6) +
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
    arrow = arrow(length = unit(0.15, "cm")),
    color = "black",
    linewidth = 0.7
  ) +
  geom_text_repel(
    data = loadings,
    aes(x = label_x, y = label_y, label = Variable),
    size = 5,
    fontface = "bold",
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "NA",
    max.overlaps = Inf
  ) +
  scale_color_brewer(palette = "Set1", name = "Severity") +
  labs(
    x = sprintf("PC2 (%.1f%%)", var_explained[2]),
    y = sprintf("PC1 (%.1f%%)", var_explained[1])
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16)
  )

p_env_biplot

ggsave("composition_results/fig2_environmental_pca_biplot.png",
       p_env_biplot, width = 10, height = 5, dpi = 300)

#not used similar to Pca by region
# env_vars <- data %>%
#   select(isolates, plot, subplot, Severity, region, county,
#          tree_ht, tree_dbh, basel_area, slope, elevation) %>%
#   distinct(plot, .keep_all = TRUE) %>%
#   drop_na(tree_ht, tree_dbh, basel_area, elevation)
# 
# env_matrix <- env_vars %>%
#   select(tree_ht, tree_dbh, basel_area, slope, elevation) %>%
#   scale()
# 
# env_pca <- prcomp(env_matrix, center = FALSE, scale. = FALSE)
# 
# env_scores <- as.data.frame(env_pca$x) %>%
#   bind_cols(env_vars %>% select(plot, Severity, region, county))
# 
# var_explained <- summary(env_pca)$importance[2, ] * 100
# var_explained

# not sure i need this plot
# p_env_pca <- ggplot(env_scores,
#                     aes(x = PC1, y = PC2, color = Severity, shape = region)) +
#   geom_point(size = 4, alpha = 0.7) +
#   scale_color_manual(values = c("Low" = "#66C2A5",
#                                 "Moderate" = "#FC8D62",
#                                 "High" = "#E78AC3")) +
#   labs(title = "Environmental PCA",
#        x = sprintf("PC1 (%.1f%%)", var_explained[1]),
#        y = sprintf("PC2 (%.1f%%)", var_explained[2])) +
#   theme_classic(base_size = 14) +
#   theme(legend.position = "right",
#         plot.title = element_text(face = "bold", size = 16))
# p_env_pca
# ggsave("ordination_results/fig_environmental_pca.png",
#        p_env_pca, width = 10, height = 7, dpi = 300)



## PERMANOVA (Bray–Curtis, controlling for region)

bray_dist_all <- vegdist(species_data_all, method = "bray")

perm_full <- adonis2(bray_dist_all ~ region + Severity,
                     data = plot_metadata_all,
                     permutations = 999,
                     method = "bray",
                     by = "margin")

print(perm_full)
region_r2   <- perm_full$R2[1]
region_p    <- perm_full$`Pr(>F)`[1]
severity_r2 <- perm_full$R2[2]
severity_p  <- perm_full$`Pr(>F)`[2]

cat("\nINTERPRETATION:\n")
cat(sprintf("  Region effect: R² = %.3f, p = %.4f\n", region_r2, region_p))
cat(sprintf("  Severity effect (controlling for region): R² = %.3f, p = %.4f\n",
            severity_r2, severity_p))

# Dispersion test
disp_severity <- betadisper(bray_dist_all, plot_metadata_all$Severity)
disp_test <- permutest(disp_severity, permutations = 999)
cat(sprintf("  Dispersion test p-value: %.4f\n",
            disp_test$tab$`Pr(>F)`[1]))
disp_test$tab$F[1]

## Region-stratified PERMANOVA (North, South)

# North
north_plots <- plot_metadata_all %>% filter(region == "North") %>% pull(plot)
if (length(north_plots) >= 5) {
  species_north <- species_data_all[north_plots, , drop = FALSE]
  species_north <- species_north[, colSums(species_north) > 0, drop = FALSE]
  metadata_north <- plot_metadata_all %>% filter(region == "North")
  bray_north <- vegdist(species_north, method = "bray")
  perm_north <- adonis2(bray_north ~ Severity,
                        data = metadata_north,
                        permutations = 999,
                        method = "bray")
  print(perm_north)
}

# South
south_plots <- plot_metadata_all %>% filter(region == "South") %>% pull(plot)
if (length(south_plots) >= 9) {
  species_south <- species_data_all[south_plots, , drop = FALSE]
  species_south <- species_south[, colSums(species_south) > 0, drop = FALSE]
  metadata_south <- plot_metadata_all %>% filter(region == "South")
  bray_south <- vegdist(species_south, method = "bray")
  perm_south <- adonis2(bray_south ~ Severity,
                        data = metadata_south,
                        permutations = 999,
                        method = "bray")
  print(perm_south)
}




## 2.3 Top 10 species by region 

species_region <- data_filtered %>%
  count(region, species_clean, name = "count") %>%
  group_by(region) %>%
  mutate(percentage = 100 * count / sum(count)) %>%
  arrange(region, desc(percentage)) %>%
  slice_head(n = 10)

print(species_region)

## 2.4 Functional composition by severity (Fig 4)

functional_composition <- data %>%
  filter(!is.na(functionalclass),
         functionalclass != "",
         functionalclass != "NA") %>%
  group_by(Severity, functionalclass) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Severity) %>%
  mutate(proportion = count / sum(count) * 100)

functional_composition

p_functional_composition <- ggplot(
  functional_composition,
  aes(x = Severity, y = proportion, fill = functionalclass)
) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(
    x = "Disease Severity",
    y = "Percentage (%)",
    fill = "Functional Class"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

ggsave("fungal_diversity_results1/fig_functional_composition.png",
       p_functional_composition, width = 10, height = 5, dpi = 300)




# 5. KEY TAXA EFFECT (indicator species, abundance, prevalence, high-severity taxa,)
## 5.1 Community matrix restricted to Endophytes + Pathogens

species_matrix_ep <- data_filtered %>%
  group_by(plot, species_clean) %>%
  summarise(abundance = n(), .groups = "drop") %>%
  pivot_wider(names_from = species_clean,
              values_from = abundance,
              values_fill = 0)

plot_metadata_ep <- data_filtered %>%
  distinct(plot, .keep_all = TRUE) %>%
  select(plot, Severity, region, county, proportion) %>%
  filter(!is.na(Severity))

species_data_ep <- species_matrix_ep %>%
  filter(plot %in% plot_metadata_ep$plot) %>%
  column_to_rownames("plot") %>%
  as.matrix()

plot_metadata_ep <- plot_metadata_ep %>%
  filter(plot %in% rownames(species_data_ep)) %>%
  arrange(match(plot, rownames(species_data_ep)))

cat(sprintf("✓ Matrix: %d plots × %d species (endophytes + pathogens)\n",
            nrow(species_data_ep), ncol(species_data_ep)))

## 5.2 Indicator species analysis
indval_result <- multipatt(species_data_ep,
                           plot_metadata_ep$Severity,
                           control = how(nperm = 999))

# All results (for inspection)
summary(indval_result, alpha = 1, indvalcomp = TRUE)

# Significant indicators
indval_sig <- indval_result$sign %>%
  rownames_to_column("Species") %>%
  filter(p.value < 0.05) %>%
  arrange(p.value)

if (nrow(indval_sig) > 0) {
  indval_sig <- indval_sig %>%
    left_join(
      data_filtered %>%
        select(species_clean, functionalclass) %>%
        distinct(),
      by = c("Species" = "species_clean")
    ) %>%
    mutate(
      Associated_Severity = case_when(
        s.Low == 1      ~ "Low",
        s.Moderate == 1 ~ "Moderate",
        s.High == 1     ~ "High",
        TRUE            ~ "Multiple"
      )
    )
  
  print(indval_sig %>%
          select(Species, functionalclass, Associated_Severity, stat, p.value))
  
  write.csv(indval_sig,
            "key_taxa_results/indicator_species_endophytes_pathogens.csv",
            row.names = FALSE)
  
  p_indval <- ggplot(indval_sig,
                     aes(x = reorder(Species, stat), y = stat,
                         fill = Associated_Severity)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Low" = "#66C2A5",
                                 "Moderate" = "#FC8D62",
                                 "High" = "#E78AC3")) +
    labs(title = "Indicator Species: Endophytes and Plant Pathogens",
         subtitle = "Species significantly associated with disease severity levels",
         x = "Species",
         y = "Indicator Value (Specificity × Fidelity)",
         fill = "Associated with") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  ggsave("key_taxa_results/fig_indicator_species.png",
         p_indval,
         width = 10,
         height = max(6, nrow(indval_sig) * 0.4),
         dpi = 300)
}

## 5.3 Relative abundance by functional class and severity

species_rel_abund_ep <- species_data_ep / rowSums(species_data_ep) * 100

abund_long_ep <- species_rel_abund_ep %>%
  as.data.frame() %>%
  rownames_to_column("plot") %>%
  pivot_longer(-plot, names_to = "Species", values_to = "Rel_Abundance") %>%
  left_join(plot_metadata_ep, by = "plot") %>%
  left_join(
    data_filtered %>%
      select(species_clean, functionalclass) %>%
      distinct(),
    by = c("Species" = "species_clean")
  )

func_summary <- abund_long_ep %>%
  group_by(functionalclass, Severity) %>%
  summarise(
    Mean_Abundance = mean(Rel_Abundance, na.rm = TRUE),
    SE = sd(Rel_Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_func_abund <- ggplot(func_summary,
                       aes(x = Severity, y = Mean_Abundance,
                           fill = functionalclass)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = Mean_Abundance - SE,
                    ymax = Mean_Abundance + SE),
                position = position_dodge(0.7),
                width = 0.2) +
  scale_fill_manual(values = c("Endophytes" = "#8DD3C7",
                               "Plant Pathogens" = "#FB8072")) +
  labs(title = "Functional Class Abundance by Disease Severity",
       x = "Disease Severity",
       y = "Mean Relative Abundance (%)",
       fill = "Functional Class") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16))

ggsave("key_taxa_results/fig_functional_class_abundance.png",
       p_func_abund, width = 10, height = 7, dpi = 300)

# ## 5.4 Top 10 species abundance by severity (Not needed )
# 
# top_species_abund <- abund_long_ep %>%
#   group_by(Species, functionalclass) %>%
#   summarise(Total_Abundance = sum(Rel_Abundance),
#             .groups = "drop") %>%
#   arrange(desc(Total_Abundance)) %>%
#   head(10)
# 
# top_species_abund_detail <- abund_long_ep %>%
#   filter(Species %in% top_species_abund$Species) %>%
#   group_by(Species, functionalclass, Severity) %>%
#   summarise(
#     Mean_Abundance = mean(Rel_Abundance),
#     SE = sd(Rel_Abundance) / sqrt(n()),
#     .groups = "drop"
#   )

# p_top_species <- ggplot(top_species_abund_detail,
#                         aes(x = Severity, y = Mean_Abundance, fill = Severity)) +
#   geom_col(width = 0.7) +
#   geom_errorbar(aes(ymin = Mean_Abundance - SE,
#                     ymax = Mean_Abundance + SE),
#                 width = 0.2) +
#   facet_wrap(~ Species + functionalclass,
#              scales = "free_y", ncol = 3) +
#   scale_fill_manual(values = c("Low" = "#66C2A5",
#                                "Moderate" = "#FC8D62",
#                                "High" = "#E78AC3")) +
#   labs(title = "Top Species Abundance by Disease Severity",
#        subtitle = "Top 10 endophytes and plant pathogens",
#        x = "Disease Severity",
#        y = "Mean Relative Abundance (%)") +
#   theme_bw(base_size = 10) +
#   theme(legend.position = "none",
#         strip.background = element_rect(fill = "white"),
#         plot.title = element_text(face = "bold", size = 14))
# 
# ggsave("key_taxa_results/fig_top_species_by_severity.png",
#        p_top_species,
#        width = 12,
#        height = max(8, ceiling(nrow(top_species_abund) / 3) * 3),
#        dpi = 300)

## 5.5 Prevalence analysis (high vs low severity) + heatmap

# species_pa_ep <- (species_data_ep > 0) * 1
# 
# prevalence_data <- species_pa_ep %>%
#   as.data.frame() %>%
#   rownames_to_column("plot") %>%
#   pivot_longer(-plot, names_to = "Species", values_to = "Present") %>%
#   left_join(plot_metadata_ep, by = "plot") %>%
#   left_join(
#     data_filtered %>%
#       select(species_clean, functionalclass) %>%
#       distinct(),
#     by = c("Species" = "species_clean")
#   )
# 
# prevalence_summary <- prevalence_data %>%
#   group_by(Species, functionalclass, Severity) %>%
#   summarise(
#     Prevalence = mean(Present) * 100,
#     n_plots    = n(),
#     n_present  = sum(Present),
#     .groups    = "drop"
#   ) %>%
#   pivot_wider(
#     id_cols    = c(Species, functionalclass),
#     names_from = Severity,
#     values_from = c(Prevalence, n_present),
#     names_sep = "_"
#   )
# 
# write.csv(prevalence_summary,
#           "key_taxa_results/prevalence_by_severity.csv",
#           row.names = FALSE)
# 
# # Clean species names and collapse duplicates for heatmap fig 5
# prev_heat <- prevalence_summary %>%
#   mutate(Species = str_trim(Species)) %>%
#   group_by(Species) %>%
#   summarise(
#     Low      = max(Prevalence_Low,      na.rm = TRUE),
#     Moderate = max(Prevalence_Moderate, na.rm = TRUE),
#     High     = max(Prevalence_High,     na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# prevalence_matrix <- prev_heat %>%
#   column_to_rownames("Species") %>%
#   as.matrix()
# 
# 
# func_annotation <- data_filtered %>%
#   select(species_clean, functionalclass) %>%
#   distinct() %>%
#   filter(species_clean %in% rownames(prevalence_matrix)) %>%
#   column_to_rownames("species_clean")
# png("fungal_diversity_results1/fig_prevalence_heatmap.png",
#     width = 10,
#     height = max(6, nrow(prevalence_matrix) * 0.3),
#     units = "in", res = 300)
# 
# pheatmap(
#   prevalence_matrix,
#   cluster_cols = FALSE,
#   cluster_rows = TRUE,
#   color = colorRampPalette(c("white", "#E78AC3"))(100),
#   display_numbers = TRUE,
#   number_format = "%.0f",
#   fontsize_number = 9,
#   angle_col = 0,
#   cellwidth = 120,
#   cellheight = 18,
#   annotation_row = func_annotation,
#   annotation_colors = list(
#     functionalclass = c("Endophytes" = "#8DD3C7",
#                         "Plant Pathogens" = "#FB8072",
#                         "Saprotrophs"       = "#635456" )
#   )
# )
# 
# dev.off()


#fig 6

#----------------------------------#
# 2. Build prevalence matrix
#----------------------------------#
prev_heat <- prevalence_summary %>%
  mutate(Species = clean_species_name(Species)) %>%
  group_by(Species) %>%
  summarise(
    Low      = if (all(is.na(Prevalence_Low))) NA_real_ else max(Prevalence_Low, na.rm = TRUE),
    Moderate = if (all(is.na(Prevalence_Moderate))) NA_real_ else max(Prevalence_Moderate, na.rm = TRUE),
    High     = if (all(is.na(Prevalence_High))) NA_real_ else max(Prevalence_High, na.rm = TRUE),
    .groups = "drop"
  )

prevalence_matrix <- prev_heat %>%
  distinct(Species, .keep_all = TRUE) %>%
  column_to_rownames("Species") %>%
  as.matrix()

#----------------------------------#
# 3. Single functional-class annotation
#    (one class per species)
#----------------------------------#
func_annotation <- data %>%
  mutate(
    Species = clean_species_name(species_clean),
    functionalclass = str_squish(functionalclass)
  ) %>%
  filter(!is.na(functionalclass)) %>%
  filter(Species %in% rownames(prevalence_matrix)) %>%
  distinct(Species, .keep_all = TRUE) %>%
  select(Species, functionalclass) %>%
  column_to_rownames("Species")

func_annotation <- func_annotation[rownames(prevalence_matrix), , drop = FALSE]

func_annotation$functionalclass <- factor(
  func_annotation$functionalclass,
  levels = c("Endophyte", "Plant Pathogen", "Saprotroph"))

#----------------------------------#
# 4. Draw heatmap
#----------------------------------#
png(
  "fungal_diversity_results1/fig_prevalence_heatmap.png",
  width = 10,
  height = max(7, nrow(prevalence_matrix) * 0.28),
  units = "in",
  res = 600
)

ph <- pheatmap(
  prevalence_matrix,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  treeheight_row = 70,
  color = colorRampPalette(c("white", "#E78AC3"))(100),
  border_color = "grey65",
  display_numbers = TRUE,
  number_format = "%.0f",
  fontsize_number = 9,
  fontsize_row = 11,
  fontsize_col = 12,
  cellwidth = 42,
  cellheight = 18,
  angle_col = 0,
  annotation_row = func_annotation,
  annotation_colors = list(
    functionalclass = c(
      "Endophyte" = "#8DD3C7",
      "Plant Pathogen" = "#FB8072",
      "Saprotroph" = "#B3B3B3"
    )
  ),
  annotation_names_row = FALSE,   # removes overlapping annotation label text
  annotation_legend = TRUE,      # removes side legend
  legend = TRUE,                 # removes heatmap scale legend
  silent = TRUE
)

# Italicize species names
row_name_idx <- which(ph$gtable$layout$name == "row_names")
ph$gtable$grobs[[row_name_idx]]$gp <- gpar(
  fontsize = 11,
  fontface = "italic"
)

grid.newpage()
grid.draw(ph$gtable)

dev.off()



# --- Key taxa summary (high severity + core taxa) -------------------------

dir.create("key_taxa_results", showWarnings = FALSE)

key_taxa_all <- tibble::tibble(
  Species          = character(),
  Functional_Class = character(),
  Evidence         = character()
)

## 1. Indicator species for HIGH severity (IndVal)
if (exists("indval_sig") && nrow(indval_sig) > 0) {
  high_indicators <- indval_sig %>%
    filter(Associated_Severity == "High", p.value < 0.05) %>%
    transmute(
      Species,
      Functional_Class = functionalclass,
      Evidence         = "Indicator (High)"
    )
  
  key_taxa_all <- bind_rows(key_taxa_all, high_indicators)
}

## 2. Taxa biased toward HIGH severity (prevalence) not needed
high_prev_taxa <- prevalence_summary %>%
  filter(
    !is.na(Prevalence_High),
    Prevalence_High > 50,
    is.na(Prevalence_Low) | Prevalence_Low < 30
  ) %>%
  transmute(
    Species,
    Functional_Class = functionalclass,
    Evidence         = "High prevalence"
  )
high_prev_taxa
key_taxa_all <- bind_rows(key_taxa_all, high_prev_taxa)

## 3. Core dominant taxa (high overall prevalence in any severity) dont need this again
core_taxa <- prevalence_summary %>%
  mutate(
    Overall_Prevalence = pmax(
      Prevalence_Low,
      Prevalence_Moderate,
      Prevalence_High,
      na.rm = TRUE
    )
  ) %>%
  filter(Overall_Prevalence > 50) %>%
  transmute(
    Species,
    Functional_Class = functionalclass,
    Evidence         = "Core dominant"
  )

key_taxa_all <- bind_rows(key_taxa_all, core_taxa)

## 4. Collapse to one row per species × functional class
key_taxa_summary <- key_taxa_all %>%
  group_by(Species, Functional_Class) %>%
  summarise(
    Evidence_Types = paste(sort(unique(Evidence)), collapse = "; "),
    Evidence_Count = n_distinct(Evidence),
    .groups = "drop"
  ) %>%
  arrange(desc(Evidence_Count), Species)

# Save table
if (nrow(key_taxa_summary) > 0) {
  write.csv(
    key_taxa_summary,
    "key_taxa_results/KEY_TAXA_SUMMARY.csv",
    row.names = FALSE
  )
}




## 5.7 Sydowia polyspora – grouping isolates for FST analyses

sydowia <- data %>%
  filter(speciesidentity == "Sydowia polyspora") %>%
  mutate(
    region = str_to_title(region),
    region = if_else(region %in% c("North", "South"), region, NA_character_)
  ) %>%
  filter(!is.na(region), !is.na(Severity)) %>%
  mutate(fst_group = paste(region, Severity, sep = "_"))
sydowia

# 1. Get all isolates per region (North/South only)
region_totals <- data %>%
  mutate(region = str_to_title(region)) %>%
  filter(region %in% c("North", "South")) %>%
  count(region, name = "total_isolates")

## 2. Get Sydowia counts per region (you already have sydowia)
syd_counts <- sydowia %>%
  count(region, name = "sydowia_isolates")

## 3. Join and compute percentages
syd_region_pct <- syd_counts %>%
  left_join(region_totals, by = "region") %>%
  mutate(perc_of_region = 100 * sydowia_isolates / total_isolates)

syd_region_pct

syd_within_pct <- sydowia %>%
  count(region, name = "sydowia_isolates") %>%
  mutate(perc_within_syd = 100 * sydowia_isolates / sum(sydowia_isolates))

syd_within_pct

pop_dir <- "fungal_diversity_results1/fst_pops"
dir.create(pop_dir, showWarnings = FALSE, recursive = TRUE)

pop_list <- split(sydowia$isolates, sydowia$fst_group)
for (grp in names(pop_list)) {
  fname <- file.path(pop_dir, paste0("pop_", grp, ".txt"))
  readr::write_lines(pop_list[[grp]], fname)
}

