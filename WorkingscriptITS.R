# ============================================================================
# COMPREHENSIVE FUNGAL PATHOBIOME DIVERSITY ANALYSIS
# Testing Amplification and Composition Effects on BSNB Disease Severity
# ============================================================================
# ENHANCED FUNGAL PATHOBIOME DIVERSITY ANALYSIS
# With Correlation Analysis, Combined Plots, and Pie Charts
# ============================================================================

library(tidyverse)
library(vegan)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(RColorBrewer)
library(scales)
library(reshape2)
library(pheatmap)
library(ggpubr)  # For correlation plots
library(tidytext)
library(ggplot2)

library(dplyr)
library(stringr)
library(tidyr)

library(ggrepel)
library(indicspecies)
# Load data
data <- read.csv("data3.csv")
data <- data %>% janitor::clean_names()

# Create output folder
dir.create("fungal_diversity_results1", showWarnings = FALSE)

# ============================================================================
# DATA PREPARATION
# ============================================================================

# Create disease severity categories
data <- data %>%
  mutate(
    Severity = case_when(
      proportion >= 30 ~ "High",
      proportion >= 10 & proportion < 30 ~ "Moderate",
      proportion < 10 ~ "Low",
      TRUE ~ "Unknown"
    ),
    Severity = factor(Severity, levels = c("Low", "Moderate", "High"))
  )

# Clean up taxonomic data
data <- data %>%
  mutate(
    species_clean = if_else(is.na(speciesidentity) | speciesidentity %in% c("NA",""),
                            paste0("Unidentified_", coalesce(genus, "sp")), speciesidentity),
    order_clean   = if_else(is.na(order) | order %in% c("NA",""),
                            "Unknown_Order", order),
    lineage_clean = if_else(is.na(lineage_class) | lineage_class %in% c("NA",""),
                            "Unknown_Class", lineage_class),
    genus_clean = if_else(is.na(genus) | genus %in% c("NA",""),
                          "Unknown_Genus", genus)
  )

cat("‚úì Data loaded and prepared\n")
cat(sprintf("  Total samples: %d\n", nrow(data)))
cat(sprintf("  Unique species: %d\n", length(unique(data$species_clean))))
cat(sprintf("  Unique genera: %d\n", length(unique(data$genus_clean))))
cat(sprintf("  Unique orders: %d\n", length(unique(data$order_clean))))
cat(sprintf("  Unique lineages/classes: %d\n", length(unique(data$lineage_clean))))
cat(sprintf("  Unique functionalclass: %d\n",
            length(unique(data$functionalclass))))
unique(data$functionalclass)

#Amplification effect analysis

#Figure 1
# Lineage composition by severity
composition_severity <- data %>%
  group_by(Severity, lineage_clean) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Severity) %>%
  mutate(proportion = count / sum(count) * 100)

# Plot: Composition by severity
p9_composition <- ggplot(composition_severity, aes(x = Severity, y = proportion, fill = lineage_clean)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(composition_severity$lineage_clean)))) +
  labs(title = "Fungal Lineage/Class Composition by Disease Severity",
       subtitle = "Testing Composition Effect",
       x = "Disease Severity", 
       y = "Proportion (%)", 
       fill = "Lineage/Class") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10)
  )

ggsave("fungal_diversity_results1/fig9_composition_severity.png", 
       p9_composition, width = 12, height = 7, dpi = 300)
p9_composition
cat("‚úì Composition plots created\n")

allowed_classes <- c("Plant Pathogens", "Endophytes")

filtered <- data %>%
  mutate(functionalclass = str_trim(functionalclass)) %>%
  filter(!is.na(functionalclass),
         functionalclass %in% allowed_classes) %>%
  filter(!is.na(species_clean), species_clean != "", species_clean != "NA")

# --- pick top 5 species among those two classes --------------------------
top_species <- filtered %>%
  count(species_clean, sort = TRUE) %>%
  slice_max(n, n = 5) %>%
  pull(species_clean)
top_species

# --- composition and proportions (by Severity and functional class) -------
species_composition <- filtered %>%
  mutate(species_grouped = ifelse(species_clean %in% top_species,
                                  species_clean, "Other")) %>%
  group_by(functionalclass, Severity, species_grouped) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(proportion = 100 * count / sum(count)) %>%
  ungroup() %>%
  filter(species_grouped != "Other")     # <—— removes “Other” completely
species_composition
print(species_composition)

# --- plot ------------------------------------------------------------------
p11_species <- ggplot(species_composition,
                      aes(x = Severity, y = proportion, fill = species_grouped)) +
  geom_col(width = 0.7) +
  labs(
    title = "Top Species by Disease Severity (Endophytes & Plant Pathogens)",
    x = "Disease Severity",
    y = "Proportion (%)",
    fill = "Species"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        legend.text = element_text(size = 10))

ggsave("fungal_diversity_results1/fig12_top_species_composition.png",
       p11_species, width = 10, height = 7, dpi = 300)
p11_species

p_functional_taxa <- p11_species | p9_composition
ggsave("fungal_diversity_results1/Fig1a_funnctional_species_region.png", 
       p_functional_taxa, width = 14, height = 5, dpi = 300)
p_functional_taxa


#Figure 2
# 2. COMBINED RICHNESS PLOTS
# ============================================================================

cat("\n=== CREATING COMBINED RICHNESS PLOTS ===\n")
# 1. Global totals (unique across *all* plots / severities)
total_species  <- n_distinct(data$species_clean,  na.rm = TRUE)
total_lineages <- n_distinct(data$lineage_clean, na.rm = TRUE)

# 2. Richness by Severity, then convert to % of global richness
richness_severity <- data %>%
  # Drop NA severity if you don't want it in the table
  filter(!is.na(Severity)) %>%
  group_by(Severity) %>% 
  summarise(
    Species = n_distinct(species_clean),
    Lineage = n_distinct(lineage_clean),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Species, Lineage),
    names_to  = "Taxonomic_Level",
    values_to = "Richness"
  ) %>%
  mutate(
    Total = dplyr::case_when(
      Taxonomic_Level == "Species" ~ total_species,
      Taxonomic_Level == "Lineage" ~ total_lineages
    ),
    Percent = (Richness / Total) * 100
  ) %>%
  select(Severity, Taxonomic_Level, Richness, Percent)

richness_severity


# Combined plot: All richness levels by severity
p1_combined <- ggplot(richness_severity, aes(x = Severity, y = Percent, fill = Taxonomic_Level)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Fungal Richness by Disease Severity",
       subtitle = "Testing Amplification Effect Across Taxonomic Levels",
       x = "Disease Severity", 
       y = "Number of Taxa",
       fill = "Taxonomic Level") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

ggsave("fungal_diversity_results1/fig1_combined_richness_severity.png", 
       p1_combined, width = 12, height = 6, dpi = 300)

# Richness by region

richness_region <- data %>%
  group_by(region) %>%
  summarise(
    Species = n_distinct(species_clean),
    Lineage = n_distinct(lineage_clean),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Species, Lineage),
    names_to  = "Taxonomic_Level",
    values_to = "Richness"
  ) %>%
  mutate(
    Total = dplyr::case_when(
      Taxonomic_Level == "Species" ~ total_species,
      Taxonomic_Level == "Lineage" ~ total_lineages
    ),
    Percent = (Richness / Total) * 100
  ) %>%
  select(region, Taxonomic_Level, Richness, Percent)

richness_region


# Combined plot: All richness levels by region
p2_combined <- ggplot(richness_region, aes(x = region, y = Percent, fill = Taxonomic_Level)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Fungal Richness by Region",
       x = "Region", 
       y = "Number of Taxa",
       fill = "Taxonomic Level") +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

ggsave("fungal_diversity_results1/fig2_combined_richness_region.png", 
       p2_combined, width = 12, height = 6, dpi = 300)
p2_combined
cat("‚úì Combined richness plots created\n")

# Combined genus + species plot
p_combined_taxa <- p1_combined | p2_combined
ggsave("fungal_diversity_results1/fig5_combined_genus_species_region.png", 
       p_combined_taxa, width = 14, height = 5, dpi = 300)
p_combined_taxa
cat("‚úì Genus and species percentage plots created\n")


# compositional effect analysis

# "HYPOTHESIS: Disease severity is associated with fungal community composition\n\n"
cat("UNIT OF ANALYSIS: PLOT (each plot = one tree)\n\n")


# Species abundance matrix (plot × species)
species_matrix <- data %>%
  group_by(plot, species_clean) %>%
  summarise(abundance = n(), .groups = "drop") %>%
  pivot_wider(names_from = species_clean, 
              values_from = abundance, 
              values_fill = 0)

# Plot metadata
plot_metadata <- data %>%
  distinct(plot, .keep_all = TRUE) %>%
  select(plot, Severity, region, county, proportion) %>%
  filter(!is.na(Severity), !is.na(region))

# Prepare matrix
species_data <- species_matrix %>%
  filter(plot %in% plot_metadata$plot) %>%
  column_to_rownames("plot") %>%
  as.matrix()

# Match metadata
plot_metadata <- plot_metadata %>%
  filter(plot %in% rownames(species_data)) %>%
  arrange(match(plot, rownames(species_data)))

cat(sprintf("✓ Community matrix created: %d plots × %d species\n", 
            nrow(species_data), ncol(species_data)))
cat(sprintf("  Total plots analyzed: %d\n", nrow(plot_metadata)))
cat(sprintf("  North region: %d plots\n", sum(plot_metadata$region == "North")))
cat(sprintf("  South region: %d plots\n", sum(plot_metadata$region == "South")))
cat("\n")

# ============================================================================
# 3. PCA ORDINATION
# ============================================================================

cat("=== PCA ORDINATION ===\n")

# Hellinger transformation (recommended for species data)
species_hell <- decostand(species_data, method = "hellinger")

# Run PCA
comm_pca <- rda(species_hell)

# Extract scores
pca_scores <- as.data.frame(scores(comm_pca, display = "sites")) %>%
  rownames_to_column("plot") %>%
  left_join(plot_metadata, by = "plot")

# Variance explained
pca_var <- eigenvals(comm_pca) / sum(eigenvals(comm_pca)) * 100

cat(sprintf("✓ PCA completed\n"))
cat(sprintf("  PC1 explains %.1f%% of variance\n", pca_var[1]))
cat(sprintf("  PC2 explains %.1f%% of variance\n", pca_var[2]))
cat("\n")

# Save PCA scores
write.csv(pca_scores, "composition_results/pca_scores.csv", row.names = FALSE)

# ============================================================================
# 4. PCA PLOTS
# ============================================================================

cat("=== CREATING PCA PLOTS ===\n")

# Plot 1: PCA by Severity
p_pca_severity <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Severity, shape = region)) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(aes(color = Severity), level = 0.95, linewidth = 1) +
  scale_color_manual(values = c("Low" = "#66C2A5", 
                                "Moderate" = "#FC8D62", 
                                "High" = "#E78AC3")) +
  labs(title = "PCA: Fungal Community Composition by Disease Severity",
       subtitle = "95% confidence ellipses | Shape indicates region",
       x = sprintf("PC1 (%.1f%%)", pca_var[1]),
       y = sprintf("PC2 (%.1f%%)", pca_var[2])) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 16))
p_pca_severity

ggsave("composition_results/fig1_pca_by_severity.png", 
       p_pca_severity, width = 12, height = 8, dpi = 300)

# Plot 2: PCA by Region
p_pca_region <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = region, shape = Severity)) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(aes(color = region), level = 0.95, linewidth = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PCA: Fungal Community Composition by Region",
       subtitle = "95% confidence ellipses | Shape indicates severity",
       x = sprintf("PC1 (%.1f%%)", pca_var[1]),
       y = sprintf("PC2 (%.1f%%)", pca_var[2])) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 16))
p_pca_region
ggsave("composition_results/fig2_pca_by_region.png", 
       p_pca_region, width = 12, height = 8, dpi = 300)

# Combined plot
p_pca_combined <- p_pca_severity / p_pca_region
ggsave("composition_results/fig3_pca_combined.png", 
       p_pca_combined, width = 12, height = 14, dpi = 300)

cat("✓ PCA plots created\n\n")

# ============================================================================
# 5. PERMANOVA: FULL DATASET WITH REGION CONTROL
# ============================================================================

cat(strrep("=", 80) %+% "\n")
cat("PERMANOVA: TESTING COMPOSITION DIFFERENCES\n")
cat(strrep("=", 80) %+% "\n\n")

# Calculate Bray-Curtis distance
bray_dist <- vegdist(species_data, method = "bray")

# PERMANOVA controlling for region
cat("--- PERMANOVA: Region + Severity (Full Dataset) ---\n\n")

perm_full <- adonis2(bray_dist ~ region + Severity, 
                     data = plot_metadata,
                     permutations = 999,
                     method = "bray",
                     by = "margin")  # Type III SS

print(perm_full)
cat("\n")

# Extract results
region_r2 <- perm_full$R2[1]
region_p <- perm_full$`Pr(>F)`[1]
severity_r2 <- perm_full$R2[2]
severity_p <- perm_full$`Pr(>F)`[2]

cat("INTERPRETATION:\n")
cat(sprintf("  Region effect: R² = %.3f, p = %.4f %s\n", 
            region_r2, region_p,
            ifelse(region_p < 0.05, "✓ SIGNIFICANT", "✗ NOT SIGNIFICANT")))
cat(sprintf("  Severity effect (controlling for region): R² = %.3f, p = %.4f %s\n", 
            severity_r2, severity_p,
            ifelse(severity_p < 0.05, "✓ SIGNIFICANT", "✗ NOT SIGNIFICANT")))
cat("\n")

# Test dispersion homogeneity
cat("--- Testing Dispersion Homogeneity ---\n")
disp_severity <- betadisper(bray_dist, plot_metadata$Severity)
disp_test <- permutest(disp_severity, permutations = 999)

cat(sprintf("  Dispersion test p-value: %.4f\n", disp_test$tab$`Pr(>F)`[1]))
cat(sprintf("  Dispersion: %s\n", 
            ifelse(disp_test$tab$`Pr(>F)`[1] < 0.05, 
                   "⚠ Heterogeneous (interpret with caution)", 
                   "✓ Homogeneous")))
cat("\n")

# ============================================================================
# 6. REGION-STRATIFIED PERMANOVA
# ============================================================================

cat(strrep("=", 80) %+% "\n")
cat("REGION-STRATIFIED ANALYSIS\n")
cat(strrep("=", 80) %+% "\n\n")

# North region
north_plots <- plot_metadata %>% filter(region == "North") %>% pull(plot)
if(length(north_plots) >= 10) {
  species_north <- species_data[north_plots, , drop = FALSE]
  species_north <- species_north[, colSums(species_north) > 0, drop = FALSE]
  metadata_north <- plot_metadata %>% filter(region == "North")
  
  bray_north <- vegdist(species_north, method = "bray")
  
  cat("--- NORTH REGION ONLY ---\n")
  cat(sprintf("  Plots: %d\n\n", nrow(species_north)))
  
  perm_north <- adonis2(bray_north ~ Severity, 
                        data = metadata_north,
                        permutations = 999,
                        method = "bray")
  
  print(perm_north)
  cat("\n")
  
  north_r2 <- perm_north$R2[1]
  north_p <- perm_north$`Pr(>F)`[1]
  
  cat(sprintf("  Severity effect: R² = %.3f, p = %.4f %s\n", 
              north_r2, north_p,
              ifelse(north_p < 0.05, "✓ SIGNIFICANT", "✗ NOT SIGNIFICANT")))
  cat("\n")
} else {
  cat("--- NORTH REGION ---\n")
  cat("  Insufficient sample size (< 10 plots)\n\n")
  north_r2 <- NA
  north_p <- NA
}

# South region
south_plots <- plot_metadata %>% filter(region == "South") %>% pull(plot)
if(length(south_plots) >= 10) {
  species_south <- species_data[south_plots, , drop = FALSE]
  species_south <- species_south[, colSums(species_south) > 0, drop = FALSE]
  metadata_south <- plot_metadata %>% filter(region == "South")
  
  bray_south <- vegdist(species_south, method = "bray")
  
  cat("--- SOUTH REGION ONLY ---\n")
  cat(sprintf("  Plots: %d\n\n", nrow(species_south)))
  
  perm_south <- adonis2(bray_south ~ Severity, 
                        data = metadata_south,
                        permutations = 999,
                        method = "bray")
  
  print(perm_south)
  cat("\n")
  
  south_r2 <- perm_south$R2[1]
  south_p <- perm_south$`Pr(>F)`[1]
  
  cat(sprintf("  Severity effect: R² = %.3f, p = %.4f %s\n", 
              south_r2, south_p,
              ifelse(south_p < 0.05, "✓ SIGNIFICANT", "✗ NOT SIGNIFICANT")))
  cat("\n")
} else {
  cat("--- SOUTH REGION ---\n")
  cat("  Insufficient sample size (< 10 plots)\n\n")
  south_r2 <- NA
  south_p <- NA
}







#Fig 3

# Functional composition by functional class
functional_composition <- data %>%
  filter(!is.na(functionalclass) & functionalclass != "" & functionalclass != "NA") %>%
  group_by(Severity, functionalclass) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Severity) %>%
  mutate(proportion = count / sum(count) * 100)
functional_composition
# Plot: Functional composition
p11_functional <- ggplot(functional_composition, aes(x = Severity, y = proportion, fill = functionalclass)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(
    title = "Functional Composition by Disease Severity",
    subtitle = "Distribution of functional classes across disease severity levels",
    x = "Disease Severity",
    y = "Proportion (%)",
    fill = "Functional Class"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Save figure
ggsave(
  "fungal_diversity_results1/fig11_functional_composition.png",
  p11_functional, width = 10, height = 7, dpi = 300
)
p11_functional

# Top 10 species by region
species_region <- data %>%
  count(region, species_clean, name = "count") %>%      # counts per region × species
  group_by(region) %>%                                  # <- key: group by region
  mutate(percentage = 100 * count / sum(count)) %>%     # % of each region total
  arrange(region, desc(percentage)) %>%
  slice_head(n = 10)

species_region


# Key taxa analysis

# 3. CREATE COMMUNITY MATRIX (ENDOPHYTES + PATHOGENS ONLY)
# ============================================================================

cat("=== CREATING COMMUNITY MATRIX ===\n")

# Species matrix
species_matrix <- data_filtered %>%
  group_by(plot, species_clean) %>%
  summarise(abundance = n(), .groups = "drop") %>%
  pivot_wider(names_from = species_clean, 
              values_from = abundance, 
              values_fill = 0)

# Plot metadata
plot_metadata <- data_filtered %>%
  distinct(plot, .keep_all = TRUE) %>%
  select(plot, Severity, region, county, proportion) %>%
  filter(!is.na(Severity))

# Prepare matrix
species_data <- species_matrix %>%
  filter(plot %in% plot_metadata$plot) %>%
  column_to_rownames("plot") %>%
  as.matrix()

# Match metadata
plot_metadata <- plot_metadata %>%
  filter(plot %in% rownames(species_data)) %>%
  arrange(match(plot, rownames(species_data)))

cat(sprintf("✓ Matrix: %d plots × %d species (endophytes + pathogens)\n", 
            nrow(species_data), ncol(species_data)))
cat("\n")

# ============================================================================
# 4. INDICATOR SPECIES ANALYSIS
# ============================================================================

cat(strrep("=", 80) %+% "\n")
cat("INDICATOR SPECIES ANALYSIS\n")
cat(strrep("=", 80) %+% "\n\n")

cat("Identifying endophytes and pathogens associated with severity levels\n\n")

# Run indicator species analysis
indval_result <- multipatt(species_data, 
                           plot_metadata$Severity,
                           control = how(nperm = 999))

cat("--- Indicator Species Results ---\n\n")
summary(indval_result, alpha = 1,indvalcomp = TRUE)
cat("\n")
which(colnames(species_data) == "Sydowia polyspora")
indval_result$sign[ "Sydowia polyspora", ]
# Extract significant indicators
indval_sig <- indval_result$sign %>%
  rownames_to_column("Species") %>%
  filter(p.value < 0.05) %>%
  arrange(p.value)

if(nrow(indval_sig) > 0) {
  cat(sprintf("✓ Found %d significant indicator species (p < 0.05)\n\n", 
              nrow(indval_sig)))
  
  # Add functional class info
  indval_sig <- indval_sig %>%
    left_join(
      data_filtered %>% 
        select(species_clean, functionalclass) %>% 
        distinct(),
      by = c("Species" = "species_clean")
    ) %>%
    mutate(
      Associated_Severity = case_when(
        s.Low == 1 ~ "Low",
        s.Moderate == 1 ~ "Moderate",
        s.High == 1 ~ "High",
        TRUE ~ "Multiple"
      )
    )
  
  print(indval_sig %>% select(Species, functionalclass, Associated_Severity, stat, p.value))
  
  write.csv(indval_sig, 
            "key_taxa_results/indicator_species_endophytes_pathogens.csv", 
            row.names = FALSE)
  
  # Plot indicator species
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
  
  ggsave("key_taxa_results/fig1_indicator_species.png", 
         p_indval, 
         width = 10, 
         height = max(6, nrow(indval_sig) * 0.4), 
         dpi = 300)
  
  cat("\n✓ Indicator species plot created\n")
  
} else {
  cat("⚠ No significant indicator species found\n")
}

cat("\n")

# ============================================================================
# 5. ABUNDANCE BY SEVERITY AND FUNCTIONAL CLASS
# ============================================================================

cat(strrep("=", 80) %+% "\n")
cat("ABUNDANCE ANALYSIS BY FUNCTIONAL CLASS\n")
cat(strrep("=", 80) %+% "\n\n")

# Calculate relative abundance per plot
species_rel_abund <- species_data / rowSums(species_data) * 100

# Convert to long format
abund_long <- species_rel_abund %>%
  as.data.frame() %>%
  rownames_to_column("plot") %>%
  pivot_longer(-plot, names_to = "Species", values_to = "Rel_Abundance") %>%
  left_join(plot_metadata, by = "plot") %>%
  left_join(
    data_filtered %>% 
      select(species_clean, functionalclass) %>% 
      distinct(),
    by = c("Species" = "species_clean")
  )

# Summary by functional class and severity
func_summary <- abund_long %>%
  group_by(functionalclass, Severity) %>%
  summarise(
    Mean_Abundance = mean(Rel_Abundance, na.rm = TRUE),
    SE = sd(Rel_Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

cat("Mean relative abundance by functional class and severity:\n")
print(func_summary)
cat("\n")

# Plot functional class abundance
p_func_abund <- ggplot(func_summary, 
                       aes(x = Severity, y = Mean_Abundance, fill = functionalclass)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = Mean_Abundance - SE, ymax = Mean_Abundance + SE),
                position = position_dodge(0.7), width = 0.2) +
  scale_fill_manual(values = c("Endophytes" = "#8DD3C7", 
                               "Plant Pathogens" = "#FB8072")) +
  labs(title = "Functional Class Abundance by Disease Severity",
       x = "Disease Severity", 
       y = "Mean Relative Abundance (%)",
       fill = "Functional Class") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16))
p_func_abund
ggsave("key_taxa_results/fig2_functional_class_abundance.png", 
       p_func_abund, width = 10, height = 7, dpi = 300)

cat("✓ Functional class abundance plot created\n\n")

# ============================================================================
# 6. TOP SPECIES ANALYSIS
# ============================================================================

cat(strrep("=", 80) %+% "\n")
cat("TOP SPECIES ANALYSIS\n")
cat(strrep("=", 80) %+% "\n\n")

# Identify top species by abundance
top_species <- abund_long %>%
  group_by(Species, functionalclass) %>%
  summarise(Total_Abundance = sum(Rel_Abundance), .groups = "drop") %>%
  arrange(desc(Total_Abundance)) %>%
  head(10)

cat("Top 10 species (endophytes + plant pathogens):\n")
print(top_species)
cat("\n")

# Abundance of top species by severity
top_species_abund <- abund_long %>%
  filter(Species %in% top_species$Species) %>%
  group_by(Species, functionalclass, Severity) %>%
  summarise(
    Mean_Abundance = mean(Rel_Abundance),
    SE = sd(Rel_Abundance) / sqrt(n()),
    .groups = "drop"
  )

# Plot top species
p_top_species <- ggplot(top_species_abund, 
                        aes(x = Severity, y = Mean_Abundance, fill = Severity)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = Mean_Abundance - SE, ymax = Mean_Abundance + SE),
                width = 0.2) +
  facet_wrap(~ Species + functionalclass, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Low" = "#66C2A5", 
                               "Moderate" = "#FC8D62", 
                               "High" = "#E78AC3")) +
  labs(title = "Top Species Abundance by Disease Severity",
       subtitle = "Top 10 endophytes and plant pathogens",
       x = "Disease Severity", 
       y = "Mean Relative Abundance (%)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"),
        plot.title = element_text(face = "bold", size = 14))
p_top_species
ggsave("key_taxa_results/fig3_top_species_by_severity.png", 
       p_top_species, 
       width = 12, 
       height = max(8, ceiling(nrow(top_species)/3) * 3), 
       dpi = 300)

# ============================================================================
# 7. PREVALENCE ANALYSIS (HIGH VS LOW SEVERITY)
# ============================================================================

cat(strrep("=", 80) %+% "\n")
cat("PREVALENCE ANALYSIS: HIGH vs LOW SEVERITY\n")
cat(strrep("=", 80) %+% "\n\n")

# Presence/absence matrix
species_pa <- (species_data > 0) * 1

# Calculate prevalence
prevalence_data <- species_pa %>%
  as.data.frame() %>%
  rownames_to_column("plot") %>%
  pivot_longer(-plot, names_to = "Species", values_to = "Present") %>%
  left_join(plot_metadata, by = "plot") %>%
  left_join(
    data_filtered %>% 
      select(species_clean, functionalclass) %>% 
      distinct(),
    by = c("Species" = "species_clean")
  )

# Prevalence by severity
prevalence_summary <- prevalence_data %>%
  group_by(Species, functionalclass, Severity) %>%
  summarise(
    Prevalence = mean(Present) * 100,
    n_plots = n(),
    n_present = sum(Present),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Severity, 
              values_from = c(Prevalence, n_present),
              names_sep = "_")

prevalence_summary <- prevalence_data %>%
  group_by(Species, functionalclass, Severity) %>%
  summarise(
    Prevalence = mean(Present) * 100,
    n_plots    = n(),           # or n_distinct(plot) if appropriate
    n_present  = sum(Present),
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols    = c(Species, functionalclass),   # <– important
    names_from = Severity,
    values_from = c(Prevalence, n_present),
    names_sep = "_"
  )


cat("Prevalence by severity level:\n")
print(prevalence_summary %>% 
        select(Species, functionalclass, Prevalence_Low, Prevalence_Moderate, Prevalence_High))
cat("\n")

write.csv(prevalence_summary, 
          "key_taxa_results/prevalence_by_severity.csv", 
          row.names = FALSE)

# Create heatmap of prevalence
prevalence_matrix <- prevalence_summary %>%
  select(Species, Prevalence_Low, Prevalence_Moderate, Prevalence_High) %>%
  column_to_rownames("Species") %>%
  as.matrix()
library(stringr)

# clean and collapse to 1 row/species
prev_heat <- prevalence_summary %>%
  mutate(Species = str_trim(Species)) %>%          # remove trailing spaces
  group_by(Species) %>%
  summarise(
    Prevalence_Low      = max(Prevalence_Low,      na.rm = TRUE),
    Prevalence_Moderate = max(Prevalence_Moderate, na.rm = TRUE),
    Prevalence_High     = max(Prevalence_High,     na.rm = TRUE),
    .groups = "drop"
  )

prevalence_matrix <- prev_heat %>%
  column_to_rownames("Species") %>%
  as.matrix()

# Add functional class annotation
func_annotation <- data_filtered %>%
  select(species_clean, functionalclass) %>%
  distinct() %>%
  filter(species_clean %in% rownames(prevalence_matrix)) %>%
  column_to_rownames("species_clean")

png("key_taxa_results/fig4_prevalence_heatmap.png", 
    width = 10, height = max(8, nrow(prevalence_matrix) * 0.3), 
    units = "in", res = 300)

pheatmap(prevalence_matrix,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("white", "#E78AC3"))(100),
         main = "Species Prevalence by Disease Severity",
         display_numbers = TRUE,
         number_format = "%.0f",
         fontsize_number = 9,
         angle_col = 0,
         cellwidth = 50,
         cellheight = 15,
         annotation_row = func_annotation,
         annotation_colors = list(
           functionalclass = c("Endophytes" = "#8DD3C7", 
                               "Plant Pathogens" = "#FB8072")
         ))

dev.off()

cat("✓ Prevalence heatmap created\n\n")

# ============================================================================
# 8. KEY TAXA SUMMARY
# ============================================================================

cat(strrep("=", 80) %+% "\n")
cat("KEY TAXA SUMMARY\n")
cat(strrep("=", 80) %+% "\n\n")

# Identify taxa associated with HIGH severity
high_severity_taxa <- data.frame(
  Species = character(),
  Functional_Class = character(),
  Evidence = character(),
  stringsAsFactors = FALSE
)

# From indicator species
if(exists("indval_sig") && nrow(indval_sig) > 0) {
  high_indicators <- indval_sig %>%
    filter(Associated_Severity == "High") %>%
    select(Species, functionalclass) %>%
    mutate(Evidence = "Indicator Species")
  
  high_severity_taxa <- bind_rows(high_severity_taxa, 
                                  high_indicators %>% 
                                    rename(Functional_Class = functionalclass))
}

# From prevalence (>50% in High, <30% in Low)
high_prev_taxa <- prevalence_summary %>%
  filter(Prevalence_High > 50, Prevalence_Low < 30) %>%
  select(Species, functionalclass) %>%
  mutate(Evidence = "High Prevalence")
high_prev_taxa

high_severity_taxa <- bind_rows(high_severity_taxa, 
                                high_prev_taxa %>% 
                                  rename(Functional_Class = functionalclass))

high_severity_taxa 
# Remove duplicates and count evidence
key_taxa_summary <- high_severity_taxa %>%
  group_by(Species, Functional_Class) %>%
  summarise(
    Evidence_Types = paste(unique(Evidence), collapse = "; "),
    Evidence_Count = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(Evidence_Count), Species)

if(nrow(key_taxa_summary) > 0) {
  cat("Taxa associated with HIGH disease severity:\n\n")
  print(key_taxa_summary)
  cat("\n")
  
  write.csv(key_taxa_summary, 
            "key_taxa_results/KEY_TAXA_HIGH_SEVERITY.csv", 
            row.names = FALSE)
  
  cat(sprintf("✓ Identified %d key taxa associated with high severity\n", 
              nrow(key_taxa_summary)))
  
  # Breakdown by functional class
  cat("\nBreakdown by functional class:\n")
  func_breakdown <- key_taxa_summary %>%
    group_by(Functional_Class) %>%
    summarise(Count = n(), .groups = "drop")
  print(func_breakdown)
  
} else {
  cat("⚠ No taxa strongly associated with high severity identified\n")
}

cat("\n" %+% strrep("=", 80) %+% "\n")
cat("Analysis complete! Check 'key_taxa_results/' folder for outputs\n")
cat(strrep("=", 80) %+% "\n")






# CORRELATION ANALYSIS - DIVERSITY vs SEVERITY

cat("\n=== CORRELATION ANALYSIS: DIVERSITY vs SEVERITY ===\n")

# Calculate diversity metrics per plot
diversity_per_plot <- data %>%
  group_by(plot, proportion, Severity, region, county) %>%
  summarise(
    Species_Richness = n_distinct(species_clean),
    Genus_Richness = n_distinct(genus_clean),
    Order_Richness = n_distinct(order_clean),
    Lineage_Richness = n_distinct(lineage_clean),
    Total_Isolates = n(),
    .groups = "drop"
  )

# Calculate Shannon diversity per plot
shannon_per_plot <- data %>%
  group_by(plot, lineage_clean) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(plot) %>%
  summarise(
    Shannon_Diversity = diversity(count, index = "shannon"),
    Simpson_Diversity = diversity(count, index = "simpson"),
    .groups = "drop"
  )

# Merge diversity metrics
diversity_per_plot <- diversity_per_plot %>%
  left_join(shannon_per_plot, by = "plot")

# Save diversity per plot table
write.csv(diversity_per_plot, 
          "fungal_diversity_results1/diversity_metrics_per_plot.csv", 
          row.names = FALSE)

#Fig 4a
#OPTION 1: BOXPLOT BY SEVERITY CATEGORY
# Best for: Showing distribution without implying linear relationship
# ============================================================================

# Plot 1: Species Richness by Severity Category
p1_box <- ggplot(diversity_per_plot, aes(x = Severity, y = Species_Richness, fill = Severity)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  scale_fill_manual(values = c("Low" = "#66C2A5", "Moderate" = "#FC8D62", "High" = "#E78AC3")) +
  labs(title = "Species Richness ",
       subtitle = sprintf("Spearman rho = %.3f, p = %.4f ", 
                          cor_species$estimate, cor_species$p.value),
       x = "Disease Severity Category", 
       y = "Species Richness") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 16))
p1_box
ggsave("correlation_plots_alternative/fig1_species_richness_boxplot.png", 
       p1_box, width = 8, height = 6, dpi = 300)

#fig4b
# Plot 2: Shannon Diversity by Severity Category
p2_box <- ggplot(diversity_per_plot, aes(x = Severity, y = Shannon_Diversity, fill = Severity)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  scale_fill_manual(values = c("Low" = "#66C2A5", "Moderate" = "#FC8D62", "High" = "#E78AC3")) +
  labs(title = "Shannon Diversity",
       subtitle = sprintf("Spearman rho = %.3f, p = %.4f ", 
                          cor_shannon$estimate, cor_shannon$p.value),
       x = "Disease Severity Category", 
       y = "Shannon Diversity Index") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 16))
p2_box
ggsave("correlation_plots_alternative/fig2_shannon_diversity_boxplot.png", 
       p2_box, width = 8, height = 6, dpi = 300)

# Combined genus + species plot
p_correlation_disease <- p1_box | p2_box | p7_ridge
ggsave("correlation_plots_alternative/fig5_combined_genus_species_region.png", 
       p_correlation_disease, width = 14, height = 5, dpi = 300)
p_correlation_disease

#Fig4c
#  RIDGELINE PLOT
# Best for: Comparing distributions across severity levels
# ============================================================================

p7_ridge <- ggplot(diversity_per_plot, aes(x = Species_Richness, y = Severity, fill = Severity)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5) +
  scale_fill_manual(values = c("Low" = "#66C2A5", "Moderate" = "#FC8D62", "High" = "#E78AC3")) +
  labs(title = "Species Richness",
       x = "Species Richness", 
       y = "Disease Severity Category") +
  #theme_ridges(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 16))
p7_ridge
ggsave("correlation_plots_alternative/fig7_species_richness_ridgeline.png", 
       p7_ridge, width = 10, height = 6, dpi = 300)


#Barplot Alternative plot(not included yet)

# Calculate summary statistics
summary_stats <- diversity_per_plot %>%
  group_by(Severity) %>%
  summarise(
    Mean_Richness = mean(Species_Richness, na.rm = TRUE),
    SE_Richness = sd(Species_Richness, na.rm = TRUE) / sqrt(n()),
    Mean_Shannon = mean(Shannon_Diversity, na.rm = TRUE),
    SE_Shannon = sd(Shannon_Diversity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p8_bar <- ggplot(summary_stats, aes(x = Severity, y = Mean_Richness, fill = Severity)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_Richness - SE_Richness, 
                    ymax = Mean_Richness + SE_Richness), 
                width = 0.2, size = 0.8) +
  scale_fill_manual(values = c("Low" = "#66C2A5", "Moderate" = "#FC8D62", "High" = "#E78AC3")) +
  labs(title = "Mean Species Richness by Disease Severity",
       subtitle = "Error bars show standard error",
       x = "Disease Severity Category", 
       y = "Mean Species Richness") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 16))
p8_bar
ggsave("correlation_plots_alternative/fig8_species_richness_barplot.png", 
       p8_bar, width = 8, height = 6, dpi = 300)

p9_bar <- ggplot(summary_stats, aes(x = Severity, y = Mean_Shannon, fill = Severity)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_Shannon - SE_Shannon, 
                    ymax = Mean_Shannon + SE_Shannon), 
                width = 0.2, size = 0.8) +
  scale_fill_manual(values = c("Low" = "#66C2A5", "Moderate" = "#FC8D62", "High" = "#E78AC3")) +
  labs(title = "Mean Shannon Diversity by Disease Severity",
       subtitle = "Error bars show standard error",
       x = "Disease Severity Category", 
       y = "Mean Shannon Diversity Index") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 16))
p9_bar
ggsave("correlation_plots_alternative/fig9_shannon_diversity_barplot.png", 
       p9_bar, width = 8, height = 6, dpi = 300)

cat("✓ Bar plots with error bars created\n")

#fig 5

cat("\n=== ENVIRONMENTAL PCA ===\n")

# Select environmental variables (remove NAs)
env_vars <- data %>%
  select(isolates, plot, subplot, Severity, region, county,
         tree_ht, tree_dbh, basel_area, slope, elevation) %>%
  distinct(plot, .keep_all = TRUE) %>%  # One row per plot
  drop_na(tree_ht, tree_dbh, basel_area, elevation)

# Standardize environmental variables
env_matrix <- env_vars %>%
  select(tree_ht, tree_dbh, basel_area, slope, elevation) %>%
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
  geom_point(size = 4, alpha = 0.7) +
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

# PCA Biplot with loadings
loadings <- as.data.frame(env_pca$rotation[,1:2]) %>%
  mutate(Variable = rownames(.))

p_env_biplot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Severity), size = 3, alpha = 0.6) +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = PC1*5, yend = PC2*5),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black", size = 0.8) +
  geom_text_repel(data = loadings, 
                  aes(x = PC1*5, y = PC2*5, label = Variable),
                  size = 4, fontface = "bold") +
  scale_color_manual(values = c("Low" = "#66C2A5", "Moderate" = "#FC8D62", "High" = "#E78AC3")) +
  labs(title = "Environmental PCA Biplot",
       subtitle = "Arrows show environmental variable loadings",
       x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2])) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 16))
p_env_biplot
ggsave("ordination_results/fig2_environmental_pca_biplot.png", 
       p_env_biplot, width = 12, height = 8, dpi = 300)







#4. REGION COMPARISON
# ============================================================================
# Create species × plot matrix
species_matrix <- data %>%
  group_by(plot, species_clean) %>%
  summarise(abundance = n(), .groups = "drop") %>%
  pivot_wider(names_from = species_clean, values_from = abundance, values_fill = 0)

# Extract plot info
plot_info <- data %>%
  distinct(plot, .keep_all = TRUE) %>%
  select(plot, Severity, region, county, proportion)

# Prepare matrix for PCA
species_data <- species_matrix %>%
  select(-plot) %>%
  as.matrix()

rownames(species_data) <- species_matrix$plot

# Hellinger transformation (recommended for species data)
species_hell <- decostand(species_data, method = "hellinger")

# Run PCA
comm_pca <- rda(species_hell)

# Extract scores
comm_scores <- as.data.frame(scores(comm_pca, display = "sites")) %>%  
  mutate(plot = rownames(.)) %>%          # keep as character
  left_join(plot_info, by = "plot")

# Variance explained
comm_var <- eigenvals(comm_pca) / sum(eigenvals(comm_pca)) * 100

cat(sprintf("  PC1 explains %.1f%% of variance\n", comm_var[1]))
cat(sprintf("  PC2 explains %.1f%% of variance\n", comm_var[2]))
cat("\n=== REGIONAL COMPARISON ===\n")
# PCA by region
p_pca_region <- ggplot(comm_scores, aes(x = PC1, y = PC2, color = region)) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(aes(color = region), level = 0.95, size = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Fungal Community PCA by Region",
       subtitle = "95% confidence ellipses shown",
       x = sprintf("PC1 (%.1f%%)", comm_var[1]),
       y = sprintf("PC2 (%.1f%%)", comm_var[2])) +
  theme_classic (base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 16))
p_pca_region
ggsave("ordination_results/fig7_pca_by_region.png", 
       p_pca_region, width = 10, height = 7, dpi = 300)

# Combined genus + species plot
pCA_combined <- p_env_pca | p_pca_region 
ggsave("correlation_plots_alternative/fig5_PCA_species_region.png", 
       pCA_combined, width = 14, height = 5, dpi = 300)
pCA_combined






sydowia <- data %>%
  filter(speciesidentity == "Sydowia polyspora")

sydowia <- sydowia %>%
  mutate(
    region = str_to_title(region),          # "north" → "North"
    region = if_else(region %in% c("North", "South"), region, NA_character_)
  )


sydowia <- sydowia %>%
  filter(!is.na(region), !is.na(Severity)) %>%
  mutate(
    fst_group = paste(region, Severity, sep = "_")
  )

sydowia %>%
  count(fst_group)
# Make sure this really is the VCF sample name column:
sydowia %>% select(isolates, fst_group) %>% head()

# Directory to store the pop files
pop_dir <- "fungal_diversity_results1/fst_pops"
dir.create(pop_dir, showWarnings = FALSE, recursive = TRUE)

sydowia %>%
  group_by(fst_group) %>%
  group_walk(~ {
    fname <- file.path(pop_dir, paste0("pop_", unique(.x$fst_group), ".txt"))
    readr::write_lines(.x$isolates, fname)
  })


id_col <- "isolates"   # <-- CHANGE THIS to your real column name
pop_dir <- "fungal_diversity_results1/fst_pops"
dir.create(pop_dir, showWarnings = FALSE, recursive = TRUE)

sydowia %>%
  group_by(fst_group) %>%
  group_walk(~ {
    fname <- file.path(pop_dir, paste0("pop_", unique(.x$fst_group), ".txt"))
    # use [[ ]] so we can pass the column name as a string
    readr::write_lines(.x[[id_col]], fname)
  })

table(sydowia$fst_group)
head(sydowia[, c("isolates", "fst_group")])

library(readr)

pop_dir <- "fungal_diversity_results1/fst_pops"
dir.create(pop_dir, showWarnings = FALSE, recursive = TRUE)

## Split isolates into a list by fst_group
pop_list <- split(sydowia$isolates, sydowia$fst_group)

## Write one file per group
for (grp in names(pop_list)) {
  fname <- file.path(pop_dir, paste0("pop_", grp, ".txt"))
  write_lines(pop_list[[grp]], fname)
}

