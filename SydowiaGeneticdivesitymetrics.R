
setwd("/Users/temitopefolorunso/Desktop/Folorunso dissertation proposal April 2024/BrownSpotProject")
# Create output folder
dir.create("genetic_diversityresults", showWarnings = FALSE)

library(tidyverse)

# Load and reshape
assembly_file <- "Sydowiaassnew_stats_combined.txt"
assembly_data <- read.delim(assembly_file, header = TRUE)

assembly_long <- assembly_data %>%
  pivot_longer(cols = c(Total_Length, Scaffolds,Contigs,Contig_N50, GC_Content),
               names_to = "Metric", values_to = "Value")

# Compute mean and max
summary_stats <- assembly_long %>%
  group_by(Metric) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    Max = max(Value, na.rm = TRUE),
    SD = sd(Value, na.rm =TRUE)
  ) %>%
  pivot_longer(cols = c(Mean, Max, SD), names_to = "Statistic", values_to = "Value")

# Separate plots
plot1_data <- summary_stats %>%
  filter(Metric %in% c("Total_Length","Contig_N50"))

plot2_data <- summary_stats %>%
  filter(Metric %in% c( "Contigs", "GC_Content"))

# ---- Plot 1: Total_Length & N50 with labels ----
p1 <- ggplot(plot1_data, aes(x = Metric, y = Value, fill = Statistic)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  geom_text(aes(label = round(Value, 0)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Mean" = "skyblue", "Max" = "tomato")) +
  labs(title = "Assembly Stats:Total_Length & GC_Content",
       x = "Metric", y = "Value") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

# ---- Plot 2: GC_Content & Contigs with labels ----
p2 <- ggplot(plot2_data, aes(x = Metric, y = Value, fill = Statistic)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  geom_text(aes(label = round(Value, 0)), 
           position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Mean" = "skyblue", "Max" = "tomato")) +
  labs(title = "Assembly Stats:Contigs&Contig_N50",
       x = "Metric", y = "Value") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

# Display
print(p1)
print(p2)

Sydowiasssta_combined <- p1 | p2
ggsave("Sydowiasssta_combined.png", 
       Sydowiasssta_combined, width = 14, height = 6, dpi = 300)
Sydowiasssta_combined

view (summary_stats)


# 
# library(tidyverse)
# library(patchwork)
# 
# # ---------------------------
# # 1) Where your FST files are
# # ---------------------------
# FST_DIR <- "/Users/temitopefolorunso/Desktop/Folorunso dissertation proposal April 2024/BrownSpotProject/diversityfst1"
# 
# # Read and combine all FST files in a folder
# fst_all <- list.files(pattern = "\\.windowed\\.weir\\.fst$", full.names = TRUE) %>%
#   purrr::map_dfr(function(f) {
#     x <- readr::read_tsv(f, comment = "#", show_col_types = FALSE)
#     
#     # pick an Fst column robustly
#     x <- x %>%
#       mutate(
#         Fst = case_when(
#           "WEIGHTED_FST" %in% names(x) ~ WEIGHTED_FST,
#           "MEAN_FST"     %in% names(x) ~ MEAN_FST,
#           TRUE                         ~ NA_real_
#         )
#       ) %>%
#       filter(!is.na(Fst))
#     
#     bn <- basename(f)
#     
#     # extract "popA" and "popB" from filenames like fst_<A>_vs_<B>_5000bp.windowed.weir.fst
#     m <- stringr::str_match(bn, "^fst_(.+)_vs_(.+)_\\d+bp\\.windowed\\.weir\\.fst$")
#     popA <- m[,2]; popB <- m[,3]
#     
#     # severity and region parsing
#     sevA <- case_when(
#       str_ends(popA, "_High") ~ "High",
#       str_ends(popA, "_LowModerate") ~ "LowModerate",
#       TRUE ~ "All"
#     )
#     sevB <- case_when(
#       str_ends(popB, "_High") ~ "High",
#       str_ends(popB, "_LowModerate") ~ "LowModerate",
#       TRUE ~ "All"
#     )
#     
#     regA <- str_replace(popA, "_(High|LowModerate)$", "")
#     regB <- str_replace(popB, "_(High|LowModerate)$", "")
#     
#     # classify comparison type
#     ctype <- case_when(
#       regA == regB && sevA != sevB ~ "Disease",  # within-region severity comparison
#       regA != regB && sevA == sevB ~ "Region",   # between-region at same severity (or All)
#       TRUE ~ "Mixed"               # both differ (optional bucket)
#     )
#     
#     clabel <- case_when(
#       ctype == "Disease" ~ paste0(regA, ": ", sevA, " vs ", sevB),
#       ctype == "Region"  ~ paste0(sevA, ": ", regA, " vs ", regB),
#       TRUE               ~ paste0(regA, "(", sevA, ") vs ", regB, "(", sevB, ")")
#     )
#     
#     x %>%
#       mutate(
#         file = bn,
#         popA = popA, popB = popB,
#         regionA = regA, regionB = regB,
#         severityA = sevA, severityB = sevB,
#         contrast_type = ctype,
#         contrast = clabel,
#         # IMPORTANT: x-axis scaffold should be scaffold only (no sample prefix)
#         SCAFFOLD = stringr::str_extract(CHROM, "scaffold\\d+")
#       )
#   })
# fst_plot_df <- fst_all %>% filter(contrast_type %in% c("Disease","Region"))
# 
# make_manhattan <- function(df, out_png) {
#   
#   df <- df %>%
#     group_by(contrast, SCAFFOLD) %>%
#     mutate(max_pos = max(BIN_END, na.rm = TRUE)) %>%
#     ungroup()
#   
#   scaffold_sizes <- df %>%
#     distinct(contrast, SCAFFOLD, max_pos) %>%
#     arrange(contrast, SCAFFOLD) %>%
#     group_by(contrast) %>%
#     mutate(
#       offset = lag(cumsum(max_pos), default = 0),
#       mid_pos = offset + max_pos / 2
#     ) %>%
#     ungroup()
#   
#   plot_df <- df %>%
#     left_join(scaffold_sizes, by = c("contrast","SCAFFOLD","max_pos")) %>%
#     mutate(
#       genome_pos = BIN_START + offset,
#       col_group  = as.integer(factor(SCAFFOLD)) %% 2
#     )
#   
#   p <- ggplot(plot_df, aes(genome_pos, Fst, colour = factor(col_group))) +
#     geom_point(alpha = 0.8, size = 0.7) +
#     facet_wrap(~ contrast, ncol = 1, scales = "free_y") +
#     scale_x_continuous(
#       breaks = scaffold_sizes$mid_pos,
#       labels = scaffold_sizes$SCAFFOLD,
#       expand = expansion(mult = c(0.001, 0.001))
#     ) +
#     labs(x = "Scaffold", y = "Fst",
#          title = paste0("Genome-wide Fst (5 kb windows) — ", unique(df$contrast_type))) +
#     theme_bw(base_size = 12) +
#     theme(
#       panel.grid.major.x = element_blank(),
#       panel.grid.minor.x = element_blank(),
#       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
#       legend.position = "none"
#     )
#   
#   ggsave(out_png, p, width = 10, height = 10, dpi = 300)
#   p
# }
# 
# p_region  <- make_manhattan(fst_plot_df %>% filter(contrast_type=="Region"),
#                             "fst_manhattan_region.png")
# p_disease <- make_manhattan(fst_plot_df %>% filter(contrast_type=="Disease"),
#                             "fst_manhattan_disease.png")
# 
# # ---------------------------
# # 5) Manhattan-style plots
# # ---------------------------
# fst_all <- fst_all %>%
#   group_by(contrast, SCAFFOLD) %>%
#   mutate(max_pos = max(BIN_END, na.rm = TRUE)) %>%
#   ungroup()
# 
# scaffold_sizes_all <- fst_all %>%
#   distinct(contrast, SCAFFOLD, max_pos) %>%
#   arrange(contrast, SCAFFOLD) %>%
#   group_by(contrast) %>%
#   mutate(
#     offset = lag(cumsum(max_pos), default = 0),
#     mid_pos = offset + max_pos / 2
#   ) %>%
#   ungroup()
# 
# fst_plot_all <- fst_all %>%
#   left_join(scaffold_sizes_all, by = c("contrast", "SCAFFOLD", "max_pos")) %>%
#   mutate(
#     genome_pos = BIN_START + offset,
#     col_group  = as.integer(factor(SCAFFOLD)) %% 2
#   )
# 
# theme_pub <- theme_bw(base_size = 12) +
#   theme(
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
#     legend.position = "none"
#   )
# 
# p_fst_all <- ggplot(fst_plot_all, aes(x = genome_pos, y = Fst, colour = factor(col_group))) +
#   geom_point(alpha = 0.8, size = 0.7) +
#   facet_wrap(~ contrast, ncol = 1, scales = "free_y") +
#   scale_x_continuous(
#     breaks = scaffold_sizes_all$mid_pos,
#     labels = scaffold_sizes_all$SCAFFOLD,
#     expand = expansion(mult = c(0.001, 0.001))
#   ) +
#   labs(
#     x = "Scaffold",
#     y = "Fst",
#     title = "Genome-wide Fst for pairwise contrasts (5 kb windows)"
#   ) +
#   theme_pub
# 
# ggsave(file.path(FST_DIR, "fst_all_contrasts_manhattan.png"), p_fst_all,
#        width = 10, height = 10, dpi = 300)
# print(p_fst_all)
# 
# # ---------------------------
# # 6) Boxplot grouped as requested
# #    - Disease vs Region separate
# #    - Region labels have NO High/LowModerate
# # ---------------------------
# fst_all <- fst_all %>%
#   mutate(
#     # Make a stable plotting order
#     contrast = fct_inorder(contrast),
#     group_type = factor(group_type, levels = c("Disease", "Region"))
#   )
# 
# p_fst_box <- ggplot(fst_all, aes(x = contrast, y = Fst)) +
#   geom_boxplot(outlier.size = 0.4, alpha = 0.9) +
#   facet_wrap(~ group_type, scales = "free_x", ncol = 1) +
#   labs(
#     x = NULL,
#     y = "Fst",
#     title = "Fst distribution by contrast (Disease vs Region)"
#   ) +
#   theme_bw(base_size = 12) +
#   theme(
#     axis.text.x = element_text(angle = 20, hjust = 1),
#     legend.position = "none"
#   )
# 
# ggsave(file.path(FST_DIR, "fst_group_boxplot_grouped.png"), p_fst_box,
#        width = 10, height = 6, dpi = 300)
# print(p_fst_box)
# 
# # Optional combined figure
# Fst_combined <- p_fst_all / p_fst_box
# ggsave(file.path(FST_DIR, "Fst_combined.png"), Fst_combined,
#        width = 12, height = 14, dpi = 300)
# 



suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# =============================================================================
# CONFIG
# =============================================================================
FST_DIR   <- "/Users/temitopefolorunso/Desktop/Folorunso dissertation proposal April 2024/BrownSpotProject/diversityfst1"         # folder containing fst_*.windowed.weir.fst
WINDOW_BP <- 5000

# --- What you want plotted (EXACT sets) ---
# Disease (4)
DISEASE_PATTERNS <- c(
  "North_LowModerate_vs_North_High",
  "North_LowModerate_vs_South_High",
  "South_High_vs_North_High",
  "South_High_vs_South_LowModerate"
)

# Region (3)
REGION_PATTERNS <- c(
  "South_vs_North",
  "North_vs_OAWEST",
  "South_vs_OAWEST"
)

# Output files
OUT_SUMMARY_CSV     <- file.path(FST_DIR, "fst_summary_selected.csv")
OUT_MH_DISEASE_PNG  <- file.path(FST_DIR, "fst_manhattan_DISEASE.png")
OUT_MH_REGION_PNG   <- file.path(FST_DIR, "fst_manhattan_REGION.png")
OUT_BOX_DISEASE_PNG <- file.path(FST_DIR, "fst_box_DISEASE.png")
OUT_BOX_REGION_PNG  <- file.path(FST_DIR, "fst_box_REGION.png")
OUT_COMBINED_PNG    <- file.path(FST_DIR, "fst_selected_combined.png")

# =============================================================================
# HELPERS
# =============================================================================
# Parse popA/popB from standard filename:
# fst_<A>_vs_<B>_<N>bp.windowed.weir.fst
parse_filename_pops <- function(path) {
  bn <- basename(path)
  m <- stringr::str_match(bn, "^fst_(.+)_vs_(.+)_\\d+bp\\.windowed\\.weir\\.fst$")
  if (any(is.na(m))) {
    stop("Bad filename (doesn't match fst_<A>_vs_<B>_<N>bp.windowed.weir.fst): ", bn)
  }
  list(popA = m[,2], popB = m[,3], basename = bn)
}

# Normalize contrast key so "A_vs_B" equals "B_vs_A"
norm_pair_key <- function(a, b) {
  x <- sort(c(a, b))
  paste0(x[1], "_vs_", x[2])
}

# Extract scaffold number (as numeric string) from CHROM.
# Works for:
#   scaffold1, scaffold_1, BSxx...scaffold1, BSxx...scaffold_12, etc.
extract_scaffold_num <- function(chrom) {
  num <- stringr::str_extract(chrom, "(?i)scaffold[_-]?\\d+")
  num <- stringr::str_extract(num, "\\d+")
  # If there is no scaffold pattern, keep original (but this will make ugly labels)
  ifelse(is.na(num), chrom, num)
}

# For plotting, we want scaffolds ordered numerically (1,2,3...)
scaffold_order_numeric <- function(scaffold_vec) {
  suppressWarnings(as.integer(scaffold_vec))
}

# Pretty labels for facets (exactly as you described)
facet_label_map <- function(key) {
  dplyr::case_when(
    key == norm_pair_key("North_LowModerate","North_High") ~ "North Low/Moderate vs North High",
    key == norm_pair_key("North_LowModerate","South_High") ~  "North Low/Moderate vs South High",
    key == norm_pair_key("South_High","North_High")        ~ "South High vs North High",
    key == norm_pair_key("South_High","South_LowModerate") ~ "South High vs South Low/Moderate",
    key == norm_pair_key("South","North")                  ~ "South vs North",
    key == norm_pair_key("North","OAWEST")                 ~ "North vs OAWEST",
    key == norm_pair_key("South","OAWEST")                 ~ "South vs OAWEST",
    TRUE ~ key
  )
}

make_manhattan <- function(df, out_png, title) {
  if (nrow(df) == 0) {
    message("No data for: ", title, " (skipping)")
    return(NULL)
  }
  
  # per-contrast scaffold sizes + offsets
  df2 <- df %>%
    group_by(contrast, SCAFFOLD) %>%
    mutate(max_pos = max(BIN_END, na.rm = TRUE)) %>%
    ungroup()
  
  scaffold_sizes <- df2 %>%
    distinct(contrast, SCAFFOLD, max_pos, SCAFFOLD_INT) %>%
    arrange(contrast, SCAFFOLD_INT) %>%
    group_by(contrast) %>%
    mutate(
      offset = lag(cumsum(max_pos), default = 0),
      mid_pos = offset + max_pos / 2
    ) %>%
    ungroup()
  
  plot_df <- df2 %>%
    left_join(scaffold_sizes, by = c("contrast", "SCAFFOLD", "max_pos", "SCAFFOLD_INT")) %>%
    mutate(
      genome_pos = BIN_START + offset,
      col_group  = as.integer(factor(SCAFFOLD_INT)) %% 2
    )
  
  p <- ggplot(plot_df, aes(x = genome_pos, y = Fst, colour = factor(col_group))) +
    geom_point(alpha = 0.8, size = 0.7) +
    facet_wrap(~ contrast, ncol = 1, scales = "free_y") +
    scale_x_continuous(
      breaks = scaffold_sizes$mid_pos,
      labels = scaffold_sizes$SCAFFOLD,   # <- ONLY scaffold number like 1,2,3...
      expand = expansion(mult = c(0.001, 0.001)),
      guide  = guide_axis(angle = 90, check.overlap = TRUE)
    ) +
    labs(x = "Scaffold", y = "Fst", title = title) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      legend.position = "none"
    )
  
  ggsave(out_png, p, width = 20, height = 15, dpi = 600)
  p
}

make_box <- function(df, out_png, title) {
  if (nrow(df) == 0) {
    message("No data for: ", title, " (skipping)")
    return(NULL)
  }
  p <- ggplot(df, aes(x = contrast, y = Fst, fill = contrast)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.9) +
    labs(x = NULL, y = "Fst", title = title) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1),
          legend.position = "none")
  
  ggsave(out_png, p, width = 10, height = 5, dpi = 300)
  p
}

# =============================================================================
# 1) DISCOVER FILES
# =============================================================================
fst_files <- list.files(
  path = FST_DIR,
  pattern = "\\.windowed\\.weir\\.fst$",
  full.names = TRUE
)
if (length(fst_files) == 0) stop("No *.windowed.weir.fst found in: ", normalizePath(FST_DIR))

# =============================================================================
# 2) READ + FILTER TO YOUR EXACT CONTRAST SETS
# =============================================================================
# Build lookup of file -> normalized key
file_index <- purrr::map_dfr(fst_files, function(f) {
  pops <- parse_filename_pops(f)
  key  <- norm_pair_key(pops$popA, pops$popB)
  tibble(file = f, popA = pops$popA, popB = pops$popB, key = key)
})

# Normalize your requested pattern lists into normalized keys
disease_keys <- purrr::map_chr(DISEASE_PATTERNS, function(pat) {
  parts <- strsplit(pat, "_vs_")[[1]]
  norm_pair_key(parts[1], parts[2])
})
region_keys <- purrr::map_chr(REGION_PATTERNS, function(pat) {
  parts <- strsplit(pat, "_vs_")[[1]]
  norm_pair_key(parts[1], parts[2])
})

wanted_keys <- unique(c(disease_keys, region_keys))

# Subset to files matching any wanted key
wanted_files <- file_index %>% filter(key %in% wanted_keys)
if (nrow(wanted_files) == 0) {
  stop("None of your requested contrasts were found.\nRequested keys:\n  ",
       paste(wanted_keys, collapse = "\n  "),
       "\n\nAvailable keys:\n  ",
       paste(sort(unique(file_index$key)), collapse = "\n  "))
}

# Read data
fst_all <- purrr::map_dfr(wanted_files$file, function(f) {
  pops <- parse_filename_pops(f)
  key  <- norm_pair_key(pops$popA, pops$popB)
  
  x <- readr::read_tsv(f, comment = "#", show_col_types = FALSE)
  
  x <- x %>%
    mutate(
      Fst = case_when(
        "WEIGHTED_FST" %in% names(x) ~ WEIGHTED_FST,
        "MEAN_FST"     %in% names(x) ~ MEAN_FST,
        TRUE                         ~ NA_real_
      ),
      key = key,
      contrast = facet_label_map(key),
      contrast_type = case_when(
        key %in% disease_keys ~ "Disease",
        key %in% region_keys  ~ "Region",
        TRUE ~ "Other"
      ),
      SCAFFOLD = extract_scaffold_num(CHROM),
      SCAFFOLD_INT = scaffold_order_numeric(SCAFFOLD)
    ) %>%
    filter(!is.na(Fst))
  
  x
})

# Keep only Disease + Region and keep your desired ordering
fst_all <- fst_all %>%
  filter(contrast_type %in% c("Disease","Region")) %>%
  mutate(
    contrast = factor(
      contrast,
      levels = c(
        facet_label_map(disease_keys),
        facet_label_map(region_keys)
      )
    )
  )

# Split sets
fst_disease <- fst_all %>% filter(contrast_type == "Disease")
fst_region  <- fst_all %>% filter(contrast_type == "Region")

# =============================================================================
# 3) SUMMARY
# =============================================================================
fst_summary <- fst_all %>%
  group_by(contrast_type, contrast) %>%
  summarise(
    n_windows  = n(),
    mean_Fst   = mean(Fst, na.rm = TRUE),
    sd_Fst     = sd(Fst, na.rm = TRUE),
    median_Fst = median(Fst, na.rm = TRUE),
    Q1_Fst     = quantile(Fst, 0.25, na.rm = TRUE),
    Q3_Fst     = quantile(Fst, 0.75, na.rm = TRUE),
    min_Fst    = min(Fst, na.rm = TRUE),
    max_Fst    = max(Fst, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(contrast_type, contrast)

print(fst_summary)
readr::write_csv(fst_summary, OUT_SUMMARY_CSV)
message("Wrote: ", OUT_SUMMARY_CSV)

# =============================================================================
# 4) PLOTS (EXACTLY AS REQUESTED)
# =============================================================================
p_mh_dis <- make_manhattan(
  fst_disease,
  OUT_MH_DISEASE_PNG,
  paste0("Genome-wide Fst — DISEASE (", WINDOW_BP, " bp windows)")
)

p_mh_reg <- make_manhattan(
  fst_region,
  OUT_MH_REGION_PNG,
  paste0("Genome-wide Fst — REGION (", WINDOW_BP, " bp windows)")
)

p_box_dis <- make_box(
fst_disease,
OUT_BOX_DISEASE_PNG,
"Fst distribution — DISEASE"
)
# 
p_box_reg <- make_box(
fst_region,
OUT_BOX_REGION_PNG,
"Fst distribution — REGION"
)

# Combined figure (optional)
plots <- list()
if (!is.null(p_mh_dis)) plots <- c(plots, list(p_mh_dis))
if (!is.null(p_mh_reg)) plots <- c(plots, list(p_mh_reg))
if (!is.null(p_box_dis)) plots <- c(plots, list(p_box_dis))
if (!is.null(p_box_reg)) plots <- c(plots, list(p_box_reg))

if (length(plots) >= 2) {
  combined <- wrap_plots(plots, ncol = 2)
  ggsave(OUT_COMBINED_PNG, combined, width = 25, height = 10, dpi = 600)
  message("Wrote: ", OUT_COMBINED_PNG)
}




library(patchwork)

#-----------------------------#
# Manhattan-only figure
#-----------------------------#
mh_plots <- list()
if (!is.null(p_mh_dis)) mh_plots <- c(mh_plots, list(p_mh_dis))
if (!is.null(p_mh_reg)) mh_plots <- c(mh_plots, list(p_mh_reg))

if (length(mh_plots) >= 1) {
  combined_manhattan <- wrap_plots(mh_plots, ncol = 1)
  ggsave(
    "manhattan_combined.png",
    combined_manhattan,
    width = 18,
    height = 12,
    dpi = 600
  )
  message("Wrote: manhattan_combined.png")
}

#-----------------------------#
# Boxplot-only figure
#-----------------------------#
box_plots <- list()
if (!is.null(p_box_dis)) box_plots <- c(box_plots, list(p_box_dis))
if (!is.null(p_box_reg)) box_plots <- c(box_plots, list(p_box_reg))

if (length(box_plots) >= 1) {
  combined_box <- wrap_plots(box_plots, ncol = 1)
  ggsave(
    "boxplot_combined.png",
    combined_box,
    width = 14,
    height = 10,
    dpi = 600
  )
  message("Wrote: boxplot_combined.png")
}
message("Done.")










#Genome-wide nucleotide diversity (π) and Tajima’s


folder1 <- "/Users/temitopefolorunso/Desktop/Folorunso dissertation proposal April 2024/BrownSpotProject/diversity_stats"

pi_file1      <- file.path(folder1, "pi_5000bp.windowed.pi")
tajima_file1  <- file.path(folder1,"tajimaD_5000bp.Tajima.D")
window_size   <- 5000    # for titles
species_name  <- "Sydowia sp."

## --------- READ DATA -------------
message("Reading π and Tajima's D files...")

pi_df1     <- read_tsv(pi_file1, comment = "#", show_col_types = FALSE)
tajima_df1 <- read_tsv(tajima_file1, comment = "#", show_col_types = FALSE)

## Expected columns:
## π:      CHROM, BIN_START, BIN_END, N_VARIANTS, PI
## Tajima: CHROM, BIN_START, BIN_END (optional), N_SNPS, TajimaD

## --------- CLEAN SCAFFOLD NAMES & NUMERIC ORDER ---------
pi_df1 <- pi_df1 %>%
  mutate(
    SCAFFOLD   = sub(".*(scaffold[0-9]+)$", "\\1", CHROM),
    scaffold_n = as.integer(sub("scaffold", "", SCAFFOLD))
  )

tajima_df1 <- tajima_df1 %>%
  mutate(
    SCAFFOLD   = sub(".*(scaffold[0-9]+)$", "\\1", CHROM),
    scaffold_n = as.integer(sub("scaffold", "", SCAFFOLD))
  )

## --------- BUILD GENOME-WIDE COORDINATES (ALL SCAFFOLDS) ---------
## Use π table to define scaffold lengths
# scaffold_sizes1 <- pi_df1 %>%
#   group_by(SCAFFOLD, scaffold_n) %>%
#   summarise(max_pos = max(BIN_END, na.rm = TRUE), .groups = "drop") %>%
#   arrange(scaffold_n) %>%                      # numeric order
#   mutate(
#     offset = lag(cumsum(max_pos), default = 0),
#     mid_pos = offset + max_pos / 2
#   )
# 
# ## Add genome-wide position to π and Tajima tables
# pi_plot_df1 <- pi_df1 %>%
#   left_join(scaffold_sizes1, by = c("SCAFFOLD", "scaffold_n")) %>%
#   mutate(genome_pos = BIN_START + offset)
# 
# tajima_plot_df1 <- tajima_df1 %>%
#   left_join(scaffold_sizes1, by = c("SCAFFOLD", "scaffold_n")) %>%
#   mutate(genome_pos = BIN_START + offset)


# Build ordered scaffold table
scaffold_sizes1 <- pi_df1 %>%
  group_by(SCAFFOLD, scaffold_n) %>%
  summarise(max_pos = max(BIN_END, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(scaffold_n)) %>%
  arrange(scaffold_n, SCAFFOLD) %>%
  mutate(
    offset = cumsum(lag(max_pos, default = 0)),
    mid_pos = offset + max_pos / 2,
    scaffold_group = rep(c("A", "B"), length.out = n())   # force true alternation
  )

# Add genome-wide position and ordered color group
pi_plot_df1 <- pi_df1 %>%
  left_join(
    scaffold_sizes1 %>% select(SCAFFOLD, scaffold_n, offset, mid_pos, scaffold_group),
    by = c("SCAFFOLD", "scaffold_n")
  ) %>%
  mutate(genome_pos = BIN_START + offset)

tajima_plot_df1 <- tajima_df1 %>%
  left_join(
    scaffold_sizes1 %>% select(SCAFFOLD, scaffold_n, offset, mid_pos, scaffold_group),
    by = c("SCAFFOLD", "scaffold_n")
  ) %>%
  mutate(genome_pos = BIN_START + offset)

## --------- SUMMARY VALUES FOR RESULTS ---------
pi_summary <- pi_plot_df1 %>%
  summarise(
    n_windows   = n(),
    mean_pi     = mean(PI, na.rm = TRUE),
    sd_pi       = sd(PI, na.rm = TRUE),
    median_pi   = median(PI, na.rm = TRUE),
    Q1_pi       = quantile(PI, 0.25, na.rm = TRUE),
    Q3_pi       = quantile(PI, 0.75, na.rm = TRUE),
    min_pi      = min(PI, na.rm = TRUE),
    max_pi      = max(PI, na.rm = TRUE)
  )

tajima_summary <- tajima_plot_df1 %>%
  summarise(
    n_windows      = n(),
    mean_TajimaD   = mean(TajimaD, na.rm = TRUE),
    sd_TajimaD     = sd(TajimaD, na.rm = TRUE),
    median_TajimaD = median(TajimaD, na.rm = TRUE),
    Q1_TajimaD     = quantile(TajimaD, 0.25, na.rm = TRUE),
    Q3_TajimaD     = quantile(TajimaD, 0.75, na.rm = TRUE),
    min_TajimaD    = min(TajimaD, na.rm = TRUE),
    max_TajimaD    = max(TajimaD, na.rm = TRUE)
  )

## Print to console so you can copy into Results
message("π summary:")
print(pi_summary)
message("\nTajima's D summary:")
print(tajima_summary)

## Save as CSV for record
write_csv(pi_summary,     "pi_summary_global.csv")
write_csv(tajima_summary, "tajimaD_summary_global.csv")

## --------- LONG FORMAT FOR COMBINED MANHATTAN PANEL ---------
long_df1 <- bind_rows(
  pi_plot_df1 %>%
    transmute(SCAFFOLD, scaffold_n, genome_pos, value = PI, stat = "Pi"),
  tajima_plot_df1 %>%
    transmute(SCAFFOLD, scaffold_n, genome_pos, value = TajimaD, stat = "TajimaD")
)

long_df1$stat <- factor(long_df1$stat, levels = c("Pi", "TajimaD"))

## Alternating scaffold colours
long_df1 <- long_df1 %>%
  left_join(
    scaffold_sizes1 %>% transmute(SCAFFOLD, scaffold_n, scaffold_index = row_number()),
    by = c("SCAFFOLD", "scaffold_n")
  ) %>%
  mutate(col_group = scaffold_index %% 2)  # 0/1

## Base theme
theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    strip.background = element_rect(fill = "grey90", colour = NA),
    legend.position = "none"
  )

## --------- COMBINED MANHATTAN PANEL: Pi + TajimaD ---------
p_panel1 <- ggplot(
  pi_plot_df1,
  aes(x = genome_pos, y = value, color = scaffold_group)
) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_smooth(
    aes(group = stat),
    method = "loess", se = FALSE, span = 0.3,
    linewidth = 0.5, colour = "black"
  ) +
  facet_grid(rows = vars(stat), scales = "free_y", switch = "y") +
  scale_x_continuous(
    breaks = scaffold_sizes1$mid_pos,
    labels = scaffold_sizes1$SCAFFOLD,
    expand = expansion(mult = c(0.001, 0.001))
  ) +
  labs(
    x = "Scaffold",
    y = "Statistic value",
    title = paste0(
      species_name,
      " – genome-wide diversity (", window_size / 1000, " kb windows)"
    )
  ) +
  theme_pub +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0)
  ) +
  scale_colour_manual(values = c("#1f77b4", "#2ca02c"))

ggsave("diversity_panel_pi_tajimaD_allScaffolds.png",
       p_panel1, width = 10, height = 7, dpi = 300)
print(p_panel1)



scaffold_sizes1 %>%
  count(scaffold_n) %>%
  filter(n > 1)

## --------- OPTIONAL: separate Manhattan plots ---------
p_pi1 <- ggplot(pi_plot_df1,
                aes(x = genome_pos, y = PI,
                    colour = factor(as.integer(factor(SCAFFOLD)) %% 2))) +
  geom_point(alpha = 0.7, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3,
              linewidth = 0.5, colour = "black") +
  scale_x_continuous(
    breaks = scaffold_sizes1$mid_pos,
    labels = scaffold_sizes1$SCAFFOLD,
    expand = expansion(mult = c(0.001, 0.001))
  ) +
  labs(
    x = "Scaffold",
    y = "Nucleotide diversity (π)",
    #title = paste0("nucleotide diversity )
  ) +
  theme_pub +
  scale_colour_manual(values = c("#1f77b4", "#2ca02c"))

ggsave("pi_manhattan_allScaffolds.png", p_pi1,
       width = 12, height = 4, dpi = 300)

p_td1 <- ggplot(tajima_plot_df1,
                aes(x = genome_pos, y = TajimaD,
                    colour = factor(as.integer(factor(SCAFFOLD)) %% 2))) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(alpha = 0.7, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3,
              linewidth = 0.5, colour = "black") +
  scale_x_continuous(
    breaks = scaffold_sizes1$mid_pos,
    labels = scaffold_sizes1$SCAFFOLD,
    expand = expansion(mult = c(0.001, 0.001))
  ) +
  labs(
    x = "Scaffold",
    y = "Tajima's D",
    #title = paste0("Tajima's D)
  ) +
  theme_pub +
  scale_colour_manual(values = c("#1f77b4", "#2ca02c"))

ggsave("tajimaD_manhattan_allScaffolds.png", p_td1,
       width = 12, height = 4, dpi = 300)





library(dplyr)
library(ggplot2)
library(grid)

## --------- COMBINE Pi + TajimaD INTO ONE LONG TABLE ---------
panel_df1 <- bind_rows(
  pi_plot_df1 %>%
    transmute(
      SCAFFOLD,
      scaffold_n,
      scaffold_group,
      genome_pos,
      value = PI,                        # change PI if your column name differs
      stat  = "Nucleotide diversity (π)"
    ),
  tajima_plot_df1 %>%
    transmute(
      SCAFFOLD,
      scaffold_n,
      scaffold_group,
      genome_pos,
      value = TajimaD,                   # change TajimaD if your column name differs
      stat  = "Tajima's D"
    )
) %>%
  mutate(
    stat = factor(stat, levels = c("Nucleotide diversity (π)", "Tajima's D"))
  )

## --------- COMBINED MANHATTAN PANEL: Pi + TajimaD ---------
p_panel1 <- ggplot(
  panel_df1,
  aes(x = genome_pos, y = value, colour = scaffold_group)
) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_smooth(
    aes(group = 1),
    method = "loess",
    se = FALSE,
    span = 0.25,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_hline(
    data = data.frame(stat = "Nucleotide diversity (π)", yint = -Inf),
    aes(yintercept = yint),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.5
  ) +
  facet_grid(
    rows = vars(stat),
    scales = "free_y",
    switch = "y"
  ) +
  scale_x_continuous(
    breaks = scaffold_sizes1$mid_pos,
    labels = scaffold_sizes1$SCAFFOLD,
    expand = expansion(mult = c(0.001, 0.001))
  ) +
  scale_colour_manual(values = c("A" = "#2ca25f", "B" = "#3182bd")) +
  labs(
    x = "Scaffold",
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(
      angle = 90,
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 7
    ),
    panel.spacing.y = unit(0.3, "lines"),
    panel.border = element_blank(),
    axis.line.x.top = element_blank(),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.5),
    axis.line.y.left = element_line(colour = "black", linewidth = 0.5),
    axis.line.y.right = element_blank()
  )

ggsave("diversity_panel_pi_tajimaD_allScaffolds.png",
       p_panel1, width = 15, height = 8, dpi = 600)
print(p_panel1)





#group analysis
suppressPackageStartupMessages({
  library(tidyverse)
  library(purrr)
  library(patchwork)
})

folder <- "/Users/temitopefolorunso/Desktop/Folorunso dissertation proposal April 2024/BrownSpotProject/diversity_statsnew1"

pi_files <- c(
  "South"              = file.path(folder, "pi_South5000bp.windowed.pi"),
  "South Low/Moderate"  = file.path(folder, "pi_South_LowModerate_5000bp.windowed.pi"),
  "South High"         = file.path(folder, "pi_South_High_5000bp.windowed.pi"),
  "North"              = file.path(folder, "pi_North_5000bp.windowed.pi"),
  "North Low/Moderate"  = file.path(folder, "pi_North_LowModerate_5000bp.windowed.pi"),
  "North High"         = file.path(folder, "pi_North_High_5000bp.windowed.pi"),
  "OAWEST"             = file.path(folder, "pi_OAWEST_5000bp.windowed.pi")
)

taj_files <- c(
  "South"              = file.path(folder, "tajima_South5000bp.Tajima.D"),
  "South Low/Moderate"  = file.path(folder, "tajima_South_LowModerate_5000bp.Tajima.D"),
  "South High"         = file.path(folder, "tajima_South_High_5000bp.Tajima.D"),
  "North"              = file.path(folder, "tajima_North_5000bp.Tajima.D"),
  "North Low/Moderate"  = file.path(folder, "tajima_North_LowModerate_5000bp.Tajima.D"),
  "North High"         = file.path(folder, "tajima_North_High_5000bp.Tajima.D"),
  "OAWEST"             = file.path(folder, "tajima_OAWEST_5000bp.Tajima.D")
)

# =============================================================================
# 2) READERS (with a simple file-exists check)
# =============================================================================
read_pi <- function(file, group) {
  if (!file.exists(file)) {
    warning("Missing π file for ", group, ": ", file, call. = FALSE)
    return(tibble())
  }
  readr::read_tsv(file, comment = "#", show_col_types = FALSE) %>%
    mutate(Group = group)
}

read_taj <- function(file, group) {
  if (!file.exists(file)) {
    warning("Missing TajimaD file for ", group, ": ", file, call. = FALSE)
    return(tibble())
  }
  readr::read_tsv(file, comment = "#", show_col_types = FALSE) %>%
    mutate(Group = group)
}

# =============================================================================
# 3) LOAD + CLEAN
# =============================================================================
pi_all  <- purrr::map2_dfr(unname(pi_files),  names(pi_files),  read_pi)
taj_all <- purrr::map2_dfr(unname(taj_files), names(taj_files), read_taj)

# Standard vcftools columns are usually: CHROM BIN_START BIN_END N_VARIANTS PI
# and for Tajima.D: CHROM BIN_START BIN_END N_SNPS TajimaD
# Keep it flexible but enforce key columns exist:
if (!"PI" %in% names(pi_all) && nrow(pi_all) > 0) {
  stop("π files loaded, but column 'PI' not found. Check file format.")
}
if (!"TajimaD" %in% names(taj_all) && nrow(taj_all) > 0) {
  stop("Tajima files loaded, but column 'TajimaD' not found. Check file format.")
}

# Remove NA and impossible values
pi_all  <- pi_all  %>% filter(!is.na(PI), PI > 0)
taj_all <- taj_all %>% filter(!is.na(TajimaD))

# =============================================================================
# 4) GROUP ORDER + COLORS (expanded)
# =============================================================================
group_order <- c(
  "South",
  "South Low/Moderate",
  "South High",
  "North",
  "North Low/Moderate",
  "North High",
  "OAWEST"
)

pi_all$Group  <- factor(pi_all$Group,  levels = group_order)
taj_all$Group <- factor(taj_all$Group, levels = group_order)

# Optional: keep your original palette idea, add new groups
group_colors <- c(
  "South"              = "#1d91c0",
  "South Low/Moderate"       = "#253494",
  "South High"         = "#4c7fb9",
  "North"              = "#feb24c",
  "North Low/Moderate"       = "#f03b20",
  "North High"         = "#bd0026",
  "OAWEST"             = "#31a354"
)

# =============================================================================
# 5) π BOXPLOT
# =============================================================================
p_pi_box <- ggplot(pi_all, aes(x = Group, y = PI, fill = Group)) +
  geom_boxplot(outlier.size = 0.4, width = 0.3) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_y_log10() +
  labs(x = NULL, y = "Nucleotide diversity (π)") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave("pi_group_boxplot.png", p_pi_box, width = 5, height = 4, dpi = 300)

# =============================================================================
# 6) Tajima's D BOXPLOT
# =============================================================================
p_taj_box <- ggplot(taj_all, aes(x = Group, y = TajimaD, fill = Group)) +
  geom_boxplot(outlier.size = 0.4, width = 0.3) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(x = NULL, y = "Tajima's D") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave("tajimaD_group_boxplot.png", p_taj_box, width = 5, height = 4, dpi = 300)

# =============================================================================
# 7) COMBINE
# =============================================================================
group_combined <- p_pi1 / p_td1   | p_pi_box  / p_taj_box  
ggsave("group_combined.png", group_combined, width = 16, height = 6, dpi = 300)

print(group_combined)




group_combined <- p_pi_box | p_taj_box
ggsave("group_combined.png", 
       group_combined, width = 20, height = 10, dpi = 300)
group_combined 


# ROH analysis

library(tidyverse)
library(janitor)

## ==============================
## 1. FILE PATHS
## ==============================

folder2 <- "/Users/temitopefolorunso/Desktop/Folorunso dissertation proposal April 2024/BrownSpotProject/hetero_homozy1"

het_file       <- file.path(folder2, "het_individuals.het")
roh_indiv_file <- file.path(folder2, "sydowia_roh.hom.indiv")
roh_seg_file   <- file.path(folder2, "sydowia_roh.hom")

metadata_file  <- "Sydowiaregion.txt"
id_col         <- "isolates"
group_col      <- "region_disease"

## ==============================
## 2. READ HETEROZYGOSITY
## ==============================

het_df <- read_table(het_file, col_names = TRUE, show_col_types = FALSE) %>%
  # vcftools --het columns: INDV, O(HOM), E(HOM), N(NM), F
  rename(
    O_HOM   = `O(HOM)`,
    E_HOM   = `E(HOM)`,
    N_SITES = N_SITES
  ) %>%
  mutate(
    Ho  = (N_SITES - O_HOM) / N_SITES,   # observed heterozygosity
    INDV = as.character(INDV)
  )

## ==============================
## 3. READ PLINK ROH SUMMARY (per individual)
## ==============================

roh_indiv <- read_table(roh_indiv_file,
                        col_names = TRUE,
                        show_col_types = FALSE) %>%
  mutate(
    IID        = as.character(IID),
    total_Mb_ROH = KB / 1000          # convert kb → Mb
  )

## ==============================
## 4. READ ROH SEGMENTS TABLE
## ==============================

roh_seg <- read_table(roh_seg_file,
                      col_names = TRUE,
                      show_col_types = FALSE) %>%
  mutate(
    IID    = as.character(IID),
    seg_Mb = KB / 1000
  )

## ==============================
## 5. ADD GROUP INFO FROM Sydowiaregion.txt
## ==============================

meta <- read_delim(metadata_file,
                   delim = "\t",
                   col_names = TRUE,
                   show_col_types = FALSE) %>%
  clean_names()

meta_small <- meta %>%
  select(all_of(c(id_col, group_col))) %>%
  rename(
    IID   = !!id_col,
    Group = !!group_col
  ) %>%
  mutate(
    IID   = as.character(IID),
    Group = as.factor(Group)
  )

het_df <- het_df %>%
  left_join(meta_small, by = c("INDV" = "IID"))

roh_indiv <- roh_indiv %>%
  left_join(meta_small, by = "IID")

roh_seg <- roh_seg %>%
  left_join(meta_small, by = "IID")

## =================================================
## 6. MERGE HET + ROH PER INDIVIDUAL (for joint plots)
## =================================================

het_roh <- het_df %>%
  rename(IID = INDV) %>%
  left_join(
    roh_indiv %>% select(IID, NSEG, total_Mb_ROH),
    by = "IID"
  )

## =================================================
## 6b. SUMMARY TABLES FOR RESULTS (NEW)
## =================================================

## 6b.1 Heterozygosity by group
het_group_summary <- het_df %>%
  group_by(Group) %>%
  summarise(
    n_individuals = n(),
    mean_Ho       = mean(Ho, na.rm = TRUE),
    sd_Ho         = sd(Ho, na.rm = TRUE),
    median_Ho     = median(Ho, na.rm = TRUE),
    Q1_Ho         = quantile(Ho, 0.25, na.rm = TRUE),
    Q3_Ho         = quantile(Ho, 0.75, na.rm = TRUE),
    min_Ho        = min(Ho, na.rm = TRUE),
    max_Ho        = max(Ho, na.rm = TRUE),
    .groups       = "drop"
  )

message("\nHeterozygosity summary by group:")
print(het_group_summary)
write_csv(het_group_summary, "heterozygosity_summary_by_group.csv")

## 6b.2 Total ROH length (Mb) by group
roh_group_summary <- roh_indiv %>%
  group_by(Group) %>%
  summarise(
    n_individuals  = n(),
    mean_ROH_Mb    = mean(total_Mb_ROH, na.rm = TRUE),
    sd_ROH_Mb      = sd(total_Mb_ROH, na.rm = TRUE),
    median_ROH_Mb  = median(total_Mb_ROH, na.rm = TRUE),
    Q1_ROH_Mb      = quantile(total_Mb_ROH, 0.25, na.rm = TRUE),
    Q3_ROH_Mb      = quantile(total_Mb_ROH, 0.75, na.rm = TRUE),
    min_ROH_Mb     = min(total_Mb_ROH, na.rm = TRUE),
    max_ROH_Mb     = max(total_Mb_ROH, na.rm = TRUE),
    .groups        = "drop"
  )

message("\nTotal ROH (Mb) summary by group:")
print(roh_group_summary)
write_csv(roh_group_summary, "roh_total_summary_by_group.csv")

## 6b.3 ROH segment length distribution by group
roh_seg_summary <- roh_seg %>%
  group_by(Group) %>%
  summarise(
    n_segments   = n(),
    mean_seg_Mb  = mean(seg_Mb, na.rm = TRUE),
    sd_seg_Mb    = sd(seg_Mb, na.rm = TRUE),
    median_seg_Mb = median(seg_Mb, na.rm = TRUE),
    Q1_seg_Mb    = quantile(seg_Mb, 0.25, na.rm = TRUE),
    Q3_seg_Mb    = quantile(seg_Mb, 0.75, na.rm = TRUE),
    min_seg_Mb   = min(seg_Mb, na.rm = TRUE),
    max_seg_Mb   = max(seg_Mb, na.rm = TRUE),
    .groups      = "drop"
  )

message("\nROH segment length summary by group:")
print(roh_seg_summary)
write_csv(roh_seg_summary, "roh_segment_summary_by_group.csv")

## 6b.4 Correlation: ROH vs Heterozygosity
cor_overall <- cor(het_roh$total_Mb_ROH,
                   het_roh$Ho,
                   use = "complete.obs")

cor_by_group <- het_roh %>%
  group_by(Group) %>%
  summarise(
    n_individuals = sum(!is.na(total_Mb_ROH) & !is.na(Ho)),
    cor_ROH_Ho    = cor(total_Mb_ROH, Ho, use = "complete.obs"),
    .groups       = "drop"
  )

message("\nCorrelation between total ROH (Mb) and Ho – overall:")
print(cor_overall)
message("\nCorrelation between total ROH (Mb) and Ho – by group:")
print(cor_by_group)
write_csv(cor_by_group, "correlation_ROH_Ho_by_group.csv")

## ==============================
## 7. PLOTS  (unchanged)
## ==============================

## 7.1 Observed heterozygosity per individual
p_het <- ggplot(het_df,
                aes(x = reorder(INDV, Ho),
                    y = Ho,
                    fill = Group)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Isolates",
    y = "Observed heterozygosity (H\u2092)",
    title = "Per-individual heterozygosity"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 6)
  )

ggsave("het_per_individual.png", p_het, width = 6, height = 6, dpi = 300)
p_het




## 7.2 Total ROH length per individual (Mb)
p_roh_tot <- ggplot(roh_indiv,
                    aes(x = reorder(IID, total_Mb_ROH),
                        y = total_Mb_ROH,
                        fill = Group)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Isolates",
    y = "Total ROH length (Mb)",
    title = "Total ROH per individual"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 6)
  )

ggsave("roh_total_per_individual.png", p_roh_tot,
       width = 10, height = 6, dpi = 300)
p_roh_tot






## 7.3 Relationship: ROH vs heterozygosity
p_scatter <- ggplot(het_roh,
                    aes(x = total_Mb_ROH,
                        y = Ho,
                        colour = Group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE,
              colour = "black", linewidth = 0.6) +
  labs(
    x = "Total ROH length (Mb)",
    y = "Observed heterozygosity (H\u2092)",
    #title = "ROH vs heterozygosity across individuals"
  ) +
  theme_bw(base_size = 13)

ggsave("roh_vs_heterozygosity_scatter.png", p_scatter,
       width = 7, height = 4, dpi = 300)
p_scatter

## 7.4 Distribution of ROH segment lengths (all isolates) # plot not used.
p_seg_hist <- ggplot(roh_seg, aes(x = seg_Mb)) +
  geom_histogram(bins = 40, fill = "blue2") +
  labs(
    x = "ROH segment length (Mb)",
    y = "Count",
    title = "Distribution of ROH segment lengths"
  ) +
  theme_bw(base_size = 12)

ggsave("roh_segment_length_hist.png", p_seg_hist,
       width = 5.5, height = 4.5, dpi = 300)
p_seg_hist

## 7.5 Boxplots of total ROH length by group
p_roh_box <- ggplot(roh_indiv,
                    aes(x = Group,
                        y = total_Mb_ROH,
                        fill = Group)) +
  geom_boxplot(outlier.size = 0.7, width = 0.6) +
  labs(
    x = NULL,
    y = "Total ROH length (Mb)",
    #title = "Total ROH by group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave("roh_total_by_group_boxplot.png", p_roh_box,
       width = 5.5, height = 4.5, dpi = 300)
p_roh_box

Foh_combined <- p_het | p_roh_tot
ggsave("Foh_combined.png",
       Foh_combined, width = 14, height = 5, dpi = 300)

Foh_combined1 <- p_roh_box | p_het_box | p_scatter
ggsave("Foh_combined1.png",
       Foh_combined1, width = 14, height = 5, dpi = 300)


## 7.5 Boxplots of hetero  by group
library(ggplot2)
library(dplyr)

# Make sure Group exists (from your left_join) and Ho is numeric
het_plot_df <- het_df %>%
  filter(!is.na(Group), !is.na(Ho)) %>%
  mutate(
    Group = droplevels(as.factor(Group)),
    Ho = as.numeric(Ho)
  )

p_het_box <- ggplot(het_plot_df, aes(x = Group, y = Ho, fill = Group)) +
  geom_boxplot(outlier.size = 0.7, width = 0.6) +
  labs(
    x = NULL,
    y = "Heterozygosity (Ho)",
    #title = "Heterozygosity by group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

ggsave("heterozygosity_by_group_boxplot.png", p_het_box,
       width = 5.5, height = 4.5, dpi = 300)

p_het_box

Fhet_combined <- p_het | p_het_box
ggsave("Fhet_combined.png",
       Fhet_combined, width = 14, height = 5, dpi = 300)





























pi_df <- read_tsv("pi_5000bp.windowed.pi", comment="#", show_col_types=FALSE)

pi_df2 <- pi_df %>%
  mutate(
    Total_SNPs = sum(N_VARIANTS),
    Percent = 100 * N_VARIANTS / Total_SNPs
  )

head(pi_df2)






## ================================
##  Genome-wide π and Tajima's D plots
##  using vcftools window outputs
## ================================

library(tidyverse)

## --------- USER SETTINGS ---------
# Change these paths / names as needed
pi_file      <- "pi_10000bp.windowed.pi"
tajima_file  <- "tajimaD_10000bp.Tajima.D"
window_size  <- 10000    # for axis labels
species_name <- "Sydowia sp."  # used in titles (optional)

## --------- READ DATA -------------
message("Reading π and Tajima's D files...")

pi_df <- read_tsv(pi_file, comment = "#", show_col_types = FALSE)
tajima_df <- read_tsv(tajima_file, comment = "#", show_col_types = FALSE)

## Expected columns:
## π:      CHROM, BIN_START, BIN_END, N_VARIANTS, PI
## Tajima: CHROM, BIN_START, N_SNPS, TajimaD

## --------- CLEAN SCAFFOLD NAMES ---------
## Extract just "scaffoldX" from things like "BS06_240_394-4_scaffold3"

pi_df <- pi_df %>%
  mutate(
    SCAFFOLD = sub(".*(scaffold[0-9]+)$", "\\1", CHROM)
  )

tajima_df <- tajima_df %>%
  mutate(
    SCAFFOLD = sub(".*(scaffold[0-9]+)$", "\\1", CHROM)
  )

## --------- BUILD GENOME-WIDE COORDINATES ---------
## 1) Get scaffold lengths from π table (BIN_END max per scaffold)
scaffold_sizes <- pi_df %>%
  group_by(SCAFFOLD) %>%
  summarise(max_pos = max(BIN_END, na.rm = TRUE), .groups = "drop") %>%
  arrange(SCAFFOLD) %>%
  mutate(
    offset = lag(cumsum(max_pos), default = 0),
    mid_pos = offset + max_pos / 2    # for axis label placement
  )


## 2) Add genome-wide position to π and Tajima tables
pi_plot_df <- pi_df %>%
  left_join(scaffold_sizes, by = "SCAFFOLD") %>%
  mutate(genome_pos = BIN_START + offset)

tajima_plot_df <- tajima_df %>%
  left_join(scaffold_sizes, by = "SCAFFOLD") %>%
  mutate(genome_pos = BIN_START + offset)




#keeping the biggest scallfold only ## --------- BUILD GENOME-WIDE COORDINATES (largest scaffolds only) ---------

# How many / how big?
n_scaffolds_keep <- 30        # e.g., top 30 by length (change as you like)
# OR use a length threshold instead:
# min_len_bp <- 500000        # only scaffolds ≥ 500 kb

# 1) Get all scaffold sizes
scaffold_sizes_all <- pi_df %>%
  group_by(SCAFFOLD) %>%
  summarise(max_pos = max(BIN_END, na.rm = TRUE), .groups = "drop")

# 2A) keep the N largest scaffolds
scaffold_sizes <- scaffold_sizes_all %>%
  arrange(desc(max_pos)) %>%           # largest first
  slice_head(n = n_scaffolds_keep)

# If you prefer a length cutoff instead, use this instead of 2A:
# scaffold_sizes <- scaffold_sizes_all %>%
#   filter(max_pos >= min_len_bp)

# 3) add offsets and midpoints *only for chosen scaffolds*
scaffold_sizes <- scaffold_sizes %>%
  arrange(SCAFFOLD) %>%                # order on x-axis (or arrange(desc(max_pos)))
  mutate(
    offset = lag(cumsum(max_pos), default = 0),
    mid_pos = offset + max_pos / 2
  )

# 4) Add genome-wide position to π and Tajima tables, subset to top scaffolds
pi_plot_df <- pi_df %>%
  filter(SCAFFOLD %in% scaffold_sizes$SCAFFOLD) %>%      # keep only big ones
  left_join(scaffold_sizes, by = "SCAFFOLD") %>%
  mutate(genome_pos = BIN_START + offset)

tajima_plot_df <- tajima_df %>%
  filter(SCAFFOLD %in% scaffold_sizes$SCAFFOLD) %>%      # keep only big ones
  left_join(scaffold_sizes, by = "SCAFFOLD") %>%
  mutate(genome_pos = BIN_START + offset)







## --------- LONG FORMAT FOR COMBINED PANEL ---------
long_df <- bind_rows(
  pi_plot_df  %>% transmute(SCAFFOLD, genome_pos, value = PI,     stat = "Pi"),
  tajima_plot_df %>% transmute(SCAFFOLD, genome_pos, value = TajimaD, stat = "TajimaD")
)

long_df$stat <- factor(long_df$stat, levels = c("Pi", "TajimaD"))

## Alternating scaffold colours for Manhattan-style look
long_df <- long_df %>%
  left_join(
    scaffold_sizes %>% transmute(SCAFFOLD, scaffold_index = row_number()),
    by = "SCAFFOLD"
  ) %>%
  mutate(
    col_group = scaffold_index %% 2  # 0 / 1 for alternating colour
  )

## A simple base theme
theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    strip.background = element_rect(fill = "grey90", colour = NA),
    legend.position = "none"
  )

## --------- 2. Combined Manhattan-style panel (Pi + TajimaD) ---------
p_panel <- ggplot(long_df, aes(x = genome_pos, y = value,
                               colour = factor(col_group))) +
  geom_point(alpha = 0.6, size = 0.8) +
  ## optional loess smoothing per scaffold & stat
  geom_smooth(aes(group = stat), method = "loess", se = FALSE, span = 0.3,
              linewidth = 0.5, colour = "black") +
  facet_grid(rows = vars(stat), scales = "free_y", switch = "y") +
  scale_x_continuous(
    breaks = scaffold_sizes$mid_pos,
    labels = scaffold_sizes$SCAFFOLD,
    expand = expansion(mult = c(0.001, 0.001))
  ) +
  labs(
    x = "Scaffold",
    y = "Statistic value",
    title = paste0(species_name, " – genome-wide diversity (", window_size / 1000, " kb windows)")
  ) +
  theme_pub +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0)
  ) +
  scale_colour_manual(values = c("#1f77b4", "#2ca02c"))  # two alternating colours

## Save and/or show
ggsave("diversity_panel_pi_tajimaD.png", p_panel,
       width = 10, height = 7, dpi = 300)
print(p_panel)

## --------- 3. Pi-only Manhattan plot ---------
p_pi <- ggplot(pi_plot_df, aes(x = genome_pos, y = PI,
                               colour = factor(as.integer(factor(SCAFFOLD)) %% 2))) +
  geom_point(alpha = 0.7, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3,
              linewidth = 0.5, colour = "black") +
  scale_x_continuous(
    breaks = scaffold_sizes$mid_pos,
    labels = scaffold_sizes$SCAFFOLD,
    expand = expansion(mult = c(0.001, 0.001))
  ) +
  labs(
    x = "Scaffold",
    y = "Nucleotide diversity (π)",
    title = paste0("Genome-wide nucleotide diversity (", window_size / 1000, " kb windows)")
  ) +
  theme_pub +
  scale_colour_manual(values = c("#1f77b4", "#2ca02c"))

ggsave("pi_manhattan.png", p_pi, width = 10, height = 4, dpi = 300)
print(p_pi)

## --------- 4. Tajima's D Manhattan plot ---------
p_td <- ggplot(tajima_plot_df, aes(x = genome_pos, y = TajimaD,
                                   colour = factor(as.integer(factor(SCAFFOLD)) %% 2))) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(alpha = 0.7, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3,
              linewidth = 0.5, colour = "black") +
  scale_x_continuous(
    breaks = scaffold_sizes$mid_pos,
    labels = scaffold_sizes$SCAFFOLD,
    expand = expansion(mult = c(0.001, 0.001))
  ) +
  labs(
    x = "Scaffold",
    y = "Tajima's D",
    title = paste0("Genome-wide Tajima's D (", window_size / 1000, " kb windows)")
  ) +
  theme_pub +
  scale_colour_manual(values = c("#1f77b4", "#2ca02c"))

ggsave("tajimaD_manhattan.png", p_td, width = 10, height = 4, dpi = 300)
print(p_td)

## --------- 5. Tajima's D distribution (hist + density) ---------
p_td_hist <- ggplot(tajima_df, aes(x = TajimaD)) +
  geom_histogram(bins = 40, fill = "skyblue", color = "black") +
  geom_vline(xintercept = mean(tajima_df$TajimaD, na.rm = TRUE),
             linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Tajima's D",
    y = "Window count",
    title = "Genome-wide Tajima's D distribution 10kb"
  )

ggsave("tajimaD_hist.png", p_td_hist, width = 5, height = 4, dpi = 300)
print(p_td_hist)

p_td_density <- ggplot(tajima_df, aes(x = TajimaD)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = mean(tajima_df$TajimaD, na.rm = TRUE),
             linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Tajima's D",
    y = "Density",
    title = "Tajima's D density"
  )

ggsave("tajimaD_density.png", p_td_density, width = 5, height = 4, dpi = 300)
print(p_td_density)











# 5kb plot

## ================================
##  Genome-wide π and Tajima's D plots
##  using vcftools window outputs
## ================================

library(tidyverse)

## --------- USER SETTINGS ---------
# Change these paths / names as needed
pi_file1      <- "Sydowiapi_5000bp.windowed.pi"
tajima_file1  <- "SydowiatajimaD_5000bp.Tajima.D"
window_size  <- 5000    # for axis labels
species_name <- "Sydowia sp."  # used in titles (optional)

## --------- READ DATA -------------
message("Reading π and Tajima's D files...")

pi_df1 <- read_tsv(pi_file1, comment = "#", show_col_types = FALSE)
tajima_df1 <- read_tsv(tajima_file1, comment = "#", show_col_types = FALSE)

## Expected columns:
## π:      CHROM, BIN_START, BIN_END, N_VARIANTS, PI
## Tajima: CHROM, BIN_START, N_SNPS, TajimaD

## --------- CLEAN SCAFFOLD NAMES ---------
## Extract just "scaffoldX" from things like "BS06_240_394-4_scaffold3"

pi_df1 <- pi_df1 %>%
  mutate(
    SCAFFOLD = sub(".*(scaffold[0-9]+)$", "\\1", CHROM)
  )

tajima_df1 <- tajima_df1 %>%
  mutate(
    SCAFFOLD = sub(".*(scaffold[0-9]+)$", "\\1", CHROM)
  )

## --------- BUILD GENOME-WIDE COORDINATES ---------
## 1) Get scaffold lengths from π table (BIN_END max per scaffold)
scaffold_sizes1 <- pi_df1 %>%
  group_by(SCAFFOLD) %>%
  summarise(max_pos = max(BIN_END, na.rm = TRUE), .groups = "drop") %>%
  arrange(SCAFFOLD) %>%
  mutate(
    offset = lag(cumsum(max_pos), default = 0),
    mid_pos = offset + max_pos / 2    # for axis label placement
  )

## 2) Add genome-wide position to π and Tajima tables
pi_plot_df1 <- pi_df1 %>%
  left_join(scaffold_sizes1, by = "SCAFFOLD") %>%
  mutate(genome_pos = BIN_START + offset)

tajima_plot_df1 <- tajima_df1 %>%
  left_join(scaffold_sizes1, by = "SCAFFOLD") %>%
  mutate(genome_pos = BIN_START + offset)



#keeping 30 largest scaffold
#keeping the biggest scallfold only ## --------- BUILD GENOME-WIDE COORDINATES (largest scaffolds only) ---------

# How many / how big?
n_scaffolds_keep <- 30        # e.g., top 30 by length (change as you like)
# OR use a length threshold instead:
# min_len_bp <- 500000        # only scaffolds ≥ 500 kb

# 1) Get all scaffold sizes
scaffold_sizes_all1 <- pi_df1 %>%
  group_by(SCAFFOLD) %>%
  summarise(max_pos = max(BIN_END, na.rm = TRUE), .groups = "drop")

# 2A) keep the N largest scaffolds
scaffold_sizes1 <- scaffold_sizes_all1 %>%
  arrange(desc(max_pos)) %>%           # largest first
  slice_head(n = n_scaffolds_keep)

# If you prefer a length cutoff instead, use this instead of 2A:
# scaffold_sizes <- scaffold_sizes_all %>%
#   filter(max_pos >= min_len_bp)

# 3) add offsets and midpoints *only for chosen scaffolds*
scaffold_sizes1 <- scaffold_sizes1 %>%
  arrange(SCAFFOLD) %>%                # order on x-axis (or arrange(desc(max_pos)))
  mutate(
    offset = lag(cumsum(max_pos), default = 0),
    mid_pos = offset + max_pos / 2
  )

# 4) Add genome-wide position to π and Tajima tables, subset to top scaffolds
pi_plot_df1 <- pi_df1 %>%
  filter(SCAFFOLD %in% scaffold_sizes1$SCAFFOLD) %>%      # keep only big ones
  left_join(scaffold_sizes1, by = "SCAFFOLD") %>%
  mutate(genome_pos = BIN_START + offset)

tajima_plot_df1 <- tajima_df1 %>%
  filter(SCAFFOLD %in% scaffold_sizes1$SCAFFOLD) %>%      # keep only big ones
  left_join(scaffold_sizes1, by = "SCAFFOLD") %>%
  mutate(genome_pos = BIN_START + offset)


## --------- LONG FORMAT FOR COMBINED PANEL ---------
long_df1 <- bind_rows(
  pi_plot_df1  %>% transmute(SCAFFOLD, genome_pos, value = PI,     stat = "Pi"),
  tajima_plot_df1 %>% transmute(SCAFFOLD, genome_pos, value = TajimaD, stat = "TajimaD")
)

long_df1$stat <- factor(long_df1$stat, levels = c("Pi", "TajimaD"))

## Alternating scaffold colours for Manhattan-style look
long_df1 <- long_df1 %>%
  left_join(
    scaffold_sizes1 %>% transmute(SCAFFOLD, scaffold_index = row_number()),
    by = "SCAFFOLD"
  ) %>%
  mutate(
    col_group = scaffold_index %% 2  # 0 / 1 for alternating colour
  )

## A simple base theme
theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    strip.background = element_rect(fill = "grey90", colour = NA),
    legend.position = "none"
  )

## --------- 2. Combined Manhattan-style panel (Pi + TajimaD) ---------
p_panel1 <- ggplot(long_df1, aes(x = genome_pos, y = value,
                               colour = factor(col_group))) +
  geom_point(alpha = 0.6, size = 0.8) +
  ## optional loess smoothing per scaffold & stat
  geom_smooth(aes(group = stat), method = "loess", se = FALSE, span = 0.3,
              linewidth = 0.5, colour = "black") +
  facet_grid(rows = vars(stat), scales = "free_y", switch = "y") +
  scale_x_continuous(
    breaks = scaffold_sizes1$mid_pos,
    labels = scaffold_sizes1$SCAFFOLD,
    expand = expansion(mult = c(0.001, 0.001))
  ) +
  labs(
    x = "Scaffold",
    y = "Statistic value",
    title = paste0(species_name, " – genome-wide diversity (", window_size / 1000, " kb windows)")
  ) +
  theme_pub +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0)
  ) +
  scale_colour_manual(values = c("#1f77b4", "#2ca02c"))  # two alternating colours

## Save and/or show
ggsave("diversity_panel_pi_tajimaD1.png", p_panel1,
       width = 10, height = 7, dpi = 300)
print(p_panel1)

tajimaall_combined <- p_panel1 | 
ggsave("Tajimaallkb.png", 
       tajimaall_combined, width = 14, height = 5, dpi = 300)
tajimaall_combined 

## --------- 3. Pi-only Manhattan plot ---------
p_pi1 <- ggplot(pi_plot_df1, aes(x = genome_pos, y = PI,
                               colour = factor(as.integer(factor(SCAFFOLD)) %% 2))) +
  geom_point(alpha = 0.7, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3,
              linewidth = 0.5, colour = "black") +
  scale_x_continuous(
    breaks = scaffold_sizes1$mid_pos,
    labels = scaffold_sizes1$SCAFFOLD,
    expand = expansion(mult = c(0.001, 0.001))
  ) +
  labs(
    x = "Scaffold",
    y = "Nucleotide diversity (π)",
    title = paste0("Genome-wide nucleotide diversity (", window_size / 1000, " kb windows)")
  ) +
  theme_pub +
  scale_colour_manual(values = c("#1f77b4", "#2ca02c"))

ggsave("pi_manhattan1.png", p_pi1, width = 10, height = 4, dpi = 300)
print(p_pi1)

## --------- 4. Tajima's D Manhattan plot ---------
p_td1 <- ggplot(tajima_plot_df1, aes(x = genome_pos, y = TajimaD,
                                   colour = factor(as.integer(factor(SCAFFOLD)) %% 2))) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(alpha = 0.7, size = 0.8) +
  geom_smooth(method = "loess", se = FALSE, span = 0.3,
              linewidth = 0.5, colour = "black") +
  scale_x_continuous(
    breaks = scaffold_sizes1$mid_pos,
    labels = scaffold_sizes1$SCAFFOLD,
    expand = expansion(mult = c(0.001, 0.001))
  ) +
  labs(
    x = "Scaffold",
    y = "Tajima's D",
    title = paste0("Genome-wide Tajima's D (", window_size / 1000, " kb windows)")
  ) +
  theme_pub +
  scale_colour_manual(values = c("#1f77b4", "#2ca02c"))

ggsave("tajimaD_manhattan1.png", p_td1, width = 10, height = 4, dpi = 300)
print(p_td1)

## --------- 5. Tajima's D distribution (hist + density) ---------
p_td_hist1 <- ggplot(tajima_df1, aes(x = TajimaD)) +
  geom_histogram(bins = 40, fill = "skyblue", color = "black") +
  geom_vline(xintercept = mean(tajima_df1$TajimaD, na.rm = TRUE),
             linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Tajima's D",
    y = "Window count",
    title = "Genome-wide Tajima's D distribution 5kb"
  )

tajimadistrib_combined <- p_td_hist1 | p_td_hist
ggsave("tajimadistrib_combined.png", 
       tajimadistrib_combined, width = 14, height = 5, dpi = 300)
tajimadistrib_combined 

ggsave("tajimaD_hist1.png", p_td_hist, width = 5, height = 4, dpi = 300)
print(p_td_hist1)

p_td_density <- ggplot(tajima_df, aes(x = TajimaD)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = mean(tajima_df$TajimaD, na.rm = TRUE),
             linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Tajima's D",
    y = "Density",
    title = "Tajima's D density"
  )

ggsave("tajimaD_density.png", p_td_density, width = 5, height = 4, dpi = 300)
print(p_td_density)

## --------- 6. π heatmap by scaffold & position ---------
## Coarsen positions to kb for prettier axes
pi_heat_df <- pi_plot_df %>%
  mutate(
    pos_kb = BIN_START / 1000
  )

p_heat <- ggplot(pi_heat_df, aes(x = pos_kb, y = SCAFFOLD, fill = PI)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Position (kb)",
    y = "Scaffold",
    fill = "π",
    title = paste0("π heatmap (", window_size / 1000, " kb windows)")
  )

ggsave("pi_heatmap.png", p_heat, width = 6, height = 8, dpi = 300)
print(p_heat)

message("All plots generated.")







#Group Analysis


library(tidyverse)

## ---- user names ----
pi_files <- c(
  GroupI  = "pi_South_LowModerate_5000bp.windowed.pi",
  GroupII = "pi_South_High_5000bp.windowed.pi",
  GroupIII  = "pi_North_LowModerate_5000bp.windowed.pi",
  GroupIV = "pi_North_High_5000bp.windowed.pi"
)

taj_files <- c(
  GroupI  = "tajima_South_LowModerate_5000bp.Tajima.D",
  GroupII = "tajima_South_High_5000bp.Tajima.D",
  GroupIII  = "tajima_North_LowModerate_5000bp.Tajima.D",
  GroupIV = "tajima_North_High_5000bp.Tajima.D")

## ---- helper to read + label ----
read_pi <- function(file, group){
  read_tsv(file, comment = "#", show_col_types = FALSE) %>%
    mutate(Group = group)
}

read_taj <- function(file, group){
  read_tsv(file, comment = "#", show_col_types = FALSE) %>%
    mutate(Group = group)
}

pi_all <- map2_dfr(pi_files, names(pi_files), read_pi)
taj_all <- map2_dfr(taj_files, names(taj_files), read_taj)

## optional: remove crazy outliers / NAs
pi_all  <- pi_all  %>% filter(!is.na(PI), PI > 0)
taj_all <- taj_all %>% filter(!is.na(TajimaD))

## ---- Boxplot for π ----
p_pi_box <- ggplot(pi_all, aes(x = Group, y = PI, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_log10() +                     # like Fig A with 1e-x axis
  labs(
    x = NULL,
    y = "Nucleotide diversity (π)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave("pi_group_boxplot.png", p_pi_box, width = 4, height = 4, dpi = 300)
print(p_pi_box)

## ---- Boxplot for Tajima's D ----
p_taj_box <- ggplot(taj_all, aes(x = Group, y = TajimaD, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(
    x = NULL,
    y = "Tajima's D"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave("tajimaD_group_boxplot.png", p_taj_box, width = 4, height = 4, dpi = 300)
print(p_taj_box)




library(tidyverse)
library(patchwork)

## All pairwise FST files (6 combinations)
fst_files <- c(
  "fst_North_High_vs_North_LowModerate_5000bp.windowed.weir.fst",
  "fst_North_High_vs_South_High_5000bp.windowed.weir.fst",
  "fst_North_High_vs_South_LowModerate_5000bp.windowed.weir.fst",
  "fst_North_LowModerate_vs_South_High_5000bp.windowed.weir.fst",
  "fst_North_LowModerate_vs_South_LowModerate_5000bp.windowed.weir.fst",
  "fst_South_High_vs_South_LowModerate_5000bp.windowed.weir.fst"
)

## Read and combine, add Fst and contrast labels
fst_all <- purrr::map_dfr(fst_files, function(f) {
  x <- readr::read_tsv(f, comment = "#", show_col_types = FALSE)
  
  x %>%
    mutate(
      Fst = dplyr::case_when(
        "WEIGHTED_FST" %in% names(x) ~ WEIGHTED_FST,
        "MEAN_FST"     %in% names(x) ~ MEAN_FST,
        TRUE                          ~ NA_real_
      ),
      contrast = dplyr::case_when(
        grepl("North_High_vs_North_LowModerate",  f) ~ "North High vs North LowMod",
        grepl("North_High_vs_South_High",         f) ~ "North High vs South High",
        grepl("North_High_vs_South_LowModerate",  f) ~ "North High vs South LowMod",
        grepl("North_LowModerate_vs_South_High",  f) ~ "North LowMod vs South High",
        grepl("North_LowModerate_vs_South_LowModerate", f) ~ "North LowMod vs South LowMod",
        grepl("South_High_vs_South_LowModerate",  f) ~ "South High vs South LowMod",
        TRUE ~ basename(f)
      ),
      SCAFFOLD = sub(".*(scaffold[0-9]+)$", "\\1", CHROM)
    ) %>%
    filter(!is.na(Fst))
})

## Per-scaffold sizes & offsets within each contrast
fst_all <- fst_all %>%
  group_by(contrast, SCAFFOLD) %>%
  mutate(max_pos = max(BIN_END, na.rm = TRUE)) %>%
  ungroup()

scaffold_sizes_all <- fst_all %>%
  distinct(contrast, SCAFFOLD, max_pos) %>%
  arrange(contrast, SCAFFOLD) %>%
  group_by(contrast) %>%
  mutate(
    offset = dplyr::lag(cumsum(max_pos), default = 0),
    mid_pos = offset + max_pos / 2
  ) %>%
  ungroup()



fst_plot_all <- fst_all %>%
  left_join(scaffold_sizes_all,
            by = c("contrast", "SCAFFOLD", "max_pos")) %>%
  mutate(
    genome_pos = BIN_START + offset,
    col_group  = as.integer(factor(SCAFFOLD)) %% 2
  )

fst_plot_all


## Common theme
theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    legend.position = "none"
  )

## Manhattan-style plot for all contrasts
p_fst_all <- ggplot(
  fst_plot_all,
  aes(x = genome_pos, y = Fst, colour = factor(col_group))
) +
  geom_point(alpha = 0.8, size = 0.8) +
  facet_wrap(~ contrast, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = scaffold_sizes_all$mid_pos,
    labels = scaffold_sizes_all$SCAFFOLD,
    expand = expansion(mult = c(0.001, 0.001))
  ) +
  labs(
    x = "Scaffold",
    y = "Fst",
    title = "Genome-wide Fst for pairwise population contrasts (5 kb windows)"
  ) +
  scale_colour_manual(values = c("grey40", "blue")) +
  theme_pub

ggsave("fst_all_contrasts_manhattan.png", p_fst_all,
       width = 10, height = 7, dpi = 300)
print(p_fst_all)

## Boxplot of Fst by contrast
p_fst_box <- ggplot(fst_all, aes(x = contrast, y = Fst, fill = contrast)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.9) +
  labs(
    x = NULL,
    y = "Fst",
    title = "Fst distribution by contrast"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "none")

ggsave("fst_group_boxplot.png", p_fst_box,
       width = 6, height = 4, dpi = 300)
print(p_fst_box)

## Optional combined panels
Fst_combined <- p_fst_all | p_fst_box
ggsave("Fst_combined.png", Fst_combined, width = 14, height = 5, dpi = 300)
Fst_combined










#Group Analysis


library(tidyverse)

## ---- user names ----
pi_files <- c(
  GroupI  = "pi_South_LowModerate_5000bp.windowed.pi",
  GroupII = "pi_South_High_5000bp.windowed.pi",
  GroupIII  = "pi_North_LowModerate_5000bp.windowed.pi",
  GroupIV = "pi_North_High_5000bp.windowed.pi"
)

taj_files <- c(
  GroupI  = "tajima_South_LowModerate_5000bp.Tajima.D",
  GroupII = "tajima_South_High_5000bp.Tajima.D",
  GroupIII  = "tajima_North_LowModerate_5000bp.Tajima.D",
  GroupIV = "tajima_North_High_5000bp.Tajima.D")

## ---- helper to read + label ----
read_pi <- function(file, group){
  read_tsv(file, comment = "#", show_col_types = FALSE) %>%
    mutate(Group = group)
}

read_taj <- function(file, group){
  read_tsv(file, comment = "#", show_col_types = FALSE) %>%
    mutate(Group = group)
}

pi_all <- map2_dfr(pi_files, names(pi_files), read_pi)
taj_all <- map2_dfr(taj_files, names(taj_files), read_taj)

## optional: remove crazy outliers / NAs
pi_all  <- pi_all  %>% filter(!is.na(PI), PI > 0)
taj_all <- taj_all %>% filter(!is.na(TajimaD))

## ---- Boxplot for π ----
p_pi_box <- ggplot(pi_all, aes(x = Group, y = PI, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_log10() +                     # like Fig A with 1e-x axis
  labs(
    x = NULL,
    y = "Nucleotide diversity (π)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave("pi_group_boxplot.png", p_pi_box, width = 4, height = 4, dpi = 300)
print(p_pi_box)

## ---- Boxplot for Tajima's D ----
p_taj_box <- ggplot(taj_all, aes(x = Group, y = TajimaD, fill = Group)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(
    x = NULL,
    y = "Tajima's D"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave("tajimaD_group_boxplot.png", p_taj_box, width = 4, height = 4, dpi = 300)
print(p_taj_box)






library(tidyverse)







