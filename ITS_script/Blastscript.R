


# ----- packages -----
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# ===== 1) Read BLAST TSV =====
# make sure the file name matches yours
in_tsv  <- "blast_combined1.tsv"
all_csv <- "blast_combined1.csv"
top_csv <- "blast_top_hits1.csv"
summary_csv <- "blast_species_summary1.csv"

# read as tab-separated; don't convert strings to factors
blast <- read.delim(in_tsv, sep = "\t", header = TRUE, check.names = FALSE)

# (Optional) check required columns
required_cols <- c("qseqid","sseqid","pident","length","qlen","slen",
                   "qstart","qend","sstart","send","evalue","bitscore","stitle")
missing_cols <- setdiff(required_cols, names(blast))
if (length(missing_cols) > 0) {
  stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
}

# Coerce numeric columns defensively (in case they were read as character)
num_cols <- c("pident","length","qlen","slen","qstart","qend","sstart","send","evalue","bitscore")
blast <- blast %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.))))

# ===== 2) Write the full table as CSV =====
write.csv(blast, all_csv, row.names = FALSE)

# ===== 3) Pick top hit per query =====
# Tie-breakers: highest identity, then highest bitscore, then lowest evalue, then longest alignment
blast_top <- blast %>%
  group_by(qseqid) %>%
  arrange(desc(pident), desc(bitscore), evalue, desc(length), .by_group = TRUE) %>%
  slice(1L) %>%
  ungroup()

# ===== 4) Parse species and genus from stitle =====
# We'll try to capture a binomial like "Genus species" (optionally with subspecies/strain info).
# Strategy:
#  - Prefer the *first* good-looking "Genus species" (Capitalized Genus + lowercase epithet)
#  - Ignore words like "sp.", "cf.", numbers, etc.
#  - If no clean binomial is found, leave species as NA

extract_species <- function(title) {
  if (is.na(title) || title == "") return(NA_character_)
  # Try a fairly robust regex: Capitalized genus + space + lowercase epithet (letters only).
  # We also allow hyphens in epithet because some fungi include hyphens.
  m <- str_match(title, "\\b([A-Z][a-z]+)\\s+([a-z][a-z-]+)\\b")
  if (!is.na(m[1,1])) {
    species <- paste(m[1,2], m[1,3])
    # Filter out cases like "Genus sp", "Genus cf", etc.
    if (grepl("\\b(sp|cf|aff)\\b", species)) return(NA_character_)
    return(species)
  }
  NA_character_
}

blast_top <- blast_top %>%
  mutate(
    species = vapply(stitle, extract_species, character(1)),
    genus   = ifelse(!is.na(species), word(species, 1), NA_character_)
  )

# ===== 5) Write top hits CSV =====
write.csv(blast_top, top_csv, row.names = FALSE)

# ===== (Optional) species summary =====
species_summary <- blast_top %>%
  count(species, sort = TRUE) %>%
  rename(n_queries_with_top_hit = n)

write.csv(species_summary, summary_csv, row.names = FALSE)

# ===== Console messages =====
cat("\nWrote:\n",
    " - All hits CSV:    ", all_csv, "\n",
    " - Top hits CSV:    ", top_csv, "\n",
    " - Species summary: ", summary_csv, "\n\n", sep = "")




suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr)
})

# ===== 1) Read BLAST TSV =====
in_tsv      <- "blast1030_combined.tsv"
all_csv     <- "blast1030_combined.csv"
top_csv     <- "blast1030_top_hits.csv"
summary_csv <- "blast1030_species_summary.csv"

blast <- read.delim(in_tsv, sep = "\t", header = TRUE, check.names = FALSE)

required_cols <- c("qseqid","sseqid","pident","length","qlen","slen",
                   "qstart","qend","sstart","send","evalue","bitscore","stitle")
missing_cols <- setdiff(required_cols, names(blast))
if (length(missing_cols) > 0) stop("Missing expected columns: ", paste(missing_cols, collapse=", "))

num_cols <- c("pident","length","qlen","slen","qstart","qend","sstart","send","evalue","bitscore")
blast <- blast %>% mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.))))

# ===== 2) Write the full table as CSV =====
write.csv(blast, all_csv, row.names = FALSE)

# ===== species extractor (used for both top hits and global summary) =====
extract_species <- function(title){
  if (is.na(title) || title == "") return(NA_character_)
  m <- str_match(title, "\\b([A-Z][a-z]+)\\s+([a-z][a-z-]+)\\b")
  if (!is.na(m[1,1])) {
    sp <- paste(m[1,2], m[1,3])
    if (grepl("\\b(sp|cf|aff)\\b", sp)) return(NA_character_)
    return(sp)
  }
  NA_character_
}

# ===== 3) Top hit per query (unchanged) =====
blast_top <- blast %>%
  group_by(qseqid) %>%
  arrange(desc(pident), desc(bitscore), evalue, desc(length), .by_group = TRUE) %>%
  slice(1L) %>%
  ungroup() %>%
  mutate(
    species = vapply(stitle, extract_species, character(1)),
    genus   = ifelse(!is.na(species), word(species, 1), NA_character_)
  )

write.csv(blast_top, top_csv, row.names = FALSE)

# ===== 4) GLOBAL species summary across ALL hits (this is the new behavior) =====
species_summary <- blast %>%
  mutate(species = vapply(stitle, extract_species, character(1))) %>%
  filter(!is.na(species)) %>%
  summarise(.by = species,
            n_hits = n(),
            n_unique_queries = n_distinct(qseqid)) %>%
  arrange(desc(n_hits))

write.csv(species_summary, summary_csv, row.names = FALSE)

cat("\nWrote:\n",
    " - All hits CSV:    ", all_csv, "\n",
    " - Top hits CSV:    ", top_csv, "\n",
    " - Species summary: ", summary_csv, " (ranked by total hits)\n\n", sep = "")

