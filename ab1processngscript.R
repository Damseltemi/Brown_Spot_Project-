# R Script for Processing Sanger Sequencing (AB1) Files

# --- 1. Installation of Required Packages ---
# This section checks for and installs any missing packages.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c("sangerseqR", "Biostrings", "DECIPHER", "dplyr", "ggplot2", "stringr")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# --- 2. Load Libraries ---

suppressPackageStartupMessages({
  library(sangerseqR)
  library(Biostrings)
  library(DECIPHER)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(tools)
})

# --- 3. User Settings ---
# Adjust these parameters to fit your project needs.

# ---- 1. User Settings ----
in_dir       <- "ITSNew1023"           # Folder containing your .ab1 files
out_dir      <- "New1023_ITS_out"      # Main output folder for results

# Primer trimming (optional)
trim_primer  <- TRUE
its5_forward <- "GGAAGTAAAAGTCGTAACAAGG"  # ITS5 primer sequence

# Minimum sequence length (after optional trimming)
min_len      <- 250

# ---- 2. Setup Output Folders ----
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fa_dir <- file.path(out_dir, "per_read_fasta")
dir.create(fa_dir, showWarnings = FALSE)

# Clean old FASTAs so ONLY current PASS reads appear
unlink(list.files(fa_dir, pattern = "\\.fasta$", full.names = TRUE), force = TRUE)

# ---- 3. Find AB1 files ----
ab1_files <- list.files(in_dir, pattern = "\\.ab1$", full.names = TRUE, ignore.case = TRUE)
if (length(ab1_files) == 0) {
  stop(paste("No .ab1 files found in the directory:", in_dir))
}
message(sprintf("Found %d AB1 files to process.", length(ab1_files)))

# ---- 4. Helper: read simple primary sequence ----
process_ab1_simple <- function(file_path) {
  ab1_data   <- read.abif(file_path)
  sanger_obj <- sangerseq(ab1_data)
  seq_char   <- as.character(primarySeq(sanger_obj))
  qual_scores <- NULL
  if (!is.null(ab1_data@data$PCON.2)) qual_scores <- ab1_data@data$PCON.2
  list(sequence_char = seq_char,
       quality       = qual_scores,
       raw_len       = nchar(seq_char))
}

# ---- 5. Processing Loop ----
qc_results    <- list()
all_sequences <- character()
all_names     <- character()

for (f in ab1_files) {
  sample_id <- tools::file_path_sans_ext(basename(f))
  msg_head  <- sprintf("[%s] ", sample_id)
  
  tryCatch({
    pr       <- process_ab1_simple(f)
    seq_char <- pr$sequence_char
    
    # Optional primer trimming (simple exact match)
    if (trim_primer && nchar(its5_forward) >= 10) {
      primer_pos <- regexpr(its5_forward, seq_char, ignore.case = TRUE)
      if (primer_pos[1] > 0) {
        start_pos <- primer_pos[1] + attr(primer_pos, "match.length")
        if (start_pos <= nchar(seq_char)) {
          seq_char <- substr(seq_char, start_pos, nchar(seq_char))
        }
      }
    }
    
    # Metrics
    len      <- nchar(seq_char)
    n_counts <- str_count(seq_char, "[Nn]")
    n_frac   <- ifelse(len > 0, n_counts / len, 1)
    
    # (Optional) mean quality if trace has PCON.2; not used for filtering here
    mean_qual <- NA_real_
    if (!is.null(pr$quality) && length(pr$quality) > 0) {
      qual_len <- min(length(pr$quality), len)
      if (qual_len > 0) mean_qual <- mean(pr$quality[1:qual_len], na.rm = TRUE)
    }
    
    # ---- Length filter (hard) ----
    if (len < min_len) {
      qc_results[[sample_id]] <- tibble(
        sample       = sample_id,
        raw_len      = pr$raw_len,
        final_len    = len,
        mean_quality = ifelse(is.na(mean_qual), NA, round(mean_qual, 2)),
        N_count      = n_counts,
        N_frac       = round(n_frac, 4),
        status       = "TooShort"
      )
      message(msg_head, "SKIPPED (", len, " bp < ", min_len, ")")
      next  # IMPORTANT: do not write any FASTA for short reads
    }
    
    # Passed → record QC and write per-read FASTA
    qc_results[[sample_id]] <- tibble(
      sample       = sample_id,
      raw_len      = pr$raw_len,
      final_len    = len,
      mean_quality = ifelse(is.na(mean_qual), NA, round(mean_qual, 2)),
      N_count      = n_counts,
      N_frac       = round(n_frac, 4),
      status       = "Processed"
    )
    
    # Write individual FASTA (PASS only)
    individual_seq <- DNAStringSet(seq_char)
    names(individual_seq) <- sample_id
    writeXStringSet(individual_seq, file.path(fa_dir, paste0(sample_id, ".fasta")))
    
    # Add to combined (PASS only)
    all_sequences <- c(all_sequences, seq_char)
    all_names     <- c(all_names, sample_id)
    
    qual_msg <- ifelse(is.na(mean_qual), "No qual", paste("AvgQual:", round(mean_qual, 1)))
    message(msg_head, "OK (", len, " bp, ", qual_msg, ")")
    
  }, error = function(e) {
    message(msg_head, "SKIPPED (", conditionMessage(e), ")")
    qc_results[[sample_id]] <<- tibble(
      sample        = sample_id,
      raw_len       = NA, final_len = NA, mean_quality = NA,
      N_count       = NA, N_frac = NA,
      status        = "Failed",
      error_message = str_trunc(conditionMessage(e), 100)
    )
  })
}

# ---- 6. Output & Summary ----
final_qc_table <- bind_rows(qc_results)

# Combined FASTA (PASS only)
if (length(all_sequences) > 0) {
  combined <- DNAStringSet(all_sequences); names(combined) <- all_names
  writeXStringSet(combined, file.path(out_dir, "all_sequences.fasta"))
  message(sprintf("\nSuccessfully processed %d sequences (≥ %d bp).", length(all_sequences), min_len))
} else {
  message("\nNo sequences passed the length filter.")
}

# CSV report
report_path <- file.path(out_dir, "processing_report.csv")
write.csv(final_qc_table, report_path, row.names = FALSE)
message(sprintf("Processing report saved to: %s", report_path))

# ---- 7. Plots (PASS only) ----
processed_data <- final_qc_table %>%
  filter(status == "Processed", !is.na(final_len))

if (nrow(processed_data) > 0) {
  # Length histogram
  p_len <- ggplot(processed_data, aes(x = final_len)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
    theme_minimal() +
    labs(
      title = "Sequence Length Distribution (PASS)",
      subtitle = paste("Reads ≥", min_len, "bp"),
      x = "Sequence Length (bp)", y = "Count"
    )
  ggsave(file.path(out_dir, "length_distribution.png"), p_len, width = 8, height = 6)
  
  # Length vs mean quality (if available)
  qual_data <- processed_data[!is.na(processed_data$mean_quality), ]
  if (nrow(qual_data) > 5) {
    p_q <- ggplot(qual_data, aes(x = final_len, y = mean_quality)) +
      geom_point(alpha = 0.6, color = "darkgreen") +
      geom_smooth(method = "loess", se = FALSE, color = "red") +
      theme_minimal() +
      labs(
        title = "Sequence Length vs Mean Quality (PASS)",
        x = "Length (bp)", y = "Mean Phred Quality"
      )
    ggsave(file.path(out_dir, "length_vs_quality.png"), p_q, width = 8, height = 6)
  }
}

# ---- 8. Console Summary ----
cat("\n=== SUMMARY ===\n")
cat("Total files:", nrow(final_qc_table), "\n")
cat("Processed (≥", min_len, "bp):", sum(final_qc_table$status == "Processed", na.rm = TRUE), "\n")
cat("Too short (<", min_len, "bp):", sum(final_qc_table$status == "TooShort", na.rm = TRUE), "\n")
cat("Failed:", sum(final_qc_table$status == "Failed", na.rm = TRUE), "\n")








# R Script for Processing Sanger Sequencing (AB1) Files - NO FILTERING VERSION

# --- 1. Installation of Required Packages ---
# This section checks for and installs any missing packages.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c("sangerseqR", "Biostrings", "dplyr", "ggplot2", "stringr")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# --- 2. Load Libraries ---

suppressPackageStartupMessages({
  library(sangerseqR)
  library(Biostrings)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# --- 3. User Settings ---
# Adjust these parameters to fit your project needs.

in_dir       <- "merged_samples"           # Folder containing your .ab1 files
out_dir      <- "sanger_ITS_out_no_filter"  # Main output folder for results

# Primer Trimming Parameters (optional)
trim_primer  <- FALSE  # Set to FALSE to skip primer trimming
its5_forward <- "GGAAGTAAAAGTCGTAACAAGG" # ITS5 primer sequence

# --- 4. Setup Environment ---

# Create output directories
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fa_dir <- file.path(out_dir, "per_read_fasta")
dir.create(fa_dir, showWarnings = FALSE)

# Find all .ab1 files
ab1_files <- list.files(in_dir, pattern = "\\.ab1$", full.names = TRUE, ignore.case = TRUE)

if (length(ab1_files) == 0) {
  stop(paste("No .ab1 files found in the directory:", in_dir))
}

message(sprintf("Found %d AB1 files to process.", length(ab1_files)))

# --- 5. Simple Processing Function ---

process_ab1_simple <- function(file_path) {
  
  # Read the ABIF file and create sangerseq object
  ab1_data <- read.abif(file_path)
  sanger_obj <- sangerseq(ab1_data)
  
  # Get the primary sequence directly
  primary_seq <- primarySeq(sanger_obj)
  
  # Convert to DNAString
  dna <- DNAString(as.character(primary_seq))
  
  # Get quality scores if available
  qual_scores <- NULL
  if (!is.null(ab1_data@data$PCON.2)) {
    qual_scores <- ab1_data@data$PCON.2
  }
  
  return(list(sequence = dna, quality = qual_scores, raw_len = length(dna)))
}

# --- 6. Processing Loop ---

# Initialize containers for results
qc_results <- list()
all_clean_seqs <- DNAStringSet()

for (f in ab1_files) {
  sample_id <- tools::file_path_sans_ext(basename(f))
  msg_head  <- sprintf("[%s] ", sample_id)
  
  tryCatch({
    # 1. Read and process the AB1 file
    processed_read <- process_ab1_simple(f)
    dna <- processed_read$sequence
    
    # 2. Optional primer trimming (simple string matching)
    if (trim_primer && nchar(its5_forward) >= 10) {
      # Simple primer removal from 5' end
      seq_char <- as.character(dna)
      primer_pos <- regexpr(its5_forward, seq_char, ignore.case = TRUE)
      if (primer_pos[1] > 0) {
        # Remove primer and everything before it
        start_pos <- primer_pos[1] + attr(primer_pos, "match.length")
        if (start_pos <= nchar(seq_char)) {
          seq_char <- substr(seq_char, start_pos, nchar(seq_char))
          dna <- DNAString(seq_char)
        }
      }
    }
    
    # 3. Calculate basic metrics (no filtering applied)
    len <- length(dna)
    n_counts <- letterFrequency(dna, letters = "N")
    n_frac <- ifelse(len > 0, n_counts / len, 1)
    
    # Calculate mean quality if available
    mean_qual <- NA
    if (!is.null(processed_read$quality) && length(processed_read$quality) > 0) {
      # Adjust quality vector length to match sequence length if needed
      qual_len <- min(length(processed_read$quality), len)
      if (qual_len > 0) {
        mean_qual <- mean(processed_read$quality[1:qual_len], na.rm = TRUE)
      }
    }
    
    # Store metrics (all sequences are kept - no filtering)
    qc_results[[sample_id]] <- tibble(
      sample = sample_id,
      raw_len = processed_read$raw_len,
      final_len = len,
      mean_quality = ifelse(is.na(mean_qual), NA, round(mean_qual, 2)),
      N_count = as.numeric(n_counts),
      N_frac  = round(n_frac, 4),
      status = "Processed"
    )
    
    # 4. Save Individual FASTA File
    names(dna) <- sample_id
    out_fa <- file.path(fa_dir, paste0(sample_id, ".fasta"))
    writeXStringSet(DNAStringSet(dna), out_fa)
    
    # 5. Add to the combined set
    all_clean_seqs <- append(all_clean_seqs, DNAStringSet(dna))
    
    qual_msg <- ifelse(is.na(mean_qual), "No qual", paste("AvgQual:", round(mean_qual, 1)))
    message(msg_head, "OK (", length(dna), " bp, ", qual_msg, ")")
    
  }, error = function(e) {
    # Error handling
    message(msg_head, "SKIPPED (", conditionMessage(e), ")")
    # Add a failed entry to QC results
    qc_results[[sample_id]] <<- tibble(
      sample = sample_id,
      raw_len = NA, final_len = NA, mean_quality = NA,
      N_count = NA, N_frac = NA, status = "Failed",
      error_message = str_trunc(conditionMessage(e), 100)
    )
  })
}

# --- 7. Final Output Generation ---

# Combine all QC results into a single data frame
final_qc_table <- bind_rows(qc_results)

# Write the combined FASTA file
if (length(all_clean_seqs) > 0) {
  writeXStringSet(all_clean_seqs, file.path(out_dir, "all_sequences.fasta"))
  message(sprintf("\nSuccessfully processed %d sequences. Combined FASTA file saved.", length(all_clean_seqs)))
} else {
  message("\nNo sequences were successfully processed.")
}

# Write the processing report
report_path <- file.path(out_dir, "processing_report.csv")
write.csv(final_qc_table, report_path, row.names = FALSE)
message(sprintf("Processing report saved to: %s", report_path))

# Generate and save a summary plot
processed_data <- final_qc_table[final_qc_table$status == "Processed" & !is.na(final_qc_table$final_len), ]

if (nrow(processed_data) > 0) {
  summary_plot <- ggplot(processed_data, aes(x = final_len)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
    theme_minimal() +
    labs(
      title = "Sanger Sequencing Processing Summary",
      subtitle = paste("Successfully processed", nrow(processed_data), "of", nrow(final_qc_table), "sequences"),
      x = "Sequence Length (bp)",
      y = "Count"
    )
  
  # Add quality information if available
  if (any(!is.na(processed_data$mean_quality))) {
    summary_plot <- summary_plot +
      labs(caption = paste("Mean quality range:", 
                           round(min(processed_data$mean_quality, na.rm = TRUE), 1), "-",
                           round(max(processed_data$mean_quality, na.rm = TRUE), 1)))
  }
  
  plot_path <- file.path(out_dir, "processing_summary.png")
  ggsave(plot_path, summary_plot, width = 8, height = 6)
  message(sprintf("Summary plot saved to: %s", plot_path))
}

# Print summary statistics
cat("\n=== PROCESSING SUMMARY ===\n")
cat("Total files found:", nrow(final_qc_table), "\n")
cat("Successfully processed:", sum(final_qc_table$status == "Processed", na.rm = TRUE), "\n")
cat("Failed:", sum(final_qc_table$status == "Failed", na.rm = TRUE), "\n")

if (nrow(processed_data) > 0) {
  cat("Length range:", min(processed_data$final_len), "-", max(processed_data$final_len), "bp\n")
  cat("Mean length:", round(mean(processed_data$final_len), 1), "bp\n")
  if (any(!is.na(processed_data$mean_quality))) {
    cat("Quality range:", round(min(processed_data$mean_quality, na.rm = TRUE), 1), "-", 
        round(max(processed_data$mean_quality, na.rm = TRUE), 1), "\n")
  }
}

message("\nProcessing complete - no filtering applied.")










# --- 3. User Settings ---
in_dir       <- "Newplate10302025"           # Folder containing your .ab1 files
out_dir      <- "New10302025_ITS_out"      # Main output folder for results

# Primer Trimming Parameters (optional)
trim_primer  <- TRUE
its5_forward <- "GGAAGTAAAAGTCGTAACAAGG"

# Minimum length filter (after optional trimming)
min_len      <- 250

# --- 4. Setup Environment ---
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fa_dir <- file.path(out_dir, "per_read_fasta")
dir.create(fa_dir, showWarnings = FALSE)

ab1_files <- list.files(in_dir, pattern = "\\.ab1$", full.names = TRUE, ignore.case = TRUE)

# FIX: stop only if none are found
if (length(ab1_files) == 0) {
  stop(paste("No .ab1 files found in the directory:", in_dir))
}

message(sprintf("Found %d AB1 files to process.", length(ab1_files)))

# --- 5. Simple Processing Function ---
process_ab1_simple <- function(file_path) {
  ab1_data <- read.abif(file_path)
  sanger_obj <- sangerseq(ab1_data)
  primary_seq <- primarySeq(sanger_obj)
  seq_char <- as.character(primary_seq)
  dna <- DNAString(seq_char)
  
  qual_scores <- NULL
  if (!is.null(ab1_data@data$PCON.2)) qual_scores <- ab1_data@data$PCON.2
  
  list(sequence = dna, sequence_char = seq_char, quality = qual_scores, raw_len = nchar(seq_char))
}

# --- 6. Processing Loop ---
qc_results <- list()
all_sequences <- character()
all_names <- character()

for (f in ab1_files) {
  sample_id <- tools::file_path_sans_ext(basename(f))
  msg_head  <- sprintf("[%s] ", sample_id)
  
  tryCatch({
    pr <- process_ab1_simple(f)
    seq_char <- pr$sequence_char
    dna <- pr$sequence
    
    # Optional primer trimming
    if (trim_primer && nchar(its5_forward) >= 10) {
      primer_pos <- regexpr(its5_forward, seq_char, ignore.case = TRUE)
      if (primer_pos[1] > 0) {
        start_pos <- primer_pos[1] + attr(primer_pos, "match.length")
        if (start_pos <= nchar(seq_char)) {
          seq_char <- substr(seq_char, start_pos, nchar(seq_char))
          dna <- DNAString(seq_char)
        }
      }
    }
    
    # Metrics
    len <- nchar(seq_char)
    n_counts <- stringr::str_count(seq_char, "[Nn]")
    n_frac <- ifelse(len > 0, n_counts / len, 1)
    
    mean_qual <- NA
    if (!is.null(pr$quality) && length(pr$quality) > 0) {
      qual_len <- min(length(pr$quality), len)
      if (qual_len > 0) mean_qual <- mean(pr$quality[1:qual_len], na.rm = TRUE)
    }
    
    # --- LENGTH FILTER ---
    if (len < min_len) {
      qc_results[[sample_id]] <- tibble::tibble(
        sample = sample_id,
        raw_len = pr$raw_len,
        final_len = len,
        mean_quality = ifelse(is.na(mean_qual), NA, round(mean_qual, 2)),
        N_count = n_counts, N_frac = round(n_frac, 4),
        status = "TooShort"
      )
      message(msg_head, "SKIPPED (", len, " bp < ", min_len, ")")
      next
    }
    
    # Save QC for kept reads
    qc_results[[sample_id]] <- tibble::tibble(
      sample = sample_id,
      raw_len = pr$raw_len,
      final_len = len,
      mean_quality = ifelse(is.na(mean_qual), NA, round(mean_qual, 2)),
      N_count = n_counts, N_frac = round(n_frac, 4),
      status = "Processed"
    )
    
    # Write individual FASTA
    individual_seq <- Biostrings::DNAStringSet(seq_char)
    names(individual_seq) <- sample_id
    Biostrings::writeXStringSet(individual_seq, file.path(fa_dir, paste0(sample_id, ".fasta")))
    
    # Add to combined
    all_sequences <- c(all_sequences, seq_char)
    all_names <- c(all_names, sample_id)
    
    qual_msg <- ifelse(is.na(mean_qual), "No qual", paste("AvgQual:", round(mean_qual, 1)))
    message(msg_head, "OK (", len, " bp, ", qual_msg, ")")
    
  }, error = function(e) {
    message(msg_head, "SKIPPED (", conditionMessage(e), ")")
    qc_results[[sample_id]] <<- tibble::tibble(
      sample = sample_id,
      raw_len = NA, final_len = NA, mean_quality = NA,
      N_count = NA, N_frac = NA, status = "Failed",
      error_message = stringr::str_trunc(conditionMessage(e), 100)
    )
  })
}

# --- 7. Final Output Generation ---
final_qc_table <- dplyr::bind_rows(qc_results)

if (length(all_sequences) > 0) {
  combined <- Biostrings::DNAStringSet(all_sequences); names(combined) <- all_names
  Biostrings::writeXStringSet(combined, file.path(out_dir, "all_sequences.fasta"))
  message(sprintf("\nSuccessfully processed %d sequences (>= %d bp).", length(all_sequences), min_len))
} else {
  message("\nNo sequences passed the length filter.")
}

# Report
write.csv(final_qc_table, file.path(out_dir, "processing_report.csv"), row.names = FALSE)

# Summary plots use only kept reads
processed_data <- final_qc_table %>%
  dplyr::filter(status == "Processed", !is.na(final_qc_table$final_len))

# (Plot code unchanged…)


