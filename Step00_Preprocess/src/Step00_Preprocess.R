# =============================================================================
# Step00_Preprocess.R — Proteomics data preprocessing
# Column cleaning → quantile normalization → batch correction →
# column dedup → outlier removal → missing-rate filtering → KNN imputation
#
# Input:  ../../20250610_PV_report.pg_matrix.csv
#         If CSV is absent, auto-converts from:
#         ../../20250610_wsy_PV_report.pg_matrix.tsv (reviewer-provided raw DIA export)
# Output: imputed_temporal_matrix.tsv
# =============================================================================

rm(list = ls())
setwd(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))))

library(tidyverse)
library(preprocessCore)

cat("Step00 — Preprocessing...\n")

# ---- Reviewer alternative: convert raw DIA-NN TSV to pipeline CSV ----
CSV_FILE <- "../../20250610_PV_report.pg_matrix.csv"
if (!file.exists(CSV_FILE)) {
  TSV_FILE <- "../../20250610_wsy_PV_report.pg_matrix.tsv"
  if (file.exists(TSV_FILE)) {
    cat("  Converting reviewer TSV to pipeline CSV...\n")
    raw <- read_tsv(TSV_FILE, show_col_types = FALSE,
                    na = character(),        # preserve empty strings
                    guess_max = 10000)
    # Drop DIA-NN metadata columns (Protein.Group, Protein.Names,
    # First.Protein.Description), keep Genes + sample intensities
    drop_cols <- intersect(c("Protein.Group", "Protein.Names",
                             "First.Protein.Description"), names(raw))
    raw <- raw[, setdiff(names(raw), drop_cols)]
    # Reverse Excel date corruption: OCT1 was auto-formatted to 1-Oct
    raw$Genes[raw$Genes == "OCT1"] <- "1-Oct"
    write_csv(raw, CSV_FILE, na = "")
    cat(sprintf("  -> %s (%d x %d)\n", CSV_FILE, nrow(raw), ncol(raw)))
  } else {
    stop(sprintf("Neither %s nor %s found.", CSV_FILE, TSV_FILE))
  }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Column cleaning and metadata
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pg_matrix_file <- "../../20250610_PV_report.pg_matrix.csv"
if (!file.exists(pg_matrix_file)) stop(sprintf("Input '%s' not found.", pg_matrix_file))

pg_matrix_raw <- read_csv(pg_matrix_file, show_col_types = FALSE)

# Clean column names
original_colnames <- colnames(pg_matrix_raw)
cleaned_colnames <- original_colnames
cleaned_colnames[2:length(cleaned_colnames)] <- str_extract(
  original_colnames[2:length(original_colnames)],
  "(?<=_)[bB]\\d+.*(?=\\.mzML)"
)
pg_matrix_clean <- pg_matrix_raw
colnames(pg_matrix_clean) <- cleaned_colnames
pg_matrix_clean <- pg_matrix_clean %>% rename(Genes = 1)

write_csv(pg_matrix_clean, "../output/cleaned_pg_matrix.csv")

# Re-read from disk for exact floating-point reproducibility
pg_matrix_clean <- read_csv("../output/cleaned_pg_matrix.csv", show_col_types = FALSE)
sample_ids <- colnames(pg_matrix_clean)[2:ncol(pg_matrix_clean)]

# Build sample metadata
run_order_inst_A <- c("b1", "b4", "b7", "b10", "b13", "b16",
                       "b11", "b12", "b17", "b18", "b21")
run_order_inst_B <- c("b2", "b3", "b5", "b6", "b8", "b9",
                       "b14", "b15", "b19", "b20", "b22")

metadata <- tibble(SampleID = sample_ids) %>%
  mutate(
    RunBatch   = str_extract(SampleID, "b\\d+"),
    Instrument = if_else(RunBatch %in% run_order_inst_A, "Instrument_A", "Instrument_B"),
    RunOrder   = case_when(
      Instrument == "Instrument_A" ~ match(RunBatch, run_order_inst_A),
      Instrument == "Instrument_B" ~ match(RunBatch, run_order_inst_B)
    ),
    PrepBatch  = {
      n <- as.numeric(str_extract(RunBatch, "\\d+"))
      case_when(
        n %in% 1:6   ~ "Prep_1",  n %in% 7:12  ~ "Prep_2",
        n %in% 13:18 ~ "Prep_3",  n %in% 19:22 ~ "Prep_4", TRUE ~ "Unknown")
    },
    Group = case_when(
      str_detect(SampleID, "pool") ~ "Pool", str_detect(SampleID, "QC") ~ "QC",
      str_detect(SampleID, "_C")   ~ "C",    str_detect(SampleID, "_H") ~ "H",
      str_detect(SampleID, "_A")   ~ "A",    TRUE ~ "Unknown"),
    Time = case_when(
      Group %in% c("C", "H", "A") ~ as.numeric(str_extract(SampleID, "(?<=_[CHA])\\d+")),
      TRUE ~ NA_real_),
    SampleType = case_when(
      Group == "Pool" ~ "Pool", Group == "QC" ~ "Quality_Control",
      TRUE ~ "Biological_Sample"),
    IsTechRep = str_ends(SampleID, "_rep"),
    IsBioRep  = str_detect(SampleID, "(_|-)\\d$")
  ) %>%
  arrange(Instrument, RunOrder) %>%
  select(SampleID, RunBatch, Instrument, RunOrder, PrepBatch,
         Group, Time, SampleType, IsTechRep, IsBioRep)

write_csv(metadata, "../output/experiment_metadata.csv")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Quantile normalization
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pg_matrix  <- read_csv("../output/cleaned_pg_matrix.csv", show_col_types = FALSE)
metadata   <- read_csv("../output/experiment_metadata.csv", show_col_types = FALSE)

metadata <- metadata %>% mutate(Sample_Order = 1:nrow(.))

# Gene deduplication
num_cols <- pg_matrix %>% select(where(is.numeric)) %>% names()
pg_mr <- pg_matrix %>%
  mutate(.mr = rowSums(is.na(select(., all_of(num_cols)))) / length(num_cols),
         .idx = seq_len(n()))

dup_genes_b <- unique(pg_mr$Genes[duplicated(pg_mr$Genes)])

if (length(dup_genes_b) > 0) {
  keep_idx <- pg_mr %>%
    group_by(Genes) %>%
    summarise(.idx = .idx[which.min(.mr)], .groups = "drop") %>%
    pull(.idx)
  pg_matrix <- pg_mr[sort(keep_idx), ] %>% select(-.mr, -.idx)
  cat(sprintf("  Gene dedup: %d removed\n", nrow(pg_mr) - nrow(pg_matrix)))
} else {
  pg_matrix <- pg_mr %>% select(-.mr, -.idx)
}

expr_data <- pg_matrix %>%
  column_to_rownames("Genes") %>%
  select(all_of(metadata$SampleID)) %>%
  as.matrix()

expr_log2 <- log2(expr_data)
expr_log2[is.infinite(expr_log2)] <- NA

# Quantile normalization
expr_normalized <- normalize.quantiles(expr_log2)
rownames(expr_normalized) <- rownames(expr_log2)
colnames(expr_normalized) <- colnames(expr_log2)

normalized_df <- as.data.frame(expr_normalized) %>%
  rownames_to_column("Genes")

write_csv(normalized_df, "../output/normalized_pg_matrix.csv")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Batch correction
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

expr_normalized_df <- read_csv("../output/normalized_pg_matrix.csv", show_col_types = FALSE)
metadata            <- read_csv("../output/experiment_metadata.csv", show_col_types = FALSE)

# Gene deduplication
numeric_cols <- expr_normalized_df %>% select(where(is.numeric)) %>% names()
if (length(numeric_cols) == 0) {
  sample_cols <- setdiff(names(expr_normalized_df),
                         c("Genes", "Protein.Group", "Precursor.Id", "Genes"))
  numeric_cols <- sample_cols
}

df_with_missing_rate <- expr_normalized_df %>%
  mutate(missing_rate = rowSums(is.na(select(., all_of(numeric_cols)))) /
           length(numeric_cols))

duplicate_genes <- df_with_missing_rate %>%
  group_by(Genes) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(count > 1) %>%
  pull(Genes)

if (length(duplicate_genes) > 0) {
  cat(sprintf("  Stage-C gene dedup: %d removed\n", length(duplicate_genes)))

  deduplicated_rows <- df_with_missing_rate %>%
    filter(Genes %in% duplicate_genes) %>%
    group_by(Genes) %>%
    arrange(missing_rate, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()

  non_duplicate_rows <- df_with_missing_rate %>%
    filter(!(Genes %in% duplicate_genes))

  expr_normalized_df <- bind_rows(non_duplicate_rows, deduplicated_rows) %>%
    arrange(Genes) %>%
    select(-missing_rate)
} else {
  expr_normalized_df <- df_with_missing_rate %>% select(-missing_rate)
}

# Build expression matrix
expr_matrix <- expr_normalized_df %>%
  column_to_rownames("Genes") %>%
  select(all_of(metadata$SampleID)) %>%
  as.matrix()

# --- Step 1: Pool-based RunBatch correction ---

pool_samples_all <- metadata %>%
  filter(SampleType == "Pool") %>%
  pull(SampleID)
global_pool_anchor <- rowMeans(expr_matrix[, pool_samples_all], na.rm = TRUE)

unique_run_batches <- unique(metadata$RunBatch)
data_step1 <- matrix(NA, nrow = nrow(expr_matrix), ncol = ncol(expr_matrix),
                     dimnames = list(rownames(expr_matrix), colnames(expr_matrix)))

for (rb in unique_run_batches) {
  pool_in_batch <- metadata %>%
    filter(RunBatch == rb, SampleType == "Pool") %>%
    pull(SampleID)

  samples_in_batch <- metadata %>% filter(RunBatch == rb) %>% pull(SampleID)

  if (length(pool_in_batch) == 0) {
    cat(sprintf("    Warning: No Pool in %s, skipping.\n", rb))
    data_step1[, samples_in_batch] <- expr_matrix[, samples_in_batch]
    next
  }

  local_pool <- expr_matrix[, pool_in_batch[1]]
  correction <- local_pool - global_pool_anchor

  for (sid in samples_in_batch) {
    data_step1[, sid] <- expr_matrix[, sid] - correction
  }
}
# --- Step 2: RunBatch-level QC correction ---

data_step2 <- data_step1

for (rb in unique_run_batches) {
  qc_in_rb <- metadata %>%
    filter(RunBatch == rb, SampleType == "Quality_Control") %>%
    pull(SampleID)

  if (length(qc_in_rb) == 0) {
    cat(sprintf("    Warning: No QC in %s, skipping.\n", rb))
    next
  }

  qc_ref <- rowMeans(data_step1[, qc_in_rb, drop = FALSE], na.rm = TRUE)
  bio_in_rb <- metadata %>%
    filter(RunBatch == rb, SampleType == "Biological_Sample") %>%
    pull(SampleID)

  if (length(bio_in_rb) > 0) {
    for (sid in bio_in_rb) {
      data_step2[, sid] <- data_step1[, sid] - qc_ref
    }
  }
}
# --- Save batch-corrected matrix ---

non_pool_samples <- metadata %>% filter(SampleType != "Pool") %>% pull(SampleID)
expr_final_biological <- data_step2[, non_pool_samples]

final_corrected_df <- as.data.frame(expr_final_biological) %>%
  rownames_to_column("Genes") %>%
  left_join(select(expr_normalized_df, Genes), by = "Genes") %>%
  select(Genes, everything())

write_tsv(final_corrected_df, "../output/batch_corrected_matrix.tsv")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Column dedup and filtering
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data <- read_tsv("../output/batch_corrected_matrix.tsv", show_col_types = FALSE)

# Remove QC columns
data_no_qc <- data %>% select(-contains("QC"))
write_tsv(data_no_qc, "../output/batch_corrected_no_qc.tsv")

# Column dedup
sample_cols <- setdiff(colnames(data_no_qc), "Genes")

extract_core_name <- function(col_name) {
  core <- str_remove(col_name, "^b\\d+_")
  core <- str_remove(core, "_(rep|\\d+)$")
  core <- str_remove(core, "-\\d+$")
  core <- str_remove(core, "_\\d+$")
  return(core)
}

sample_mapping <- data.frame(
  original = sample_cols,
  core     = sapply(sample_cols, extract_core_name, USE.NAMES = FALSE)
)

duplicate_cores <- sample_mapping %>%
  group_by(core) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1) %>%
  pull(core)


final_cols  <- character()
cols_to_keep <- character()

for (cn in duplicate_cores) {
  group_cols <- sample_mapping %>% filter(core == cn) %>% pull(original)
  na_rates <- sapply(group_cols, function(col) {
    sum(is.na(data_no_qc[[col]])) / nrow(data_no_qc)
  })
  best <- group_cols[which.min(na_rates)]
  cols_to_keep <- c(cols_to_keep, best)
  final_cols   <- c(final_cols, cn)
}

non_duplicate_cols <- sample_mapping %>%
  filter(!core %in% duplicate_cores) %>%
  pull(original)

for (col in non_duplicate_cols) {
  cn <- extract_core_name(col)
  cols_to_keep <- c(cols_to_keep, col)
  final_cols   <- c(final_cols, cn)
}

data_final <- data_no_qc %>% select(Genes, all_of(cols_to_keep))
names(data_final)[-1] <- final_cols

# Sort columns by group and time
sample_names <- setdiff(colnames(data_final), "Genes")

sort_samples <- function(sample_names) {
  data.frame(name = sample_names) %>%
    mutate(
      group  = str_extract(name, "^[CHA]"),
      number = as.numeric(str_extract(name, "\\d+"))
    ) %>%
    arrange(factor(group, levels = c("C", "H", "A")), number) %>%
    pull(name)
}

sorted_samples <- sort_samples(sample_names)
data_sorted <- data_final %>% select(Genes, all_of(sorted_samples))

write_tsv(data_sorted, "../output/dedup_sorted_matrix.tsv")

# Remove outliers
outliers_to_remove <- c("H53", "C80", "A48", "A44")
E02_data <- data_sorted %>% select(-any_of(outliers_to_remove))
write_tsv(E02_data, "../output/outliers_removed_matrix.tsv")
cat(sprintf("  Outliers removed: %s\n", paste(intersect(outliers_to_remove, colnames(data_sorted)), collapse = ", ")))

# Filter by missing rate
missing_rates <- E02_data %>%
  select(-Genes) %>%
  rowwise() %>%
  summarise(missing_rate = sum(is.na(c_across(everything()))) /
              (ncol(E02_data) - 1)) %>%
  pull(missing_rate)

keep_rows <- missing_rates <= 0.40
E03_final_matrix <- E02_data[keep_rows, ]

write_tsv(E03_final_matrix, "../output/filtered_matrix.tsv")
cat(sprintf("  Missing-rate filter: %d proteins removed\n", sum(!keep_rows)))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. K-means imputation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

E03_data <- read_tsv("../output/filtered_matrix.tsv", show_col_types = FALSE)

expr_matrix <- E03_data %>%
  column_to_rownames("Genes") %>%
  as.matrix()

# K-means imputation
K <- 5
complete_rows <- complete.cases(expr_matrix)
complete_data <- expr_matrix[complete_rows, , drop = FALSE]

if (nrow(complete_data) < K) {
  K <- nrow(complete_data)
  cat(sprintf("    K reduced to %d.\n", K))
}

set.seed(123)
km <- kmeans(complete_data, centers = K, nstart = 25)
cluster_centers <- km$centers

imputed_matrix <- expr_matrix
missing_count <- 0

for (i in 1:nrow(expr_matrix)) {
  if (any(is.na(expr_matrix[i, ]))) {
    row_data <- expr_matrix[i, ]
    available_cols <- !is.na(row_data)

    if (sum(available_cols) > 0) {
      distances <- numeric(K)
      for (j in 1:K) {
        distances[j] <- sqrt(sum((row_data[available_cols] -
                                   cluster_centers[j, available_cols])^2,
                                 na.rm = TRUE))
      }
      closest_cluster <- which.min(distances)
      missing_positions <- is.na(row_data)
      imputed_matrix[i, missing_positions] <-
        cluster_centers[closest_cluster, missing_positions]
      missing_count <- missing_count + sum(missing_positions)
    }
  }
}

# Save imputed matrix
imputed_temporal <- imputed_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Genes")

write_tsv(imputed_temporal, "../output/imputed_temporal_matrix.tsv")

cat(sprintf("  K-means imputation: %d values filled, %d proteins x %d samples\n",
            missing_count, nrow(imputed_temporal), ncol(imputed_temporal) - 1))
cat("Step00 — Done.\n")
