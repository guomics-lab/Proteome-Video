# =============================================================================
# Step01_Overview.R — Fig.1 and Extended Fig.1: QC and global overview
#
# Extended Fig.1:
#   ef1a  Sample-pair Pearson correlation violin
#   ef1b  Protein identification count by treatment
#   ef1c  t-SNE before batch correction
#   ef1d  PVCA before batch correction
#   ef1e  t-SNE after batch correction
#   ef1f  PVCA after batch correction
#   ef1g  Stage-wise centroid distance permutation test
#
# Fig.1:
#   fig1a Flow schematic (Adobe Illustrator, no code)
#   fig1b Temporal evolution of centroid distances (LOESS)
#   fig1c PCoA trajectories across 5 stages
#
# Prerequisite: Run Step00_Preprocess.R first.
# =============================================================================

rm(list = ls())
setwd(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))))

library(tidyverse)
library(Rtsne)
library(pvca)
library(Biobase)
library(vegan)
library(patchwork)

cat("Step01 — Overview...\n")

# ---- Paths to Step00 outputs ----
METADATA       <- "../../Step00_Preprocess/output/experiment_metadata.csv"
NORMALIZED     <- "../../Step00_Preprocess/output/normalized_pg_matrix.csv"
OUTLIERS_REM   <- "../../Step00_Preprocess/output/outliers_removed_matrix.tsv"
CORRECTED_NOQC <- "../../Step00_Preprocess/output/batch_corrected_no_qc.tsv"
CORRECTED_FULL <- "../../Step00_Preprocess/output/batch_corrected_matrix.tsv"
IMPUTED_FINAL  <- "../../Step00_Preprocess/output/imputed_temporal_matrix.tsv"


# =============================================================================
# Extended Fig.1a — Sample-pair Pearson correlation violin
# =============================================================================

metadata <- read_csv(METADATA, show_col_types = FALSE)
expr_matrix <- read_csv(NORMALIZED, show_col_types = FALSE) %>%
  column_to_rownames("Genes") %>%
  as.matrix()

common_samples <- intersect(colnames(expr_matrix), metadata$SampleID)
expr_matrix <- expr_matrix[, common_samples]
metadata <- metadata %>% filter(SampleID %in% common_samples)

completeness_rate <- rowSums(!is.na(expr_matrix)) / ncol(expr_matrix)
expr_filtered <- expr_matrix[completeness_rate > 0.5, ]
global_min <- min(expr_filtered, na.rm = TRUE)
imputation_value <- if (global_min > 0) global_min * 0.8 else global_min * 1.1
expr_filtered[is.na(expr_filtered)] <- imputation_value

all_correlations <- tibble()

pool_samples <- metadata %>% filter(SampleType == "Pool") %>% pull(SampleID)
if (length(pool_samples) >= 2) {
  pool_pairs <- combn(pool_samples, 2, simplify = FALSE)
  pool_correlations <- map_dbl(pool_pairs, ~cor(expr_filtered[, .x[1]], expr_filtered[, .x[2]], method = "pearson"))
  all_correlations <- bind_rows(all_correlations, tibble(Correlation = pool_correlations, Category = "Pool"))
}

qc_samples <- metadata %>% filter(SampleType == "Quality_Control") %>% pull(SampleID)
if (length(qc_samples) >= 2) {
  qc_pairs <- combn(qc_samples, 2, simplify = FALSE)
  qc_correlations <- map_dbl(qc_pairs, ~cor(expr_filtered[, .x[1]], expr_filtered[, .x[2]], method = "pearson"))
  all_correlations <- bind_rows(all_correlations, tibble(Correlation = qc_correlations, Category = "Yeast QC Samples"))
}

potential_replicates_meta <- metadata %>% filter(IsTechRep == TRUE | IsBioRep == TRUE)
potential_originals_ids <- metadata %>% filter(IsTechRep == FALSE & IsBioRep == FALSE) %>% pull(SampleID)
tech_rep_pairs <- list()
bio_rep_pairs <- list()

for (i in seq_len(nrow(potential_replicates_meta))) {
  rep_row <- potential_replicates_meta[i, ]
  rep_id <- rep_row$SampleID
  potential_matches <- potential_originals_ids[startsWith(rep_id, potential_originals_ids)]
  if (length(potential_matches) > 0) {
    best_match_original <- potential_matches[which.max(nchar(potential_matches))]
    if (rep_row$IsTechRep == TRUE) tech_rep_pairs <- append(tech_rep_pairs, list(c(best_match_original, rep_id)))
    if (rep_row$IsBioRep == TRUE) bio_rep_pairs <- append(bio_rep_pairs, list(c(best_match_original, rep_id)))
  }
}

if (length(tech_rep_pairs) > 0) {
  tech_rep_correlations <- map_dbl(tech_rep_pairs, ~cor(expr_filtered[, .x[1]], expr_filtered[, .x[2]], method = "pearson"))
  all_correlations <- bind_rows(all_correlations, tibble(Correlation = tech_rep_correlations, Category = "Technical Replicates"))
}
if (length(bio_rep_pairs) > 0) {
  bio_rep_correlations <- map_dbl(bio_rep_pairs, ~cor(expr_filtered[, .x[1]], expr_filtered[, .x[2]], method = "pearson"))
  all_correlations <- bind_rows(all_correlations, tibble(Correlation = bio_rep_correlations, Category = "Biological Replicates"))
}

category_counts <- tibble(
  Category_Full = c("Pool", "Yeast QC Samples", "Technical Replicates", "Biological Replicates"),
  n = c(length(pool_samples), length(qc_samples), length(tech_rep_pairs), length(bio_rep_pairs))
)
plot_labels <- category_counts %>%
  mutate(label = paste0(Category_Full, "\nn=", n)) %>%
  pull(label, name = Category_Full)

category_order <- c("Pool", "Yeast QC Samples", "Technical Replicates", "Biological Replicates")
existing_categories <- intersect(category_order, unique(all_correlations$Category))
all_correlations$Category <- factor(all_correlations$Category, levels = existing_categories)

p_ef1a <- ggplot(all_correlations, aes(x = Category, y = Correlation)) +
  geom_violin(aes(fill = Category), trim = TRUE, color = "black", alpha = 1,
              draw_quantiles = c(0.25, 0.75)) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "red", color = "black") +
  scale_x_discrete(labels = plot_labels[existing_categories]) +
  scale_y_continuous(limits = c(0.65, 0.95), breaks = seq(0.65, 0.95, by = 0.1)) +
  labs(x = "", y = "Pearson Correlation") +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(legend.position = "none",
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", size = 14))

ggsave("../figure/ef1a.pdf", p_ef1a, width = 9, height = 6, device = cairo_pdf)


# =============================================================================
# Extended Fig.1b — Protein identification count by treatment
# =============================================================================

protein_data <- read_tsv(OUTLIERS_REM, show_col_types = FALSE)
sample_cols <- setdiff(colnames(protein_data), "Genes")
expression_matrix_b <- protein_data %>%
  column_to_rownames("Genes") %>%
  as.matrix()

sample_protein_counts <- map_dfr(sample_cols, function(sample_id) {
  protein_count <- sum(!is.na(expression_matrix_b[, sample_id]))
  treatment_char <- str_extract(sample_id, "^[CHA]")
  tibble(SampleID = sample_id, treatment = treatment_char, protein_count = protein_count)
}) %>%
  mutate(
    treatment_full = recode(treatment,
      "A" = "Acetic acid (A)",
      "C" = "Control (C)",
      "H" = "Hydrogen peroxide (H)"
    ),
    treatment_full = factor(treatment_full, levels = c("Acetic acid (A)", "Control (C)", "Hydrogen peroxide (H)"))
  )

sample_size_labels <- sample_protein_counts %>%
  group_by(treatment_full) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0("n=", n))

top_y_limit <- max(sample_protein_counts$protein_count, na.rm = TRUE) * 1.1

p_ef1b <- ggplot(sample_protein_counts, aes(x = treatment_full, y = protein_count, color = treatment)) +
  geom_jitter(width = 0.25, alpha = 0.4, size = 3) +
  geom_boxplot(width = 0.6, alpha = 0.3, lwd = 1.2, outlier.shape = NA) +
  geom_text(data = sample_size_labels,
            aes(x = treatment_full, y = top_y_limit, label = label),
            size = 10, fontface = "bold", color = "black") +
  scale_color_manual(values = c("A" = "#FF6347", "C" = "#2E8B57", "H" = "#1E90FF")) +
  scale_y_continuous(limits = c(0, top_y_limit * 1.05)) +
  labs(x = "Treatment group", y = "# Proteins identified") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 24, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, face = "bold"))

ggsave("../figure/ef1b.pdf", p_ef1b, width = 10, height = 8, device = "pdf")


# =============================================================================
# Shared t-SNE function (used by ef1c, ef1e)
# =============================================================================
create_tsne_plot <- function(expr_matrix, metadata_df, file_prefix) {
  completeness_rate <- rowSums(!is.na(expr_matrix)) / ncol(expr_matrix)
  expr_filtered <- expr_matrix[completeness_rate > 0.5, ]
  global_min <- min(expr_filtered, na.rm = TRUE)
  imputation_value <- if (global_min > 0) global_min * 0.8 else global_min * 1.1
  expr_filtered[is.na(expr_filtered)] <- imputation_value

  data_for_tsne <- t(expr_filtered)
  unique_rows_mask <- !duplicated(data_for_tsne)
  data_for_tsne_unique <- data_for_tsne[unique_rows_mask, ]
  metadata_subset_unique <- metadata_df %>% filter(SampleID %in% rownames(data_for_tsne_unique))

  adaptive_perplexity <- min(10, floor((nrow(data_for_tsne_unique) - 1) / 3))
  set.seed(42)
  tsne_result <- Rtsne(data_for_tsne_unique, perplexity = adaptive_perplexity)

  plot_data <- tibble(
    SampleID = rownames(data_for_tsne_unique),
    tSNE1 = tsne_result$Y[, 1],
    tSNE2 = tsne_result$Y[, 2]
  ) %>% left_join(metadata_subset_unique, by = "SampleID")

  p_inst <- ggplot(plot_data, aes(x = tSNE1, y = tSNE2, color = Instrument)) +
    geom_point(alpha = 0.95, size = 0.5) +
    stat_ellipse(linewidth = 1.0, level = 0.95) +
    labs(title = "t-SNE by Instrument", subtitle = "95% Confidence Ellipses") +
    theme_bw() + theme(legend.position = "bottom")

  p_prep <- ggplot(plot_data, aes(x = tSNE1, y = tSNE2, color = as.factor(PrepBatch))) +
    geom_point(alpha = 0.95, size = 0.5) +
    stat_ellipse(linewidth = 1.0, level = 0.95) +
    labs(title = "t-SNE by PrepBatch", subtitle = "95% Confidence Ellipses", color = "PrepBatch") +
    theme_bw() + theme(legend.position = "bottom")

  (p_inst | p_prep) +
    plot_annotation(title = paste("t-SNE Batch Effect:", file_prefix),
                    theme = theme(plot.title = element_text(face = "bold", size = 16)))
}


# =============================================================================
# Extended Fig.1c — t-SNE before batch correction
# =============================================================================

expr_before_tsne <- read_csv(NORMALIZED, show_col_types = FALSE) %>%
  column_to_rownames("Genes") %>% as.matrix()

p_ef1c <- create_tsne_plot(expr_before_tsne, metadata, "Before Correction")
ggsave("../figure/ef1c.pdf", p_ef1c, width = 12, height = 7, device = "pdf")


# =============================================================================
# Shared PVCA function (used by ef1d, ef1f)
# =============================================================================
run_pvca <- function(expr_matrix, metadata_df, biological_groups, file_prefix) {
  bio_metadata <- metadata_df %>% filter(Group %in% biological_groups)
  common_samples <- intersect(colnames(expr_matrix), bio_metadata$SampleID)
  if (length(common_samples) == 0) stop("No biological samples found for PVCA.")

  expr_bio_only <- expr_matrix[, common_samples]
  metadata_bio_only <- bio_metadata %>% filter(SampleID %in% common_samples)

  completeness_rate <- rowSums(!is.na(expr_bio_only)) / ncol(expr_bio_only)
  expr_filtered <- expr_bio_only[completeness_rate >= 0.5, ]
  gene_means <- rowMeans(expr_filtered, na.rm = TRUE)
  for (i in seq_len(nrow(expr_filtered))) {
    expr_filtered[i, is.na(expr_filtered[i, ])] <- gene_means[i]
  }

  pData_df <- data.frame(
    PrepBatch  = as.factor(metadata_bio_only$PrepBatch),
    Instrument = as.factor(metadata_bio_only$Instrument),
    Group      = as.factor(metadata_bio_only$Group),
    Time       = as.factor(metadata_bio_only$Time),
    row.names  = colnames(expr_filtered)
  )
  eset <- ExpressionSet(assayData = as.matrix(expr_filtered),
                         phenoData = AnnotatedDataFrame(pData_df))

  pvca_result <- pvcaBatchAssess(abatch = eset,
                                  batch.factors = c("PrepBatch", "Instrument", "Group", "Time"),
                                  threshold = 0.6)

  variance_data <- data.frame(Factor = pvca_result$label,
                               Variance = as.numeric(pvca_result$dat))
  variance_data$Category <- case_when(
    variance_data$Factor %in% c("PrepBatch", "Instrument") ~ "Technical",
    variance_data$Factor %in% c("Group", "Time")          ~ "Biological",
    grepl(":", variance_data$Factor)                       ~ "Interaction",
    variance_data$Factor == "resid"                        ~ "Residual",
    TRUE                                                   ~ "Other"
  )

  cat("\nVariance Proportion Table (", file_prefix, "):\n")
  print_table <- variance_data %>%
    mutate(Percentage = paste0(round(Variance * 100, 2), "%")) %>%
    arrange(desc(Variance)) %>%
    select(Factor, Category, Percentage)
  print(as.data.frame(print_table), row.names = FALSE)

  ggplot(variance_data, aes(x = reorder(Factor, Variance), y = Variance, fill = Category)) +
    geom_bar(stat = "identity", color = "white") +
    coord_flip() +
    scale_fill_manual(values = c("Technical" = "#E31A1C", "Biological" = "#2E8B57",
                                  "Interaction" = "#1E1E3F", "Residual" = "#808080")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = paste("PVCA:", file_prefix),
         subtitle = paste("Based on", nrow(expr_filtered), "features from C/H/A biological samples"),
         x = "Factors", y = "Weighted Proportion of Variance") +
    theme_bw() + theme(plot.title = element_text(face = "bold"))
}


# =============================================================================
# Extended Fig.1d — PVCA before batch correction
# =============================================================================

expr_before_pvca <- read_csv(NORMALIZED, show_col_types = FALSE) %>%
  column_to_rownames("Genes") %>% as.matrix()

p_ef1d <- run_pvca(expr_before_pvca, metadata, c("C", "H", "A"), "Before Correction")
ggsave("../figure/ef1d.pdf", p_ef1d, width = 9, height = 6, device = "pdf")


# =============================================================================
# Extended Fig.1e — t-SNE after batch correction
# =============================================================================

expr_after_tsne <- read_tsv(CORRECTED_NOQC, show_col_types = FALSE) %>%
  column_to_rownames("Genes") %>% as.matrix()
common_after <- intersect(colnames(expr_after_tsne), metadata$SampleID)
expr_after_tsne <- expr_after_tsne[, common_after]

p_ef1e <- create_tsne_plot(expr_after_tsne, metadata, "After Correction")
ggsave("../figure/ef1e.pdf", p_ef1e, width = 12, height = 7, device = "pdf")


# =============================================================================
# Extended Fig.1f — PVCA after batch correction
# =============================================================================

expr_after_pvca <- read_tsv(CORRECTED_FULL, show_col_types = FALSE) %>%
  column_to_rownames("Genes") %>% as.matrix()

p_ef1f <- run_pvca(expr_after_pvca, metadata, c("C", "H", "A"), "After Correction")
ggsave("../figure/ef1f.pdf", p_ef1f, width = 9, height = 6, device = "pdf")


# =============================================================================
# Shared PCoA and permutation setup (used by ef1g, fig1b, fig1c)
# =============================================================================

k_stages <- 5
n_perms <- 999

protein_pcoa <- read_tsv(IMPUTED_FINAL, show_col_types = FALSE)
sample_ids <- setdiff(colnames(protein_pcoa), "Genes")

pcoa_metadata <- tibble(SampleID = sample_ids) %>%
  mutate(
    treatment  = factor(str_extract(SampleID, "^[CHA]"), levels = c("C", "A", "H")),
    time_point = as.numeric(str_extract(SampleID, "\\d+$"))
  ) %>%
  mutate(
    stage_num = cut(time_point, breaks = k_stages, labels = FALSE, include.lowest = TRUE)
  ) %>%
  mutate(time_stage = factor(paste0("Stage ", stage_num)))

expr_pcoa <- protein_pcoa %>% column_to_rownames("Genes") %>% t()
expr_pcoa_scaled <- scale(expr_pcoa[pcoa_metadata$SampleID, ])
dist_matrix <- dist(expr_pcoa_scaled)
pcoa_full <- cmdscale(dist_matrix, k = 10, eig = TRUE)
pc_coords <- as.data.frame(pcoa_full$points)
colnames(pc_coords) <- paste0("PC", 1:10)
pcoa_data <- bind_cols(pcoa_metadata, pc_coords)

calc_centroid_significance <- function(df, group_col, pc_cols, perms = 999) {
  groups <- levels(df[[group_col]])
  pairs <- combn(groups, 2, simplify = FALSE)
  results <- list()
  for (p in pairs) {
    sub_df <- df[df[[group_col]] %in% p, ]
    c1 <- colMeans(sub_df[sub_df[[group_col]] == p[1], pc_cols])
    c2 <- colMeans(sub_df[sub_df[[group_col]] == p[2], pc_cols])
    obs_dist <- sqrt(sum((c1 - c2)^2))
    null_dists <- replicate(perms, {
      shuffled_labels <- sample(sub_df[[group_col]])
      nc1 <- colMeans(sub_df[shuffled_labels == p[1], pc_cols])
      nc2 <- colMeans(sub_df[shuffled_labels == p[2], pc_cols])
      sqrt(sum((nc1 - nc2)^2))
    })
    p_val <- (sum(null_dists >= obs_dist) + 1) / (perms + 1)
    results[[paste(sort(p), collapse = " vs ")]] <- tibble(obs_dist = obs_dist, p_val = p_val)
  }
  bind_rows(results, .id = "Comparison")
}

treat_colors <- c("C" = "#2E8B57", "A" = "#FF6347", "H" = "#1E90FF")
comp_colors <- c("A vs C" = "#e31a1c", "C vs H" = "#1f78b4", "A vs H" = "#6a3d9a")


# =============================================================================
# Extended Fig.1g — Stage-wise centroid distance permutation test
# =============================================================================

pc_cols <- paste0("PC", 1:10)

stage_stats <- pcoa_data %>%
  group_by(time_stage) %>%
  group_modify(~ calc_centroid_significance(.x, "treatment", pc_cols, perms = n_perms)) %>%
  ungroup()

stage_stats <- stage_stats %>%
  mutate(sig_label = case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE          ~ "ns"
  ))

p_ef1g <- ggplot(stage_stats, aes(x = time_stage, y = obs_dist, fill = Comparison)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
  geom_text(aes(label = sig_label), position = position_dodge(0.8),
            vjust = -0.5, fontface = "bold") +
  scale_fill_manual(values = comp_colors) +
  labs(title = "Centroid Distances by Stage (PC1-10 Space)", y = "Distance", x = "") +
  theme_bw(base_size = 14) + theme(legend.position = "top")

ggsave("../figure/ef1g.pdf", p_ef1g, width = 8, height = 6, device = "pdf")
write_csv(stage_stats, "../output/ef1g_SourceData.csv")


# =============================================================================
# Fig.1b — Temporal evolution of centroid distances (LOESS)
# =============================================================================

time_points <- sort(unique(pcoa_metadata$time_point))
dist_dynamic <- tibble()
for (tp in time_points) {
  tp_data <- pcoa_data %>% filter(time_point == tp)
  if (length(unique(tp_data$treatment)) > 1) {
    res <- calc_centroid_significance(tp_data, "treatment", pc_cols, perms = 1)
    dist_dynamic <- bind_rows(dist_dynamic, res %>% mutate(time_point = tp))
  }
}

p_fig1b <- ggplot(dist_dynamic, aes(x = time_point, y = obs_dist, color = Comparison, fill = Comparison)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "loess", span = 0.3, alpha = 0.2, size = 1.2) +
  scale_color_manual(values = comp_colors) +
  scale_fill_manual(values = comp_colors) +
  labs(title = "Temporal Evolution of Centroid Distances (PC1-10)",
       x = "Time Point", y = "Euclidean Distance (PC1-10)") +
  theme_bw(base_size = 14) + theme(legend.position = "bottom")

ggsave("../figure/fig1b.pdf", p_fig1b, width = 10, height = 6, device = "pdf")


# =============================================================================
# Fig.1c — PCoA trajectories across 5 stages
# =============================================================================

var_expl <- pcoa_full$eig / sum(pcoa_full$eig)

p_fig1c <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = treatment, shape = treatment)) +
  geom_point(alpha = 0.6, size = 2) +
  stat_ellipse(size = 0.8) +
  facet_wrap(~ time_stage, nrow = 1) +
  scale_color_manual(values = treat_colors) +
  labs(title = "PCoA Trajectories (5 Stages)",
       x = paste0("PC1 (", round(var_expl[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(var_expl[2] * 100, 1), "%)"),
       color = "Treatment") +
  theme_bw(base_size = 14) + theme(legend.position = "bottom")

ggsave("../figure/fig1c.pdf", p_fig1c, width = 15, height = 4.5, device = "pdf")

cat("Step01 — Done.\n")
