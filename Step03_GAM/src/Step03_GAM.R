# =============================================================================
# Step03_GAM.R — Fig.3 and EFig.4: GAM differential analysis and GSEA
#
# NOTE: The original study did not save the random seed for gseKEGG, so exact
# GSEA reproduction is not possible. Two paths are provided:
#
#   Algorithm-level: This script runs the full GAM + GSEA pipeline with a
#   fixed seed (20250707L). Figures may differ slightly from the publication.
#
#   Source-data: original/plot_from_paper.R reproduces the published figures
#   directly from the original GSEA intermediates saved in original/data/.
#
# GAM fitting (mgcv::gam, REML) is deterministic and unaffected by the seed;
# fig3c, ef4a, ef4d, ef4e are identical to the published version.
#
# Prerequisites: Step02_Clustering completed
# =============================================================================

rm(list = ls())
setwd(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))))

suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(progress)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(scales)
})

# ---- Paths ----
F_SMOOTHED  <- "../../Step02_Clustering/output/smoothed_average.tsv"
ID_MAPPING  <- "../../gene_id_mapping.tsv"

# Output file names (semantic, no K_advanced prefixes)
GAM_META       <- "../output/GAM_metadata.tsv"
GAM_INTEGRATED <- "../output/GAM_integrated_data.tsv"
GAM_RAW        <- "../output/GAM_raw_results.tsv"
GAM_CORRECTED  <- "../output/GAM_results_corrected.tsv"
GAM_EFFECT     <- "../output/GAM_effect_sizes.tsv"
GAM_RANKING    <- "../output/GAM_maxdiff_ranking.csv"
GAM_TEMP_DIR   <- "../output/GAM_temporal_GSEA"
GAM_STATIC_DIR <- "../output/GAM_static_GSEA"

cat("Step03 — GAM...\n")

dir.create("../output", showWarnings = FALSE, recursive = TRUE)
dir.create("../figure", showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# 1. PER-PROTEIN GAM FITTING
# =============================================================================
if (all(file.exists(c(GAM_META, GAM_INTEGRATED, GAM_RAW, GAM_CORRECTED)))) {
} else {
  require_file <- function(p) if (!file.exists(p)) stop("Missing: ", p)
  require_file(F_SMOOTHED)

  sample_ids <- colnames(read_tsv(F_SMOOTHED, n_max = 0, show_col_types = FALSE))[-1]
  metadata <- tibble(SampleID = sample_ids) %>%
    mutate(Treatment = str_extract(SampleID, "^[CHA]")) %>%
    group_by(Treatment) %>%
    mutate(Position = row_number() - 1) %>% ungroup() %>%
    mutate(Time_Numeric = Position, Time = as.character(Position)) %>%
    dplyr::select(SampleID, Treatment, Position, Time_Numeric, Time)
  write_tsv(metadata, GAM_META)

  protein_matrix <- read_tsv(F_SMOOTHED, show_col_types = FALSE)
  names(protein_matrix)[1] <- "ProteinID"

  integrated_data <- protein_matrix %>%
    pivot_longer(-ProteinID, names_to = "SampleID", values_to = "Abundance") %>%
    left_join(metadata, by = "SampleID") %>%
    mutate(
      Treatment = factor(Treatment, levels = c("C", "H", "A")),
      Time = as.numeric(Time_Numeric)
    ) %>%
    dplyr::select(ProteinID, SampleID, Time, Treatment, Abundance) %>%
    arrange(ProteinID, Treatment, Time) %>%
    filter(!is.na(Abundance))
  write_tsv(integrated_data, GAM_INTEGRATED)
  analyze_protein_gam <- function(protein_id, data, k_value = 7) {
    protein_data <- data %>% filter(ProteinID == protein_id) %>% arrange(Treatment, Time)
    if (nrow(protein_data) < 20) return(tibble(ProteinID = protein_id, Status = "INSUFFICIENT_DATA"))
    tc <- protein_data %>% count(Treatment) %>% pull(n)
    if (length(tc) < 3 || any(tc < 6)) return(tibble(ProteinID = protein_id, Status = "UNBALANCED"))
    if (var(protein_data$Abundance, na.rm = TRUE) < 1e-8) return(tibble(ProteinID = protein_id, Status = "NO_VARIANCE"))
    tryCatch({
      model <- gam(Abundance ~ Treatment + s(Time, k = k_value) +
                     s(Time, by = Treatment, k = k_value),
                   data = protein_data, method = "REML")
      smry <- summary(model)
      result <- tibble(ProteinID = protein_id, Status = "SUCCESS")
      st <- smry$s.table
      for (i in seq_len(nrow(st))) {
        tn <- rownames(st)[i]
        suffix <- case_when(
          grepl("^s\\(Time\\)$", tn) ~ "_Control_Trend",
          grepl("s\\(Time\\):TreatmentH", tn) ~ "_Diff_H_vs_C",
          grepl("s\\(Time\\):TreatmentA", tn) ~ "_Diff_A_vs_C",
          TRUE ~ NA_character_)
        if (is.na(suffix)) next
        result[[paste0("EDF", suffix)]]    <- st[i, "edf"]
        result[[paste0("F_stat", suffix)]] <- st[i, "F"]
        result[[paste0("P_value", suffix)]] <- st[i, "p-value"]
      }
      pt <- smry$p.table
      if (!is.null(pt) && nrow(pt) > 0) {
        if ("TreatmentH" %in% rownames(pt)) {
          result$P_value_TreatmentH_Param  <- pt["TreatmentH", "Pr(>|t|)"]
          result$Estimate_TreatmentH_Param <- pt["TreatmentH", "Estimate"]
        }
        if ("TreatmentA" %in% rownames(pt)) {
          result$P_value_TreatmentA_Param  <- pt["TreatmentA", "Pr(>|t|)"]
          result$Estimate_TreatmentA_Param <- pt["TreatmentA", "Estimate"]
        }
      }
      result$AIC <- AIC(model); result$R_squared <- smry$r.sq
      result$Deviance_Explained <- smry$dev.expl
      result
    }, error = function(e) tibble(ProteinID = protein_id, Status = paste0("ERROR: ", e$message)))
  }

  proteins <- unique(integrated_data$ProteinID)
  cat(sprintf("  Fitting GAM for %d proteins...\n", length(proteins)))
  pb <- progress_bar$new(total = length(proteins), format = "  [:bar] :percent eta: :eta", clear = FALSE, width = 60)
  gam_results <- map_dfr(proteins, function(pid) { pb$tick(); analyze_protein_gam(pid, integrated_data) })
  write_tsv(gam_results, GAM_RAW)

  # BH correction
  gam_corrected <- gam_results %>% filter(Status == "SUCCESS")
  p_cols <- c("P_value_Control_Trend", "P_value_Diff_H_vs_C", "P_value_Diff_A_vs_C")
  adj_cols <- c("Adj_P_Control_Trend", "Adj_P_Diff_H_vs_C", "Adj_P_Diff_A_vs_C")
  sig_cols <- c("Significant_Control_Trend", "Significant_Diff_H_vs_C", "Significant_Diff_A_vs_C")
  for (i in seq_along(p_cols)) {
    if (p_cols[i] %in% names(gam_corrected)) {
      gam_corrected[[adj_cols[i]]] <- p.adjust(gam_corrected[[p_cols[i]]], method = "BH")
      gam_corrected[[sig_cols[i]]] <- gam_corrected[[adj_cols[i]]] < 0.05
    }
  }
  param_p <- c("P_value_TreatmentH_Param", "P_value_TreatmentA_Param")
  param_a <- c("Adj_P_TreatmentH_Param", "Adj_P_TreatmentA_Param")
  for (i in seq_along(param_p)) {
    if (param_p[i] %in% names(gam_corrected))
      gam_corrected[[param_a[i]]] <- p.adjust(gam_corrected[[param_p[i]]], method = "BH")
  }
  write_tsv(gam_corrected, GAM_CORRECTED)
  cat(sprintf("  GAM fitted: %d proteins\n", nrow(gam_corrected)))
}


# =============================================================================
# 2. EFFECT-SIZE COMPUTATION
# =============================================================================
effect_cols <- c("ProteinID", "Adj_P_Diff_H_vs_C", "Adj_P_Diff_A_vs_C",
                 "MaxDiff_H_vs_C", "MaxDiff_A_vs_C", "MaxAbsDiff_H_vs_C",
                 "MaxAbsDiff_A_vs_C", "RMSD_H_vs_C", "RMSD_A_vs_C",
                 "Abs_AUC_Diff_H_vs_C", "Abs_AUC_Diff_A_vs_C")

if (file.exists(GAM_EFFECT) && all(effect_cols %in% names(read_tsv(GAM_EFFECT, n_max = 0, show_col_types = FALSE)))) {
} else {
  integrated_data <- read_tsv(GAM_INTEGRATED, show_col_types = FALSE) %>%
    mutate(Treatment = factor(Treatment, levels = c("C", "H", "A")))
  gam_corrected <- read_tsv(GAM_CORRECTED, show_col_types = FALSE)

  sig_proteins <- gam_corrected %>%
    filter(Significant_Diff_H_vs_C | Significant_Diff_A_vs_C) %>% pull(ProteinID)

  compute_effect_sizes <- function(protein_id, data, k_value = 7) {
    protein_data <- data %>% filter(ProteinID == protein_id)
    if (nrow(protein_data) < 20 || n_distinct(protein_data$Treatment) < 3)
      return(tibble(ProteinID = protein_id))
    tryCatch({
      model <- gam(Abundance ~ Treatment + s(Time, k = k_value) +
                     s(Time, by = Treatment, k = k_value),
                   data = protein_data, method = "REML")
      time_grid <- seq(min(data$Time, na.rm = TRUE), max(data$Time, na.rm = TRUE), length.out = 101)
      newdata <- expand.grid(Time = time_grid, Treatment = levels(data$Treatment))
      newdata$Treatment <- factor(newdata$Treatment, levels = levels(data$Treatment))
      pred <- predict(model, newdata = newdata, se.fit = TRUE)
      newdata$Fit <- pred$fit
      pred_matrix <- matrix(newdata$Fit, nrow = length(time_grid),
                            dimnames = list(NULL, levels(data$Treatment)))
      diff_H <- pred_matrix[, "H"] - pred_matrix[, "C"]
      diff_A <- pred_matrix[, "A"] - pred_matrix[, "C"]
      max_idx_H <- which.max(abs(diff_H)); max_idx_A <- which.max(abs(diff_A))
      tibble(
        ProteinID = protein_id,
        Adj_P_Diff_H_vs_C = gam_corrected$Adj_P_Diff_H_vs_C[gam_corrected$ProteinID == protein_id][1],
        Adj_P_Diff_A_vs_C = gam_corrected$Adj_P_Diff_A_vs_C[gam_corrected$ProteinID == protein_id][1],
        MaxDiff_H_vs_C = diff_H[max_idx_H], MaxDiff_A_vs_C = diff_A[max_idx_A],
        MaxAbsDiff_H_vs_C = max(abs(diff_H)), MaxAbsDiff_A_vs_C = max(abs(diff_A)),
        RMSD_H_vs_C = sqrt(mean(diff_H^2)), RMSD_A_vs_C = sqrt(mean(diff_A^2)),
        Abs_AUC_Diff_H_vs_C = sum(abs(diff_H)) / length(time_grid),
        Abs_AUC_Diff_A_vs_C = sum(abs(diff_A)) / length(time_grid)
      )
    }, error = function(e) tibble(ProteinID = protein_id))
  }

  cat(sprintf("  Computing effect sizes for %d significant proteins...\n", length(sig_proteins)))
  pb <- progress_bar$new(total = length(sig_proteins), format = "  [:bar] :percent eta: :eta", clear = FALSE, width = 60)
  effect_results <- map_dfr(sig_proteins, function(pid) { pb$tick(); compute_effect_sizes(pid, integrated_data) })
  write_tsv(effect_results, GAM_EFFECT)

  # Full-protein MaxDiff ranking
  all_proteins <- unique(integrated_data$ProteinID)
  pb2 <- progress_bar$new(total = length(all_proteins), format = "  [:bar] :percent eta: :eta", clear = FALSE, width = 60)
  ranking <- map_dfr(all_proteins, function(pid) {
    pb2$tick()
    pd <- integrated_data %>% filter(ProteinID == pid)
    if (nrow(pd) < 20 || n_distinct(pd$Treatment) < 3) return(tibble(ProteinID = pid))
    tryCatch({
      model <- gam(Abundance ~ Treatment + s(Time, k = 7) + s(Time, by = Treatment, k = 7),
                   data = pd, method = "REML")
      tg <- seq(min(pd$Time, na.rm = TRUE), max(pd$Time, na.rm = TRUE), length.out = 101)
      nd <- expand.grid(Time = tg, Treatment = factor(levels(pd$Treatment), levels = levels(pd$Treatment)))
      pred <- predict(model, newdata = nd)
      pm <- matrix(pred, nrow = 101)
      colnames(pm) <- levels(pd$Treatment)
      dH <- pm[, "H"] - pm[, "C"]; dA <- pm[, "A"] - pm[, "C"]
      tibble(ProteinID = pid,
             MaxDiff_H_vs_C = dH[which.max(abs(dH))],
             MaxDiff_A_vs_C = dA[which.max(abs(dA))])
    }, error = function(e) tibble(ProteinID = pid))
  })
  write_csv(ranking, GAM_RANKING)
}


# =============================================================================
# 3. TEMPORAL KEGG GSEA
# =============================================================================
temporal_outputs <- c(
  file.path(GAM_TEMP_DIR, "GSEA_Filtered_NES_Acetic acid vs Control.csv"),
  file.path(GAM_TEMP_DIR, "GSEA_Filtered_NES_Hhdrogen peroxide vs Control.csv")
)

if (all(file.exists(temporal_outputs))) {
} else {
  library(clusterProfiler)
  library(org.Sc.sgd.db)
  if (!dir.exists(GAM_TEMP_DIR)) dir.create(GAM_TEMP_DIR, recursive = TRUE)

  # ---- Build ID mapper (ENTREZID, matching original) ----
  id_mapper <- read_tsv(ID_MAPPING, show_col_types = FALSE) %>%
    dplyr::select(Original_ID, GENENAME, ENTREZID) %>%
    filter(!is.na(Original_ID)) %>%
    distinct(Original_ID, .keep_all = TRUE)

  integrated_data <- read_tsv(GAM_INTEGRATED, show_col_types = FALSE) %>%
    mutate(Treatment = factor(Treatment, levels = c("C", "H", "A")))

  # ---- Per-protein GAM on 96-point grid, bin 6 → 16 ----
  time_grid <- seq(min(integrated_data$Time), max(integrated_data$Time), length.out = 16 * 6)

  calc_binned_gam <- function(pid, data, grid) {
    p_data <- data %>% filter(ProteinID == pid)
    if (nrow(p_data) < 20) return(NULL)
    model <- tryCatch(
      gam(Abundance ~ Treatment + s(Time, k = 7) + s(Time, by = Treatment, k = 7),
          data = p_data, method = "REML"),
      error = function(e) NULL)
    if (is.null(model)) return(NULL)
    new_data <- expand_grid(Time = grid, Treatment = factor(c("C", "A", "H"), levels = c("C", "H", "A")))
    new_data$Pred <- predict(model, newdata = new_data)
    pred_c <- new_data$Pred[new_data$Treatment == "C"]
    pred_a <- new_data$Pred[new_data$Treatment == "A"]
    pred_h <- new_data$Pred[new_data$Treatment == "H"]
    bin_sum <- function(x) colSums(matrix(x, nrow = 6, ncol = 16))
    list(AC = bin_sum(pred_a - pred_c), HC = bin_sum(pred_h - pred_c))
  }

  proteins <- unique(integrated_data$ProteinID)
  pb <- progress_bar$new(total = length(proteins), format = "  [:bar] :percent eta: :eta", clear = FALSE, width = 60)
  results <- list()
  for (pid in proteins) { pb$tick(); out <- calc_binned_gam(pid, integrated_data, time_grid); if (!is.null(out)) results[[pid]] <- out }

  list_to_matrix <- function(result_list, component) {
    mat <- do.call(rbind, lapply(result_list, function(x) x[[component]]))
    rownames(mat) <- names(result_list)
    colnames(mat) <- sprintf("%.0f", seq(3, 93, by = 6))
    mat
  }
  mat_ac <- list_to_matrix(results, "AC")
  mat_hc <- list_to_matrix(results, "HC")
  write.csv(mat_ac, file.path(GAM_TEMP_DIR, "Matrix_Binned_Diff_Acetic_vs_Control.csv"))
  write.csv(mat_hc, file.path(GAM_TEMP_DIR, "Matrix_Binned_Diff_H2O2_vs_Control.csv"))

  # ---- Gene name translation helper ----
  translate_gene_string <- function(id_string, lookup_vec) {
    ids <- str_split(id_string, "/")[[1]]
    names <- lookup_vec[ids]; names[is.na(names)] <- ids[is.na(names)]
    paste(names, collapse = "/")
  }

  # ---- Run GSEA per time point ----
  run_temporal_gsea <- function(rank_matrix, comparison_name) {
    rank_df <- as.data.frame(rank_matrix) %>% rownames_to_column("Original_ID")

    ranked <- rank_df %>%
      left_join(id_mapper, by = "Original_ID") %>%
      filter(!is.na(ENTREZID), ENTREZID != "") %>%
      distinct(ENTREZID, .keep_all = TRUE)

    all_results <- list()
    for (tp in colnames(rank_matrix)) {
      gene_list <- ranked[[tp]]; names(gene_list) <- ranked$ENTREZID
      gene_list <- sort(gene_list, decreasing = TRUE)
      set.seed(20250707L)
      gse <- tryCatch(
        gseKEGG(geneList = gene_list, organism = "sce", keyType = "ncbi-geneid",
                pvalueCutoff = 1, minGSSize = 10, maxGSSize = 500,
                verbose = FALSE, nPermSimple = 1000),
        error = function(e) NULL)
      if (!is.null(gse) && nrow(gse) > 0) {
        all_results[[tp]] <- as.data.frame(gse) %>%
          dplyr::select(ID, Description, NES, pvalue, core_enrichment) %>%
          mutate(TimePoint = tp)
      }
    }
    if (length(all_results) == 0) return(invisible(NULL))

    combined <- bind_rows(all_results)
    # Build NES matrix — rownames = Description (matching original format)
    hm_mat <- combined %>%
      dplyr::select(Description, TimePoint, NES) %>%
      pivot_wider(names_from = TimePoint, values_from = NES, values_fill = 0) %>%
      column_to_rownames("Description")
    pval_mat <- combined %>%
      dplyr::select(Description, TimePoint, pvalue) %>%
      pivot_wider(names_from = TimePoint, values_from = pvalue, values_fill = 1) %>%
      column_to_rownames("Description")
    valid_cols <- intersect(colnames(rank_matrix), colnames(hm_mat))
    hm_mat <- hm_mat[, valid_cols, drop = FALSE]
    pval_mat <- pval_mat[, valid_cols, drop = FALSE]

    # Blacklist + filter
    blacklist <- c("Viral","virus","HIV","Infection","Bacterial","cancer","disease",
                   "Systemic","Lupus","Huntington","Alzheimer","Parkinson","Methane",
                   "IgSF","Malaria","Metabolic pathways",
                   "Biosynthesis of secondary metabolites","Autophagy - other")
    clean_rows <- !grepl(paste(blacklist, collapse = "|"), rownames(hm_mat), ignore.case = TRUE)
    max_abs_nes <- apply(hm_mat, 1, function(x) max(abs(x), na.rm = TRUE))
    sig_counts <- rowSums(pval_mat < 0.05, na.rm = TRUE)
    keep <- clean_rows & max_abs_nes > 1.5 & sig_counts >= 2
    hm_filt <- hm_mat[keep, , drop = FALSE]
    pval_filt <- pval_mat[rownames(hm_filt), colnames(hm_filt), drop = FALSE]
    if (nrow(hm_filt) < 1) return(invisible(NULL))

    write.csv(hm_filt, file.path(GAM_TEMP_DIR, paste0("GSEA_Filtered_NES_", comparison_name, ".csv")))

    # Detailed CSV with translated gene names
    entrez_to_gene <- setNames(id_mapper$GENENAME, id_mapper$ENTREZID)
    detailed <- combined %>%
      filter(Description %in% rownames(hm_filt)) %>%
      rowwise() %>%
      mutate(Genes = translate_gene_string(core_enrichment, entrez_to_gene)) %>%
      ungroup() %>%
      dplyr::select(TimePoint, Description, NES, pvalue, Genes) %>%
      arrange(Description, TimePoint)
    write_csv(detailed, file.path(GAM_TEMP_DIR, paste0("GSEA_Detailed_", comparison_name, ".csv")))

    # Full heatmap (ef4f/ef4g — filtered matrix + significance stars, matching reference)
    sig_full <- matrix("", nrow = nrow(pval_filt), ncol = ncol(pval_filt), dimnames = dimnames(pval_filt))
    sig_full[pval_filt < 0.05] <- "*"; sig_full[pval_filt < 0.01] <- "**"
    plot_h <- max(6, nrow(hm_filt) * 0.25)
    pheatmap(hm_filt, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
             color = colorRampPalette(c("navy","white","firebrick3"))(100),
             border_color = "grey95",
             display_numbers = sig_full, number_color = "black",
             fontsize_row = 10, fontsize_col = 10,
             main = comparison_name,
             filename = file.path(GAM_TEMP_DIR, paste0("Heatmap_", comparison_name, ".pdf")),
             width = 12, height = plot_h)
  }

  run_temporal_gsea(mat_ac, "Acetic acid vs Control")
  run_temporal_gsea(mat_hc, "Hhdrogen peroxide vs Control")
}


# =============================================================================
# 4. STATIC KEGG GSEA (matching original)
# =============================================================================
static_outputs <- c(
  file.path(GAM_STATIC_DIR, "Acetic_vs_Control_KEGG_data.csv"),
  file.path(GAM_STATIC_DIR, "H2O2_vs_Control_KEGG_data.csv")
)

if (all(file.exists(static_outputs))) {
} else {
  library(clusterProfiler)
  library(org.Sc.sgd.db)
  if (!dir.exists(GAM_STATIC_DIR)) dir.create(GAM_STATIC_DIR, recursive = TRUE)

  id_mapper <- read_tsv(ID_MAPPING, show_col_types = FALSE) %>%
    dplyr::select(Original_ID, ORF, GENENAME, ENTREZID) %>%
    filter(!is.na(Original_ID)) %>% distinct(Original_ID, .keep_all = TRUE)

  ranked_base <- read_csv(GAM_RANKING, show_col_types = FALSE) %>%
    left_join(id_mapper, by = c("ProteinID" = "Original_ID"))

  prepare_gene_list <- function(df, metric_col, id_col) {
    out <- df %>% filter(!is.na(.data[[metric_col]]), !is.na(.data[[id_col]]), .data[[id_col]] != "") %>%
      distinct(.data[[id_col]], .keep_all = TRUE) %>% arrange(desc(.data[[metric_col]]))
    setNames(out[[metric_col]], out[[id_col]])
  }

  translate_core <- function(gse_result) {
    e2g <- AnnotationDbi::mapIds(org.Sc.sgd.db, keys(org.Sc.sgd.db, "ENTREZID"),
                                  column = "GENENAME", keytype = "ENTREZID")
    result_df <- as.data.frame(gse_result)
    result_df$core_enrichment <- sapply(result_df$core_enrichment, function(cs) {
      ids <- strsplit(cs, "/")[[1]]
      gn <- e2g[ids]; gn[is.na(gn)] <- ids[is.na(gn)]
      paste(gn, collapse = "/")
    })
    gse_result@result$core_enrichment <- result_df$core_enrichment
    gse_result
  }

  for (comp in c("H_vs_C", "A_vs_C")) {
    metric <- paste0("MaxDiff_", comp)
    comp_name <- if (comp == "H_vs_C") "H2O2_vs_Control" else "Acetic_vs_Control"
    gene_list <- prepare_gene_list(ranked_base, metric, "ENTREZID")
    set.seed(20250707L)
    gse <- tryCatch(
      gseKEGG(geneList = gene_list, organism = "sce", keyType = "ncbi-geneid",
              pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE,
              minGSSize = 10, maxGSSize = 500),
      error = function(e) { cat(sprintf("  Error: %s\n", e$message)); NULL })
    if (!is.null(gse) && nrow(gse) > 0) {
      gse <- translate_core(gse)
      write_csv(as.data.frame(gse), file.path(GAM_STATIC_DIR, paste0(comp_name, "_KEGG_data.csv")))
    }
  }
}


# =============================================================================
# 5. FIGURE PANELS
# =============================================================================
# ---- 5a. Fig.3A/B — Selected protein trace index ----
trace_map <- tribble(
  ~Panel, ~ProteinID, ~Source_pdf,
  "fig3a1", "RPL3",  "All_Protein_Plots_PDF_Protein_name/RPL3.pdf",
  "fig3a2", "KRE6",  "All_Protein_Plots_PDF_Protein_name/KRE6.pdf",
  "fig3a3", "LEU1",  "All_Protein_Plots_PDF_Protein_name/LEU1.pdf",
  "fig3a4", "FAS2",  "All_Protein_Plots_PDF_Protein_name/FAS2.pdf",
  "fig3b1", "HSP82", "All_Protein_Plots_PDF_Protein_name/HSP82.pdf",
  "fig3b2", "SSA4",  "All_Protein_Plots_PDF_Protein_name/SSA4.pdf"
)
write_tsv(trace_map, "../figure/fig3_trace_panel_index.tsv")


# ---- 5b. Fig.3C — Dynamic volcano ----

effect_file <- GAM_EFFECT
final_results <- read_tsv(effect_file, show_col_types = FALSE)

create_volcano <- function(results_df, title, p_col, effect_col, stability_col, y_cap = 8) {
  p_threshold <- 0.05; effect_threshold <- 1
  stability_threshold <- median(results_df[[stability_col]], na.rm = TRUE)
  plot_data <- results_df %>%
    filter(!is.na(.data[[effect_col]]), !is.na(.data[[p_col]])) %>%
    mutate(
      log_p = pmin(-log10(.data[[p_col]]), y_cap),
      Category = case_when(
        .data[[p_col]] >= p_threshold ~ "Insignificant",
        .data[[p_col]] < p_threshold & abs(.data[[effect_col]]) > effect_threshold &
          .data[[stability_col]] > stability_threshold & .data[[effect_col]] > 0 ~ "Stable Up",
        .data[[p_col]] < p_threshold & abs(.data[[effect_col]]) > effect_threshold &
          .data[[stability_col]] > stability_threshold & .data[[effect_col]] < 0 ~ "Stable Down",
        .data[[p_col]] < p_threshold ~ "Significant, other",
        TRUE ~ "Insignificant"),
      Category = factor(Category, c("Stable Up","Stable Down","Significant, other","Insignificant"))) %>%
    arrange(desc(Category))
  label_data <- plot_data %>% filter(Category %in% c("Stable Up","Stable Down")) %>%
    group_by(Category) %>% slice_max(abs(.data[[effect_col]]), n = 20) %>% ungroup()
  counts <- plot_data %>% count(Category)
  up_n   <- counts$n[match("Stable Up", counts$Category)]; up_n <- ifelse(is.na(up_n), 0, up_n)
  down_n <- counts$n[match("Stable Down", counts$Category)]; down_n <- ifelse(is.na(down_n), 0, down_n)
  volcano_colors <- c("Stable Up" = "#D73027", "Stable Down" = "#4575B4",
                      "Significant, other" = "darkgrey", "Insignificant" = "lightgrey")
  ggplot(plot_data, aes(x = .data[[effect_col]], y = log_p)) +
    geom_point(aes(color = Category), alpha = 0.7, size = 1.8) +
    geom_text_repel(data = label_data, aes(label = ProteinID), size = 3,
                    box.padding = 0.3, max.overlaps = Inf, show.legend = FALSE) +
    geom_vline(xintercept = c(-effect_threshold, effect_threshold), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = volcano_colors, drop = FALSE) +
    labs(title = title, subtitle = sprintf("Up: %d | Down: %d", up_n, down_n),
         x = "Signed MaxDiff", y = expression(-log[10] ~ "(Adjusted P-value)")) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5), legend.position = "bottom")
}

p_h <- create_volcano(final_results, "H2O2 vs Control", "Adj_P_Diff_H_vs_C", "MaxDiff_H_vs_C", "RMSD_H_vs_C")
p_a <- create_volcano(final_results, "Acetic acid vs Control", "Adj_P_Diff_A_vs_C", "MaxDiff_A_vs_C", "RMSD_A_vs_C")
p_combined <- p_h + p_a + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("../figure/fig3c.pdf", p_combined, width = 16, height = 8, device = cairo_pdf)


# ---- 5c. Fig.3D/E — Refined temporal GSEA heatmaps ----

plot_refined_heatmap <- function(comp_name, target_keywords, output_file) {
  nes_file <- file.path(GAM_TEMP_DIR, paste0("GSEA_Filtered_NES_", comp_name, ".csv"))
  detail_file <- file.path(GAM_TEMP_DIR, paste0("GSEA_Detailed_", comp_name, ".csv"))
  if (!file.exists(nes_file) || !file.exists(detail_file)) {
    cat(sprintf("  [WARN] Missing GSEA data for %s\n", comp_name)); return()
  }
  full_nes <- read.csv(nes_file, row.names = 1, check.names = FALSE)
  detail_df <- read_csv(detail_file, show_col_types = FALSE)
  keep_rows <- unique(unlist(lapply(target_keywords, function(kw) {
    grep(kw, rownames(full_nes), ignore.case = TRUE)
  })))
  if (length(keep_rows) == 0) { cat(sprintf("  No target pathways for %s\n", comp_name)); return() }
  refined_nes <- full_nes[keep_rows, , drop = FALSE]
  pval_wide <- detail_df %>% dplyr::select(Description, TimePoint, pvalue) %>%
    mutate(TimePoint = as.character(TimePoint)) %>%
    pivot_wider(names_from = TimePoint, values_from = pvalue) %>% column_to_rownames("Description")
  colnames(pval_wide) <- as.character(colnames(pval_wide))
  refined_pval <- matrix(1, nrow = nrow(refined_nes), ncol = ncol(refined_nes), dimnames = dimnames(refined_nes))
  cr <- intersect(rownames(pval_wide), rownames(refined_pval))
  cc <- intersect(colnames(pval_wide), colnames(refined_pval))
  refined_pval[cr, cc] <- as.matrix(pval_wide[cr, cc]); refined_pval[is.na(refined_pval)] <- 1
  sig_mat <- matrix("", nrow = nrow(refined_pval), ncol = ncol(refined_pval), dimnames = dimnames(refined_pval))
  sig_mat[refined_pval < 0.05] <- "*"; sig_mat[refined_pval < 0.01] <- "**"
  pheatmap(refined_nes, scale = "none", cluster_rows = TRUE, cluster_cols = FALSE,
           color = colorRampPalette(c("navy","white","firebrick3"))(100),
           border_color = "grey90", display_numbers = sig_mat, number_color = "black",
           fontsize_number = 12, fontsize_row = 11, fontsize_col = 10,
           main = comp_name, filename = output_file,
           width = 8, height = 3 + nrow(refined_nes) * 0.42)
}

plot_refined_heatmap("Acetic acid vs Control",
  c("Oxidative phosphorylation","SNARE interactions","Propanoate","Glyoxylate",
    "Autophagy","Cell cycle - yeast","Ribosome biogenesis in eukaryotes","Tyrosine"),
  "../figure/fig3d.pdf")
plot_refined_heatmap("Hhdrogen peroxide vs Control",
  c("DNA replication","Base excision repair","Mitophagy","Citrate cycle","SNARE","Ribosome biogenesis"),
  "../figure/fig3e.pdf")


# ---- 5d. EFig.4A — GAM significance overview ----

gam_corrected <- read_tsv(GAM_CORRECTED, show_col_types = FALSE)
sig_counts <- tribble(
  ~Comparison, ~Count,
  "Control\nbaseline trend", sum(gam_corrected$Significant_Control_Trend, na.rm = TRUE),
  "H2O2\nvs Control", sum(gam_corrected$Significant_Diff_H_vs_C, na.rm = TRUE),
  "Acetic acid\nvs Control", sum(gam_corrected$Significant_Diff_A_vs_C, na.rm = TRUE))

p_ef4a <- ggplot(sig_counts, aes(x = Comparison, y = Count, fill = Comparison)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("#2E8B57","#1E90FF","#FF6347")) +
  labs(x = NULL, y = "# Significant proteins") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", axis.text.x = element_text(color = "black", size = 13),
        panel.grid.minor = element_blank())
ggsave("../figure/ef4a.pdf", p_ef4a, width = 8, height = 5, device = cairo_pdf)


# ---- 5e. EFig.4B/C — Static KEGG dotplots ----

plot_static_kegg <- function(source_csv, output_file, result_name, top_n = 15, plot_width = 7, wrap_width = 15) {
  if (!file.exists(source_csv)) { cat(sprintf("  [WARN] %s not found\n", source_csv)); return() }
  gsea_df <- read_csv(source_csv, show_col_types = FALSE)
  plot_data <- gsea_df %>%
    mutate(Regulation = if_else(NES > 0, "Activated", "Suppressed")) %>%
    group_by(Regulation) %>%
    slice_min(p.adjust, n = top_n, with_ties = FALSE) %>% ungroup() %>%
    mutate(Description = str_wrap(Description, width = wrap_width)) %>%
    arrange(Regulation, NES) %>%
    mutate(Description = fct_inorder(Description))
  p <- ggplot(plot_data, aes(x = NES, y = Description)) +
    geom_point(aes(size = setSize, color = p.adjust)) +
    facet_grid(Regulation ~ ., scales = "free_y", space = "free_y") +
    scale_size_continuous(name = "Gene Count", range = c(3, 8)) +
    scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
    labs(title = paste("GSEA -", result_name, "KEGG"), x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          axis.text.y = element_text(size = 12, lineheight = 0.9),
          strip.text.y = element_text(face = "bold", size = 14),
          panel.spacing = unit(1, "lines"), legend.position = "right")
  ggsave(output_file, p, device = cairo_pdf,
         width = plot_width, height = max(6, 2 + 0.5 * nrow(plot_data)), limitsize = FALSE)
}

plot_static_kegg(file.path(GAM_STATIC_DIR, "Acetic_vs_Control_KEGG_data.csv"), "../figure/ef4b.pdf", "Acetic_vs_Control")
plot_static_kegg(file.path(GAM_STATIC_DIR, "H2O2_vs_Control_KEGG_data.csv"), "../figure/ef4c.pdf", "H2O2_vs_Control")


# ---- 5f. EFig.4D/E — GAM MaxDiff vs ESR (Python) ----
ret <- system2("python", c("Step3_efig4.py"))
if (ret != 0) cat("  [WARN] Step3_efig4.py failed\n")


# ---- 5g. EFig.4F/G — Full temporal GSEA heatmaps ----
for (nm in c("Acetic acid vs Control", "Hhdrogen peroxide vs Control")) {
  src <- file.path(GAM_TEMP_DIR, paste0("Heatmap_", nm, ".pdf"))
  dst <- if (grepl("Acetic", nm)) "../figure/ef4f.pdf" else "../figure/ef4g.pdf"
  if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
}


# ---- 5h. EFig.4 source map ----
ef4_map <- tribble(
  ~Panel, ~Output_file, ~Source_file,
  "ef4a", "ef4a.pdf", GAM_CORRECTED,
  "ef4b", "ef4b.pdf", file.path(GAM_STATIC_DIR, "Acetic_vs_Control_KEGG_data.csv"),
  "ef4c", "ef4c.pdf", file.path(GAM_STATIC_DIR, "H2O2_vs_Control_KEGG_data.csv"),
  "ef4d", "ef4d.pdf", GAM_RANKING,
  "ef4e", "ef4e.pdf", GAM_RANKING,
  "ef4f", "ef4f.pdf", file.path(GAM_TEMP_DIR, "Heatmap_Acetic acid vs Control.pdf"),
  "ef4g", "ef4g.pdf", file.path(GAM_TEMP_DIR, "Heatmap_Hhdrogen peroxide vs Control.pdf")
)
write_tsv(ef4_map, "../output/ef4_source_map.tsv")

cat("Step03 — Done.\n")
