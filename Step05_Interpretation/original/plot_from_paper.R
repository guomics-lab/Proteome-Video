# =============================================================================
# plot_from_paper.R — Reproduce published Signed Saliency figures from
# original GSEA intermediates (the original gseKEGG runs lacked a fixed seed)
# =============================================================================

script_dir <- dirname(normalizePath(sub("^--file=", "",
  commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])))
setwd(script_dir)

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

DATA_DIR <- "."

cat("Plotting from original GSEA intermediates...\n")


plot_signed_heatmap <- function(analysis_name, output_file, panel_label) {
  nes_file    <- file.path(DATA_DIR, paste0(analysis_name, "_Filtered_NES.csv"))
  detail_file <- file.path(DATA_DIR, paste0(analysis_name, "_Detailed_Info.csv"))

  if (!file.exists(nes_file) || !file.exists(detail_file)) {
    cat(sprintf("  [WARN] Missing data for %s\n", analysis_name)); return()
  }

  hm_mat  <- read.csv(nes_file, row.names = 1, check.names = FALSE)
  pval_df <- read_csv(detail_file, show_col_types = FALSE)

  pval_wide <- pval_df %>%
    dplyr::select(Description, TimePoint, pvalue) %>%
    mutate(TimePoint = as.character(TimePoint)) %>%
    pivot_wider(names_from = TimePoint, values_from = pvalue, values_fill = 1) %>%
    column_to_rownames("Description")
  colnames(pval_wide) <- as.character(colnames(pval_wide))

  pval_filt <- matrix(1, nrow = nrow(hm_mat), ncol = ncol(hm_mat),
                       dimnames = dimnames(hm_mat))
  cr <- intersect(rownames(pval_wide), rownames(pval_filt))
  cc <- intersect(colnames(pval_wide), colnames(pval_filt))
  pval_filt[cr, cc] <- as.matrix(pval_wide[cr, cc])
  pval_filt[is.na(pval_filt)] <- 1

  sig_marker <- matrix("", nrow = nrow(pval_filt), ncol = ncol(pval_filt))
  sig_marker[pval_filt < 0.05] <- "*"
  sig_marker[pval_filt < 0.01] <- "**"

  plot_height <- max(6, nrow(hm_mat) * 0.25)

  pheatmap(hm_mat, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           border_color = "grey95",
           display_numbers = sig_marker, number_color = "black",
           fontsize_number = 10, fontsize_row = 10, fontsize_col = 10,
           main = paste0(panel_label, " (Signed Saliency)"),
           filename = output_file, width = 12, height = plot_height)
}


# ---- Signed Saliency heatmaps ----
plot_signed_heatmap("Acetic_vs_Control", "Original_paper_fig5d.pdf", "Fig5D: Acetic vs Control")
plot_signed_heatmap("H2O2_vs_Control",   "Original_paper_fig5e.pdf", "Fig5E: H2O2 vs Control")

# ---- EFig5C/D copies ----
for (pair in list(c("Original_paper_fig5d.pdf", "Original_paper_efig5c.pdf"),
                   c("Original_paper_fig5e.pdf", "Original_paper_efig5d.pdf"))) {
  if (file.exists(pair[1])) file.copy(pair[1], pair[2], overwrite = TRUE)
}

cat("Done.\n")
