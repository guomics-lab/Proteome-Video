# =============================================================================
# Step05_GSEA.R — Signed saliency GSEA, PC1 KEGG, and scatter panels
# Called by Step05_interpret.py after gradient computation.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Sc.sgd.db)
  library(ggrepel)
})

SANDBOX   <- "../.."
INTERMED  <- "../output"
GENE_ID   <- file.path(SANDBOX, "gene_id_mapping.tsv")

cat("Step05 — GSEA...\n")

GAM_DIR <- file.path(SANDBOX, "Step03_GAM", "output", "GAM_temporal_GSEA")

# Ensure binned GAM diffs exist (compute from Step02 if Step03 hasn't run)
mat_ac <- file.path(GAM_DIR, "Matrix_Binned_Diff_Acetic_vs_Control.csv")
mat_hc <- file.path(GAM_DIR, "Matrix_Binned_Diff_H2O2_vs_Control.csv")
if (!all(file.exists(mat_ac, mat_hc))) {
  library(mgcv)
  library(progress)
  SMOOTHED <- file.path(SANDBOX, "Step02_Clustering", "output", "smoothed_average.tsv")
  integrated_data <- read_tsv(SMOOTHED, show_col_types = FALSE)
  names(integrated_data)[1] <- "ProteinID"
  sample_ids <- names(integrated_data)[-1]
  metadata <- tibble(SampleID = sample_ids) %>%
    mutate(Treatment = str_extract(SampleID, "^[CHA]"),
           Position = as.numeric(str_extract(SampleID, "\\d+$")),
           Time = Position)
  integrated_data <- integrated_data %>%
    pivot_longer(-ProteinID, names_to = "SampleID", values_to = "Abundance") %>%
    left_join(metadata, by = "SampleID") %>%
    mutate(Treatment = factor(Treatment, levels = c("C", "H", "A"))) %>%
    dplyr::select(ProteinID, SampleID, Time, Treatment, Abundance) %>%
    arrange(ProteinID, Treatment, Time) %>% filter(!is.na(Abundance))
  time_grid <- seq(min(integrated_data$Time), max(integrated_data$Time), length.out = 96)
  calc_binned_gam <- function(pid, data, grid) {
    p_data <- data %>% filter(ProteinID == pid)
    if (nrow(p_data) < 20) return(NULL)
    model <- tryCatch(
      gam(Abundance ~ Treatment + s(Time, k = 7) + s(Time, by = Treatment, k = 7),
          data = p_data, method = "REML"), error = function(e) NULL)
    if (is.null(model)) return(NULL)
    nd <- expand_grid(Time = grid, Treatment = factor(c("C","A","H"), levels = c("C","H","A")))
    nd$Pred <- predict(model, newdata = nd)
    pc <- nd$Pred[nd$Treatment == "C"]; pa <- nd$Pred[nd$Treatment == "A"]; ph <- nd$Pred[nd$Treatment == "H"]
    bin_sum <- function(x) colSums(matrix(x, nrow = 6, ncol = 16))
    list(AC = bin_sum(pa - pc), HC = bin_sum(ph - pc))
  }
  proteins <- unique(integrated_data$ProteinID)
  pb <- progress_bar$new(total = length(proteins), format = "  [:bar] :percent eta: :eta", clear = FALSE, width = 60)
  results <- list()
  for (pid in proteins) { pb$tick(); out <- calc_binned_gam(pid, integrated_data, time_grid); if (!is.null(out)) results[[pid]] <- out }
  list_to_matrix <- function(rl, comp) {
    m <- do.call(rbind, lapply(rl, function(x) x[[comp]]))
    rownames(m) <- names(rl); colnames(m) <- sprintf("%.0f", seq(3, 93, by = 6)); m
  }
  dir.create(GAM_DIR, showWarnings = FALSE, recursive = TRUE)
  write.csv(list_to_matrix(results, "AC"), mat_ac)
  write.csv(list_to_matrix(results, "HC"), mat_hc)
}


# ---- Section A: Signed Saliency GSEA ----

id_mapper <- read_tsv(GENE_ID, show_col_types = FALSE) %>%
  dplyr::select(1, 3, 5) %>%
  set_names(c("Original_ID", "GENENAME", "ENTREZID")) %>%
  filter(!is.na(Original_ID)) %>%
  distinct(Original_ID, .keep_all = TRUE)

blacklist <- c("Viral","virus","HIV","Infection","Bacterial","Bacteria","Staphylococcus",
               "cancer","disease","Systemic","Lupus","Huntington","Alzheimer","Parkinson","Prion",
               "Methane","IgSF","Amoebiasis","Malaria","Leishmaniasis","Graft","Allograft","Autoimmune",
               "Metabolic pathways","Biosynthesis of secondary metabolites",
               "Microbial metabolism in diverse environments",
               "Autophagy - other","Autophagy – other",
               "Longevity regulating pathway - multiple species","Efferocytosis")
bl_pattern <- paste(blacklist, collapse = "|")

translate_gene_string <- function(id_string, lookup_vec) {
  ids <- str_split(id_string, "/")[[1]]
  names <- lookup_vec[ids]
  names[is.na(names)] <- ids[is.na(names)]
  paste(names, collapse = "/")
}

analysis_pairs <- list(
  list(name = "Acetic_vs_Control",
       grad_file = "Diff_Gradients_Acetic_vs_Control_temporal_gradients.csv",
       gam_file  = "Matrix_Binned_Diff_Acetic_vs_Control.csv"),
  list(name = "H2O2_vs_Control",
       grad_file = "Diff_Gradients_H2O2_vs_Control_temporal_gradients.csv",
       gam_file  = "Matrix_Binned_Diff_H2O2_vs_Control.csv")
)

run_signed_saliency <- function(pair_info) {
  aname <- pair_info$name

  f_grad <- file.path(INTERMED, pair_info$grad_file)
  f_gam  <- file.path(GAM_DIR, pair_info$gam_file)
  if (!file.exists(f_grad) || !file.exists(f_gam)) {
    return(NULL)
  }

  df_grad <- read.csv(f_grad, row.names = 1, check.names = FALSE)
  df_gam  <- read.csv(f_gam, row.names = 1, check.names = FALSE)

  # Rename T0-T15 → absolute time (3,9,...,93) in gradient columns
  abs_times <- as.character(seq(3, 93, by = 6))
  names(abs_times) <- paste0("T", 0:15)
  for (nm in names(abs_times)) {
    if (nm %in% colnames(df_grad)) {
      colnames(df_grad)[colnames(df_grad) == nm] <- abs_times[nm]
    }
  }

  common_genes <- intersect(rownames(df_grad), rownames(df_gam))
  common_cols  <- intersect(colnames(df_grad), colnames(df_gam))
  if (length(common_cols) == 0) {
    # Try: GAM cols might be "X3","X9"... vs gradient cols "3","9"...
    cat("    [WARN] No common columns, trying name fix...\n")
    colnames(df_gam) <- gsub("^X", "", colnames(df_gam))
    common_cols <- intersect(colnames(df_grad), colnames(df_gam))
  }
  if (length(common_cols) == 0) {
    cat("    [SKIP] Still no common columns\n"); return(NULL)
  }

  mat_grad <- as.matrix(df_grad[common_genes, common_cols])
  mat_gam  <- as.matrix(df_gam[common_genes, common_cols])
  mat_score <- abs(mat_grad) * sign(mat_gam)

  mapped <- as.data.frame(mat_score) %>%
    rownames_to_column("Original_ID") %>%
    left_join(id_mapper, by = "Original_ID") %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  if (nrow(mapped) < 10) return(NULL)

  # Export ranked list
  mapped %>%
    dplyr::select(GENENAME, ENTREZID, all_of(common_cols)) %>%
    pivot_longer(cols = -c(GENENAME, ENTREZID), names_to = "Time", values_to = "Score") %>%
    arrange(Time, desc(Score)) %>%
    write_csv(file.path(INTERMED, paste0("Ranked_List_", aname, ".csv")))

  abs_labels <- sprintf("%.0f", seq(0, length(common_cols) - 1) * 6 + 3)
  all_nes <- list()

  set.seed(20250707L)
  for (i in seq_along(common_cols)) {
    col_name <- common_cols[i]
    tp_label <- abs_labels[i]
    gene_list <- mapped[[col_name]]
    names(gene_list) <- mapped$ENTREZID
    gene_list <- sort(gene_list, decreasing = TRUE)

    gse_res <- tryCatch({
      gseKEGG(geneList = gene_list, organism = "sce", keyType = 'ncbi-geneid',
              pvalueCutoff = 1.0, minGSSize = 10, maxGSSize = 500,
              verbose = FALSE, nPermSimple = 1000)
    }, error = function(e) NULL)

    if (!is.null(gse_res) && nrow(gse_res) > 0) {
      all_nes[[tp_label]] <- as.data.frame(gse_res) %>%
        dplyr::select(ID, Description, NES, pvalue, core_enrichment) %>%
        mutate(TimePoint = tp_label)
    }
  }

  if (length(all_nes) == 0) { cat("    No pathways\n"); return(NULL) }
  combined <- bind_rows(all_nes)

  hm_mat <- combined %>% dplyr::select(Description, TimePoint, NES) %>%
    pivot_wider(names_from = TimePoint, values_from = NES, values_fill = 0) %>%
    column_to_rownames("Description")
  pv_mat <- combined %>% dplyr::select(Description, TimePoint, pvalue) %>%
    pivot_wider(names_from = TimePoint, values_from = pvalue, values_fill = 1) %>%
    column_to_rownames("Description")

  ordered_cols <- abs_labels[abs_labels %in% colnames(hm_mat)]
  hm_mat <- hm_mat[, ordered_cols, drop = FALSE]
  pv_mat <- pv_mat[, ordered_cols, drop = FALSE]

  keep_clean <- !grepl(bl_pattern, rownames(hm_mat), ignore.case = TRUE)
  is_valid  <- abs(hm_mat) > 1.5 & pv_mat < 0.05
  keep_gen  <- rowSums(is_valid, na.rm = TRUE) >= 2
  early_cols <- intersect(c("3","9"), colnames(is_valid))
  keep_early <- if (length(early_cols) > 0) rowSums(is_valid[, early_cols, drop = FALSE], na.rm = TRUE) >= 1 else FALSE
  final_keep <- keep_clean & (keep_gen | keep_early)

  # Rescue specific pathways of interest
  if (aname == "Acetic_vs_Control") {
    rescue_pathways <- c("Basal transcription factors")
    rescue_idx <- which(rownames(hm_mat) %in% rescue_pathways)
    if (length(rescue_idx) > 0) {
      final_keep[rescue_idx] <- TRUE
    }
  }

  hm_filt  <- hm_mat[final_keep, , drop = FALSE]
  pv_filt  <- pv_mat[final_keep, , drop = FALSE]

  if (nrow(hm_filt) > 1) {
    write.csv(hm_filt, file.path(INTERMED, paste0(aname, "_Filtered_NES.csv")))
    entrez_map <- setNames(id_mapper$GENENAME, id_mapper$ENTREZID)
    combined %>% filter(Description %in% rownames(hm_filt)) %>%
      rowwise() %>% mutate(Genes = translate_gene_string(core_enrichment, entrez_map)) %>%
      dplyr::select(TimePoint, Description, NES, pvalue, Genes) %>%
      arrange(Description, match(TimePoint, ordered_cols)) %>%
      write_csv(file.path(INTERMED, paste0(aname, "_Detailed_Info.csv")))

    sig_marker <- matrix("", nrow = nrow(pv_filt), ncol = ncol(pv_filt))
    sig_marker[pv_filt < 0.05] <- "*"
    sig_marker[pv_filt < 0.01] <- "**"

    panel_name <- if (aname == "Acetic_vs_Control") "fig5d" else "fig5e"
    pdf_out <- paste0("../figure/", panel_name, ".pdf")
    plot_height <- max(6, nrow(hm_filt) * 0.25)

    pheatmap(hm_filt, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE,
             color = colorRampPalette(c("navy","white","firebrick3"))(100),
             border_color = "grey95",
             display_numbers = sig_marker, number_color = "black",
             fontsize_number = 10, fontsize_row = 10, fontsize_col = 10,
             main = paste0(aname, " (Signed Saliency)"),
             filename = pdf_out, width = 12, height = plot_height)
  }
}

for (pair in analysis_pairs) {
  run_signed_saliency(pair)
}


# ---- Section B: PC1 KEGG ----

pc1_csv <- file.path(INTERMED, "gsea_ranked_list_PC1.csv")
if (!file.exists(pc1_csv)) {
  cat("  [WARN] PC1 ranking not found\n")
} else {
  id_map_kegg <- read_tsv(GENE_ID, show_col_types = FALSE) %>%
    dplyr::select(1, 3, 5) %>%
    set_names(c("Original_ID", "GENENAME", "ENTREZID")) %>%
    filter(!is.na(Original_ID)) %>%
    distinct(Original_ID, .keep_all = TRUE)

  pc1_df <- read_csv(pc1_csv, show_col_types = FALSE)
  colnames(pc1_df)[1] <- "gene_name"
  metric_col <- names(pc1_df)[2]

  ranked <- pc1_df %>%
    left_join(id_map_kegg, by = c("gene_name" = "Original_ID")) %>%
    filter(!is.na(ENTREZID), ENTREZID != "") %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    arrange(desc(.data[[metric_col]]))

  kegg_list <- ranked[[metric_col]]
  names(kegg_list) <- ranked$ENTREZID

  set.seed(20250707L)
  gse_res <- tryCatch({
    gseKEGG(geneList = kegg_list, organism = "sce", keyType = 'ncbi-geneid',
            pvalueCutoff = 0.05, pAdjustMethod = "BH",
            minGSSize = 5, maxGSSize = 500, verbose = FALSE)
  }, error = function(e) NULL)

  if (!is.null(gse_res) && nrow(gse_res) > 0) {
    # Translate core enrichment
    e2g <- AnnotationDbi::mapIds(org.Sc.sgd.db, keys(org.Sc.sgd.db, "ENTREZID"),
                                  column = "GENENAME", keytype = "ENTREZID")
    df_res <- as.data.frame(gse_res)
    df_res$core_enrichment <- sapply(df_res$core_enrichment, function(cs) {
      ids <- strsplit(cs, "/")[[1]]
      gn <- e2g[ids]; gn[is.na(gn)] <- ids[is.na(gn)]
      paste(gn, collapse = "/")
    })

    write_csv(df_res, file.path(INTERMED, "PC1_KEGG_data.csv"))

    plot_data <- df_res %>%
      mutate(Regulation = if_else(NES > 0, "Activated", "Suppressed")) %>%
      group_by(Regulation) %>%
      slice_min(order_by = p.adjust, n = 15, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(Description = str_wrap(Description, width = 50)) %>%
      arrange(Regulation, NES) %>%
      mutate(Description = forcats::fct_inorder(Description))

    if (nrow(plot_data) > 0) {
      p <- ggplot(plot_data, aes(x = NES, y = Description)) +
        geom_point(aes(size = setSize, color = p.adjust)) +
        facet_grid(Regulation ~ ., scales = "free_y", space = "free_y") +
        scale_size_continuous(name = "Gene Count", range = c(3, 8)) +
        scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
        labs(title = "PC1 — KEGG GSEA", x = "Normalized Enrichment Score (NES)", y = NULL) +
        theme_bw(base_size = 14) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
              axis.text.y = element_text(size = 12, lineheight = 0.9),
              strip.text.y = element_text(face = "bold", size = 14),
              panel.spacing = unit(1, "lines"), legend.position = "right")

      plot_height <- max(6, 2 + 0.5 * nrow(plot_data))
      ggsave("../figure/fig5c.pdf", p, width = 8, height = plot_height, device = cairo_pdf, limitsize = FALSE)
    }
  } else {
    cat("  [WARN] No significant KEGG pathways for PC1\n")
  }
}


# ---- Section C: ExtFig5C/D copies ----
for (pair in list(c("fig5d.pdf", "efig5c.pdf"), c("fig5e.pdf", "efig5d.pdf"))) {
  src <- file.path("../figure", pair[1])
  dst <- file.path("../figure", pair[2])
  if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
}


# ---- Section D: ExtFig5A/B Signed Saliency Scatters ----

TIME_MAP <- c("3"=1, "9"=2, "15"=3, "21"=4, "27"=5, "33"=6, "39"=7, "45"=8,
              "51"=9, "57"=10, "63"=11, "69"=12, "75"=13, "81"=14, "87"=15, "93"=16)

minmax_normalize_signed <- function(x) {
  x / max(abs(x), na.rm = TRUE)
}

plot_signed_saliency_scatter <- function(treatment, time_stage, gene_str,
                                          top_n_labels, output_file, panel_label) {
  # File paths
  f_grad <- file.path(INTERMED, paste0("Diff_Gradients_", treatment, "_temporal_gradients.csv"))
  f_gam  <- file.path(GAM_DIR, paste0("Matrix_Binned_Diff_", treatment, ".csv"))
  if (!file.exists(f_grad) || !file.exists(f_gam)) {
    return()
  }

  df_grad <- read.csv(f_grad, row.names = 1, check.names = FALSE)
  df_gam  <- read.csv(f_gam, row.names = 1, check.names = FALSE)

  # Map time string to column names (matching ef5a.R logic)
  time_numeric <- as.character(as.numeric(gsub("[^0-9]", "", time_stage)))
  time_idx     <- TIME_MAP[time_numeric]
  col_grad     <- sprintf("T%d", time_idx - 1)
  gam_col      <- time_numeric
  if (!gam_col %in% colnames(df_gam)) {
    gam_col <- as.character(as.integer(time_numeric))
  }

  if (!col_grad %in% colnames(df_grad)) return()
  if (!gam_col %in% colnames(df_gam)) return()

  common <- intersect(rownames(df_grad), rownames(df_gam))

  v_gam  <- df_gam[common, gam_col]
  v_grad <- df_grad[common, col_grad]
  v_sal  <- abs(v_grad) * sign(v_gam)
  v_gam_norm <- minmax_normalize_signed(v_gam)
  v_sal_norm <- minmax_normalize_signed(v_sal)

  # Parse gene list (strip optional scores in parentheses, matching ef5a.R parse_gene_string)
  raw_genes <- trimws(unlist(strsplit(gene_str, ";")))
  raw_genes <- raw_genes[raw_genes != ""]
  gene_list <- toupper(gsub("\\(.*\\)", "", raw_genes))
  gene_list <- trimws(gene_list)

  hl_found  <- gene_list[gene_list %in% common]
  lbl_genes <- gene_list[1:min(top_n_labels, length(gene_list))]
  lbl_found <- lbl_genes[lbl_genes %in% common]

  # Stats
  r_pearson  <- cor(v_gam_norm, v_sal_norm, method = "pearson")
  r_spearman <- cor(v_gam_norm, v_sal_norm, method = "spearman")
  q_pp <- sum(v_gam_norm > 0 & v_sal_norm > 0)
  q_pn <- sum(v_gam_norm > 0 & v_sal_norm < 0)
  q_np <- sum(v_gam_norm < 0 & v_sal_norm > 0)
  q_nn <- sum(v_gam_norm < 0 & v_sal_norm < 0)

  # Build data frame
  dat <- data.frame(
    GAM_Norm = v_gam_norm,
    Sal_Norm = v_sal_norm,
    Gene     = common,
    Highlight = common %in% hl_found,
    Label     = common %in% lbl_found,
    Gene_Display = "",
    stringsAsFactors = FALSE
  )

  # Format gene names: capitalize first letter, lowercase rest
  for (i in which(dat$Label)) {
    g <- dat$Gene[i]
    dat$Gene_Display[i] <- paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))))
  }

  dat_regular   <- dat[!dat$Highlight, ]
  dat_highlight <- dat[dat$Highlight, ]
  dat_labeled   <- dat[dat$Label, ]

  treatment_label <- gsub("_", " ", treatment)

  p <- ggplot() +
    geom_point(data = dat_regular,
               aes(x = GAM_Norm, y = Sal_Norm),
               color = "gray40", size = 0.8, alpha = 0.3) +
    geom_point(data = dat_highlight,
               aes(x = GAM_Norm, y = Sal_Norm),
               color = "red", size = 2.5, alpha = 0.8) +
    {if (nrow(dat_labeled) > 0) {
      geom_text_repel(data = dat_labeled,
                      aes(x = GAM_Norm, y = Sal_Norm, label = Gene_Display),
                      size = 3.5, color = "darkred", fontface = "bold",
                      box.padding = 0.8, point.padding = 0.5,
                      segment.color = "red", segment.size = 0.4,
                      max.overlaps = 30, min.segment.length = 0)
    }} +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.4) +
    geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed", linewidth = 0.5, alpha = 0.6) +
    scale_x_continuous(limits = c(-1, 1), expand = c(0.02, 0), breaks = seq(-1, 1, by = 0.5)) +
    scale_y_continuous(limits = c(-1, 1), expand = c(0.02, 0), breaks = seq(-1, 1, by = 0.5)) +
    coord_fixed(ratio = 1) +
    labs(
      title = sprintf("%s: %s — %s", panel_label, treatment_label, time_stage),
      subtitle = sprintf("r=%.3f, ρ=%.3f | ++:%d, +−:%d, −+:%d, −−:%d | Highlighted: %d genes",
                         r_pearson, r_spearman, q_pp, q_pn, q_np, q_nn, length(hl_found)),
      x = "GAM (normalized)", y = "Signed Saliency (normalized)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text  = element_text(size = 9)
    )

  ggsave(file.path("../figure", output_file), p, width = 8, height = 8, device = cairo_pdf)
}


# EFig5A1: Acetic 3min, 4 genes, top4 labeled
plot_signed_saliency_scatter("Acetic_vs_Control", "3 min",
  "TAF3;SPT15;TOA2;TFB2", top_n_labels = 4,
  output_file = "efig5a1.pdf", panel_label = "EFig5A1")

# EFig5A2: Acetic 3min, 16 genes, top5 labeled
plot_signed_saliency_scatter("Acetic_vs_Control", "3 min",
  "ATP18;CYC1;RIP1;COX13;VMA7;ATP14;COX4;ATP4;ATP19;ATP7;STV1;VMA10;QCR2;TIM11;VMA6;ATP16",
  top_n_labels = 5, output_file = "efig5a2.pdf", panel_label = "EFig5A2")

# EFig5A3: Acetic 27min, 21 genes, top21 labeled
plot_signed_saliency_scatter("Acetic_vs_Control", "27 min",
  "CCZ1;ARC18;KOG1;UME6;PEP5;VPS30;ATG13;ARC35;VPS45;LCB1;VPS33;TAP42;SNF1;ARG82;PPH22;TPK3;PHO85;KCS1;RIM15;SUI2;TOR2",
  top_n_labels = 21, output_file = "efig5a3.pdf", panel_label = "EFig5A3")

# EFig5B1: H2O2 3min, 10 genes, top5 labeled
plot_signed_saliency_scatter("H2O2_vs_Control", "3 min",
  "DPB4;CDC9;RFC4;POL3;DPB2;RAD27;DPB3;POL30;POL31;UNG1",
  top_n_labels = 5, output_file = "efig5b1.pdf", panel_label = "EFig5B1")

# EFig5B2: H2O2 21min, 9 mitophagy genes, top5 labeled
plot_signed_saliency_scatter("H2O2_vs_Control", "21 min",
  "TOR2;MDM10;CKA2;CKB1;UBP3;UME6;SIN3;CKA1;FMC1",
  top_n_labels = 5, output_file = "efig5b2.pdf", panel_label = "EFig5B2")

# EFig5B3: H2O2 21min, 9 genes, top5 labeled
plot_signed_saliency_scatter("H2O2_vs_Control", "21 min",
  "YAT1;FAA3;PEX19;IDP2;PEX14;NPY1;RSM26;IDP1;FAA4",
  top_n_labels = 5, output_file = "efig5b3.pdf", panel_label = "EFig5B3")


cat("Step05_GSEA — Done.\n")
