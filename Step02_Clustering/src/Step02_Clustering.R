# =============================================================================
# Step02_Clustering.R — Fig.2, EFig.2, EFig.3
#
# maSigPro clustering → WGCNA network → GO/KEGG enrichment → ESR Fisher test
# Prerequisites: Step2_prepare.py, Step00_Preprocess completed
# =============================================================================

rm(list = ls())
setwd(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))))

suppressPackageStartupMessages({
  library(tidyverse)
  library(conflicted)
})

conflicts_prefer(base::setdiff)
conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)

# ---- Paths ----
IMPUTED    <- "../../Step00_Preprocess/output/imputed_temporal_matrix.tsv"
ID_MAPPING <- "../../gene_id_mapping.tsv"
F_SMOOTHED <- "../output/smoothed_average.tsv"
F_DOWNSAMP <- "../output/downsampled_median.tsv"

cat("Step02 — Clustering...\n")

dir.create("../output", recursive = TRUE, showWarnings = FALSE)
dir.create("../figure", recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# PART 1 — DATA PREPARATION
# =============================================================================

# ---- 1a. Smoothing & Downsampling (Python) ----
ret <- system2("python", c("Step2_prepare.py"))
if (ret != 0) stop("Step2_prepare.py failed.")


# ---- 1b. FFT Optimal Binning ----

protein_data <- read_tsv(IMPUTED, show_col_types = FALSE)
sample_cols <- setdiff(colnames(protein_data), "Genes")

metadata <- tibble(SampleID = sample_cols) %>%
  mutate(
    Treatment = str_extract(SampleID, "^[CHA]"),
    Time = as.numeric(str_extract(SampleID, "\\d+$"))
  ) %>%
  arrange(factor(Treatment, levels = c("C", "H", "A")), Time)

create_time_series_matrix <- function(treatment_name) {
  trt_samples <- metadata %>%
    filter(Treatment == treatment_name) %>% arrange(Time)
  protein_data %>%
    column_to_rownames("Genes") %>%
    dplyr::select(all_of(trt_samples$SampleID)) %>% as.matrix()
}

calculate_protein_psd <- function(protein_timeseries, fs) {
  detrended <- protein_timeseries - mean(protein_timeseries, na.rm = TRUE)
  detrended[is.na(detrended)] <- mean(detrended, na.rm = TRUE)
  n <- length(detrended)
  fft_result <- fft(detrended)
  half_n <- floor(n / 2) + 1
  psd <- (Mod(fft_result[1:half_n]))^2 / (fs * n)
  freqs <- seq(0, fs / 2, length.out = half_n)
  data.frame(frequency = freqs, psd = psd)
}

treatments_fft <- c("C", "H", "A")
fs <- 1.0

treatment_psds <- map(treatments_fft, function(tr) {
  mat <- create_time_series_matrix(tr)
  n_proteins <- nrow(mat)
  all_psds <- map(seq_len(n_proteins), ~ calculate_protein_psd(as.numeric(mat[.x, ]), fs))
  bind_rows(all_psds) %>%
    group_by(frequency) %>%
    summarise(mean_psd = mean(psd, na.rm = TRUE),
              median_psd = median(psd, na.rm = TRUE),
              sd_psd = sd(psd, na.rm = TRUE), .groups = "drop")
})
names(treatment_psds) <- treatments_fft

find_elbow_point <- function(freqs, psds) {
  valid <- freqs > 0 & psds > 0
  log_freqs <- log10(freqs[valid])
  log_psds  <- log10(psds[valid])
  x1 <- log_freqs[1]; y1 <- log_psds[1]
  x2 <- tail(log_freqs, 1); y2 <- tail(log_psds, 1)
  a <- y2 - y1; b <- -(x2 - x1); c <- x2 * y1 - y2 * x1
  distances <- abs(a * log_freqs + b * log_psds + c) / sqrt(a^2 + b^2)
  idx <- which.max(distances)
  list(frequency = 10^log_freqs[idx], log_frequency = log_freqs[idx], log_psd = log_psds[idx])
}

elbow_results <- map(treatments_fft, function(tr) {
  res <- find_elbow_point(treatment_psds[[tr]]$frequency, treatment_psds[[tr]]$mean_psd)
  res$treatment <- tr; res
})
names(elbow_results) <- treatments_fft

f_cutoffs <- map_dbl(elbow_results, "frequency")
f_cutoff_final <- max(f_cutoffs)
bin_size_optimal <- max(round(1 / (2 * f_cutoff_final)), 2)

plot_psd <- map_dfr(treatments_fft, ~ mutate(treatment_psds[[.x]], treatment = .x)) %>%
  filter(frequency > 0, mean_psd > 0) %>%
  mutate(log_f = log10(frequency), log_p = log10(mean_psd))

plot_elbow <- map_dfr(elbow_results, as_tibble)

p_fft <- ggplot(plot_psd, aes(x = log_f, y = log_p, color = treatment)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_point(data = plot_elbow, aes(x = log_frequency, y = log_psd),
             size = 4, shape = 19, color = "black") +
  geom_point(data = plot_elbow, aes(x = log_frequency, y = log_psd, color = treatment),
             size = 3, shape = 19) +
  geom_vline(xintercept = log10(f_cutoff_final), linetype = "dashed") +
  scale_color_manual(values = c("C" = "#2E8B57", "H" = "#1E90FF", "A" = "#FF6347")) +
  labs(title = "FFT-based Optimal Binning Analysis",
       subtitle = sprintf("Cutoff: %.3f Hz | Optimal Bin Size: %d", f_cutoff_final, bin_size_optimal),
       x = "log10(Frequency) [Hz]", y = "log10(Mean Power)") +
  theme_bw() + theme(legend.position = "top")

ggsave("../output/G02_FFT_Analysis.pdf", p_fft, width = 8, height = 6)
cat(sprintf("  FFT optimal bin size: %d\n", bin_size_optimal))


# ---- 1c. maSigPro Temporal Clustering ----

library(maSigPro)
library(openxlsx)

masigpro_dir <- "../output/maSigPro_Final_Results"
if (!dir.exists(masigpro_dir)) dir.create(masigpro_dir, recursive = TRUE)

bin_size <- 3; polynomial_degree <- 4; fdr_cutoff <- 0.05
r_squared_cutoff <- 0.75; num_clusters <- 6
group_colors <- c("#FF6347", "#2E8B57", "#1E90FF")

data_matrix_original <- read_tsv(F_SMOOTHED, show_col_types = FALSE) %>%
  column_to_rownames("Genes") %>% as.matrix()

original_colnames <- colnames(data_matrix_original)
original_info <- tibble(
  SampleID = original_colnames,
  Group = str_extract(SampleID, "^[A-Z]"),
  Time  = as.numeric(str_extract(SampleID, "\\d+$"))
)

binned_matrix_list <- list()
for (tg in c("C", "H", "A")) {
  group_info <- original_info %>% filter(Group == tg) %>% arrange(Time)
  num_samples <- nrow(group_info)
  num_bins <- floor(num_samples / bin_size)
  for (i in seq_len(num_bins)) {
    idx_start <- (i - 1) * bin_size + 1
    idx_end   <- i * bin_size
    cur_samples_info <- group_info[idx_start:idx_end, ]
    new_time <- median(cur_samples_info$Time)
    binned_data <- data_matrix_original[, cur_samples_info$SampleID, drop = FALSE]
    colnames(binned_data) <- paste0(tg, new_time, "-", seq_len(bin_size))
    binned_matrix_list[[length(binned_matrix_list) + 1]] <- binned_data
  }
}
data_matrix_binned <- do.call(cbind, binned_matrix_list)

saveRDS(data_matrix_binned, file.path(masigpro_dir, "data_matrix_binned.rds"))

binned_colnames <- colnames(data_matrix_binned)
edesign_info <- tibble(
  SampleID = binned_colnames,
  Time = as.numeric(str_extract(SampleID, "(?<=^[CHA])\\d+\\.?\\d*")),
  Group = str_extract(SampleID, "^[CHA]")
) %>%
  mutate(ReplicateGroup = paste(Group, Time, sep = "_"),
         Replicates = as.numeric(factor(ReplicateGroup)))

edesign <- edesign_info %>%
  mutate(
    Control    = ifelse(Group == "C", 1, 0),
    H2O2       = ifelse(Group == "H", 1, 0),
    AceticAcid = ifelse(Group == "A", 1, 0)
  ) %>%
  dplyr::select(Time, Replicates, Control, H2O2, AceticAcid) %>%
  as.data.frame()
rownames(edesign) <- binned_colnames

design <- make.design.matrix(edesign, degree = polynomial_degree)
fit <- p.vector(data_matrix_binned, design, Q = fdr_cutoff,
                min.obs = floor(ncol(data_matrix_binned) * 0.7))
tstep <- T.fit(fit, step.method = "backward", alfa = fdr_cutoff)
sigs <- get.siggenes(tstep, rsq = r_squared_cutoff, vars = "groups")

if (is.null(sigs$summary) || nrow(sigs$summary) == 0) {
  stop("No significant genes found. Try lowering r_squared_cutoff.")
}

wb <- createWorkbook()

for (comparison_name in names(sigs$sig.genes)) {
  current_sig_genes <- sigs$sig.genes[[comparison_name]]
  if (is.null(current_sig_genes) || is.null(current_sig_genes$sig.profiles) ||
      nrow(current_sig_genes$sig.profiles) == 0) next

  cat(sprintf("  Processing: %s\n", comparison_name))
  pdf_filename <- file.path(masigpro_dir,
                            paste0("Clusters_", gsub(" ", "", comparison_name), ".pdf"))
  pdf(pdf_filename, width = 12, height = 9)
  palette(group_colors)
  see.genes_result <- see.genes(current_sig_genes, k = num_clusters,
                                group.cols = c("Control", "H2O2", "AceticAcid"),
                                show.fit = TRUE, dis = design$dis)
  palette("default")
  dev.off()

  if (!is.null(see.genes_result$cut)) {
    cluster_table <- tibble(
      Gene    = names(see.genes_result$cut),
      Cluster = paste0("Cluster_", unlist(see.genes_result$cut))
    )
    sheet_name <- paste0("Clusters_", gsub(" ", "", comparison_name))
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, cluster_table)
  }
}

if (!is.null(sigs$summary)) {
  addWorksheet(wb, "Summary_Table")
  writeData(wb, "Summary_Table", sigs$summary)
}

excel_filename <- file.path(masigpro_dir, "maSigPro_Analysis_Results.xlsx")
saveWorkbook(wb, excel_filename, overwrite = TRUE)

# Global Correlation
expr_mat <- readRDS(file.path(masigpro_dir, "data_matrix_binned.rds"))

sheet_names <- getSheetNames(excel_filename)
cluster_sheets <- grep("^Clusters_", sheet_names, value = TRUE)
target_sheet <- cluster_sheets[1]

cluster_df <- read.xlsx(excel_filename, sheet = target_sheet) %>%
  dplyr::select(Gene, Cluster) %>% distinct()

cluster_df$Global_Correlation <- NA_real_
unique_clusters <- unique(cluster_df$Cluster)

for (clus in unique_clusters) {
  genes_in_clus <- cluster_df %>% filter(Cluster == clus) %>% pull(Gene)
  valid_genes <- intersect(genes_in_clus, rownames(expr_mat))
  if (length(valid_genes) > 1) {
    sub_mat <- expr_mat[valid_genes, , drop = FALSE]
    centroid_profile <- colMeans(sub_mat, na.rm = TRUE)
    cors <- apply(sub_mat, 1, function(gp) cor(gp, centroid_profile, method = "pearson"))
    idx <- which(cluster_df$Gene %in% names(cors) & cluster_df$Cluster == clus)
    cluster_df$Global_Correlation[idx] <- cors[cluster_df$Gene[idx]]
  }
}

saveRDS(cluster_df, file.path(masigpro_dir, "cluster_info_final_with_global_cor.rds"))


# ---- 1d. GO / KEGG Enrichment ----

library(clusterProfiler)
library(org.Sc.sgd.db)

top_n_terms <- 10
background_genes <- read_tsv(F_SMOOTHED, show_col_types = FALSE) %>% pull(1)
conversion_results <- read_tsv(ID_MAPPING, show_col_types = FALSE)

sheet_names <- getSheetNames(excel_filename)
cluster_sheet <- sheet_names[grepl("^Clusters_", sheet_names)][1]
cluster_data <- read.xlsx(excel_filename, sheet = cluster_sheet)

gene_lookup <- conversion_results %>%
  filter(Status %in% c("SUCCESS", "PARTIAL")) %>%
  dplyr::select(Original_ID, Final_ID) %>% deframe()

cluster_gene_lists <- list()
min_cluster_size <- 5
for (cluster_id in unique(cluster_data$Cluster)) {
  core_genes_raw <- cluster_data %>% filter(Cluster == cluster_id) %>% pull(Gene)
  core_genes_standard <- sapply(core_genes_raw,
    function(g) if (g %in% names(gene_lookup)) gene_lookup[[g]] else g,
    USE.NAMES = FALSE)
  core_genes_final <- core_genes_standard[core_genes_standard %in% background_genes]
  if (length(core_genes_final) >= min_cluster_size) {
    cluster_gene_lists[[cluster_id]] <- core_genes_final
  }
}
perform_go_enrichment <- function(gene_lists, ontology = "BP", universe = background_genes) {
  tryCatch({
    res <- compareCluster(geneClusters = gene_lists, fun = "enrichGO",
      OrgDb = org.Sc.sgd.db, keyType = "GENENAME", ont = ontology,
      universe = universe, pAdjustMethod = "BH", pvalueCutoff = 0.05,
      qvalueCutoff = 0.2, minGSSize = 3, maxGSSize = 500)
    if (!is.null(res) && nrow(res@compareClusterResult) > 0) {
      cat(sprintf("    GO %s: %d terms\n", ontology, nrow(res@compareClusterResult)))
      return(res)
    }
    cat(sprintf("    GO %s: no terms found\n", ontology)); NULL
  }, error = function(e) { cat(sprintf("    GO %s error: %s\n", ontology, e$message)); NULL })
}

go_bp_result <- perform_go_enrichment(cluster_gene_lists, "BP")

# Export raw GO BP results for EFig.3B (before refinement)
raw_bp <- as.data.frame(go_bp_result)
if (!is.null(raw_bp) && nrow(raw_bp) > 0) {
  raw_bp <- raw_bp %>%
    mutate(GeneRatio_numeric = as.numeric(sub("/.*", "", GeneRatio)) /
             as.numeric(sub(".*/", "", GeneRatio)),
           Description_wrapped = str_wrap(Description, width = 40))
  for (cl in c("Cluster_3", "Cluster_4")) {
    raw_bp %>% filter(Cluster == cl) %>%
      write_csv(paste0("../output/EF3B_GO_BP_", cl, "_Source.csv"))
  }
}

go_mf_result <- perform_go_enrichment(cluster_gene_lists, "MF")
go_cc_result <- perform_go_enrichment(cluster_gene_lists, "CC")

genename_to_orf <- conversion_results %>%
  filter(Status %in% c("SUCCESS", "PARTIAL") & !is.na(ORF) & ORF != "") %>%
  dplyr::select(Final_ID, ORF) %>% distinct() %>% deframe()

background_orfs <- na.omit(sapply(background_genes,
  function(g) if (g %in% names(genename_to_orf)) genename_to_orf[[g]] else NA,
  USE.NAMES = FALSE))

cluster_orf_lists <- map(cluster_gene_lists, ~ na.omit(sapply(.x,
  function(g) if (g %in% names(genename_to_orf)) genename_to_orf[[g]] else NA,
  USE.NAMES = FALSE))) %>% keep(~ length(.x) >= 5)

kegg_result <- NULL
if (length(cluster_orf_lists) > 0) {
  kegg_result <- tryCatch(
    compareCluster(geneClusters = cluster_orf_lists, fun = "enrichKEGG",
      organism = "sce", universe = background_orfs, pAdjustMethod = "BH",
      pvalueCutoff = 0.05, qvalueCutoff = 0.2, minGSSize = 3, maxGSSize = 500),
    error = function(e) NULL)
}

refine_enrichment_data <- function(enrich_obj, ontology_name) {
  if (is.null(enrich_obj)) return(NULL)
  df_raw <- as.data.frame(enrich_obj)
  df_refined_list <- list()
  for (clus in unique(df_raw$Cluster)) {
    df_sub <- df_raw %>% filter(Cluster == clus)
    if (clus == "Cluster_3" && ontology_name == "GO_BP") {
      df_sub <- df_sub %>%
        mutate(temp_desc = tolower(Description),
               temp_desc = gsub("^mitochondrial ", "", temp_desc),
               temp_desc = gsub("complex ii", "complex 2", temp_desc)) %>%
        group_by(temp_desc) %>%
        slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) %>%
        ungroup() %>% dplyr::select(-temp_desc) %>%
        arrange(p.adjust) %>% head(top_n_terms)
    } else if (clus == "Cluster_4" && ontology_name == "GO_BP") {
      target_keywords <- c("rna processing", "rna splicing",
                           "rrna metabolic process", "ribosome biogenesis",
                           "protein-rna complex organization")
      best_hits <- list()
      for (k in target_keywords) {
        matches <- df_sub %>% filter(grepl(k, tolower(Description))) %>% arrange(p.adjust)
        if (nrow(matches) > 0) best_hits[[k]] <- matches[1, ]
      }
      df_sub <- bind_rows(best_hits) %>% distinct(ID, .keep_all = TRUE) %>% arrange(p.adjust)
    } else {
      df_sub <- df_sub %>% head(top_n_terms)
    }
    df_refined_list[[clus]] <- df_sub
  }
  bind_rows(df_refined_list)
}

refined_go_bp_df <- refine_enrichment_data(go_bp_result, "GO_BP")
saveRDS(refined_go_bp_df, "../output/refined_go_bp_df.rds")


# ---- 1e. WGCNA Network Analysis ----

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
})

conflicts_prefer(WGCNA::cor)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(base::intersect)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
set.seed(123)

expression_data <- read_tsv(IMPUTED, show_col_types = FALSE)
gene_ids <- expression_data$Genes
expression_matrix <- expression_data %>% select(-Genes) %>% as.data.frame()
rownames(expression_matrix) <- gene_ids
datExpr <- t(expression_matrix)

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

datExpr_clean <- datExpr[!rownames(datExpr) %in% c("C0"), ]

n_top_genes <- round(ncol(datExpr_clean) * 0.5)
variances <- apply(datExpr, 2, var, na.rm = TRUE)
top_genes <- names(sort(variances, decreasing = TRUE)[1:n_top_genes])
datExpr <- datExpr_clean[, top_genes]

sample_ids <- rownames(datExpr)
treatments <- substr(sample_ids, 1, 1)
times <- as.numeric(substr(sample_ids, 2, nchar(sample_ids)))

datTraits <- data.frame(
  is_Control     = as.numeric(treatments == "C"),
  is_H2O2        = as.numeric(treatments == "H"),
  is_Acetic      = as.numeric(treatments == "A"),
  Time           = times,
  Time_squared   = times^2,
  Time_in_H2O2   = ifelse(treatments == "H", times, 0),
  Time_in_Acetic = ifelse(treatments == "A", times, 0),
  Early_timepoint = as.numeric(times <= 25),
  Late_timepoint  = as.numeric(times >= 75),
  row.names = sample_ids
)

softPower <- 30
net <- blockwiseModules(
  datExpr, power = softPower, TOMType = "signed", networkType = "signed",
  minModuleSize = 5, deepSplit = 4, pamRespectsDendro = FALSE,
  mergeCutHeight = 0.2, numericLabels = TRUE, saveTOMs = TRUE,
  saveTOMFileBase = "M02_TOM", verbose = 0, maxBlockSize = 5000, corType = "bicor"
)

moduleColors <- labels2colors(net$colors)

# Dendrogram (for ef3d)
pdf("../output/Fig4_M02_Dendrogram_and_Modules.pdf", width = 16, height = 8)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
  "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
  guideHang = 0.05, main = "Gene dendrogram and module colors (Signed Network)",
  cex.main = 1.5, cex.lab = 1.2)
dev.off()

MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

common_samples <- intersect(rownames(MEs), rownames(datTraits))
MEs <- MEs[common_samples, ]
datTraits_sub <- datTraits[common_samples, ]

# Module-trait heatmap (for ef3e)
nSamples <- nrow(MEs)
moduleTraitCor <- cor(MEs, datTraits_sub, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
getSigStr <- function(p) {
  sapply(p, function(x) if(is.na(x)) "" else if(x<0.001) "***" else if(x<0.01) "**" else if(x<0.05) "*" else "")
}
textMatrix <- paste(signif(moduleTraitCor, 2), "\n", getSigStr(moduleTraitPvalue), sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
mc_names <- gsub("ME", "", names(MEs))
mc_sizes <- table(moduleColors)
yl <- paste0(names(MEs), " (", mc_sizes[mc_names], " proteins)")

pdf("../output/Fig4_M05_Module_Trait_Heatmap_Advanced.pdf", width = 11, height = 8)
par(mar = c(6, 12, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits_sub),
  yLabels = yl, ySymbols = names(MEs), colorLabels = FALSE,
  colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE,
  cex.text = 1.2, cex.lab.x = 1.0, cex.lab.y = 1.0, zlim = c(-1, 1),
  main = "Module-trait relationships", cex.main = 1.8)
dev.off()

# Export module assignments and eigengenes
write.csv(data.frame(Gene = colnames(datExpr), ModuleColor = moduleColors),
          "../output/WGCNA_module_assignments.csv", row.names = FALSE)

write.csv(cbind(datTraits_sub, MEs), "../output/WGCNA_eigengenes_traits.csv")

# Hub genes + Cytoscape export
modules_of_interest <- c("turquoise", "pink")

for (module_color in modules_of_interest) {
  inModule <- (moduleColors == module_color)
  module_MEs <- MEs[[paste0("ME", module_color)]]
  cor_results <- cor(datExpr[, inModule], module_MEs, use = "p")

  geneModuleMembership <- data.frame(
    ModuleMembership = as.vector(cor_results),
    row.names = rownames(cor_results)
  ) %>% arrange(desc(abs(ModuleMembership)))

  n_top <- 50
  top_genes <- rownames(head(geneModuleMembership, n_top))
  adj_sub <- adjacency(datExpr[, top_genes], power = softPower, type = "signed")
  edge_threshold <- 0.01

  edge_list <- which(adj_sub > edge_threshold, arr.ind = TRUE)
  edges <- data.frame(
    Source = rownames(adj_sub)[edge_list[, 1]],
    Target = colnames(adj_sub)[edge_list[, 2]],
    Weight = adj_sub[edge_list]
  ) %>% filter(Source != Target)

  nodes <- data.frame(
    ID = top_genes,
    ModuleMembership = geneModuleMembership[top_genes, "ModuleMembership"]
  )

  write.table(edges, paste0("../output/efig3_", module_color, "_edges_SourceData.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(nodes, paste0("../output/efig3_", module_color, "_nodes_SourceData.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
# ---- 1f. EFig.3: ESR Fisher test + bar chart (Python) ----
ret <- system2("python", c("Step2_efig3.py"))
if (ret != 0) stop("Step2_efig3.py failed.")


# =============================================================================
# PART 2 — FIGURE VISUALISATION
# =============================================================================

treatment_colors <- c("C" = "#2E8B57", "H" = "#1E90FF", "A" = "#FF6347")
treatment_names  <- c("C" = "Control", "H" = "H2O2", "A" = "Acetic Acid")


# ---- 2a. Fig.2A / 2B ----
fig2_ab <- tibble(
  Panel = c("Fig2A", "Fig2B"),
  Description = c("maSigPro Cluster 4 portrait", "maSigPro Cluster 3 portrait"),
  Upstream = "EFig.3A full six-cluster maSigPro portrait (../output/maSigPro_Final_Results/Clusters_*.pdf)",
  Reproduce_here = "No"
)
write_tsv(fig2_ab, "../output/fig2_ab_SourceData.tsv")


# ---- 2b. Fig.2C — Heatmap + GO ribbon ----

library(ComplexHeatmap)
library(circlize)
library(grid)

refined_go_bp_df <- readRDS("../output/refined_go_bp_df.rds")
raw_mat <- readRDS(file.path(masigpro_dir, "data_matrix_binned.rds"))
cluster_info <- readRDS(file.path(masigpro_dir, "cluster_info_final_with_global_cor.rds"))

conversion_df <- read_tsv(ID_MAPPING, show_col_types = FALSE) %>%
  filter(Status %in% c("SUCCESS", "PARTIAL")) %>%
  dplyr::select(Original_ID, Final_ID) %>% distinct()

std_to_orig_map <- setNames(conversion_df$Original_ID, conversion_df$Final_ID)
orig_to_std_map <- setNames(conversion_df$Final_ID, conversion_df$Original_ID)

cols <- colnames(raw_mat)
col_meta <- tibble(SampleID = cols) %>%
  mutate(
    Group = str_extract(SampleID, "^[CHA]"),
    Time  = as.numeric(str_extract(SampleID, "(?<=^[CHA])\\d+")),
    UniquePoint = paste(Group, Time, sep = "_"))

median_mat <- t(apply(raw_mat, 1, function(v) tapply(v, col_meta$UniquePoint, median)))

groups <- c("C", "H", "A")
time_points <- sort(unique(col_meta$Time))
ordered_cols <- c()
for (g in groups) ordered_cols <- c(ordered_cols, paste(g, time_points, sep = "_"))
valid_cols <- intersect(ordered_cols, colnames(median_mat))
expr_mat_med <- median_mat[, valid_cols]

get_pathway_gene_ids <- function(target_cluster) {
  sub_df <- refined_go_bp_df %>% filter(Cluster == target_cluster)
  if (nrow(sub_df) == 0) return(character(0))
  all_symbols <- unique(unlist(strsplit(paste(sub_df$geneID, collapse = "/"), "/")))
  all_ids <- std_to_orig_map[trimws(all_symbols)]
  unique(all_ids[!is.na(all_ids)])
}

genes_c3 <- cluster_info %>%
  filter(Gene %in% get_pathway_gene_ids("Cluster_3"), Cluster == "Cluster_3") %>%
  arrange(desc(Global_Correlation))

genes_c4 <- cluster_info %>%
  filter(Gene %in% get_pathway_gene_ids("Cluster_4"), Cluster == "Cluster_4") %>%
  arrange(desc(Global_Correlation)) %>% head(20)

heatmap_genes_df <- bind_rows(genes_c3, genes_c4)

target_pathways <- refined_go_bp_df %>%
  filter(Cluster %in% c("Cluster_3", "Cluster_4")) %>%
  dplyr::select(Description, Cluster) %>% distinct() %>%
  arrange(Cluster, Description) %>%
  mutate(PathwayID = row_number(), Cluster = as.character(Cluster))

pathway_colors <- c(
  RColorBrewer::brewer.pal(3, "Oranges")[2:3],
  RColorBrewer::brewer.pal(5, "PuBuGn")[2:6])
if (length(pathway_colors) < nrow(target_pathways))
  pathway_colors <- RColorBrewer::brewer.pal(nrow(target_pathways), "Set3")
names(pathway_colors) <- target_pathways$PathwayID

link_data <- data.frame()
for (g_id in heatmap_genes_df$Gene) {
  sym <- orig_to_std_map[g_id]; if (is.na(sym)) sym <- g_id
  for (i in seq_len(nrow(target_pathways))) {
    p_desc <- target_pathways$Description[i]
    p_genes <- unlist(strsplit(
      paste(refined_go_bp_df %>% filter(Description == p_desc) %>% pull(geneID), collapse = "/"), "/"))
    if (toupper(sym) %in% toupper(trimws(p_genes))) {
      link_data <- rbind(link_data, data.frame(
        Gene = g_id, PathwayID = target_pathways$PathwayID[i],
        Cluster = target_pathways$Cluster[i]))
    }
  }
}

row_labels_sym <- sapply(heatmap_genes_df$Gene, function(x) {
  s <- orig_to_std_map[x]; if (is.na(s)) s <- x; stringr::str_to_title(s)})

plot_mat <- expr_mat_med[heatmap_genes_df$Gene, ]
plot_mat_z <- t(scale(t(plot_mat)))
split_vec <- factor(heatmap_genes_df$Cluster, levels = c("Cluster_3", "Cluster_4"))

meta_df <- tibble(SampleID = colnames(plot_mat_z)) %>%
  mutate(Group = str_extract(SampleID, "^[CHA]"),
         Time  = as.numeric(str_extract(SampleID, "\\d+$")),
         Treatment = factor(case_when(
           Group == "C" ~ "Control", Group == "H" ~ "H2O2", Group == "A" ~ "Acetic Acid"
         ), levels = c("Control", "H2O2", "Acetic Acid")))

x_labels <- as.character(meta_df$Time)
show_idx <- seq(1, length(x_labels), by = 9)
x_labels[-show_idx] <- ""

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))

top_anno <- HeatmapAnnotation(
  Treatment = meta_df$Treatment,
  col = list(Treatment = c("Control" = "#2E8B57", "H2O2" = "#1E90FF", "Acetic Acid" = "#FF6347")),
  show_annotation_name = TRUE, annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 11, fontfamily = "sans"),
  simple_anno_size = unit(0.4, "cm"))

ribbon_panel_fun <- function(index, levels) {
  current_cluster <- as.character(levels)
  block_genes <- rownames(plot_mat_z)[index]
  relevant_pathways <- target_pathways %>% filter(Cluster == current_cluster)
  if (nrow(relevant_pathways) == 0) return()
  path_y_coords <- seq(0.9, 0.1, length.out = nrow(relevant_pathways))
  names(path_y_coords) <- relevant_pathways$PathwayID
  pushViewport(viewport(xscale = c(0, 1), yscale = c(0, 1)))
  for (i in seq_along(index)) {
    g_id <- block_genes[i]
    g_links <- link_data %>% filter(Gene == g_id)
    if (nrow(g_links) > 0) {
      y_start <- 1 - (i - 0.5) / length(index)
      grid.points(x = 0, y = y_start, pch = 19, size = unit(1.5, "mm"), gp = gpar(col = "grey40"))
      for (r in seq_len(nrow(g_links))) {
        pid <- as.character(g_links$PathwayID[r])
        if (pid %in% names(path_y_coords)) {
          y_end <- path_y_coords[pid]
          l_col <- pathway_colors[pid]
          grid.bezier(x = c(0, 0.3, 0.3, 0.6), y = c(y_start, y_start, y_end, y_end),
                      gp = gpar(col = paste0(l_col, "99"), lwd = 1.5))
        }
      }
    }
  }
  for (pid in names(path_y_coords)) {
    py <- path_y_coords[pid]
    desc <- relevant_pathways$Description[relevant_pathways$PathwayID == as.numeric(pid)]
    l_col <- pathway_colors[pid]
    grid.points(x = 0.6, y = py, pch = 19, size = unit(2.5, "mm"), gp = gpar(col = l_col))
    grid.text(label = desc, x = 0.62, y = py, just = "left",
              gp = gpar(fontsize = 11, fontfamily = "sans", col = "black"))
  }
  popViewport()
}

right_complex_anno <- rowAnnotation(
  Symbols = anno_text(row_labels_sym[rownames(plot_mat_z)],
    gp = gpar(fontsize = 11, fontfamily = "sans"), just = "left", location = unit(0, "npc")),
  Ribbon = anno_block(panel_fun = ribbon_panel_fun, width = unit(10, "cm"), gp = gpar(col = NA)),
  gap = unit(2, "mm"))

pdf_height <- max(10, nrow(plot_mat_z) * 0.25 + 2)
pdf("../figure/fig2c.pdf", width = 26, height = pdf_height)
ht <- Heatmap(plot_mat_z,
  name = "Z-score", split = split_vec,
  column_split = meta_df$Treatment, column_gap = unit(3, "mm"),
  cluster_rows = TRUE, show_row_dend = TRUE, row_dend_width = unit(1.5, "cm"),
  cluster_columns = FALSE, show_row_names = FALSE,
  right_annotation = right_complex_anno,
  row_title_rot = 0, row_title_gp = gpar(fontsize = 14, fontfamily = "sans", fontface = "bold"),
  column_labels = x_labels, column_names_gp = gpar(fontsize = 11, fontfamily = "sans"),
  top_annotation = top_anno, col = col_fun, border = TRUE, use_raster = TRUE)
draw(ht, merge_legend = TRUE,
     column_title = "Time (min) - Aggregated Median",
     column_title_gp = gpar(fontsize = 14, fontfamily = "sans", fontface = "bold"),
     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()
# ---- 2c. Fig.2D — WGCNA module-trait heatmap ----

library(WGCNA)

data_all <- read.csv("../output/WGCNA_eigengenes_traits.csv", row.names = 1, check.names = FALSE)
is_ME <- grepl("^ME", names(data_all))
datTraits <- data_all[, !is_ME, drop = FALSE]
MEs_all <- data_all[, is_ME, drop = FALSE]

drop_cols <- intersect(c("Time_squared", "Time^2", "Time 2", "Time2"), colnames(datTraits))
if (length(drop_cols) > 0) datTraits <- datTraits[, setdiff(colnames(datTraits), drop_cols), drop = FALSE]

target_MEs <- paste0("ME", c("turquoise", "red", "pink"))
exist_MEs <- intersect(target_MEs, names(MEs_all))
if (length(exist_MEs) == 0) stop("No selected WGCNA modules found")

MEs_subset <- MEs_all[, exist_MEs, drop = FALSE]
nSamples <- nrow(MEs_subset)
modTraitCor <- cor(MEs_subset, datTraits, use = "p")
modTraitP   <- corPvalueStudent(modTraitCor, nSamples)

yLabels <- exist_MEs
if (file.exists("../output/WGCNA_module_assignments.csv")) {
  mod_assign <- read.csv("../output/WGCNA_module_assignments.csv")
  mod_sizes  <- table(mod_assign$ModuleColor)
  cc <- gsub("ME", "", exist_MEs)
  yLabels <- paste0(exist_MEs, " (", mod_sizes[cc], " proteins)")}

textMatrix <- paste(signif(modTraitCor, 2), "\n",
  ifelse(modTraitP < 0.001, "***", ifelse(modTraitP < 0.01, "**",
    ifelse(modTraitP < 0.05, "*", ""))), sep = "")
dim(textMatrix) <- dim(modTraitCor)

pdf("../figure/fig2d.pdf", width = 10, height = 5)
par(mar = c(6, 12, 3, 3))
labeledHeatmap(Matrix = modTraitCor, xLabels = names(datTraits), yLabels = yLabels,
  colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix,
  setStdMargins = FALSE, cex.text = 1.0, cex.lab.y = 1.0, zlim = c(-1, 1),
  main = "Module-trait relationships (Selected Modules)")
dev.off()
# ---- 2d. Fig.2E + Fig.2F — ME dynamics ----

plot_data_me <- data_all %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(Treatment = factor(case_when(
    is_Control == 1 ~ "Control", is_H2O2 == 1 ~ "H2O2",
    is_Acetic == 1 ~ "Acetic Acid", TRUE ~ "Unknown"),
    levels = c("Control", "H2O2", "Acetic Acid")))

mk <- function(me_col, mod_name, panel) {
  p <- ggplot(plot_data_me, aes_string(x = "Time", y = me_col, color = "Treatment")) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.15, span = 0.8) +
    scale_color_manual(values = c("Control" = "#2E8B57", "H2O2" = "#1E90FF", "Acetic Acid" = "#FF6347")) +
    labs(title = paste0("ME", mod_name, " Module Eigengene Dynamics"),
         x = "Time Points", y = paste0("ME", mod_name, " Expression")) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(paste0("../figure/fig2", panel, ".pdf"), p, width = 16, height = 8, device = "pdf")
}

mk("MEturquoise", "turquoise", "e")
mk("MEpink", "pink", "f")


# ---- 2e. Fig.2G — Hub-protein trajectories ----

raw_df <- read_tsv(F_SMOOTHED, show_col_types = FALSE)
raw_mat <- raw_df %>% column_to_rownames(colnames(raw_df)[1]) %>% as.matrix()

conv_df <- read_tsv(ID_MAPPING, show_col_types = FALSE) %>%
  filter(Status %in% c("SUCCESS", "PARTIAL")) %>%
  dplyr::select(Original_ID, Final_ID) %>% distinct()
sym_to_orig <- setNames(conv_df$Original_ID, conv_df$Final_ID)

target_genes <- bind_rows(
  tibble(Symbol = c("RLI1","PRP43","YHR020W","MDN1","ARB1"), Module = "Turquoise"),
  tibble(Symbol = c("DAD3","SLO1","YPR010C-A","HTL1","YBR085C-A","SMD3","MRPL33","RSM19"), Module = "Pink"))
target_genes$MatrixID <- sapply(target_genes$Symbol, function(s) {
  if (s %in% rownames(raw_mat)) return(s)
  id <- sym_to_orig[s]; if (!is.na(id) && id %in% rownames(raw_mat)) return(id); NA_character_})

valid <- target_genes %>% filter(!is.na(MatrixID))
sub_mat <- raw_mat[valid$MatrixID, ]; sub_mat_z <- t(scale(t(sub_mat)))
rownames(sub_mat_z) <- valid$Symbol[match(rownames(sub_mat_z), valid$MatrixID)]

meta_g <- tibble(SampleID = colnames(sub_mat_z)) %>%
  mutate(Time = as.numeric(str_extract(SampleID, "\\d+$")),
    Treatment = factor(case_when(
      str_extract(SampleID, "^[CHA]") == "C" ~ "Control",
      str_extract(SampleID, "^[CHA]") == "H" ~ "H2O2",
      TRUE ~ "Acetic Acid"), levels = c("Control", "H2O2", "Acetic Acid")))

plot_df <- as.data.frame(sub_mat_z) %>% rownames_to_column("Symbol") %>%
  pivot_longer(-Symbol, names_to = "SampleID", values_to = "Z_Score") %>%
  left_join(meta_g, by = "SampleID") %>%
  left_join(valid %>% dplyr::select(Symbol, Module), by = "Symbol") %>%
  mutate(Symbol = factor(Symbol, levels = valid %>% arrange(desc(Module), Symbol) %>% pull(Symbol)))

p_g <- ggplot(plot_df, aes(x = Time, y = Z_Score)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.3) +
  geom_point(aes(color = Treatment), size = 1, alpha = 0.6) +
  geom_smooth(aes(color = Treatment), method = "loess", span = 0.25, se = FALSE, linewidth = 0.8) +
  facet_grid(Symbol ~ ., scales = "free_y", switch = "y") +
  scale_color_manual(values = c("Control" = "#2E8B57", "H2O2" = "#1E90FF", "Acetic Acid" = "#FF6347")) +
  theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    strip.background = element_blank(), strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 10),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
    panel.spacing.y = unit(0.1, "cm"), legend.position = "top", legend.justification = "right") +
  labs(x = "Time (min)", title = "Temporal Dynamics of Hub Proteins (Z-score)")

ggsave("../figure/fig2g.pdf", p_g, width = 10, height = n_distinct(plot_df$Symbol) * 0.6 + 2, device = "pdf")
# ---- 2f. EFig.2 — Time-course overviews ----

create_ts <- function(input_file, title, filename) {
  df <- read_tsv(input_file, show_col_types = FALSE) %>%
    dplyr::rename(Genes = 1) %>%
    pivot_longer(-Genes, names_to = "SampleID", values_to = "Expr") %>%
    mutate(treatment = str_extract(SampleID, "^[CHA]"),
           time_point = as.numeric(str_extract(SampleID, "[0-9.]+$"))) %>%
    filter(!is.na(treatment), !is.na(Expr))
  p <- ggplot(df %>% sample_frac(1),
    aes(x = time_point, y = Expr, group = interaction(Genes, treatment), color = treatment)) +
    geom_line(alpha = 0.2, linewidth = 0.2) +
    scale_color_manual(name = "Treatment", values = treatment_colors, labels = treatment_names) +
    labs(title = title, x = "Time Points", y = "Protein Expression Level") +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
          axis.title = element_text(size = 18), legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 2.5)))
  ggsave(filename, p, device = cairo_pdf, width = 12, height = 8)
}

create_ts(IMPUTED, "Original Protein Expression Time Course", "../figure/ef2a.pdf")
create_ts(F_SMOOTHED, "Smoothed (Moving Average, window = 6)", "../figure/ef2c.pdf")
create_ts(F_DOWNSAMP, "Downsampled (Median, window = 6)", "../figure/ef2d.pdf")

file.copy("../output/G02_FFT_Analysis.pdf", "../figure/ef2b.pdf", overwrite = TRUE)


# ---- 3a. EFig.3A — maSigPro Acetic cluster portrait ----
# pypdf one-liner to extract page 2
ret <- system2("python", c("-c",
  shQuote("import sys; from pypdf import PdfReader, PdfWriter; src='../output/maSigPro_Final_Results/Clusters_AceticAcidvsControl.pdf'; r=PdfReader(src); w=PdfWriter(); w.add_page(r.pages[1]); w.write('../figure/ef3a.pdf')")
))
if (ret != 0) cat("  [WARN] ef3a extraction failed\n")


# ---- 3b. EFig.3B — GO BP dotplots ----

source_dir <- "../output"

plot_go_bp <- function(source_csv, output_file, title) {
  if (!file.exists(source_csv)) {
    cat(sprintf("  [WARN] %s not found, skipping\n", source_csv)); return()
  }
  df <- readr::read_csv(source_csv, show_col_types = FALSE)
  plot_data <- df %>%
    arrange(p.adjust) %>%
    mutate(Description_wrapped = forcats::fct_inorder(Description_wrapped))
  count_breaks <- unique(round(scales::pretty_breaks(n = 4)(plot_data$Count)))
  x_lim <- max(plot_data$GeneRatio_numeric, na.rm = TRUE) * 1.25
  p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description_wrapped)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_size_continuous(name = "Count", range = c(4, 10), breaks = count_breaks) +
    scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
    scale_x_continuous(labels = scales::percent, limits = c(0, x_lim),
                       expand = expansion(mult = c(0.01, 0.05))) +
    labs(title = title, x = "Gene Ratio", y = NULL) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          axis.text.y = element_text(size = 12, lineheight = 0.9), legend.position = "right")
  ggsave(output_file, p, width = 12, height = max(5, 3 + 0.6 * nrow(plot_data)),
         device = cairo_pdf, limitsize = FALSE)
}

plot_go_bp(file.path(source_dir, "EF3B_GO_BP_Cluster_3_Source.csv"),
           "../figure/ef3b_cluster3.pdf", "GO Enrichment: Biological Process - Cluster 3")
plot_go_bp(file.path(source_dir, "EF3B_GO_BP_Cluster_4_Source.csv"),
           "../figure/ef3b_cluster4.pdf", "GO Enrichment: Biological Process - Cluster 4")


# ---- 3d. EFig.3D — WGCNA dendrogram ----
if (file.exists("../output/Fig4_M02_Dendrogram_and_Modules.pdf")) {
  file.copy("../output/Fig4_M02_Dendrogram_and_Modules.pdf", "../figure/ef3d.pdf", overwrite = TRUE)
}


# ---- 3e. EFig.3E — WGCNA module-trait heatmap ----
if (file.exists("../output/Fig4_M05_Module_Trait_Heatmap_Advanced.pdf")) {
  file.copy("../output/Fig4_M05_Module_Trait_Heatmap_Advanced.pdf", "../figure/ef3e.pdf", overwrite = TRUE)
}


# ---- 3f. EFig.3 source map ----
ef3_source_map <- tibble::tribble(
  ~Panel, ~File, ~Source,
  "ef3a", "ef3a.pdf", "../output/maSigPro_Final_Results/Clusters_AceticAcidvsControl.pdf page 2",
  "ef3b_1", "ef3b_cluster3.pdf", "../output/EF3B_GO_BP_Cluster_3_Source.csv",
  "ef3b_2", "ef3b_cluster4.pdf", "../output/EF3B_GO_BP_Cluster_4_Source.csv",
  "ef3c", "ef3c.pdf", "../output/ef3c_maSigPro_SourceData.csv + ../output/ef3c_WGCNA_SourceData.csv",
  "ef3d", "ef3d.pdf", "../output/Fig4_M02_Dendrogram_and_Modules.pdf",
  "ef3e", "ef3e.pdf", "../output/Fig4_M05_Module_Trait_Heatmap_Advanced.pdf"
)
readr::write_tsv(ef3_source_map, "../output/ef3_source_map.tsv")

cat("Step02 — Done.\n")
