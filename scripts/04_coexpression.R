#!/usr/bin/env Rscript
# Phase 4 — Cross-dataset co-expression
set.seed(42)
ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(SummarizedExperiment)
})

inp <- readRDS(file.path(ROOT, "objects", "phase0_input.rds"))
seu <- readRDS(file.path(ROOT, "objects", "seurat_wt_annotated.rds"))
bulk_obj <- readRDS(file.path(ROOT, "objects", "bulk_deseq.rds"))

shared <- Reduce(intersect, list(inp$shared_genes, rownames(bulk_obj$vsd), rownames(seu)))
counts <- GetAssayData(seu, layer = "counts")
bulk_vst <- SummarizedExperiment::assay(bulk_obj$vsd)[shared, , drop = FALSE]

# Pseudobulk single-cell: sum counts per cell type
meta <- seu@meta.data
cell_types <- sort(unique(meta$cell_type))

pseudobulk_list <- list()
for (ct in cell_types) {
  cells <- rownames(meta)[meta$cell_type == ct]
  pseudobulk_list[[ct]] <- Matrix::rowSums(counts[, cells, drop = FALSE])
}
pseudobulk_list[["overall"]] <- Matrix::rowSums(counts)
pb <- do.call(cbind, pseudobulk_list)
pb <- pb[shared, , drop = FALSE]

# Log-normalize pseudobulk (log2 CPM-like)
pb_norm <- apply(pb, 2, function(x) log2(x / sum(x) * 1e6 + 1))
bulk_norm <- apply(bulk_vst, 2, function(x) x)  # already VST

# Overall bulk mean profile for correlation
bulk_mean <- rowMeans(bulk_norm)

cor_results <- data.frame(
  comparison = character(),
  pearson = numeric(),
  spearman = numeric(),
  stringsAsFactors = FALSE
)

cor_results <- rbind(cor_results, data.frame(
  comparison = "overall_pseudobulk_vs_bulk_mean",
  pearson = cor(pb_norm[, "overall"], bulk_mean, method = "pearson"),
  spearman = cor(pb_norm[, "overall"], bulk_mean, method = "spearman", use = "complete.obs")
))

for (ct in cell_types) {
  ct_cor <- data.frame(
    comparison = paste0(ct, "_pseudobulk_vs_bulk_mean"),
    pearson = cor(pb_norm[, ct], bulk_mean, method = "pearson"),
    spearman = cor(pb_norm[, ct], bulk_mean, method = "spearman", use = "complete.obs")
  )
  cor_results <- rbind(cor_results, ct_cor)
}
write.csv(cor_results, file.path(ROOT, "tables", "pseudobulk_bulk_correlations.csv"), row.names = FALSE)

# Scatter overall
scatter_df <- data.frame(sc = pb_norm[, "overall"], bulk = bulk_mean)
p_sc <- ggplot(scatter_df, aes(sc, bulk)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "SC pseudobulk (log2 CPM+1)", y = "Bulk mean VST",
       title = paste0("Overall concordance (Spearman rho = ",
                      round(cor_results$spearman[1], 3), ")")) +
  theme_bw()
ggsave(file.path(ROOT, "figures", "pseudobulk_bulk_scatter_overall.png"), p_sc, width = 6, height = 5, dpi = 150)
ggsave(file.path(ROOT, "figures", "pseudobulk_bulk_scatter_overall.pdf"), p_sc, width = 6, height = 5)

# Per cell type scatters (top 4 by correlation)
ct_cors <- cor_results[grepl("_pseudobulk", cor_results$comparison) &
                         !grepl("overall", cor_results$comparison), ]
ct_cors <- ct_cors[order(-ct_cors$spearman), ]
top_ct <- head(gsub("_pseudobulk_vs_bulk_mean", "", ct_cors$comparison), 4)
for (ct in top_ct) {
  df <- data.frame(sc = pb_norm[, ct], bulk = bulk_mean)
  rho <- cor_results$spearman[cor_results$comparison == paste0(ct, "_pseudobulk_vs_bulk_mean")]
  p <- ggplot(df, aes(sc, bulk)) + geom_point(alpha = 0.15, size = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(title = paste0(ct, " (Spearman = ", round(rho, 3), ")"),
         x = "SC pseudobulk", y = "Bulk mean VST") + theme_bw()
  ggsave(file.path(ROOT, "figures", paste0("pseudobulk_bulk_scatter_", ct, ".png")), p, width = 6, height = 5, dpi = 150)
}

# Gene-gene co-expression modules (simplified WGCNA-style via hclust on top variable genes)
n_mod <- 500
top_var <- names(sort(apply(pb_norm, 1, var), decreasing = TRUE))[1:n_mod]
bulk_top <- bulk_norm[top_var, , drop = FALSE]
sc_top <- pb_norm[top_var, , drop = FALSE]

bulk_cor <- cor(t(bulk_top), method = "spearman")
sc_cor <- cor(t(sc_top), method = "spearman")

bulk_dist <- as.dist(1 - bulk_cor)
sc_dist <- as.dist(1 - sc_cor)
bulk_cl <- cutree(hclust(bulk_dist, method = "average"), k = 6)
sc_cl <- cutree(hclust(sc_dist, method = "average"), k = 6)

module_df <- data.frame(gene = top_var, bulk_module = bulk_cl, sc_module = sc_cl)
write.csv(module_df, file.path(ROOT, "tables", "coexpression_modules.csv"), row.names = FALSE)

# Concordant genes: same module cluster assignment
concordant <- module_df[module_df$bulk_module == module_df$sc_module, ]
concordant_rank <- concordant[order(concordant$bulk_module), ]
write.csv(concordant_rank, file.path(ROOT, "tables", "cross_modality_concordant_genes.csv"), row.names = FALSE)

saveRDS(list(pseudobulk = pb_norm, bulk_vst = bulk_norm), file.path(ROOT, "objects", "phase4_coexpression.rds"))
message("Phase 4 complete. Concordant genes: ", nrow(concordant))
