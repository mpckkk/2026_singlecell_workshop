#!/usr/bin/env Rscript
# Phase 3 — Bulk RNA-seq (DESeq2)
set.seed(42)
ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

inp <- readRDS(file.path(ROOT, "objects", "phase0_input.rds"))
bulk <- inp$bulk
shared <- inp$shared_genes
bulk <- bulk[shared, , drop = FALSE]

# Sample metadata: RP = protoplasted, RU = unprotoplasted (Denyer et al. design)
sample_ids <- colnames(bulk)
condition <- ifelse(grepl("^RP", sample_ids), "protoplasted", "unprotoplasted")
coldata <- data.frame(
  sample = sample_ids,
  condition = factor(condition, levels = c("unprotoplasted", "protoplasted")),
  row.names = sample_ids
)

dds <- DESeqDataSetFromMatrix(countData = round(bulk), colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_counts), file.path(ROOT, "tables", "bulk_normalized_counts.csv"))

# DE: protoplasted vs unprotoplasted
res <- results(dds, contrast = c("condition", "protoplasted", "unprotoplasted"))
res_df <- as.data.frame(res[order(res$padj, na.last = TRUE), ])
res_df$gene <- rownames(res_df)
write.csv(res_df, file.path(ROOT, "tables", "bulk_DE_protoplasted_vs_unprotoplasted.csv"), row.names = FALSE)

# PCA
pca <- prcomp(t(assay(vsd)), scale. = FALSE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], sample = rownames(pca$x))
pca_df <- merge(pca_df, coldata, by.x = "sample", by.y = "row.names")
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 4) + geom_text(vjust = -0.8, size = 3) +
  ggtitle("Bulk RNA-seq PCA (VST)") + theme_bw()
ggsave(file.path(ROOT, "figures", "bulk_pca.png"), p_pca, width = 7, height = 5, dpi = 150)
ggsave(file.path(ROOT, "figures", "bulk_pca.pdf"), p_pca, width = 7, height = 5)

# Sample distance heatmap
sample_dist <- dist(t(assay(vsd)))
mat <- as.matrix(sample_dist)
png(file.path(ROOT, "figures", "bulk_sample_distance_heatmap.png"), width = 700, height = 600, res = 150)
pheatmap(mat, main = "Bulk sample distance (VST)")
dev.off()
pdf(file.path(ROOT, "figures", "bulk_sample_distance_heatmap.pdf"), width = 7, height = 6)
pheatmap(mat, main = "Bulk sample distance (VST)")
dev.off()

# Top variable genes heatmap
rv <- apply(assay(vsd), 1, var)
top_genes <- names(sort(rv, decreasing = TRUE))[1:50]
mat_hm <- assay(vsd)[top_genes, ]
mat_hm <- t(scale(t(mat_hm)))
png(file.path(ROOT, "figures", "bulk_top_variable_genes_heatmap.png"), width = 800, height = 900, res = 150)
pheatmap(mat_hm, show_rownames = TRUE, fontsize_row = 6, main = "Top 50 variable genes (VST z-score)")
dev.off()
pdf(file.path(ROOT, "figures", "bulk_top_variable_genes_heatmap.pdf"), width = 8, height = 9)
pheatmap(mat_hm, show_rownames = TRUE, fontsize_row = 6, main = "Top 50 variable genes (VST z-score)")
dev.off()

# Volcano
res_df$sig <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "sig", "ns")
p_vol <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c(sig = "red", ns = "grey70")) +
  ggtitle("Bulk DE: protoplasted vs unprotoplasted") + theme_bw()
ggsave(file.path(ROOT, "figures", "bulk_volcano.png"), p_vol, width = 7, height = 5, dpi = 150)
ggsave(file.path(ROOT, "figures", "bulk_volcano.pdf"), p_vol, width = 7, height = 5)

# DEG heatmap
deg <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
if (nrow(deg) > 0) {
  deg_genes <- head(deg$gene[order(deg$padj)], min(50, nrow(deg)))
  mat_deg <- assay(vsd)[deg_genes, ]
  mat_deg <- t(scale(t(mat_deg)))
  png(file.path(ROOT, "figures", "bulk_DEG_heatmap.png"), width = 800, height = 900, res = 150)
  pheatmap(mat_deg, show_rownames = TRUE, fontsize_row = 6, main = "Significant DEGs (VST z-score)")
  dev.off()
}

saveRDS(list(dds = dds, vsd = vsd, coldata = coldata), file.path(ROOT, "objects", "bulk_deseq.rds"))
message("Phase 3 complete. DEGs (padj<0.05, |LFC|>1): ", nrow(deg))
