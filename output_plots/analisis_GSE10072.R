# ============================================================
# ANALISIS DIFFERENTIALLY EXPRESSED GENES (DEG)
# Dataset: GSE10072 - Lung Adenocarcinoma vs Normal
# Tools: GEOquery, limma, hgu133a.db, clusterProfiler, ggplot2
# ============================================================

# ============================================================
# BAGIAN 0: INSTALASI PACKAGES
# Jalankan bagian ini SATU KALI jika belum punya package-nya
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery",
  "limma",
  "hgu133a.db",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "AnnotationDbi"
), ask = FALSE)

install.packages(c("ggplot2", "pheatmap", "ggrepel", "dplyr", "RColorBrewer"))


# ============================================================
# BAGIAN 1: LOAD LIBRARY
# ============================================================

library(GEOquery)
library(limma)
library(hgu133a.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(RColorBrewer)

# Buat folder output untuk menyimpan plot
dir.create("output_plots", showWarnings = FALSE)

cat("✅ Library berhasil dimuat\n")


# ============================================================
# BAGIAN 2: DOWNLOAD & LOAD DATASET GSE10072
# ============================================================

cat("📥 Mengunduh dataset GSE10072...\n")

gse <- getGEO("GSE10072", GSEMatrix = TRUE, getGPL = FALSE)
gse_data <- gse[[1]]

cat("✅ Dataset berhasil diunduh\n")
cat("📊 Dimensi data:", dim(exprs(gse_data)), "\n")

# Lihat metadata sampel
pdata <- pData(gse_data)
cat("📋 Kolom metadata:\n")
print(colnames(pdata))

# Lihat informasi karakteristik sampel
cat("\n📋 Karakteristik sampel (beberapa baris pertama):\n")
print(head(pdata[, c("title", "characteristics_ch1")]))


# ============================================================
# BAGIAN 3: PREPROCESSING
# ============================================================

# Ambil matrix ekspresi
expr_matrix <- exprs(gse_data)

# --- 3.1 Cek apakah data perlu transformasi log2 ---
max_val <- max(expr_matrix, na.rm = TRUE)
cat("\n🔍 Nilai maksimum ekspresi:", max_val, "\n")

if (max_val > 100) {
  cat("⚙️  Melakukan transformasi log2...\n")
  expr_matrix <- log2(expr_matrix + 1)
} else {
  cat("✅ Data sudah dalam skala log\n")
}

# --- 3.2 Boxplot QC sebelum normalisasi (opsional, bisa di-skip) ---
cat("📊 Membuat boxplot QC...\n")
png("output_plots/00_boxplot_QC.png", width = 1200, height = 600)
boxplot(expr_matrix,
        las = 2,
        cex.axis = 0.5,
        main = "Distribusi Ekspresi Setelah Log2",
        ylab = "log2 Intensitas",
        col = rainbow(ncol(expr_matrix)))
dev.off()

cat("✅ Preprocessing selesai\n")


# ============================================================
# BAGIAN 4: DEFINISI GRUP SAMPEL
# ============================================================

# Cek karakteristik untuk menentukan grup
cat("\n📋 Unik nilai karakteristik sampel:\n")
print(table(pdata$characteristics_ch1))

# Tentukan grup berdasarkan karakteristik
# GSE10072: "tissue: adenocarcinoma" vs "tissue: normal"
group_labels <- ifelse(
  grepl("adenocarcinoma|tumor|cancer", pdata$characteristics_ch1, ignore.case = TRUE),
  "Tumor",
  "Normal"
)

# Verifikasi
cat("\n📊 Distribusi grup:\n")
print(table(group_labels))

# Buat factor
group_factor <- factor(group_labels, levels = c("Normal", "Tumor"))


# ============================================================
# BAGIAN 5: ANALISIS DEG DENGAN LIMMA
# ============================================================

cat("\n⚙️  Menjalankan analisis limma...\n")

# Buat design matrix
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor)

cat("Design matrix:\n")
print(head(design))

# Fit model linear
fit <- lmFit(expr_matrix, design)

# Buat contrast: Tumor vs Normal
contrast_matrix <- makeContrasts(
  Tumor_vs_Normal = Tumor - Normal,
  levels = design
)

# Hitung contrast dan empirical Bayes
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

cat("✅ Model limma selesai difit\n")

# Ambil tabel hasil DEG
deg_table <- topTable(fit2,
                      coef = "Tumor_vs_Normal",
                      number = Inf,
                      adjust.method = "BH",
                      sort.by = "P")

cat("📊 Jumlah total probe:", nrow(deg_table), "\n")

# --- 5.1 Anotasi gen (mapping probe ID ke gene symbol) ---
cat("🔬 Melakukan anotasi gen...\n")

probe_ids <- rownames(deg_table)

gene_symbols <- mapIds(hgu133a.db,
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

entrez_ids <- mapIds(hgu133a.db,
                     keys = probe_ids,
                     column = "ENTREZID",
                     keytype = "PROBEID",
                     multiVals = "first")

deg_table$Gene_Symbol <- gene_symbols
deg_table$Entrez_ID   <- entrez_ids

# Hapus probe tanpa gene symbol
deg_annotated <- deg_table[!is.na(deg_table$Gene_Symbol), ]
deg_annotated <- deg_annotated[!duplicated(deg_annotated$Gene_Symbol), ]
rownames(deg_annotated) <- deg_annotated$Gene_Symbol

cat("✅ Anotasi selesai\n")
cat("📊 Jumlah gen ternotasi:", nrow(deg_annotated), "\n")

# --- 5.2 Filter DEG signifikan ---
sig_degs <- deg_annotated %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1)

cat("\n🎯 DEG Signifikan (adj.P.Val < 0.05 & |logFC| > 1):", nrow(sig_degs), "\n")
cat("   ↑ Upregulated  :", sum(sig_degs$logFC > 0), "\n")
cat("   ↓ Downregulated:", sum(sig_degs$logFC < 0), "\n")

# Simpan hasil DEG ke CSV
write.csv(deg_annotated, "output_plots/DEG_results_all.csv", row.names = FALSE)
write.csv(sig_degs,      "output_plots/DEG_results_significant.csv", row.names = FALSE)
cat("💾 Hasil DEG disimpan ke output_plots/\n")


# ============================================================
# BAGIAN 6: VOLCANO PLOT
# ============================================================

cat("\n🌋 Membuat Volcano Plot...\n")

# Siapkan data untuk plot
volcano_data <- deg_annotated %>%
  mutate(
    neg_log10_pval = -log10(adj.P.Val),
    Regulation = case_when(
      logFC >  1 & adj.P.Val < 0.05 ~ "Upregulated",
      logFC < -1 & adj.P.Val < 0.05 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

# Pilih top gen untuk label (top 15 up + top 15 down)
top_genes <- bind_rows(
  volcano_data %>% filter(Regulation == "Upregulated")   %>% arrange(adj.P.Val) %>% head(15),
  volcano_data %>% filter(Regulation == "Downregulated") %>% arrange(adj.P.Val) %>% head(15)
)

# Buat plot
volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_pval, color = Regulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated"    = "#E74C3C",
    "Downregulated"  = "#3498DB",
    "Not Significant" = "#95A5A6"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text_repel(
    data = top_genes,
    aes(label = Gene_Symbol),
    size = 2.8,
    max.overlaps = 30,
    color = "black",
    box.padding = 0.3
  ) +
  labs(
    title    = "Volcano Plot: Tumor vs Normal (GSE10072)",
    subtitle = "Lung Adenocarcinoma — Differentially Expressed Genes",
    x        = "log2 Fold Change",
    y        = "-log10 (Adjusted P-value)",
    color    = "Regulasi"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "right"
  ) +
  annotate("text", x = max(volcano_data$logFC) * 0.8, y = -log10(0.05) + 0.5,
           label = "adj.P.Val = 0.05", size = 3, color = "black") +
  annotate("text", x = 3, y = max(volcano_data$neg_log10_pval, na.rm=TRUE) * 0.95,
           label = paste0("Up: ", sum(volcano_data$Regulation == "Upregulated")),
           color = "#E74C3C", size = 4, fontface = "bold") +
  annotate("text", x = -3, y = max(volcano_data$neg_log10_pval, na.rm=TRUE) * 0.95,
           label = paste0("Down: ", sum(volcano_data$Regulation == "Downregulated")),
           color = "#3498DB", size = 4, fontface = "bold")

ggsave("output_plots/01_volcano_plot.png",
       plot   = volcano_plot,
       width  = 10,
       height = 7,
       dpi    = 300)

print(volcano_plot)
cat("✅ Volcano Plot disimpan\n")


# ============================================================
# BAGIAN 7: HEATMAP 50 DEG TERATAS
# ============================================================

cat("\n🗺️  Membuat Heatmap 50 DEG Teratas...\n")

# Ambil 50 DEG teratas (berdasarkan adj.P.Val terkecil)
top50_genes <- deg_annotated %>%
  filter(!is.na(Gene_Symbol)) %>%
  arrange(adj.P.Val) %>%
  head(50)

# Ambil data ekspresi untuk gen tersebut
top50_expr <- expr_matrix[rownames(expr_matrix) %in% rownames(top50_genes), ]

# Jika rownames masih probe ID, ganti dengan gene symbol
matched_idx <- match(rownames(top50_expr), rownames(top50_genes))
rownames(top50_expr) <- top50_genes$Gene_Symbol[matched_idx]

# Buat annotation kolom berdasarkan grup
annotation_col <- data.frame(
  Group = group_factor,
  row.names = colnames(top50_expr)
)

# Warna annotation
ann_colors <- list(
  Group = c(Normal = "#2ECC71", Tumor = "#E74C3C")
)

# Buat heatmap
png("output_plots/02_heatmap_top50.png", width = 1400, height = 1000, res = 150)

pheatmap(
  top50_expr,
  annotation_col  = annotation_col,
  annotation_colors = ann_colors,
  scale           = "row",           # Z-score normalisasi per gen
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color           = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
  show_colnames   = FALSE,
  fontsize_row    = 8,
  main            = "Heatmap 50 DEG Teratas\nGSE10072: Lung Adenocarcinoma vs Normal",
  border_color    = NA,
  treeheight_row  = 50,
  treeheight_col  = 30
)

dev.off()
cat("✅ Heatmap disimpan\n")


# ============================================================
# BAGIAN 8: ENRICHMENT ANALYSIS — GENE ONTOLOGY (GO)
# ============================================================

cat("\n🧬 Menjalankan GO Enrichment Analysis...\n")

# Ambil Entrez ID dari DEG signifikan
sig_entrez <- sig_degs$Entrez_ID
sig_entrez <- sig_entrez[!is.na(sig_entrez)]
sig_entrez <- unique(sig_entrez)

cat("🔢 Jumlah gen untuk enrichment:", length(sig_entrez), "\n")

# Universe: semua gen ternotasi
universe_entrez <- deg_annotated$Entrez_ID
universe_entrez <- universe_entrez[!is.na(universe_entrez)]
universe_entrez <- unique(universe_entrez)

# --- GO Biological Process ---
go_bp <- enrichGO(
  gene          = sig_entrez,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

cat("✅ GO BP: ditemukan", nrow(go_bp@result[go_bp@result$p.adjust < 0.05,]), "term signifikan\n")

# --- GO Molecular Function ---
go_mf <- enrichGO(
  gene          = sig_entrez,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# --- GO Cellular Component ---
go_cc <- enrichGO(
  gene          = sig_entrez,
  universe      = universe_entrez,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# --- Plot GO BP: Dotplot ---
if (nrow(go_bp) > 0) {
  go_dotplot <- dotplot(go_bp,
                        showCategory = 20,
                        title = "GO Biological Process — Enrichment\nGSE10072: Lung Adenocarcinoma vs Normal") +
    theme(
      axis.text.y   = element_text(size = 9),
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 12)
    )

  ggsave("output_plots/03_GO_BP_dotplot.png",
         plot   = go_dotplot,
         width  = 10,
         height = 9,
         dpi    = 300)
  cat("✅ GO BP dotplot disimpan\n")
}

# --- Plot GO BP: Barplot ---
if (nrow(go_bp) > 0) {
  png("output_plots/03b_GO_BP_barplot.png", width = 1200, height = 900, res = 150)
  print(barplot(go_bp, showCategory = 15,
                title = "GO Biological Process — Top 15 Term\nGSE10072"))
  dev.off()
  cat("✅ GO BP barplot disimpan\n")
}

# --- Plot GO: Emapplot (network enrichment) ---
if (nrow(go_bp) > 5) {
  go_bp2 <- pairwise_termsim(go_bp)
  emap <- emapplot(go_bp2,
                   showCategory = 20,
                   color = "p.adjust") +
    labs(title = "GO Enrichment Network — Biological Process") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave("output_plots/03c_GO_BP_network.png",
         plot   = emap,
         width  = 12,
         height = 10,
         dpi    = 300)
  cat("✅ GO network plot disimpan\n")
}

# Simpan hasil GO ke CSV
write.csv(as.data.frame(go_bp), "output_plots/GO_BP_results.csv", row.names = FALSE)
write.csv(as.data.frame(go_mf), "output_plots/GO_MF_results.csv", row.names = FALSE)
write.csv(as.data.frame(go_cc), "output_plots/GO_CC_results.csv", row.names = FALSE)


# ============================================================
# BAGIAN 9: ENRICHMENT ANALYSIS — KEGG PATHWAY
# ============================================================

cat("\n🗺️  Menjalankan KEGG Pathway Enrichment...\n")

kegg_result <- enrichKEGG(
  gene          = sig_entrez,
  organism      = "hsa",          # Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

cat("✅ KEGG: ditemukan", nrow(kegg_result@result[kegg_result@result$p.adjust < 0.05,]), "pathway signifikan\n")

# Readable gene names
kegg_readable <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# --- KEGG Dotplot ---
if (nrow(kegg_result) > 0) {
  kegg_dot <- dotplot(kegg_readable,
                      showCategory = 20,
                      title = "KEGG Pathway Enrichment\nGSE10072: Lung Adenocarcinoma vs Normal") +
    theme(
      axis.text.y = element_text(size = 9),
      plot.title  = element_text(face = "bold", hjust = 0.5, size = 12)
    )

  ggsave("output_plots/04_KEGG_dotplot.png",
         plot   = kegg_dot,
         width  = 10,
         height = 8,
         dpi    = 300)
  cat("✅ KEGG dotplot disimpan\n")
}

# --- KEGG Barplot ---
if (nrow(kegg_result) > 0) {
  png("output_plots/04b_KEGG_barplot.png", width = 1200, height = 800, res = 150)
  print(barplot(kegg_readable, showCategory = 15,
                title = "KEGG Pathway — Top 15\nGSE10072"))
  dev.off()
  cat("✅ KEGG barplot disimpan\n")
}

# --- KEGG Cnetplot (gen-pathway network) ---
if (nrow(kegg_result) > 0) {
  cnet <- cnetplot(kegg_readable,
                   showCategory = 10,
                   colorEdge    = TRUE) +
    labs(title = "KEGG — Gene-Pathway Network") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave("output_plots/04c_KEGG_cnetplot.png",
         plot   = cnet,
         width  = 12,
         height = 10,
         dpi    = 300)
  cat("✅ KEGG cnetplot disimpan\n")
}

# Simpan hasil KEGG
write.csv(as.data.frame(kegg_result), "output_plots/KEGG_results.csv", row.names = FALSE)


# ============================================================
# BAGIAN 10: RINGKASAN HASIL
# ============================================================

cat("\n")
cat("=" , rep("=", 50), "\n", sep = "")
cat("✅ ANALISIS SELESAI — RINGKASAN HASIL\n")
cat("=" , rep("=", 50), "\n", sep = "")
cat("\n📊 DEG Summary:\n")
cat("   Total gen ternotasi    :", nrow(deg_annotated), "\n")
cat("   DEG signifikan (total) :", nrow(sig_degs), "\n")
cat("   ↑ Upregulated          :", sum(sig_degs$logFC > 0), "\n")
cat("   ↓ Downregulated        :", sum(sig_degs$logFC < 0), "\n")

cat("\n📁 File output tersimpan di folder: output_plots/\n")
cat("   - 00_boxplot_QC.png\n")
cat("   - 01_volcano_plot.png\n")
cat("   - 02_heatmap_top50.png\n")
cat("   - 03_GO_BP_dotplot.png\n")
cat("   - 03b_GO_BP_barplot.png\n")
cat("   - 03c_GO_BP_network.png\n")
cat("   - 04_KEGG_dotplot.png\n")
cat("   - 04b_KEGG_barplot.png\n")
cat("   - 04c_KEGG_cnetplot.png\n")
cat("   - DEG_results_all.csv\n")
cat("   - DEG_results_significant.csv\n")
cat("   - GO_BP_results.csv\n")
cat("   - KEGG_results.csv\n")
cat("\n🎉 Semua plot dan data siap untuk laporan!\n")
