# Modul: Analisis Ekspresi Gen Respon Doxorubicin
# Dataset: GSE39870 (Doxorubicin vs Vehicle)
# Platform: Microarray (Affymetrix Human Genome U133 Plus 2.0)
# Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG)

# PART A. PERSIAPAN LINGKUNGAN KERJA

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("GEOquery", "limma", "hgu133plus2.db"), 
                     ask = FALSE, update = FALSE)

install.packages(c("pheatmap", "ggplot2", "dplyr"))

library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133plus2.db)
library(AnnotationDbi)

# PART B. PENGAMBILAN DATA

gset <- getGEO("GSE39870", GSEMatrix = TRUE)[[1]]

ex <- exprs(gset)

# PART C. LOG2 TRANSFORMASI

qx <- as.numeric(quantile(ex, c(0,0.25,0.5,0.75,0.99,1), na.rm=TRUE))

LogTransform <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)

if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

# PART D. DEFINISI KELOMPOK SAMPEL

group_info <- pData(gset)$title

group <- ifelse(grepl("doxorubicin", group_info, ignore.case=TRUE),
                "Doxorubicin",
                "Vehicle")

gset$group <- factor(group)

print(levels(gset$group))

# PART E. DESIGN MATRIX

design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)

grup_treatment <- "Doxorubicin"
grup_control   <- "Vehicle"

contrast_formula <- paste(grup_treatment, "-", grup_control)

print(paste("Kontras yang dianalisis:", contrast_formula))

# PART F. ANALISIS LIMMA

fit <- lmFit(ex, design)

contrast_matrix <- makeContrasts(contrasts = contrast_formula,
                                 levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)

fit2 <- eBayes(fit2)

topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)

head(topTableResults)

# PART G. ANOTASI NAMA GEN (Affymetrix)

probe_ids <- rownames(topTableResults)

gene_annotation <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

# PART H. TOP 5 UP-REGULATED & DOWN-REGULATED GENES

# Up-regulated
up_genes <- topTableResults %>%
  filter(adj.P.Val < 0.01, logFC > 1) %>%
  arrange(desc(logFC)) %>%
  head(5)

up_genes[, c("GENENAME", "SYMBOL", "logFC", "adj.P.Val")]

#Down-regulated
down_genes <- topTableResults %>%
  filter(adj.P.Val < 0.01, logFC < -1) %>%
  arrange(logFC) %>%
  head(5)

down_genes[, c("GENENAME", "SYMBOL", "logFC", "adj.P.Val")]

# PART I. VOLCANO PLOT

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & 
                      volcano_data$adj.P.Val < 0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & 
                      volcano_data$adj.P.Val < 0.01] <- "DOWN"

ggplot(volcano_data, 
       aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN"="blue",
                                "NO"="grey",
                                "UP"="red")) +
  geom_vline(xintercept = c(-1,1), linetype="dashed") +
  geom_hline(yintercept = -log10(0.01), linetype="dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG Doxorubicin vs Vehicle")

# PART J. HEATMAP TOP 50

topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)

mat_heatmap <- ex[top50$PROBEID, ]

gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,
  top50$SYMBOL
)

rownames(mat_heatmap) <- gene_label

mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

annotation_col <- data.frame(
  Group = gset$group
)

rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

# PART K. BOXPLOT DISTRIBUSI NILAI EKSPRESI

group_colors <- as.numeric(gset$group)
boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)
legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  
  cex = 0.8
)

# PART L. DENSITY PLOT

library(ggplot2)

expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )

# PART M. UMAP

library(umap)  


umap_input <- t(ex)
umap_input <- scale(umap_input)

umap_result <- umap(umap_input,
                    n_neighbors = 5,
                    min_dist = 0.3,
                    metric = "euclidean",
                    random_state = 42)
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)
#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot: Doxorubicin vs Vehicle",
    x = "UMAP 1",
    y = "UMAP 2"
  )

# PART N. ANALISIS ENRICHMENT - GENE ONTOLOGY (GO)

# Install jika belum ada
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"),
                     ask = FALSE, update = FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)

# Ambil gen signifikan saja (UP + DOWN)
deg_sig <- topTableResults[
  topTableResults$adj.P.Val < 0.01 &
    abs(topTableResults$logFC) > 1,
]

up_genes <- subset(topTableResults, 
                   logFC > 1 & adj.P.Val < 0.01)

down_genes <- subset(topTableResults, 
                     logFC < -1 & adj.P.Val < 0.01)

cat("Number of upregulated genes:", nrow(up_genes), "\n")
cat("Number of downregulated genes:", nrow(down_genes), "\n")

# Ambil SYMBOL dan buang NA
gene_list <- deg_sig$SYMBOL
gene_list <- gene_list[!is.na(gene_list)]

# GO enrichment (Biological Process)
ego <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# Visualisasi GO
dotplot(ego, showCategory = 15) +
  ggtitle("GO Biological Process Enrichment")

# Top 5 Biological Processes
ego_result <- as.data.frame(ego)

top_go <- ego_result %>%
  arrange(p.adjust) %>%
  head(5)

top_go[, c("Description", "p.adjust", "GeneRatio")]

# PART O. ANALISIS ENRICHMENT - KEGG PATHWAY

# Konversi SYMBOL ke ENTREZ ID
gene_entrez <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

# KEGG enrichment
kegg <- enrichKEGG(
  gene         = gene_entrez$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

# Visualisasi KEGG
dotplot(kegg, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment")

# Top 5 Pathway
kegg_result <- as.data.frame(kegg)

top_kegg <- kegg_result %>%
  arrange(p.adjust) %>%
  head(5)

top_kegg[, c("Description", "p.adjust", "GeneRatio")]

# PART P. MENYIMPAN HASIL

write.csv(topTableResults, "Hasil_GSE39870_DEG.csv")

message("Analisis selesai. File hasil telah disimpan.")






