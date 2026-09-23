# ============================================================
# GSE297652 Single-cell RNA-seq Analysis
# Prostate Cancer
# Author: Lopamudra Basu
# ============================================================

# ------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------
library(Seurat)
library(Matrix)
library(data.table)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)

base_path <- "D:/Bioinfomatics_project/SingleCell-RNAseq-Prostate-Cancer/data"

# ------------------------------------------------------------
# Read 10X expression matrices for all five samples
# ------------------------------------------------------------

m1_counts <- Read10X(
  data.dir = file.path(base_path, "mPCa_M1_filtered_feature_bc_matrix")
)

m2_counts <- Read10X(
  data.dir = file.path(base_path, "mPCa_M2_filtered_feature_bc_matrix")
)

m3_counts <- Read10X(
  data.dir = file.path(base_path, "mPCa_M3_filtered_feature_bc_matrix")
)

m4_counts <- Read10X(
  data.dir = file.path(base_path, "mPCa_M4_filtered_feature_bc_matrix")
)

m5_counts <- Read10X(
  data.dir = file.path(base_path, "mPCa_M5_filtered_feature_bc_matrix")
)
# Check matrix dimensions
dim(m1_counts)
dim(m2_counts)
dim(m3_counts)
dim(m4_counts)
dim(m5_counts)


# ------------------------------------------------------------
# Read cell metadata
# ------------------------------------------------------------


metadata_file <- "D:/Bioinfomatics_project/SingleCell-RNAseq-Prostate-Cancer/data/GSE297652_cell_metadata.csv.gz"

metadata <- fread(metadata_file)

dim(metadata)

head(metadata)

# ------------------------------------------------------------
# Extract sample-specific cell barcodes
# ------------------------------------------------------------

m1_barcodes <- sub("^M1_", "", metadata[Sample == "M1", Barcode])
m2_barcodes <- sub("^M2_", "", metadata[Sample == "M2", Barcode])
m3_barcodes <- sub("^M3_", "", metadata[Sample == "M3", Barcode])
m4_barcodes <- sub("^M4_", "", metadata[Sample == "M4", Barcode])
m5_barcodes <- sub("^M5_", "", metadata[Sample == "M5", Barcode])

# ------------------------------------------------------------
# Subset expression matrices to cells present in metadata
# ------------------------------------------------------------
m1_counts <- m1_counts[, m1_barcodes]
m2_counts <- m2_counts[, m2_barcodes]
m3_counts <- m3_counts[, m3_barcodes]
m4_counts <- m4_counts[, m4_barcodes]
m5_counts <- m5_counts[, m5_barcodes]

# Check dimensions after subsetting
dim(m1_counts)
dim(m2_counts)
dim(m3_counts)
dim(m4_counts)
dim(m5_counts)


# Check that the matrices are sparse
class(m1_counts)
class(m2_counts)
class(m3_counts)
class(m4_counts)
class(m5_counts)

# Check total counts in each sample
sum(m1_counts)
sum(m2_counts)
sum(m3_counts)
sum(m4_counts)
sum(m5_counts)

# Check whether counts are integers
all(m1_counts@x == floor(m1_counts@x))
all(m2_counts@x == floor(m2_counts@x))
all(m3_counts@x == floor(m3_counts@x))
all(m4_counts@x == floor(m4_counts@x))
all(m5_counts@x == floor(m5_counts@x))

# -----------------------------------------
# Prepare metadata for each sample
# -----------------------------------------

meta_M1 <- metadata[Sample == "M1"]
meta_M2 <- metadata[Sample == "M2"]
meta_M3 <- metadata[Sample == "M3"]
meta_M4 <- metadata[Sample == "M4"]
meta_M5 <- metadata[Sample == "M5"]

# Set cell barcodes as row names
rownames(meta_M1) <- meta_M1$Barcode
rownames(meta_M2) <- meta_M2$Barcode
rownames(meta_M3) <- meta_M3$Barcode
rownames(meta_M4) <- meta_M4$Barcode
rownames(meta_M5) <- meta_M5$Barcode

colnames(m1_counts) <- paste0("M1_", colnames(m1_counts))
colnames(m2_counts) <- paste0("M2_", colnames(m2_counts))
colnames(m3_counts) <- paste0("M3_", colnames(m3_counts))
colnames(m4_counts) <- paste0("M4_", colnames(m4_counts))
colnames(m5_counts) <- paste0("M5_", colnames(m5_counts))

all(rownames(meta_M1) == colnames(m1_counts))
all(rownames(meta_M2) == colnames(m2_counts))
all(rownames(meta_M3) == colnames(m3_counts))
all(rownames(meta_M4) == colnames(m4_counts))
all(rownames(meta_M5) == colnames(m5_counts))


# =========================================
# Create Seurat objects for each sample
# =========================================

seurat_M1 <- CreateSeuratObject(
  counts = m1_counts,
  meta.data = meta_M1,
  project = "GSE297652_M1"
)

seurat_M2 <- CreateSeuratObject(
  counts = m2_counts,
  meta.data = meta_M2,
  project = "GSE297652_M2"
)

seurat_M3 <- CreateSeuratObject(
  counts = m3_counts,
  meta.data = meta_M3,
  project = "GSE297652_M3"
)

seurat_M4 <- CreateSeuratObject(
  counts = m4_counts,
  meta.data = meta_M4,
  project = "GSE297652_M4"
)

seurat_M5 <- CreateSeuratObject(
  counts = m5_counts,
  meta.data = meta_M5,
  project = "GSE297652_M5"
)

dim(seurat_M1)
dim(seurat_M2)
dim(seurat_M3)
dim(seurat_M4)
dim(seurat_M5)

head(seurat_M1@meta.data)

table(seurat_M1$Sample)
table(seurat_M2$Sample)
table(seurat_M3$Sample)
table(seurat_M4$Sample)
table(seurat_M5$Sample)

# =========================================
# Merge all five samples
# =========================================

seurat_combined <- merge(
  x = seurat_M1,
  y = list(seurat_M2, seurat_M3, seurat_M4, seurat_M5),
  project = "GSE297652_mPCa"
)

dim(seurat_combined)
table(seurat_combined$Sample)
seurat_combined

# =========================================
# QC verification
# =========================================

# Recalculate RNA counts from the raw count layers
seurat_combined <- JoinLayers(seurat_combined)

# Calculate QC directly from raw counts
seurat_combined[["percent.mt"]] <- PercentageFeatureSet(
  seurat_combined,
  pattern = "^MT-"
)

head(
  seurat_combined@meta.data[
    ,
    c("nCount_RNA", "nFeature_RNA", "mitoPercent", "percent.mt")
  ]
)

summary(seurat_combined$nCount_RNA)
summary(seurat_combined$nFeature_RNA)
summary(seurat_combined$mitoPercent)
summary(seurat_combined$percent.mt)

# QC plots
VlnPlot(
  seurat_combined,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

FeatureScatter(
  seurat_combined,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

FeatureScatter(
  seurat_combined,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)


# =========================================
# QC percentile analysis
# =========================================

quantile(
  seurat_combined$nCount_RNA,
  probs = c(0.90, 0.95, 0.975, 0.99, 0.995, 0.999)
)

quantile(
  seurat_combined$nFeature_RNA,
  probs = c(0.90, 0.95, 0.975, 0.99, 0.995, 0.999)
)

quantile(
  seurat_combined$percent.mt,
  probs = c(0.90, 0.95, 0.975, 0.99, 0.995, 0.999)
)

head(
  seurat_combined@meta.data[
    order(seurat_combined$nCount_RNA, decreasing = TRUE),
    c("Sample", "nCount_RNA", "nFeature_RNA", "percent.mt")
  ],
  20
)

head(
  seurat_combined@meta.data[
    order(seurat_combined$nFeature_RNA, decreasing = TRUE),
    c("Sample", "nCount_RNA", "nFeature_RNA", "percent.mt")
  ],
  20
)

# =========================================
# Count cells at possible QC cutoffs
# =========================================

cat("nFeature_RNA < 500:",
    sum(seurat_combined$nFeature_RNA < 500), "\n")

cat("nFeature_RNA > 10000:",
    sum(seurat_combined$nFeature_RNA > 10000), "\n")

cat("nCount_RNA > 50000:",
    sum(seurat_combined$nCount_RNA > 50000), "\n")

cat("nCount_RNA > 75000:",
    sum(seurat_combined$nCount_RNA > 75000), "\n")

cat("nCount_RNA > 100000:",
    sum(seurat_combined$nCount_RNA > 100000), "\n")

cat("percent.mt >= 10:",
    sum(seurat_combined$percent.mt >= 10), "\n")

table(
  seurat_combined$Sample[
    seurat_combined$nCount_RNA > 100000
  ]
)

table(
  seurat_combined$Sample[
    seurat_combined$nFeature_RNA > 10000
  ]
)


# =========================================
# Normalization
# =========================================

seurat_combined <- NormalizeData(
  seurat_combined,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

seurat_combined
Layers(seurat_combined[["RNA"]])

# =========================================
# Identify highly variable genes
# =========================================

seurat_combined <- FindVariableFeatures(
  seurat_combined,
  selection.method = "vst",
  nfeatures = 2000
)

length(VariableFeatures(seurat_combined))
head(VariableFeatures(seurat_combined), 20)
VariableFeaturePlot(seurat_combined)


# =========================================
# Check composition of variable genes
# =========================================

variable_genes <- VariableFeatures(seurat_combined)

sum(grepl("^MT-", variable_genes))

sum(grepl("^RPL|^RPS", variable_genes))

sum(variable_genes == "MALAT1")

variable_genes[
  grepl("^MT-|^RPL|^RPS|^MALAT1$", variable_genes)
]


# =========================================
# Scale highly variable genes
# =========================================

seurat_combined <- ScaleData(
  seurat_combined,
  features = VariableFeatures(seurat_combined)
)

seurat_combined


# =========================================
# PCA
# =========================================

seurat_combined <- RunPCA(
  seurat_combined,
  features = VariableFeatures(seurat_combined),
  npcs = 50
)

print(seurat_combined[["pca"]], dims = 1:10, nfeatures = 10)

ElbowPlot(
  seurat_combined,
  ndims = 50
)

# =========================================
# Construct the KNN/SNN graph
# =========================================

seurat_combined <- FindNeighbors(
  seurat_combined,
  dims = 1:20
)

seurat_combined


# =========================================
# Clustering
# =========================================

seurat_combined <- FindClusters(
  seurat_combined,
  resolution = 0.5
)

table(Idents(seurat_combined))

table(
  seurat_combined$Sample,
  Idents(seurat_combined)
)

# =========================================
# UMAP
# =========================================

seurat_combined <- RunUMAP(
  seurat_combined,
  dims = 1:20
)

DimPlot(
  seurat_combined,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

DimPlot(
  seurat_combined,
  reduction = "umap",
  group.by = "Sample"
)

# =========================================
# Find cluster marker genes
# =========================================

markers <- FindAllMarkers(
  seurat_combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

head(markers)
colnames(markers)

library(dplyr)

top10_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  )

top10_markers

write.csv(
  markers,
  "GSE297652_cluster_markers.csv",
  row.names = FALSE
)

# =========================================
# Top markers for every cluster
# =========================================

top10_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  arrange(cluster, desc(avg_log2FC))

top10_markers

marker_summary <- top10_markers %>%
  group_by(cluster) %>%
  summarise(
    top_genes = paste(gene, collapse = ", ")
  )

marker_summary

top5_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 5,
    with_ties = FALSE
  ) %>%
  arrange(cluster, desc(avg_log2FC))

top5_markers

write.csv(
  top10_markers,
  "GSE297652_top10_cluster_markers.csv",
  row.names = FALSE
)

write.csv(
  marker_summary,
  "GSE297652_cluster_marker_summary.csv",
  row.names = FALSE
)

print(marker_summary, n = 27, width = Inf)

top15_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 15,
    with_ties = FALSE
  ) %>%
  arrange(cluster, desc(avg_log2FC))

print(top15_markers, n = 405, width = Inf)


# =========================================
# Canonical marker validation
# =========================================

canonical_markers <- c(
  # T cells
  "CD3D", "CD3E", "TRBC1", "TRBC2",
  
  # CD8 / cytotoxic
  "CD8A", "CD8B", "CCL5", "NKG7", "GNLY",
  
  # B cells
  "CD19", "CD79A", "MS4A1", "CD37", "CD74",
  
  # Plasma cells
  "JCHAIN", "MZB1", "SDC1", "TNFRSF17",
  
  # Myeloid/macrophage
  "LYZ", "LST1", "TYROBP", "FCER1G",
  "C1QA", "C1QB", "C1QC", "APOC1",
  
  # Mast
  "TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT",
  
  # Fibroblast
  "COL1A1", "COL1A2", "DCN", "LUM", "COL3A1",
  
  # Pericyte/smooth muscle
  "RGS5", "CSPG4", "MCAM", "PDGFRB", "ACTA2", "MYH11",
  
  # Endothelial
  "PECAM1", "VWF", "EMCN", "KDR", "ESM1",
  
  # ACKR1 endothelial
  "ACKR1", "SELE", "SELP",
  
  # Epithelial
  "EPCAM", "KRT8", "KRT18", "KRT19",
  
  # Prostate epithelial
  "KLK2", "KLK3", "KLK4", "MSMB", "ACPP", "AMACR",
  
  # Neuroendocrine
  "CHGA", "CHGB", "SCG3", "SCGN", "SYP", "VGF",
  
  # Ciliated
  "FOXJ1", "TPPP3", "DNAI1", "CCNO",
  
  # Cycling
  "MKI67", "TOP2A", "CDK1", "AURKB"
)

canonical_markers <- canonical_markers[
  canonical_markers %in% rownames(seurat_combined)
]

length(canonical_markers)

DotPlot(
  seurat_combined,
  features = canonical_markers,
  group.by = "seurat_clusters"
) +
  RotatedAxis()

DotPlot(
  seurat_combined,
  features = c(
    "CD3D", "CD3E", "CD8A", "CCL5",
    "NKG7", "GNLY",
    "CD19", "CD79A", "MS4A1",
    "JCHAIN", "MZB1", "TNFRSF17",
    "LYZ", "LST1", "TYROBP",
    "C1QA", "C1QB", "C1QC",
    "TPSAB1", "TPSB2", "CPA3", "MS4A2"
  )
) +
  RotatedAxis()

DotPlot(
  seurat_combined,
  features = c(
    "COL1A1", "COL1A2", "DCN", "LUM", "COL3A1",
    "RGS5", "CSPG4", "MCAM", "PDGFRB", "MYH11",
    "PECAM1", "VWF", "EMCN", "KDR", "ESM1",
    "ACKR1", "SELE", "SELP"
  )
) +
  RotatedAxis()

DotPlot(
  seurat_combined,
  features = c(
    "EPCAM", "KRT8", "KRT18", "KRT19",
    "KLK2", "KLK3", "KLK4", "MSMB", "ACPP", "AMACR",
    "KRT17", "KRT15",
    "FOXJ1", "TPPP3", "DNAI1",
    "CHGA", "CHGB", "SCG3", "SCGN", "SYP",
    "MKI67", "TOP2A", "CDK1", "AURKB"
  )
) +
  RotatedAxis()


# =========================================
# Investigate ambiguous clusters
# =========================================

table(Idents(seurat_combined))

# Number of cells in ambiguous clusters
table(Idents(seurat_combined))[
  names(table(Idents(seurat_combined))) %in% c("17", "19", "26")
]

seurat_combined@meta.data %>%
  dplyr::filter(seurat_clusters %in% c("17", "19", "26")) %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(
    cells = dplyr::n(),
    median_nCount = median(nCount_RNA),
    median_nFeature = median(nFeature_RNA),
    median_mt = median(percent.mt)
  )

FeaturePlot(
  seurat_combined,
  features = c(
    "TPSAB1",
    "TPSB2",
    "KIT",
    "LYZ",
    "TYROBP",
    "KRT4",
    "KRT7",
    "DNAI1",
    "CCNO",
    "JCHAIN",
    "IGHG1"
  ),
  ncol = 3
)

# =========================================
# Cell-level investigation of clusters
# 17, 19 and 26
# =========================================

ambiguous_cells <- WhichCells(
  seurat_combined,
  idents = c("17", "19", "26")
)

ambiguous_data <- FetchData(
  seurat_combined,
  vars = c(
    "seurat_clusters",
    "TPSAB1",
    "TPSB2",
    "KIT",
    "LYZ",
    "TYROBP",
    "KRT4",
    "KRT7",
    "DNAI1",
    "CCNO",
    "JCHAIN",
    "IGHG1"
  ),
  cells = ambiguous_cells
)

library(dplyr)

marker_positive_summary <- ambiguous_data %>%
  group_by(seurat_clusters) %>%
  summarise(
    cells = n(),
    
    TPSAB1_pos = sum(TPSAB1 > 0),
    TPSB2_pos = sum(TPSB2 > 0),
    KIT_pos = sum(KIT > 0),
    
    LYZ_pos = sum(LYZ > 0),
    TYROBP_pos = sum(TYROBP > 0),
    
    KRT4_pos = sum(KRT4 > 0),
    KRT7_pos = sum(KRT7 > 0),
    
    DNAI1_pos = sum(DNAI1 > 0),
    CCNO_pos = sum(CCNO > 0),
    
    JCHAIN_pos = sum(JCHAIN > 0),
    IGHG1_pos = sum(IGHG1 > 0)
  )

marker_positive_summary

marker_positive_percent <- marker_positive_summary %>%
  mutate(
    across(
      -c(seurat_clusters, cells),
      ~ .x / cells * 100
    )
  )

marker_positive_percent

marker_positive_percent %>%
  select(
    seurat_clusters,
    cells,
    KRT4_pos,
    KRT7_pos,
    DNAI1_pos,
    CCNO_pos,
    JCHAIN_pos,
    IGHG1_pos
  )

coexpression_summary <- ambiguous_data %>%
  group_by(seurat_clusters) %>%
  summarise(
    mast_epithelial =
      sum(
        TPSAB1 > 0 &
          (KRT4 > 0 | KRT7 > 0)
      ),
    
    mast_myeloid =
      sum(
        TPSAB1 > 0 &
          (LYZ > 0 | TYROBP > 0)
      ),
    
    ciliated_immune =
      sum(
        (DNAI1 > 0 | CCNO > 0) &
          (JCHAIN > 0 | IGHG1 > 0)
      ),
    
    cells = n()
  )

coexpression_summary

coexpression_summary %>%
  mutate(
    mast_epithelial_pct = mast_epithelial / cells * 100,
    mast_myeloid_pct = mast_myeloid / cells * 100,
    ciliated_immune_pct = ciliated_immune / cells * 100
  )

# =========================================
# Examine ambiguous clusters individually
# =========================================

VlnPlot(
  seurat_combined,
  features = c(
    "TPSAB1",
    "TPSB2",
    "KIT",
    "KRT4",
    "KRT7",
    "LYZ",
    "TYROBP"
  ),
  group.by = "seurat_clusters",
  idents = c("14", "17", "26"),
  pt.size = 0.1,
  ncol = 3
)

VlnPlot(
  seurat_combined,
  features = c(
    "DNAI1",
    "CCNO",
    "JCHAIN",
    "IGHG1",
    "FOXJ1",
    "TPPP3"
  ),
  group.by = "seurat_clusters",
  idents = c("19", "20", "22"),
  pt.size = 0.1,
  ncol = 3
)

# =========================================
# Final cell-type annotation
# =========================================

cluster_annotations <- c(
  "0"  = "Activated CD8 T cells",
  "1"  = "Prostate epithelial",
  "2"  = "Naive/memory T cells",
  "3"  = "Prostate epithelial",
  "4"  = "Pericyte/smooth muscle",
  "5"  = "Macrophages",
  "6"  = "Basal epithelial",
  "7"  = "Prostate epithelial",
  "8"  = "Neuronal-associated prostate epithelial",
  "9"  = "Prostate epithelial",
  "10" = "Endothelial",
  "11" = "Fibroblast/stromal",
  "12" = "ACKR1+ endothelial",
  "13" = "B cells",
  "14" = "Mast cells",
  "15" = "Endothelial",
  "16" = "Neuroendocrine",
  "17" = "Mast cells - myeloid-associated",
  "18" = "Activated fibroblast/stromal",
  "19" = "Ciliated epithelial",
  "20" = "Ciliated/secretory epithelial",
  "21" = "Activated macrophage/myeloid",
  "22" = "Plasma cells",
  "23" = "Cycling/proliferating epithelial",
  "24" = "Secretory/progenitor-like epithelial",
  "25" = "Prostate epithelial",
  "26" = "Mixed mast-epithelial / doublet-like"
)

seurat_combined$cell_type <- unname(
  cluster_annotations[
    as.character(seurat_combined$seurat_clusters)
  ]
)

head(
  seurat_combined@meta.data[
    ,
    c("seurat_clusters", "cell_type")
  ]
)
table(seurat_combined$cell_type)
table(
  seurat_combined$Sample,
  seurat_combined$cell_type
)

DimPlot(
  seurat_combined,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
)


# =========================================
# Short labels for clean UMAP visualization
# =========================================

cell_type_short <- c(
  "Activated CD8 T cells" = "CD8 T",
  "Prostate epithelial" = "Prostate epithelial",
  "Naive/memory T cells" = "Naive/memory T",
  "Pericyte/smooth muscle" = "Pericyte/SMC",
  "Macrophages" = "Macrophage",
  "Basal epithelial" = "Basal epithelial",
  "Neuronal-associated prostate epithelial" = "Neuronal-associated",
  "Endothelial" = "Endothelial",
  "Fibroblast/stromal" = "Fibroblast",
  "ACKR1+ endothelial" = "ACKR1+ endothelial",
  "B cells" = "B cell",
  "Mast cells" = "Mast cell",
  "Neuroendocrine" = "Neuroendocrine",
  "Mast cells - myeloid-associated" = "Mast/myeloid",
  "Activated fibroblast/stromal" = "Activated fibroblast",
  "Ciliated epithelial" = "Ciliated epithelial",
  "Ciliated/secretory epithelial" = "Ciliated/secretory",
  "Activated macrophage/myeloid" = "Activated macrophage",
  "Plasma cells" = "Plasma cell",
  "Cycling/proliferating epithelial" = "Cycling epithelial",
  "Secretory/progenitor-like epithelial" = "Secretory/progenitor",
  "Mixed mast-epithelial / doublet-like" = "Mixed/doublet-like"
)

seurat_combined$cell_type_short <- unname(
  cell_type_short[seurat_combined$cell_type]
)

DimPlot(
  seurat_combined,
  reduction = "umap",
  group.by = "cell_type_short",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
) +
  ggtitle("GSE297652 Single-cell RNA-seq: Cell-Type Annotation"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  )

ggsave(
  filename = file.path(
    fig_dir,
    "GSE297652_final_annotated_UMAP.png"
  ),
  width = 12,
  height = 9,
  dpi = 300
)

# =========================================
# UMAP BY SAMPLE
# =========================================

DimPlot(
  seurat_combined,
  reduction = "umap",
  group.by = "Sample",
  pt.size = 0.4
) +
  ggtitle(
    "GSE297652 UMAP by Sample"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    )
  )
ggsave(
  filename = file.path(
    fig_dir,
    "GSE297652_UMAP_by_sample.png"
  ),
  width = 10,
  height = 8,
  dpi = 300
)

# =========================================
# Cell-type composition by sample
# =========================================

celltype_by_sample <- table(
  seurat_combined$Sample,
  seurat_combined$cell_type
)

celltype_by_sample

celltype_percent <- prop.table(
  celltype_by_sample,
  margin = 1
) * 100

round(celltype_percent, 2)

celltype_total <- sort(
  table(seurat_combined$cell_type),
  decreasing = TRUE
)

celltype_total

round(
  prop.table(celltype_total) * 100,
  2
)

# =========================================
# Cell-type composition by sample
# =========================================

library(ggplot2)

celltype_percent_df <- as.data.frame(celltype_percent)

colnames(celltype_percent_df) <- c(
  "Sample",
  "CellType",
  "Percent"
)

ggplot(
  celltype_percent_df,
  aes(
    x = Sample,
    y = Percent,
    fill = CellType
  )
) +
  geom_bar(
    stat = "identity",
    position = "stack"
  ) +
  labs(
    title = "Cell-Type Composition Across Samples",
    x = "Sample",
    y = "Percentage of Cells",
    fill = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

celltype_percent_clean <- celltype_percent_df %>%
  filter(
    CellType != "Mixed mast-epithelial / doublet-like"
  ) %>%
  group_by(Sample) %>%
  mutate(
    Percent = Percent / sum(Percent) * 100
  )

ggplot(
  celltype_percent_clean,
  aes(
    x = Sample,
    y = Percent,
    fill = CellType
  )
) +
  geom_bar(
    stat = "identity",
    position = "stack"
  ) +
  labs(
    title = "Cell-Type Composition Across Samples",
    x = "Sample",
    y = "Percentage of Cells",
    fill = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

write.csv(
  celltype_by_sample,
  "GSE297652_celltype_counts_by_sample.csv"
)

write.csv(
  celltype_percent,
  "GSE297652_celltype_percent_by_sample.csv"
)

# =========================================
# Marker validation heatmap
# =========================================

marker_genes <- c(
  # T cells
  "CD3D", "CD3E", "CD8A", "CCL5", "IL7R", "CCR7",
  
  # B / plasma
  "CD79A", "MS4A1", "JCHAIN", "MZB1", "IGHG1",
  
  # Myeloid
  "LYZ", "TYROBP", "FCER1G", "C1QA", "C1QC",
  
  # Mast
  "TPSAB1", "TPSB2", "KIT",
  
  # Endothelial
  "PECAM1", "VWF", "EMCN", "KDR", "ACKR1",
  
  # Fibroblast / stromal
  "COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA",
  
  # Pericyte / smooth muscle
  "RGS5", "MCAM", "ACTA2", "MYH11",
  
  # Epithelial
  "EPCAM", "KRT8", "KRT18", "KRT19", "KRT7",
  
  # Basal epithelial
  "KRT5", "KRT14", "KRT17",
  
  # Ciliated
  "FOXJ1", "TPPP3", "DNAI1", "CCNO",
  
  # Prostate epithelial
  "KLK2", "KLK3", "MSMB", "ACPP",
  
  # Neuroendocrine
  "CHGA", "CHGB", "SYP",
  
  # Cycling
  "MKI67", "TOP2A"
)

# Keep only genes present in the dataset
marker_genes_present <- marker_genes[
  marker_genes %in% rownames(seurat_combined)
]

length(marker_genes_present)
marker_genes_present

DotPlot(
  seurat_combined,
  features = marker_genes_present,
  group.by = "cell_type_short"
) +
  RotatedAxis()


# =========================================
# Marker heatmap for annotated cell types
# =========================================

DoHeatmap(
  seurat_combined,
  features = marker_genes_present,
  group.by = "cell_type_short",
  size = 3
) +
  ggtitle("Marker Gene Expression Across Cell Types")

heatmap_markers <- c(
  "CD3D", "CD8A", "IL7R",
  "CD79A", "JCHAIN", "MZB1",
  "LYZ", "C1QA", "TYROBP",
  "TPSAB1", "TPSB2", "KIT",
  "PECAM1", "VWF", "ACKR1",
  "COL1A1", "DCN", "LUM",
  "RGS5", "ACTA2", "MYH11",
  "EPCAM", "KRT8", "KRT19",
  "KRT5", "KRT14", "KRT17",
  "FOXJ1", "TPPP3", "DNAI1",
  "KLK2", "KLK3", "MSMB",
  "CHGA", "CHGB",
  "MKI67", "TOP2A"
)

DoHeatmap(
  seurat_combined,
  features = heatmap_markers,
  group.by = "cell_type_short",
  size = 3
) +
  ggtitle("Canonical Marker Expression by Cell Type")

# =========================================
# Average marker expression by cell type
# =========================================

avg_marker_expression <- AverageExpression(
  seurat_combined,
  assays = "RNA",
  features = heatmap_markers,
  group.by = "cell_type_short",
  slot = "data"
)$RNA

pheatmap::pheatmap(
  avg_marker_expression,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 9,
  fontsize_col = 8,
  angle_col = 45,
  main = "Canonical Marker Expression by Cell Type"
)

fig_dir <- "D:/Bioinfomatics_project/SingleCell-RNAseq-Prostate-Cancer/results/figures"

pheatmap::pheatmap(
  avg_marker_expression,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 9,
  fontsize_col = 8,
  angle_col = 45,
  main = "Canonical Marker Expression by Cell Type",
  filename = file.path(fig_dir, "GSE297652_final_marker_heatmap.png"),
  width = 12,
  height = 8
)
# =========================================
# Save completed Seurat object
# =========================================

saveRDS(
  seurat_combined,
  file = "GSE297652_single_cell_annotated.rds"
)


# =========================================
# Cell-type-specific marker summary
# =========================================

celltype_markers <- markers %>%
  mutate(
    cell_type = unname(
      cluster_annotations[as.character(cluster)]
    )
  ) %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.5
  ) %>%
  arrange(cell_type, desc(avg_log2FC))

head(celltype_markers)

table(celltype_markers$cell_type)

top10_celltype_markers <- celltype_markers %>%
  group_by(cell_type) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  ) %>%
  arrange(cell_type, desc(avg_log2FC))

top10_celltype_markers

write.csv(
  top10_celltype_markers,
  "GSE297652_top10_celltype_markers.csv",
  row.names = FALSE
)


# =========================================
# Prepare marker genes for enrichment
# =========================================

celltype_gene_lists <- celltype_markers %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.5
  ) %>%
  group_by(cell_type) %>%
  summarise(
    genes = list(unique(gene)),
    .groups = "drop"
  )

celltype_gene_lists

# =========================================
# Convert gene symbols → Entrez IDs
# =========================================

gene_conversion <- bitr(
  unique(celltype_markers$gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(gene_conversion)

cat(
  "Total marker genes:",
  length(unique(celltype_markers$gene)),
  "\n"
)

cat(
  "Successfully mapped:",
  length(unique(gene_conversion$SYMBOL)),
  "\n"
)


# =========================================
# Define GO enrichment background
# =========================================

background_genes <- rownames(seurat_combined)

background_conversion <- bitr(
  background_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_entrez <- unique(
  background_conversion$ENTREZID
)

cat(
  "Background genes tested:",
  length(background_genes),
  "\n"
)

cat(
  "Background genes mapped to Entrez:",
  length(background_entrez),
  "\n"
)

# =========================================
# GO Biological Process enrichment
# =========================================

go_results <- list()

for (ct in unique(celltype_markers$cell_type)) {
  
  genes <- celltype_markers %>%
    filter(cell_type == ct) %>%
    pull(gene) %>%
    unique()
  
  gene_ids <- gene_conversion %>%
    filter(SYMBOL %in% genes) %>%
    pull(ENTREZID) %>%
    unique()
  
  if (length(gene_ids) >= 5) {
    
    ego <- enrichGO(
      gene = gene_ids,
      universe = background_entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    
    go_results[[ct]] <- ego
  }
}

go_results_df <- bind_rows(
  lapply(
    names(go_results),
    function(ct) {
      
      if (
        !is.null(go_results[[ct]]) &&
        nrow(as.data.frame(go_results[[ct]])) > 0
      ) {
        
        df <- as.data.frame(go_results[[ct]])
        df$cell_type <- ct
        df
        
      } else {
        NULL
      }
    }
  )
)

dim(go_results_df)

head(go_results_df)

table(go_results_df$cell_type)

write.csv(
  go_results_df,
  "GSE297652_GO_Biological_Process_Enrichment.csv",
  row.names = FALSE
)

# =========================================
# Top GO Biological Processes
# =========================================

go_top10 <- go_results_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(cell_type) %>%
  slice_min(
    order_by = p.adjust,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(cell_type, p.adjust)

go_top10

write.csv(
  go_top10,
  "GSE297652_top10_GO_terms_by_cell_type.csv",
  row.names = FALSE
)

# =========================================
# Top 5 GO Biological Processes
# =========================================

go_top5 <- go_results_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(cell_type) %>%
  slice_min(
    order_by = p.adjust,
    n = 5,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(cell_type, p.adjust)

go_top5

write.csv(
  go_top5,
  "GSE297652_top5_GO_terms_by_cell_type.csv",
  row.names = FALSE
)


go_main <- go_top5 %>%
  mutate(
    term_label = str_wrap(
      Description,
      width = 35
    ),
    significance = -log10(p.adjust)
  )

go_final_plot <- ggplot(
  go_main,
  aes(
    x = significance,
    y = term_label,
    size = Count
  )
) +
  geom_point() +
  facet_wrap(
    ~cell_type,
    ncol = 3,
    scales = "free_y"
  ) +
  labs(
    title = "Cell-Type-Specific GO Biological Processes",
    x = "-log10 Adjusted P-value",
    y = "Biological Process",
    size = "Gene Count"
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(
      size = 10,
      face = "bold"
    ),
    strip.background = element_rect(
      fill = "grey95",
      color = "grey70"
    ),
    axis.text.y = element_text(
      size = 7
    ),
    axis.text.x = element_text(
      size = 8
    ),
    axis.title = element_text(
      size = 10,
      face = "bold"
    ),
    plot.title = element_text(
      size = 17,
      face = "bold",
      hjust = 0.5
    ),
    panel.spacing = unit(
      1,
      "lines"
    )
  )

go_final_plot

fig_dir <- "D:/Bioinfomatics_project/SingleCell-RNAseq-Prostate-Cancer/results/figures"

ggsave(
  filename = file.path(
    fig_dir,
    "GSE297652_celltype_GO_enrichment_FINAL.png"
  ),
  plot = go_final_plot,
  width = 14,
  height = 15,
  dpi = 300
)


# -----------------------------------------
# Prepare KEGG background
# -----------------------------------------

background_conversion <- bitr(
  rownames(seurat_combined),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

kegg_background <- unique(
  background_conversion$ENTREZID
)

length(kegg_background)

# =========================================
# KEGG enrichment for each cell type
# =========================================

kegg_results <- lapply(
  unique(celltype_markers$cell_type),
  function(ct) {
    
    genes <- celltype_markers %>%
      filter(cell_type == ct) %>%
      pull(gene) %>%
      unique()
    
    conversion <- bitr(
      genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    entrez_genes <- unique(
      conversion$ENTREZID
    )
    
    # Need enough genes for enrichment
    if (length(entrez_genes) < 5) {
      return(NULL)
    }
    
    ek <- enrichKEGG(
      gene = entrez_genes,
      universe = kegg_background,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.20
    )
    
    if (is.null(ek)) {
      return(NULL)
    }
    
    result <- as.data.frame(ek)
    
    if (nrow(result) == 0) {
      return(NULL)
    }
    
    result$cell_type <- ct
    
    return(result)
  }
)

# Combine all cell types
kegg_results_df <- bind_rows(kegg_results)

dim(kegg_results_df)

head(kegg_results_df)

# =========================================
# Significant KEGG pathways by cell type
# =========================================

kegg_significant_counts <- kegg_results_df %>%
  filter(p.adjust < 0.05) %>%
  count(cell_type, sort = TRUE)

kegg_significant_counts

# =========================================
# Top 5 KEGG pathways per cell type
# =========================================

kegg_top5 <- kegg_results_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(cell_type) %>%
  slice_min(
    order_by = p.adjust,
    n = 5,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(cell_type, p.adjust)

kegg_top5 %>%
  dplyr::select(
    cell_type,
    Description,
    GeneRatio,
    FoldEnrichment,
    p.adjust,
    Count
  )

write.csv(
  kegg_top5,
  "GSE297652_top5_KEGG_pathways_by_cell_type.csv",
  row.names = FALSE
)

write.csv(
  kegg_results_df,
  "GSE297652_KEGG_enrichment_by_cell_type.csv",
  row.names = FALSE
)


# =========================================
# TOP 3 KEGG PATHWAYS PER CELL TYPE
# =========================================

kegg_top3 <- kegg_results_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(cell_type) %>%
  slice_min(
    order_by = p.adjust,
    n = 3,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(cell_type, p.adjust) %>%
  mutate(
    pathway_label = str_wrap(
      Description,
      width = 28
    )
  )

kegg_top3

write.csv(
  kegg_top3,
  "GSE297652_top3_KEGG_pathways_by_cell_type.csv",
  row.names = FALSE
)



# =========================================
# KEGG — ONE REPRESENTATIVE PATHWAY
# PER CELL TYPE
# =========================================

kegg_main <- kegg_results_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(cell_type) %>%
  slice_min(
    order_by = p.adjust,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    pathway_label = str_wrap(
      Description,
      width = 25
    ),
    significance = -log10(p.adjust)
  )

kegg_main %>%
  dplyr::select(
    cell_type,
    Description,
    Count,
    FoldEnrichment,
    p.adjust
  )

write.csv(
  kegg_main,
  "GSE297652_representative_KEGG_pathway_by_cell_type.csv",
  row.names = FALSE
)


# =========================================
# Final KEGG plot for GitHub
# =========================================

kegg_main <- kegg_results_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(cell_type) %>%
  slice_min(
    order_by = p.adjust,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    pathway_label = str_wrap(Description, width = 32),
    significance = -log10(p.adjust)
  )

kegg_final_plot <- ggplot(
  kegg_main,
  aes(
    x = significance,
    y = pathway_label,
    size = Count
  )
) +
  geom_point() +
  facet_wrap(
    ~cell_type,
    ncol = 3,
    scales = "free_y"
  ) +
  labs(
    title = "Representative KEGG Pathway by Cell Type",
    x = "-log10 Adjusted P-value",
    y = "KEGG Pathway",
    size = "Gene Count"
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(
      size = 10,
      face = "bold"
    ),
    strip.background = element_rect(
      fill = "grey95",
      color = "grey70"
    ),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title = element_text(
      size = 10,
      face = "bold"
    ),
    plot.title = element_text(
      size = 17,
      face = "bold",
      hjust = 0.5
    ),
    panel.spacing = unit(1, "lines")
  )

kegg_final_plot

fig_dir <- "D:/Bioinfomatics_project/SingleCell-RNAseq-Prostate-Cancer/results/figures"

ggsave(
  file.path(
    fig_dir,
    "GSE297652_representative_KEGG_pathway_FINAL.png"
  ),
  plot = kegg_final_plot,
  width = 14,
  height = 15,
  dpi = 300
)


# =========================================
# INTEGRATED CELL-TYPE SUMMARY
# Markers + GO + KEGG
# =========================================

# Top 5 markers
marker_summary <- celltype_markers %>%
  group_by(cell_type) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 5,
    with_ties = FALSE
  ) %>%
  summarise(
    top_markers = paste(gene, collapse = ", "),
    .groups = "drop"
  )

# Top 3 GO terms
go_summary <- go_results_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(cell_type) %>%
  slice_min(
    order_by = p.adjust,
    n = 3,
    with_ties = FALSE
  ) %>%
  summarise(
    top_GO = paste(Description, collapse = " | "),
    .groups = "drop"
  )

# Top 3 KEGG pathways
kegg_summary <- kegg_results_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(cell_type) %>%
  slice_min(
    order_by = p.adjust,
    n = 3,
    with_ties = FALSE
  ) %>%
  summarise(
    top_KEGG = paste(Description, collapse = " | "),
    .groups = "drop"
  )

# Combine
integrated_summary <- marker_summary %>%
  left_join(
    go_summary,
    by = "cell_type"
  ) %>%
  left_join(
    kegg_summary,
    by = "cell_type"
  )

integrated_summary

write.csv(
  integrated_summary,
  "GSE297652_integrated_celltype_summary.csv",
  row.names = FALSE
)

# =========================================
# FINAL INTEGRATED RESULTS TABLE
# =========================================

celltype_counts <- seurat_combined@meta.data %>%
  count(cell_type, name = "Cell_Count") %>%
  mutate(
    Percentage = round(
      100 * Cell_Count / sum(Cell_Count),
      2
    )
  )

final_results_table <- celltype_counts %>%
  left_join(
    integrated_summary,
    by = "cell_type"
  ) %>%
  arrange(desc(Cell_Count))

final_results_table

write.csv(
  final_results_table,
  "GSE297652_FINAL_integrated_results_table.csv",
  row.names = FALSE
)
