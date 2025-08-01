library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(SeuratData)
# Load UMI sparse matrix (triplet format)
umi_table <- read.csv("data/GSE125527_UMI_cell_table_sparse.csv")
# Load row and column names
gene_names <- read.csv("data/GSE125527_gene_id_rownames.csv", header = FALSE)[[1]]
cell_ids <- read.csv("data/GSE125527_cell_id_colnames.csv", header = FALSE)[[1]]
# Check: column = cell, row = gene
head(umi_table)
# Rename columns for clarity
colnames(umi_table) <- c("gene", "cell", "count")
# Construct sparse matrix
sparse_counts <- sparseMatrix(
  i = umi_table$gene,
  j = umi_table$cell,
  x = umi_table$count
)
# Assign row and column names
rownames(sparse_counts) <- gene_names
colnames(sparse_counts) <- cell_ids
# Check matrix shape
dim(sparse_counts)
# metdata 
# Load metadata
metadata <- read.csv("data/GSE125527_cell_metadata.csv")
rownames(metadata) <- metadata$cell_id
# Ensure metadata matches cell order
metadata <- metadata[colnames(sparse_counts), ]
# Create Seurat object
seurat <- CreateSeuratObject(counts = sparse_counts, meta.data = metadata, project = "GSE125527")
head(seurat@meta.data)
# qc
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
seurat <- subset(seurat,
                 subset = nFeature_RNA > 200 &
                   nFeature_RNA < 4000 &
                   nCount_RNA < 20000 &
                   percent.mt < 10)
#norm
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
# Examine and visualize PCA results a few different ways
print(seurat[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
DimPlot(seurat, reduction = "pca") + NoLegend()
DimHeatmap(seurat, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(seurat)
seurat <- FindNeighbors(seurat, dims = 1:10)
seurat <- FindClusters(seurat, resolution = 0.2)
# non linear 
seurat <- RunUMAP(seurat, dims = 1:10)
DimPlot(seurat, reduction = "umap", label = TRUE)
saveRDS(seurat, file = "output/uc.rds")
# marker 
# find all markers of cluster 4
cluster4.markers <- FindMarkers(seurat, ident.1 = 4)
head(cluster4.markers, n = 5)
# find all markers distinguishing cluster 3 from clusters 0 and 1,2,7
cluster3.markers <- FindMarkers(seurat, ident.1 = 3, ident.2 = c(0, 1,2,7))
head(cluster3.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE)
seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
VlnPlot(seurat, features = c("IL7R", "CCR7", "S100A4", "CD8A"))
# plot raw counts as well
VlnPlot(seurat, features = c("MS4A1", "NKG7", "GNLY"), slot = "counts", log = TRUE)
FeaturePlot(seurat, features = c("IL7R", "CCR7", "S100A4", "CD8A","MS4A1", "NKG7", "GNLY"))
seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(seurat, features = top10$gene) + NoLegend()
# assign
new.cluster.ids <- c("CD8 T", "Memory CD4 T", "Naive CD4 T", "B", "NK", "CD14+ Mono","FCGR3A+ Mono", "Memory CD4 T")
names(new.cluster.ids) <- levels(seurat)
seurat <- RenameIdents(seurat, new.cluster.ids)
seurat$seurat_annotations <- Idents(seurat)
DimPlot(seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
library(ggplot2)
final_clustering <- DimPlot(seurat, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "output/final_clustering.jpg", height = 7, width = 12, plot = final_clustering, quality = 50)
saveRDS(seurat, file = "output/seurat_final.rds")

#now we will compare the cells across conditions
seurat$new2_idents <- paste(Idents(seurat), seurat$disease_assignment, sep = "_")
Idents(seurat) <- "new2_idents"
# Set the future global size limit to 8 GB
options(future.globals.maxSize = 8000 * 1024^2)  # 8 GB , u need it due to the internal function of paralellization, otherwise U will get error as it will exceed the default of the memory size
tcell <- FindMarkers(seurat, ident.1 = "CD8 T_diseased", ident.2 = "CD8 T_healthy", verbose = FALSE)
head(tcell, n = 10)
# aggregate 
seurat$annotations <- Idents(seurat)
# pseudobulk the counts based on donor-condition-celltype
pseudo_seurat <- AggregateExpression(seurat, assays = "RNA", return.seurat = T, group.by = c("disease_assignment","patient_assignment", "seurat_annotations"))
pseudo_seurat$new_idents <- paste(pseudo_seurat$seurat_annotations, pseudo_seurat$disease_assignment, sep = "_")
# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_seurat))
Idents(pseudo_seurat) <- "new_idents"
bulk.tcr <- FindMarkers(object = pseudo_seurat, 
                            ident.1 = "CD8 T_diseased", 
                            ident.2 = "CD8 T_healthy",
                            test.use = "DESeq2")
head(bulk.tcr, n = 15)


# compare the DE P-values between the single-cell level and the pseudobulk level results
names(bulk.tcr) <- paste0(names(bulk.tcr), ".bulk")
bulk.tcr$gene <- rownames(bulk.tcr)

names(tcell) <- paste0(names(tcell), ".sc")
tcell$gene <- rownames(tcell)

merge_dat <- merge(tcell, bulk.tcr, by = "gene")
merge_dat <- merge_dat[order(merge_dat$p_val.bulk), ]

# Number of genes that are marginally significant in both; marginally significant only in bulk; and marginally significant only in single-cell
common <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                 merge_dat$p_val.sc < 0.05)]
only_sc <- merge_dat$gene[which(merge_dat$p_val.bulk > 0.05 & 
                                  merge_dat$p_val.sc < 0.05)]
only_bulk <- merge_dat$gene[which(merge_dat$p_val.bulk < 0.05 & 
                                    merge_dat$p_val.sc > 0.05)]
print(paste0('# Common: ',length(common)))
print(paste0('# Only in bulk: ',length(only_bulk)))
print(paste0('# Only in single-cell: ',length(only_sc)))
#
# create a new column to annotate sample-condition-celltype in the single-cell dataset
seurat$id.disease <- paste0(seurat$patient_assignment, "-", seurat$disease_assignment)
# generate violin plot 
seurat$new_idents <- paste(seurat$seurat_annotations, seurat$disease_assignment, sep = "_")
Idents(seurat) <- "new_idents"
print(merge_dat[merge_dat$gene%in%common[1:2],c('gene','p_val.sc','p_val.bulk')])
VlnPlot(seurat, features = common[1:10], idents = c("CD8 T_diseased", "CD8 T_healthy"), group.by = "disease_assignment") 
VlnPlot(seurat, features = common[1:2], idents = c("CD8 T_diseased", "CD8 T_healthy"), group.by = "patient_assignment", ncol = 1) # to ensure robustness

print(merge_dat[merge_dat$gene%in%c('TMCO3','MFAP1'),c('gene','p_val.sc','p_val.bulk')]) # DE only under sc
VlnPlot(seurat, features <- c('PPIB','POMC'), idents = c("CD8 T_diseased", "CD8 T_healthy"), group.by = "disease_assignment") 
VlnPlot(seurat, features <- c('PPIB','POMC'), idents = c("CD8 T_diseased", "CD8 T_healthy"), group.by = "patient_assignment", ncol = 1) # to ensure robustness

# automatic annotation
library(SingleR)
library(celldex)
# Load reference dataset (e.g., Human Primary Cell Atlas)
ref <- celldex::HumanPrimaryCellAtlasData()
# Extract normalized data matrix from Seurat object
seurat.data <- GetAssayData(seurat, assay = "RNA", slot = "data")
# Run SingleR annotation
pred <- SingleR(test = seurat.data, ref = ref, labels = ref$label.main)
# Add labels to Seurat metadata
seurat$SingleR.labels <- pred$pruned.labels
# Visualize automatic annotation on UMAP
DimPlot(seurat, group.by = "SingleR.labels", label = TRUE, repel = TRUE)
# search for suitable one
for (res in c(0.2, 0.5, 0.8, 1.0)) {
  seurat <- FindClusters(seurat, resolution = res)
  print(DimPlot(seurat, reduction = "umap", label = TRUE) + ggtitle(paste("Resolution:", res)))
}
# Set active identities to final clusters
Idents(seurat) <- seurat$RNA_snn_res.0.2  # or whatever resolution you chose
# markers again
markers <- FindAllMarkers(seurat, only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)
head(markers)


# Between conditions (e.g., control vs IBD)
conditions <- seurat@meta.data$disease_assignment
conditions <- as.factor(conditions)
Idents(seurat) <- conditions # Set metadata column
FindMarkers(seurat, ident.1 = "diseased", ident.2 = "healthy")
uc_healthy <- FindMarkers(seurat, ident.1 = "diseased", ident.2 = "healthy")

# pathway 
# Optional: convert gene symbols to Entrez IDs
library(clusterProfiler)
library(org.Hs.eg.db)
genes <- rownames(uc_healthy)
entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# GO enrichment
ego <- enrichGO(entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
dotplot(ego)
# WE CAN proceed to marker genes (all), compare  

b_clusters <- read.csv("data/GSE125527_Bcell_cluster.csv")
t_clusters <- read.csv("data/GSE125527_Tcell_cluster.csv")

combined_clusters <- rbind(b_clusters, t_clusters)
combined_clusters <- combined_clusters %>% column_to_rownames("cell_id")
seurat <- AddMetaData(seurat, metadata = combined_clusters)

DimPlot(seurat, group.by = "Cluster_id")  # or whatever column holds cluster labels




