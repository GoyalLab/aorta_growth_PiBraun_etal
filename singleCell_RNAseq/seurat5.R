# Load required libraries
library(Seurat)      # For single-cell RNA sequencing data analysis
library(harmony)     # For integrating scRNA-seq data across conditions
library(dplyr)       # For data manipulation
library(Matrix)      # For working with sparse matrices
require("ggrepel")   # For decorative labels in ggplot

# Create a directory for plot outputs
dir.create("plots")

# Read in the aggregated expression data
all_expr = readRDS("../aggr_matrix/aggr_matrix.rds")

# Create a Seurat object using the expression data
sample = CreateSeuratObject(counts=all_expr, min.cells=1, min.features=100)
# Inspect dimensions of the Seurat object
dim(sample)
# Count how many samples are present based on original identifiers
table(sample$orig.ident)

# Identify mitochondrial genes and calculate their contribution
mito.features = grep(pattern="^mt-", x=rownames(x=sample), value=TRUE) # Identify mitochondrial genes
percent.mito = Matrix::colSums(x=GetAssayData(object=sample, slot="counts")[mito.features,]) /
  Matrix::colSums(x=GetAssayData(object=sample, slot="counts"))  # Calculate the percentage of counts from mitochondrial genes
# Store percentages in the Seurat object metadata
sample[["percent.mito"]] = percent.mito 

# Generate violin plots for quality control metrics
VlnPlot(object=sample, features=c("nFeature_RNA","nCount_RNA","percent.mito"),
        group.by="orig.ident", ncol=3, pt.size=0) # Violin plot of features
# Save the quality plots to a PDF
ggsave("plots/quality_by_origin.pdf", width=16, height=9)

# Generate scatter plots to assess relationships between feature counts
FeatureScatter(object=sample, feature1="nCount_RNA", feature2="percent.mito")
FeatureScatter(object=sample, feature1="nCount_RNA", feature2="nFeature_RNA")

# Filter out low-quality cells based on counts and mitochondrial percentage
sample = subset(x=sample, subset=nCount_RNA <= 2e4 & percent.mito <= 0.075)
# Normalize the data
sample = NormalizeData(object = sample, normalization.method = "LogNormalize", 
                       scale.factor = 10000)

# Identify variable genes and scale the data
sample = FindVariableFeatures(object=sample, selection.method="vst")
sample = ScaleData(object=sample)

# Run PCA to reduce dimensionality
sample = RunPCA(object=sample, verbose=FALSE)
# Plot PCA results, colored by identified clusters
DimPlot(sample, reduction="pca", group.by="ident")
ElbowPlot(object=sample, ndims=50) # Identify the optimal number of dimensions to use

# Run Harmony to correct for batch effects
pca_dims = 1:30
sample = RunHarmony(sample, group.by.vars="orig.ident", dims.use=pca_dims)
ElbowPlot(object=sample, ndims=50, reduction="harmony")

# Run UMAP for visualization
dims_use = 1:30
sample = RunUMAP(object=sample, reduction="harmony", dims=dims_use, verbose=FALSE)
# Build a neighbor graph for clustering
sample = FindNeighbors(object=sample, reduction="harmony", dims=dims_use, verbose=FALSE)
# Identify clusters of cells
sample = FindClusters(object=sample, resolution=1, verbose=FALSE)

# Create dimension plots to visualize clustering results
set.seed(1)
DimPlot(sample, label=TRUE, reduction="umap", group.by="ident", pt.size=0.1) + labs(title=NULL)
ggsave("plots/umap_cluster.pdf", width=7, height=6)
DimPlot(sample, label=FALSE, reduction="umap", group.by="orig.ident", pt.size=0.1, shuffle=TRUE) + labs(title=NULL)
ggsave("plots/umap_origin.pdf", width=6.6, height=6)

# Save the processed Seurat object for future use
saveRDS(file="sample_seurat5.rds", sample)
# Load the Seurat object if needed
sample = readRDS("scRNA/sample_seurat5.rds")

# Generate quality control plots by cluster
VlnPlot(object=sample, features=c("nFeature_RNA","nCount_RNA","percent.mito"), 
        ncol=3, pt.size=0)
ggsave("plots/quality_by_cluster.pdf", width=30, height=5)

# Extract metadata from the Seurat object
metadata <- sample@meta.data

# Add cluster information to the metadata if it is not already included
metadata$cluster <- Idents(sample)  # Ensure we have cluster identities

# Compute summary statistics for each cluster
summary_stats <- metadata %>%
  group_by(cluster) %>%
  summarise(
    avg_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    Q25_nFeature_RNA = quantile(nFeature_RNA, 0.25, na.rm = TRUE),
    Q75_nFeature_RNA = quantile(nFeature_RNA, 0.75, na.rm = TRUE),
    
    avg_nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    Q25_nCount_RNA = quantile(nCount_RNA, 0.25, na.rm = TRUE),
    Q75_nCount_RNA = quantile(nCount_RNA, 0.75, na.rm = TRUE),
    
    avg_percent_mito = mean(percent.mito, na.rm = TRUE),
    median_percent_mito = median(percent.mito, na.rm = TRUE),
    Q25_percent_mito = quantile(percent.mito, 0.25, na.rm = TRUE),
    Q75_percent_mito = quantile(percent.mito, 0.75, na.rm = TRUE)
  )

### NOTE
# Identified low quality cluster: Remove specific clusters based on calculated statistics:
# clusters 29, 10, 30, 14, 13, 7 

# Begin analysis section
library(ggplot2)

# Define cluster and cell type mappings for visualization
sample = readRDS("scRNA/sample_seurat5.rds")
cluster = 0:30
celltype = c("Endothelial", "Smooth muscle", "Endothelial", "Smooth muscle", 
             "Fibroblast", "Smooth muscle", "Smooth muscle", "Remove", 
             "Proliferating", "Fibroblast", "Remove", "Fibroblast", 
             "Fibroblast", "Remove", "Remove", "Myeloid", "Smooth muscle", 
             "T cell", "Endothelial", "Endothelial", "Proliferating", 
             "Proliferating", "Smooth muscle", "Myeloid", "Proliferating", 
             "Endothelial", "Fibroblast", "Endothelial", "B cell", 
             "Remove", "Remove")

# Map cell types to the Seurat object based on their cluster IDs
sample$celltype = plyr::mapvalues(x=Idents(sample), from=cluster, to=celltype)

# Define clusters to keep based on biological relevance
clusters_to_keep <- c(0:6, 8, 9, 11, 12, 15:28)

# Subset the Seurat object to retain only specified clusters
filtered_sample <- subset(sample, idents = clusters_to_keep)

# Generate UMAP plots color-coded by cell type
DimPlot(filtered_sample, label=TRUE, reduction="umap", group.by="celltype", pt.size=0.1) + labs(title=NULL)
ggsave("/scRNA/plotFigures/supplementary/umap_celltype_beginUMAP.svg", width=7.5, height=6)

# Extract expression data for specific genes of interest
pecam1_data <- FetchData(filtered_sample, vars = "Pecam1")
# Plot a histogram of gene expression
ggplot(pecam1_data, aes(x = Pecam1)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
  geom_vline(xintercept = 2, color = "red", linetype = "dashed", size = 1) + 
  labs(x = "Expression Level", y = "Frequency") +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("scRNA/plotFigures/supplementary/PecamHisto.svg", width=6.5, height=6)

# Repeat the process for additional genes
pecam1_data <- FetchData(filtered_sample, vars = "Dcn")
ggplot(pecam1_data, aes(x = Dcn)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
  geom_vline(xintercept = 2, color = "red", linetype = "dashed", size = 1) + 
  labs(x = "Expression Level", y = "Frequency") +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("scRNA/plotFigures/supplementary/DcnHisto.svg", width=6.5, height=6)

# Continue with additional genes
pecam1_data <- FetchData(sample, vars = "Prox1")
ggplot(pecam1_data, aes(x = Prox1)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", size = 1) + 
  labs(title = "Histogram of Prox1 Expression", x = "Expression Level", y = "Frequency") +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("scRNA/plotFigures/supplementary/Prox1Histo.svg", width=6.5, height=6)

# Repeat for additional genes
pecam1_data <- FetchData(filtered_sample, vars = "Tagln")
ggplot(pecam1_data, aes(x = Tagln)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
  geom_vline(xintercept = 2, color = "red", linetype = "dashed", size = 1) + 
  labs(title = "Histogram of Tagln Expression", x = "Expression Level", y = "Frequency") +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), panel.background=element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("scRNA/plotFigures/supplementary/TaglnHisto.svg", width=6.5, height=6)

# Continue for additional genes
pecam1_data <- FetchData(sample, vars = "Cd3e")
ggplot(pecam1_data, aes(x = Cd3e)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", size = 1) + 
  labs(title = "Histogram of Cd3e Expression", x = "Expression Level", y = "Frequency") +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(),
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("scRNA/plotFigures/supplementary/Cd3eHisto.svg", width=6.5, height=6)

# Add binary columns based on gene expression thresholds for specific genes
filtered_sample$Pecam_pos <- ifelse(FetchData(filtered_sample, vars = "Pecam1") > 2, 1, 0)
filtered_sample$Dcn_pos <- ifelse(FetchData(filtered_sample, vars = "Dcn") > 2, 1, 0)
filtered_sample$Prox1_pos <- ifelse(FetchData(filtered_sample, vars = "Prox1") > 0.5, 1, 0)
filtered_sample$Tagln_pos <- ifelse(FetchData(filtered_sample, vars = "Tagln") > 2, 1, 0)
filtered_sample$Cd3e_pos <- ifelse(FetchData(filtered_sample, vars = "Cd3e") > 0.5, 1, 0)

# Calculate the total counts of positive features
filtered_sample$Positive_Sum <- filtered_sample$Pecam_pos + filtered_sample$Dcn_pos + 
  filtered_sample$Prox1_pos + filtered_sample$Tagln_pos + 
  filtered_sample$Cd3e_pos

# Plot the summed positive features on UMAP
FeaturePlot(filtered_sample, features = "Positive_Sum", reduction = "umap") +
  scale_color_viridis_c(option = "magma")
ggsave("scRNA/plotFigures/supplementary/UMAP_filtering_doublets.svg", width=6.5, height=6)

# Create a new metadata column based on the sum threshold
filtered_sample$Sum_Status_2 <- ifelse(filtered_sample$Positive_Sum >= 2, "Remove", "Keep")

# Generate a UMAP plot colored by the sum status
DimPlot(
  filtered_sample, 
  group.by = "Sum_Status_2", 
  reduction = "umap", 
  cols = c("Keep" = "blue", "Remove" = "red")
) +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(),
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("scRNA/plotFigures/supplementary/umap_doublets_filterinv.svg", width=6.5, height=6)

# Subset the samples to retain only those labeled as "Keep"
filtered_sample2 <- subset(filtered_sample, subset = Sum_Status_2 == "Keep")

# Plot the UMAP for filtered sample colored by cell type
DimPlot(filtered_sample2, label=TRUE, reduction="umap", group.by="celltype", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + NoLegend()
ggsave("scRNA/plotFigures/supplementary/umap_cellType_afterQualityCheck_removalAndDoubletDetection.svg", width=7.5, height=6)

# Plot UMAP for filtered sample colored by original identification
DimPlot(filtered_sample2, reduction="umap", group.by="orig.ident", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
ggsave("scRNA/plotFigures/supplementary/umap_age_afterQualityCheck_removalAndDoubletDetection.svg", width=7.5, height=6)

# Define genes of interest for dot plots
gene_list_test <- c("Mki67","Pcna","Pclaf", "Cdk1","Cks2",
                    "Cdh5","Pecam1","Erg", "Fli1", "Etv2",
                     "Dcn")

genes_flipped <- rev(gene_list_test)
# Generate dot plots for specified genes
DotPlot(filtered_sample2, features = genes_flipped, cols = c("grey90", "red")) +
  labs(x = "Gene", y = "Cluster") +
  coord_flip()
ggsave("scRNA/plotFigures/main/panel_G/proliferating_dotplot_specific_genes_by_cluster.svg", width=8, height=8)

# Assign cell types to clusters based on a provided mapping
cluster = 0:30
celltypes = c("Endothelial","Smooth muscle","Endothelial","Smooth muscle","Fibroblast",
              "Smooth muscle","Smooth muscle","Remove","Proliferating","Fibroblast","Remove",
              "Fibroblast","Fibroblast","Remove","Remove","Myeloid","Smooth muscle","T cell",
              "Endothelial","Endothelial","Proliferating","Proliferating","Smooth muscle",
              "Myeloid","Proliferating","Endothelial","Fibroblast","Endothelial","B cell",
              "Remove","Remove")

proliferating$celltype = plyr::mapvalues(x=Idents(proliferating), from=cluster, to=celltypes)

# Define the clusters to keep based on biological relevance
clusters_to_keep <- c(0:6, 8, 9, 11, 12, 15:28)

# Subset the proliferating sample to include only these clusters
filtered_sample <- subset(proliferating, idents = clusters_to_keep)

# Plot the UMAP data highlighting the cell types
DimPlot(filtered_sample, label=TRUE, reduction="umap", group.by="celltype", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18), axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + NoLegend()
ggsave("scRNA/plotFigures/supplementary/umap_cellType_afterQualityCheck_removalAndDoubletDetection.svg", width=7.5, height=6)


# Run PCA and UMAP for the endothelial cells
endothelial <- subset(filtered_sample2, subset = celltype == "Endothelial")

# Prepare cell counts for differential expression analysis
endothelial = NormalizeData(object = endothelial, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
endothelial = FindVariableFeatures(object=endothelial, selection.method="vst")
endothelial = ScaleData(object=endothelial)

# Run PCA
endothelial = RunPCA(object=endothelial, verbose=FALSE)
DimPlot(endothelial, reduction="pca", group.by="ident")
ElbowPlot(object=endothelial, ndims=50)
ggsave("scRNA/plotFigures/supplementary/ElbowPCA_endothelialCells.svg", width=6.5, height=6)

# Run Harmony to correct for batch effects
pca_dims = 1:30
endothelial = RunHarmony(endothelial, group.by.vars="orig.ident", dims.use=pca_dims)
ElbowPlot(object=endothelial, ndims=50, reduction="harmony")
ggsave("scRNA/plotFigures/supplementary/ElbowHarmony_endothelialCells.svg", width=6.5, height=6)

# Run UMAP and perform clustering
dims_use = 1:30
endothelial = RunUMAP(object=endothelial, reduction="harmony", dims=dims_use, verbose=FALSE)
endothelial = FindNeighbors(object=endothelial, reduction="harmony", dims=dims_use, verbose=FALSE)
endothelial = FindClusters(object=endothelial, resolution=0.2, verbose=FALSE)

# Plot UMAP results grouped by original identifiers
DimPlot(endothelial, reduction="umap", group.by="orig.ident", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
ggsave("scRNA/plotFigures/supplementary/umap_age_afterQualityCheck_removalAndDoubletDetection_endothelial.svg", width=7.5, height=6)

DimPlot(endothelial, label=TRUE, reduction="umap", group.by="ident", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + NoLegend()
ggsave("scRNA/plotFigures/supplementary/umap_ident_afterQualityCheck_removalAndDoubletDetection_endothelial.svg", width=7.5, height=6)

# Define genes of interest
genes = c("Pecam1", "Sox17", "Prox1", "Cd36")
for (i in 1:length(genes)){
  FeaturePlot(endothelial, features=genes, reduction = "umap") +
    scale_colour_gradientn(colours = c("#ededed", pal_material("pink", n = 30, alpha = 1, reverse = FALSE)(30))) +
                    theme(axis.title.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.x = element_blank(),
                          axis.text.y = element_blank(),
                          plot.title = element_blank(),
                          axis.line = element_blank(),
                          axis.ticks = element_blank())
                    ggsave(sprintf("scRNA/plotFigures/main/panel_K/umap_endothelial_%s.svg", genes[i]), width=6.5, height=6)
                    }

# Identify markers for each age group
Idents(endothelial) <- endothelial$orig.ident
markers <- FindAllMarkers(
  object = endothelial,
  group.by = "orig.ident",
  only.pos = TRUE,       # Only consider upregulated genes
  min.pct = 0.1,        # Filter genes expressed in at least 10% of cells
  logfc.threshold = 0.75 # Minimum log fold change
)

# Select the top genes for each age group
top_genes <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 15) # Top 15 genes per group

# Plot heatmap of the differentially expressed genes
DoHeatmap(
  object = endothelial,
  features = top_genes$gene,    # Use the top differentially expressed genes
  group.by = "orig.ident"       # Group cells by age (orig.ident)
) +
  theme_minimal()
ggsave("scRNA/plotFigures/main/panel_L/DE_genes_endothelial.png", width=6, height=12)

# Extract total counts per age and Ki67 status
data <- FetchData(filtered_sample2, vars = c("orig.ident", "celltype"))

# Calculate percentages of each cell type per age group
composition <- data %>%
  group_by(orig.ident, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(percentage = count / sum(count) * 100)

# Create a stacked bar plot to show cell type composition by age
ggplot(composition, aes(x = orig.ident, y = percentage, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Age", y = "Percentage", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 22),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18)) +
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank(), 
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("scRNA/plotFigures/main/panel_A/cell_type_composition_byAge.svg", width = 8, height = 6)

# Save total counts for later analysis
write.csv(total_counts, file = "scRNA/plotFigures/main/total_counts_ECs.csv", row.names = FALSE)

# Save the percentage data for proliferating cells
write.csv(data_long, file = "scRNA/plotFigures/main/percentageProliferating_endothelialCells.csv", row.names = FALSE)

# Save summary statistics including cell counts
write.csv(data, file = "scRNA/plotFigures/main/endothelialCells_withProliferationNumbers.csv", row.names = FALSE)

# Extract proliferative endothelial cells
endothelial <- subset(filtered_sample2, subset = celltype == "Endothelial")
table(endothelial$orig.ident)

prolif_ECs_4 <- subset(proliferating, idents = 4)
prolif_ECs_4$Pecam_pos <- ifelse(FetchData(prolif_ECs_4, vars = "Pecam1") > 2, 1, 0)  # Binary classification based on expression
prolif_ECs_pos_4 <- subset(prolif_ECs_4, subset = Pecam_pos == 1)
table(prolif_ECs_pos_4$orig.ident)  # Summarize the counts of original identifiers

# Continue the analysis for further genes of interest and visualization...