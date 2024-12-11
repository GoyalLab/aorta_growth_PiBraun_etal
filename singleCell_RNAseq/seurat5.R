library(Seurat)
library(harmony)
library(dplyr)
library(Matrix)
require("ggrepel")

# Create the plots director
dir.create("plots")

# Read in the expression data
all_expr = readRDS("../aggr_matrix/aggr_matrix.rds")

# Make Seurat object
sample = CreateSeuratObject(counts=all_expr, min.cells=1, min.features=100)
dim(sample)
table(sample$orig.ident)

# Count Mito genes and make plots
mito.features = grep(pattern="^mt-", x=rownames(x=sample), value=TRUE)
percent.mito = Matrix::colSums(x=GetAssayData(object=sample, slot="counts")[mito.features,]) /
  Matrix::colSums(x=GetAssayData(object=sample, slot="counts"))
sample[["percent.mito"]] = percent.mito
VlnPlot(object=sample, features=c("nFeature_RNA","nCount_RNA","percent.mito"),
        group.by="orig.ident", ncol=3, pt.size=0)
ggsave("plots/quality_by_origin.pdf", width=16, height=9)
FeatureScatter(object=sample, feature1="nCount_RNA", feature2="percent.mito")
FeatureScatter(object=sample, feature1="nCount_RNA", feature2="nFeature_RNA")

# Filter cells and Normalize data
sample = subset(x=sample, subset=nCount_RNA <= 2e4 & percent.mito <= 0.075)
sample = NormalizeData(object = sample, normalization.method = "LogNormalize", 
                       scale.factor = 10000)

# Find variable Genes and scale data by number of UMIs and Mito gene percentage
sample = FindVariableFeatures(object=sample, selection.method="vst")
sample = ScaleData(object=sample)

# Run PCA
sample = RunPCA(object=sample, verbose=FALSE)
DimPlot(sample, reduction="pca", group.by="ident")
ElbowPlot(object=sample, ndims=50)

# Run Harmony
pca_dims = 1:30
sample = RunHarmony(sample, group.by.vars="orig.ident", dims.use=pca_dims)
ElbowPlot(object=sample, ndims=50, reduction="harmony")

# Run UMAP and clustering
dims_use = 1:30
sample = RunUMAP(object=sample, reduction="harmony", dims=dims_use, verbose=FALSE)
sample = FindNeighbors(object=sample, reduction="harmony", dims=dims_use, verbose=FALSE)
sample = FindClusters(object=sample, resolution=1, verbose=FALSE)

# Make DM plots
set.seed(1)
DimPlot(sample, label=TRUE, reduction="umap", group.by="ident", pt.size=0.1) + labs(title=NULL)
ggsave("plots/umap_cluster.pdf", width=7, height=6)
DimPlot(sample, label=FALSE, reduction="umap", group.by="orig.ident", pt.size=0.1, shuffle=TRUE) + labs(title=NULL)
ggsave("plots/umap_origin.pdf", width=6.6, height=6)

# Save the Seurat object
saveRDS(file="sample_seurat5.rds", sample)
sample = readRDS("/Users/jones/Downloads/scRNA/sample_seurat5.rds")

# Quality by cluster
VlnPlot(object=sample, features=c("nFeature_RNA","nCount_RNA","percent.mito"), 
        ncol=3, pt.size=0)
ggsave("plots/quality_by_cluster.pdf", width=30, height=5)



library(dplyr)
sample = readRDS("/Users/jones/Downloads/scRNA/sample_seurat5.rds")
# Extract metadata
metadata <- sample@meta.data

# Add cluster information if not already in the metadata
metadata$cluster <- Idents(sample)  # Replace 'ident' with the cluster identity if not set

# Compute summary statistics
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


###new
# remove following clusters based on statistics:
#29,10,30,14,13,7


########
#####this is where the new magic starts
library(ggplot2)

sample = readRDS("/Users/jones/Downloads/scRNA/sample_seurat5.rds")
cluster = 0:30
celltype = c("Endothelial","Smooth muscle","Endothelial","Smooth muscle","Fibroblast",
             "Smooth muscle","Smooth muscle","Remove","Proliferating","Fibroblast","Remove",
             "Fibroblast","Fibroblast","Remove","Remove","Myeloid","Smooth muscle","T cell",
             "Endothelial","Endothelial","Proliferating","Proliferating","Smooth muscle",
             "Myeloid","Proliferating","Endothelial","Fibroblast","Endothelial","B cell",
             "Remove","Remove")


sample$celltype = plyr::mapvalues(x=Idents(sample), from=cluster, to=celltype)

# Define the clusters to keep
clusters_to_keep <- c(0:6, 8, 9, 11, 12, 15:28)

# Subset the Seurat object to include only these clusters
filtered_sample <- subset(sample, idents = clusters_to_keep)

DimPlot(filtered_sample, label=TRUE, reduction="umap", group.by="celltype", pt.size=0.1) + labs(title=NULL)
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_celltype_beginUMAP.svg", width=7.5, height=6)

# Extract expression data for PECAM1
pecam1_data <- FetchData(filtered_sample, vars = "Pecam1")
# Plot a histogram
ggplot(pecam1_data, aes(x = Pecam1)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
  geom_vline(xintercept = 2, color = "red", linetype = "dashed", size = 1) + 
  labs(x = "Expression Level", y = "Frequency") +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/PecamHisto.svg", width=6.5, height=6)


pecam1_data <- FetchData(filtered_sample, vars = "Dcn")
ggplot(pecam1_data, aes(x = Dcn)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
  geom_vline(xintercept = 2, color = "red", linetype = "dashed", size = 1) + 
  labs(x = "Expression Level", y = "Frequency") +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/DcnHisto.svg", width=6.5, height=6)


pecam1_data <- FetchData(sample, vars = "Prox1")
ggplot(pecam1_data, aes(x = Prox1)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", size = 1) + 
  labs(title = "Histogram of Prox1 Expression", x = "Expression Level", y = "Frequency") +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/Prox1Histo.svg", width=6.5, height=6)

"Tagln"
pecam1_data <- FetchData(filtered_sample, vars = "Tagln")
ggplot(pecam1_data, aes(x = Tagln)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
  geom_vline(xintercept = 2, color = "red", linetype = "dashed", size = 1) + 
  labs(title = "Histogram of Tagln Expression", x = "Expression Level", y = "Frequency") +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/TaglnHisto.svg", width=6.5, height=6)


"Cd3e"
pecam1_data <- FetchData(sample, vars = "Cd3e")
ggplot(pecam1_data, aes(x = Cd3e)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dashed", size = 1) + 
  labs(title = "Histogram of Cd3e Expression", x = "Expression Level", y = "Frequency") +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/Cd3eHisto.svg", width=6.5, height=6)



# Add a new metadata column based on Pecam1 expression
filtered_sample$Pecam_pos <- ifelse(FetchData(filtered_sample, vars = "Pecam1") > 2, 1, 0)
filtered_sample$Dcn_pos <- ifelse(FetchData(filtered_sample, vars = "Dcn") > 2, 1, 0)
filtered_sample$Prox1_pos <- ifelse(FetchData(filtered_sample, vars = "Prox1") > 0.5, 1, 0)
filtered_sample$Tagln_pos <- ifelse(FetchData(filtered_sample, vars = "Tagln") > 2, 1, 0)
filtered_sample$Cd3e_pos <- ifelse(FetchData(filtered_sample, vars = "Cd3e") > 0.5, 1, 0)


# Calculate the sum of positive features and add as a new metadata column
filtered_sample$Positive_Sum <- filtered_sample$Pecam_pos + filtered_sample$Dcn_pos + filtered_sample$Prox1_pos + filtered_sample$Tagln_pos + filtered_sample$Cd3e_pos

# Plot the summed feature on UMAP
FeaturePlot(filtered_sample, features = "Positive_Sum", reduction = "umap") +
  scale_color_viridis_c(option = "magma")
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/UMAP_filtering_doublets.svg", width=6.5, height=6)


# Add a new metadata column based on the sum threshold (2)
filtered_sample$Sum_Status_2 <- ifelse(filtered_sample$Positive_Sum >= 2, "Remove", "Keep")

# Plot the UMAP
DimPlot(
  filtered_sample, 
  group.by = "Sum_Status_2", 
  reduction = "umap", 
  cols = c("Keep" = "blue", "Remove" = "red")
)  +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/UMAP_doublets_filterinv.svg", width=6.5, height=6)

filtered_sample2 <- subset(filtered_sample, subset = Sum_Status_2 == "Keep")


DimPlot(filtered_sample2, label=TRUE, reduction="umap", group.by="celltype", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())  + NoLegend()
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_celltype_afterQualityCheck_removalAndDoubletDetection.svg", width=7.5, height=6)

DimPlot(filtered_sample2, reduction="umap", group.by="orig.ident", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_age_afterQualityCheck_removalAndDoubletDetection.svg", width=7.5, height=6)


gene_list_test <- c("Mki67","Pcna","Pclaf", "Cdk1","Cks2",
                      "Cdh5","Pecam1","Erg", "Fli1", "Etv2",
                    "Dcn")

genes_flipped <- rev(gene_list_test)
DotPlot(filtered_sample2, features = genes_flipped, cols = c("grey90", "red")) +
  labs(x = "Gene", y = "Cluster") +
  coord_flip()
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/define_undefinedCLusters_withMarkers.svg", width=12, height=8)


filtered_sample2 = NormalizeData(object = filtered_sample2, normalization.method = "LogNormalize",
                          scale.factor = 10000)
# Find variable Genes and scale data by number of UMIs and Mito gene percentage
filtered_sample2 = FindVariableFeatures(object=filtered_sample2, selection.method="vst")
filtered_sample2 = ScaleData(object=filtered_sample2)

# Run PCA
filtered_sample2 = RunPCA(object=filtered_sample2, verbose=FALSE)
DimPlot(filtered_sample2, reduction="pca", group.by="ident")
ElbowPlot(object=filtered_sample2, ndims=50)
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/ELbowPlt_cleanedData_PCA.svg", width=6.5, height=6)

# Run Harmony
pca_dims = 1:30
filtered_sample2 = RunHarmony(filtered_sample2, group.by.vars="orig.ident", dims.use=pca_dims)
ElbowPlot(object=filtered_sample2, ndims=50, reduction="harmony")
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/ELbowPlt_cleanedData_Harmony.svg", width=6.5, height=6)

# Run UMAP and clustering
dims_use = 1:30
filtered_sample2 = RunUMAP(object=filtered_sample2, reduction="harmony", dims=dims_use, verbose=FALSE)
filtered_sample2 = FindNeighbors(object=filtered_sample2, reduction="harmony", dims=dims_use, verbose=FALSE)
filtered_sample2 = FindClusters(object=filtered_sample2, resolution=0.2, verbose=FALSE)



# Make DM plots
DimPlot(filtered_sample2, reduction="umap", group.by="orig.ident", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_age_afterQualityCheck_removalAndDoubletDetection_NEWPCAHARMONY.svg", width=7.5, height=6)

DimPlot(filtered_sample2, label = TRUE, reduction="umap", group.by="ident", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + NoLegend()
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_ident_afterQualityCheck_removalAndDoubletDetection_NEWPCAHARMONY.svg", width=7.5, height=6)

DimPlot(filtered_sample2, label = TRUE, reduction="umap", group.by="celltype", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + NoLegend()
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_cellType_afterQualityCheck_removalAndDoubletDetection_NEWPCAHARMONY.svg", width=7.5, height=6)


library("ggsci")
library("scales")

genes = c("Pecam1", "Cdh5", "Prox1", "Dcn", "Tagln", "Cd3e", "Mki67", "Ptprc")
for (i in 1:length(genes)){
  FeaturePlot(filtered_sample2, features=genes[i],reduction = "umap") +
    scale_colour_gradientn(colours = c("#ededed",pal_material("pink", n = 30, alpha = 1, reverse = FALSE)(30))) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          plot.title=element_blank(),
          axis.line = element_blank(),
          axis.ticks=element_blank())
  ggsave(sprintf("/Users/jones/Downloads/scRNA/plotFigures/main/umap_allFilteredCells_%s.svg", genes[i]), width=6.5, height=6)
}

endothelial <- subset(filtered_sample2, subset = celltype == "Endothelial")

# Extract expression data for PECAM1
pecam1_data <- FetchData(endothelial, vars = "Mki67")
# Plot a histogram
ggplot(pecam1_data, aes(x = Mki67)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Mki671 Expression", x = "Expression Level", y = "Frequency") +
  theme_minimal()
ggsave("/Users/jones/Downloads/scRNA/plots_new/Mki67Histo.png", width=6.5, height=6)

endothelial <- subset(filtered_sample2, subset = celltype == "Endothelial")
# Add a new metadata column based on Pecam1 expression
endothelial$ki67_pos <- ifelse(FetchData(endothelial, vars = "Mki67") > 0.25, 1, 0)


# Calculate total counts for each age and ki67 status
age_ki67_counts <- table(endothelial$orig.ident, endothelial$ki67_pos)

# Convert to a data frame for easy handling
age_ki67_df <- as.data.frame(age_ki67_counts)

# Rename columns
colnames(age_ki67_df) <- c("Age", "ki67_pos", "Count")

# View the result
age_ki67_df

table(pecam$orig.ident)


# add prolif ECs from proliferating cluster
sub_sample_prol = readRDS("/Users/jones/Downloads/scRNA/proliferating/subsample1/subsample1_seurat5.rds")
clusters = 0:4
celltypes = c("Smooth muscle","Smooth muscle","Endothelial","Smooth muscle","Fibroblast")
sub_sample_prol$celltype = plyr::mapvalues(from=clusters, to=celltypes, Idents(sub_sample_prol))


cell_mapping <- data.frame(
  Cell_ID = colnames(sub_sample_prol),       # Unique cell identifiers
  Cell_Type = sub_sample_prol$celltype      # Cell types from the subset
)
head(cell_mapping)

proliferating <- subset(filtered_sample2, subset = celltype == "Proliferating")

# Create a default column with "undefined"
proliferating$celltype_mapped <- "undefined"

# Assign cell types based on the mapping
common_cells <- intersect(cell_mapping$Cell_ID, colnames(proliferating))  # Match cells
proliferating$celltype_mapped[common_cells] <- as.character(cell_mapping$Cell_Type[match(common_cells, cell_mapping$Cell_ID)])

table(proliferating$celltype_mapped)



#RUN UMAP & PCA
# Find variable Genes and scale data by number of UMIs and Mito gene percentage
proliferating = NormalizeData(object = proliferating, normalization.method = "LogNormalize", 
                           scale.factor = 10000)
proliferating = FindVariableFeatures(object=proliferating, selection.method="vst")
proliferating = ScaleData(object=proliferating)

# Run PCA
proliferating = RunPCA(object=proliferating, verbose=FALSE)
DimPlot(proliferating, reduction="pca", group.by="ident")
ElbowPlot(object=proliferating, ndims=50)
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/ElbowPCA_proliferativeCells.svg", width=6.5, height=6)

# Run Harmony
pca_dims = 1:30
proliferating = RunHarmony(proliferating, group.by.vars="orig.ident", dims.use=pca_dims)
ElbowPlot(object=proliferating, ndims=50, reduction="harmony")
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/ElbowHarmony_proliferativeCells.svg", width=6.5, height=6)

# Run UMAP and clustering
dims_use = 1:15
proliferating = RunUMAP(object=proliferating, reduction="harmony", dims=dims_use, verbose=FALSE)
proliferating = FindNeighbors(object=proliferating, reduction="harmony", dims=dims_use, verbose=FALSE)
proliferating = FindClusters(object=proliferating, resolution=0.1, verbose=FALSE)

# Make DM plots
DimPlot(proliferating, reduction="umap", group.by="orig.ident", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_age_afterQualityCheck_removalAndDoubletDetection_proliferating.svg", width=7.5, height=6)

DimPlot(proliferating, label = TRUE, reduction="umap", group.by="ident", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + NoLegend()
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_ident_afterQualityCheck_removalAndDoubletDetection_proliferating.svg", width=7.5, height=6)

DimPlot(proliferating, reduction="umap", group.by="celltype_mapped", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_ident_afterQualityCheck_removalAndDoubletDetection_proliferating_CellTypeWithUndefined.svg", width=7.5, height=6)


genes = c("Pecam1", "Sox17", "Mki67", "Pcna", "Pclaf", "Cdk1")
for (i in 1:length(genes)){
  FeaturePlot(proliferating, features=genes[i],reduction = "umap") +
    scale_colour_gradientn(colours = c("#ededed",pal_material("pink", n = 30, alpha = 1, reverse = FALSE)(30))) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          plot.title=element_blank(),
          axis.line = element_blank(),
          axis.ticks=element_blank())
  ggsave(sprintf("/Users/jones/Downloads/scRNA/plotFigures/main/umap_prolif_%s.svg", genes[i]), width=6.5, height=6)
}


Idents(proliferating) <- proliferating$seurat_clusters

Idents(proliferating) <- proliferating$seurat_clusters
genes=c("Cdh5","Pecam1","Erg", "Kdr", "Mki67","Pcna", "Hmgb2", "Pclaf","Ube2c","Top2a","Birc5","Rrm2","Cks1b","Tuba1b","Cdk1","Cks2",
        "Notch1", "Gja4","Gja5", "Gja1", "Sox17", "Acvrl1","Nrp1","Mgp", "Fbln5", "Clec14a", "Cst3",
        "Slc38a5", "Slc16a1", "Apln", "Esm1", "Dll4", "Lyve1", "Prox1", "Dcn", "Col1a1", "Acta2","Tagln")


proliferating$clor = paste(Idents(proliferating), proliferating$orig.ident)
DotPlot(proliferating, features=genes, cols=c("grey90","red"), group.by="clor") +
  labs(x=NULL, y=NULL) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("/Users/jones/Downloads/scRNA/plots_new/dotplot_specific_genes_by_cluster_origin.png", width=12, height=4)

genes_flipped <- rev(genes)
DotPlot(proliferating, features = genes_flipped, cols = c("grey90", "red")) +
  labs(x = "Gene", y = "Cluster") +
  coord_flip()
ggsave("/Users/jones/Downloads/scRNA/plots_new/proliferating_dotplot_specific_genes_by_cluster.png", width=8, height=12)


gene_list_Luisa <- c("Mki67","Pcna","Pclaf", "Cdk1","Cks2",
                      "Cdh5","Pecam1","Erg",
                      "Dcn", "Col1a1",
                      "Acta2","Tagln",
                      "Notch1", "Gja4","Gja5", "Sox17", "Acvrl1","Nrp1")

genes_flipped <- rev(gene_list_Luisa )
DotPlot(proliferating, features = genes_flipped) +
  scale_colour_gradientn(colours = c("#ededed",pal_material("pink", n = 30, alpha = 1, reverse = FALSE)(30))) +
  labs(x = "Gene", y = "Cluster") +
  coord_flip()
ggsave("/Users/jones/Downloads/scRNA/plotFigures/main/panel_G/proliferating_dotplot_specific_genes_by_cluster.svg", width=8, height=8)






gene_list_Luisa_I <- c("Mki67","Pcna","Pclaf", "Cdk1","Cks2",
                       "Cdh5","Pecam1","Erg", "Fli1", "Etv2",
                       "Dcn", "Col1a1", "Col3a1",
                       "Acta2","Tagln", "Cnn1", "Myl6",
                       "Notch1", "Gja4","Gja5", "Sox17", "Acvrl1","Nrp1",
                       "Prox1", "Lyve1", "Esm1", "Apln", 
                       "Dll4", "Nrp2", "Nr2f2")



proliferating$clor = paste(Idents(proliferating), proliferating$orig.ident)
DotPlot(proliferating, features=gene_list_Luisa_I,group.by="clor") +
  scale_colour_gradientn(colours = c("#ededed",pal_material("pink", n = 30, alpha = 1, reverse = FALSE)(30))) +
  labs(x = "Gene", y = "Cluster") +
  coord_flip() + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("/Users/jones/Downloads/scRNA/plotFigures/main/panel_I/dotplot_specific_genes_by_clusterAndAGE_FigI.svg", width=12, height=12)


#####endothelial cells

endothelial <- subset(filtered_sample2, subset = celltype == "Endothelial")

#RUN UMAP & PCA
# Find variable Genes and scale data by number of UMIs and Mito gene percentage
endothelial = NormalizeData(object = endothelial, normalization.method = "LogNormalize", 
                              scale.factor = 10000)
endothelial = FindVariableFeatures(object=endothelial, selection.method="vst")
endothelial = ScaleData(object=endothelial)

# Run PCA
endothelial = RunPCA(object=endothelial, verbose=FALSE)
DimPlot(endothelial, reduction="pca", group.by="ident")
ElbowPlot(object=endothelial, ndims=50)
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/ElbowPCA_endothelialCells.svg", width=6.5, height=6)

# Run Harmony
pca_dims = 1:30
endothelial = RunHarmony(endothelial, group.by.vars="orig.ident", dims.use=pca_dims)
ElbowPlot(object=endothelial, ndims=50, reduction="harmony")
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/ElbowHarmony_endothelialCells.svg", width=6.5, height=6)

# Run UMAP and clustering
dims_use = 1:30
endothelial = RunUMAP(object=endothelial, reduction="harmony", dims=dims_use, verbose=FALSE)
endothelial = FindNeighbors(object=endothelial, reduction="harmony", dims=dims_use, verbose=FALSE)
endothelial = FindClusters(object=endothelial, resolution=0.2, verbose=FALSE)


# Make DM plots
DimPlot(endothelial, reduction="umap", group.by="orig.ident", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_age_afterQualityCheck_removalAndDoubletDetection_endothelial.svg", width=7.5, height=6)

DimPlot(endothelial, label = TRUE, reduction="umap", group.by="ident", pt.size=0.1) + labs(title=NULL) + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + NoLegend()
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_ident_afterQualityCheck_removalAndDoubletDetection_endothelial.svg", width=7.5, height=6)



genes = c("Pecam1", "Sox17", "Prox1", "Cd36")
for (i in 1:length(genes)){
  FeaturePlot(endothelial, features=genes[i],reduction = "umap") +
    scale_colour_gradientn(colours = c("#ededed",pal_material("pink", n = 30, alpha = 1, reverse = FALSE)(30))) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          plot.title=element_blank(),
          axis.line = element_blank(),
          axis.ticks=element_blank())
  ggsave(sprintf("/Users/jones/Downloads/scRNA/plotFigures/main/panel_K/umap_endothelial_%s.svg", genes[i]), width=6.5, height=6)
}


# Module scores on GO pathways
genes1 = read.csv("/Users/jones/Downloads/scRNA/endothelial/subsample1/GO_pathways/GO_acetyl-CoA_metabolic_process.txt", sep="\t", header=FALSE)
genes = intersect(genes1[,1], rownames(endothelial))
endothelial = AddModuleScore(endothelial, features = list(genes), name="Module_score1")
FeaturePlot(endothelial, features="Module_score11", cols=c("grey90","red"), 
            min.cutoff=0) + labs(title="Acetyl-CoA Metabolic Process")
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_Acetyl_CoA_score.svg", width=4.5, height=4)
FeaturePlot(endothelial, features="Module_score11", cols=c("grey90","red"), 
            min.cutoff=0, split.by="orig.ident", pt.size=0.1, keep.scale="all")
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_Acetyl_CoA_score_split_by_origin.svg", width=12, height=3)

data <- FetchData(endothelial, vars = c("Module_score11", "orig.ident"))
# Calculate the 75% quantile for each group
quantiles <- data %>%
  group_by(orig.ident) %>%
  summarize(quantile_90 = quantile(Module_score11, 0.9))
# View the quantiles
print(quantiles)

VlnPlot(endothelial, "Module_score11", group.by="orig.ident", pt.size=0) + 
  guides(fill = "none") +
  labs(x = NULL, title = "Extracellular Matrix") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none", plot.title=element_blank() ) +
  # Add horizontal lines for each age
  geom_hline(data = quantiles, aes(yintercept = quantile_90, color = orig.ident), size = 0.5) 
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/vln_Acetyl_CoA_score.svg", width=4, height=4)

genes1 = read.csv("/Users/jones/Downloads/scRNA/endothelial/subsample1/GO_pathways/GO_ECM.txt", sep="\t", header=FALSE)
genes = intersect(genes1[,1], rownames(endothelial))
endothelial = AddModuleScore(endothelial, features = list(genes), name="Module_score1")
FeaturePlot(endothelial, features="Module_score11", cols=c("grey90","red"), 
            min.cutoff=0) + labs(title="Extracellular Matrix")
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_ECM_score.png", width=4.5, height=4)
FeaturePlot(endothelial, features="Module_score11", cols=c("grey90","red"), 
            min.cutoff=0, split.by="orig.ident", pt.size=0.1, keep.scale="all")
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_ECM_score_split_by_origin.png", width=12, height=3)


data <- FetchData(endothelial, vars = c("Module_score11", "orig.ident"))
# Calculate the 75% quantile for each group
quantiles <- data %>%
  group_by(orig.ident) %>%
  summarize(quantile_90 = quantile(Module_score11, 0.9))
# View the quantiles
print(quantiles)

VlnPlot(endothelial, "Module_score11", group.by = "orig.ident", pt.size = 0) + 
  guides(fill = "none") +
  labs(x = NULL, title = "Extracellular Matrix") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none", plot.title=element_blank() ) +
  # Add horizontal lines for each age
  geom_hline(data = quantiles, aes(yintercept = quantile_90, color = orig.ident), size = 0.5) 
ggsave("/Users/jones/Downloads/scRNA/plotFigures/main/panel_M/vln_ECM_score_90Quant.svg", width=4, height=4)

genes1 = read.csv("/Users/jones/Downloads/scRNA/endothelial/subsample1/GO_pathways/GO_fatty_acid_metabolic_process.txt", sep="\t", header=FALSE)
genes = intersect(genes1[,1], rownames(endothelial))
endothelial = AddModuleScore(endothelial, features = list(genes), name="Module_score1")
FeaturePlot(endothelial, features="Module_score11", cols=c("grey90","red"), 
            min.cutoff=0) + labs(title="Fatty Acid Metabolic Process")
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_Fatty_Acid_score.png", width=4.5, height=4)
FeaturePlot(endothelial, features="Module_score11", cols=c("grey90","red"), 
            min.cutoff=0, split.by="orig.ident", pt.size=0.1, keep.scale="all")
ggsave("/Users/jones/Downloads/scRNA/plotFigures/supplementary/umap_Fatty_Acid_score_split_by_origin.png", width=12, height=3)


data <- FetchData(endothelial, vars = c("Module_score11", "orig.ident"))
# Calculate the 75% quantile for each group
quantiles <- data %>%
  group_by(orig.ident) %>%
  summarize(quantile_90 = quantile(Module_score11, 0.90))
# View the quantiles
print(quantiles)

VlnPlot(endothelial, "Module_score11", group.by = "orig.ident", pt.size = 0) + 
  guides(fill = "none") +
  labs(x = NULL, title = "Extracellular Matrix") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none", plot.title=element_blank() ) +
  # Add horizontal lines for each age
  geom_hline(data = quantiles, aes(yintercept = quantile_90, color = orig.ident), size = 0.5) 
ggsave("/Users/jones/Downloads/scRNA/plotFigures/main/panel_M/vln_Fatty_Acid_score_90Quant.svg", width=4, height=4)


####DE gene analysis for endos


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
  slice_max(order_by = avg_log2FC, n = 15) # Top 10 genes per group


DoHeatmap(
  object = endothelial,
  features = top_genes$gene,    # Use the top differentially expressed genes
  group.by = "orig.ident"       # Group cells by age (orig.ident)
) +
  theme_minimal()
ggsave("/Users/jones/Downloads/scRNA/plotFigures/main/panel_L/DE_genes_endothelial.png", width=6, height=12)



### cell type compostion all cells by age
# Extract relevant data
data <- FetchData(filtered_sample2, vars = c("orig.ident", "celltype"))

# Calculate percentages
composition <- data %>%
  group_by(orig.ident, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(percentage = count / sum(count) * 100)

ggplot(composition, aes(x = orig.ident, y = percentage, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    x = "Age",
    y = "Percentage",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22),
        legend.text = element_text(size=16), legend.title = element_text(size=18)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
ggsave("/Users/jones/Downloads/scRNA/plotFigures/main/panel_A/cell_type_composition_byAge.svg", width=8, height=6)


####get proliferative ECs vs non proliferative
c(1, 4)


prolif_ECs <- subset(proliferating, idents = 1)
cell_counts_by_age <- table(prolif_ECs$orig.ident)
table(endothelial$orig.ident)

library(dplyr)
library(ggplot2)
library(tidyr) 

# Create data frames from your tables
prolif_counts <- as.data.frame(table(prolif_ECs$orig.ident)) %>%
  rename(Age = Var1, Proliferating = Freq)

total_counts <- as.data.frame(table(endothelial$orig.ident) + table(prolif_ECs$orig.ident)) %>%
  rename(Age = Var1, Total = Freq)

# Merge the tables
data <- left_join(prolif_counts, total_counts, by = "Age")

# Calculate non-proliferating counts and percentages
data <- data %>%
  mutate(
    Non_Proliferating = Total - Proliferating,
    Proliferating_Percentage = Proliferating / Total * 100,
    Non_Proliferating_Percentage = Non_Proliferating / Total * 100
  )

# Reshape the data for ggplot
data_long <- data %>%
  select(Age, Proliferating_Percentage, Non_Proliferating_Percentage) %>%
  pivot_longer(cols = c(Proliferating_Percentage, Non_Proliferating_Percentage),
               names_to = "Type",
               values_to = "Percentage")

ggplot(data_long, aes(x = Percentage, y = Age, fill = Type)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = sprintf("%.1f%%", Percentage), 
        x = Percentage + 5),  # Slightly offset the text to the right
    position = position_stack(vjust = 0.5),
    size = 6
  ) +
  labs(
    x = "Percentage",
    y = "Age"
  ) +
  theme_minimal() +
  theme(axis.text=element_text(size=18),axis.title=element_text(size=22)) +
  theme(axis.line = element_line(colour = "black"), panel.border=element_blank(), legend.position = "none", 
        panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  scale_fill_manual(values = c("Proliferating_Percentage" = "#d6265e", 
                               "Non_Proliferating_Percentage" = "#ededed")) + 
  scale_x_continuous(
    breaks = seq(0, 100, by = 25)  # Set breaks every 10
  )
ggsave("/Users/jones/Downloads/scRNA/plotFigures/main/panel_H/proliferative_endothelialCells_updated.svg", width=8, height=4)



# Save total_counts to a CSV file
write.csv(total_counts, file = "/Users/jones/Downloads/scRNA/plotFigures/main/total_counts_ECs.csv", row.names = FALSE)

# Save data_long to a CSV file
write.csv(data_long, file = "/Users/jones/Downloads/scRNA/plotFigures/main/percentageProliferating_endothelialCells.csv", row.names = FALSE)


# Save data_long to a CSV file
write.csv(data, file = "/Users/jones/Downloads/scRNA/plotFigures/main/endothelialCells_withProliferationNumbers.csv", row.names = FALSE)




endothelial <- subset(filtered_sample2, subset = celltype == "Endothelial")
table(endothelial$orig.ident)

prolif_ECs_4 <- subset(proliferating, idents = 4)
prolif_ECs_4$Pecam_pos <- ifelse(FetchData(prolif_ECs_4, vars = "Pecam1") > 2, 1, 0)
prolif_ECs_pos_4 <- subset(prolif_ECs_4, subset = Pecam_pos == 1)
table(prolif_ECs_pos_4$orig.ident)


