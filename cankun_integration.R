library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)
library(qs)
library(Polychrome)
#setwd("D:/Project/mingtao")

here::set_here()
message(paste("Current working directory:", here::here()))

# Load data and create objects for each dataset
A1.data <-
  Read10X_h5("data/Con0_CKDL200167803-1a-SI_GA_D3_HNNKFDSXY/filtered_feature_bc_matrix.h5")
A2.data <-
  Read10X_h5("data/Con2_CKDL210000544-1a-SI_GA_B6_HNNKFDSXY/filtered_feature_bc_matrix.h5")
A3.data <-
  Read10X_h5("data/Con5_CKDL200167807-1a-SI_GA_A3_HNNKFDSXY/filtered_feature_bc_matrix.h5")
A4.data <-
  Read10X_h5("data/Con10_CKDL200167809-1a-SI_GA_A7_HNNKFDSXY/filtered_feature_bc_matrix.h5")
A5.data <-
  Read10X_h5("data/Con14_CKDL200167811-1a-SI_GA_A9_HNNKFDSXY/filtered_feature_bc_matrix.h5")
A6.data <-
  Read10X_h5(
    "data/Con30_CKDL200167813-1a-SI_GA_A11_HNNKFDSXY/filtered_feature_bc_matrix.h5"
  )
A7.data <-
  Read10X_h5("data/N1KO0_CKDL200167804-1a-SI_GA_E3_HNNKFDSXY/filtered_feature_bc_matrix.h5")
A8.data <-
  Read10X_h5("data/N1KO2_CKDL210000545-1a-SI_GA_B7_HNNKFDSXY/filtered_feature_bc_matrix.h5")
A9.data <-
  Read10X_h5("data/N1KO5_CKDL200167808-1a-SI_GA_A4_HNNKFDSXY/filtered_feature_bc_matrix.h5")
A10.data <-
  Read10X_h5(
    "data/N1KO10_CKDL200167810-1a-SI_GA_A8_HNNKFDSXY/filtered_feature_bc_matrix.h5"
  )
A11.data <-
  Read10X_h5(
    "data/N1KO14_CKDL200167812-1a-SI_GA_A10_HNNKFDSXY/filtered_feature_bc_matrix.h5"
  )
A12.data <-
  Read10X_h5(
    "data/N1KO30_CKDL200167814-1a-SI_GA_A12_HNNKFDSXY/filtered_feature_bc_matrix.h5"
  )


A1 <-
  CreateSeuratObject(
    A1.data,
    project = "Con0",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )
A2 <-
  CreateSeuratObject(
    A2.data,
    project = "Con2",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )
A3 <-
  CreateSeuratObject(
    A3.data,
    project = "Con5",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )
A4 <-
  CreateSeuratObject(
    A4.data,
    project = "Con10",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )
A5 <-
  CreateSeuratObject(
    A5.data,
    project = "Con14",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )
A6 <-
  CreateSeuratObject(
    A6.data,
    project = "Con30",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )
A7 <-
  CreateSeuratObject(
    A7.data,
    project = "N1KO0",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )
A8 <-
  CreateSeuratObject(
    A8.data,
    project = "N1KO2",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )
A9 <-
  CreateSeuratObject(
    A9.data,
    project = "N1KO5",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )
A10 <-
  CreateSeuratObject(
    A10.data,
    project = "N1KO10",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )
A11 <-
  CreateSeuratObject(
    A11.data,
    project = "N1KO14",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )
A12 <-
  CreateSeuratObject(
    A12.data,
    project = "N1KO30",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )



# Merge all 4 datasets
cord.big <- merge(
  A1,
  y = c(A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12),
  add.cell.ids = c(
    "Con0",
    "Con2",
    "Con5",
    "Con10",
    "Con14",
    "Con30",
    "N1KO0",
    "N1KO2",
    "N1KO5",
    "N1KO10",
    "N1KO14",
    "N1KO30"
  ),
  project = "Zhao_NOTCH1"
)


head(colnames(cord.big))

combine.list <- SplitObject(cord.big, split.by = "orig.ident")
combine.list <- lapply(
  X = combine.list,
  FUN = function(x) {
    x <- NormalizeData(x)
    x <-
      FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  }
)
combine.anchors <-
  FindIntegrationAnchors(object.list = combine.list, dims = 1:30)
combined <- IntegrateData(anchorset = combine.anchors, dims = 1:30)

custom_color <-
  as.character(palette36.colors(36)[-2])[1:length(levels(Idents(combined)))]

Idents(combined) <- combined$orig.ident
VlnPlot(combined, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA"), group.by = "orig.ident", 
        ncol = 2, pt.size = 0.1)

DefaultAssay(combined) <- "RNA"
rb.genes <- rownames(combined)[grep("^RP[SL][[:digit:]]",rownames(combined))]
percent.ribo <- colSums(combined[rb.genes,])/Matrix::colSums(combined)*100
combined <- AddMetaData(combined, percent.ribo, col.name = "percent_ribo")
VlnPlot(combined, features = "percent.ribo", pt.size = 0.1) + NoLegend()

combined <- PercentageFeatureSet(combined, "^MT-", col.name = "percent_mito")
VlnPlot(combined, features = "percent_mito", pt.size = 0.1) + NoLegend()

total_counts_per_cell <- colSums(combined@assays$RNA@counts)
mito_genes <- rownames(combined)[grep("^MT-", rownames(combined))]
combined$percent_mito <- colSums(combined@assays$RNA@counts[mito_genes, ])/total_counts_per_cell

DefaultAssay(combined) <- "integrated"

FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


Idents(combined) <- combined$orig.ident

condition_idents <- as.factor(combined$orig.ident)
levels(condition_idents) <- c(rep("ctrl",6), rep("N1KO",6))
combined <- AddMetaData(combined, condition_idents, col.name = "condition")


combined <- CellCycleScoring(object = combined, g2m.features = cc.genes$g2m.genes, 
                              s.features = cc.genes$s.genes)
  
VlnPlot(combined, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
        ncol = 2, pt.size = 0.1)

# After integration analysis

DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30, k.param = 35)
combined <- FindClusters(combined, resolution = 0.3)


Idents(combined) <- combined$seurat_clusters
p1 <- DimPlot(combined, reduction = "umap", label = TRUE)
p1

Idents(combined) <- combined$orig.ident
p2 <- DimPlot(combined, reduction = "umap", label = TRUE)
p2

combined2 <- subset(combined, subset = nCount_RNA < 200000 & percent_mito < 0.25)

VlnPlot(combined2, features = c("nFeature_RNA", "nCount_RNA","percent_mito","percent.ribo"), group.by = "orig.ident", 
        ncol = 2, pt.size = 0.1)


DefaultAssay(combined) <- "integrated"
qs::qsave(combined, "combined.qsave") 

VlnPlot(combined, features = c("NOTCH1"), group.by = "orig.ident", 
        ncol = 2, pt.size = 0.1)


combined <- qs::qread('combined.qsave')
library(plotly)

umap_embed <- 1
fig <- plot_ly(mtcars, x = ~wt, y = ~hp, z = ~qsec, color = ~am, colors = c('#BF382A', '#0C4B8E'))


######### 3D UMAP

# Interacive multimodal 3D UMAP plotting of scRNA sequencing datasets
# The following is a length of code generated to create nice 
# 3D UMAP plots of seurat v3.0.0-v3.1.1 objects utilizing the visualization 
# package plot_ly

# R v3.5.3 (x64 bit) and RStudio v1.2.1335 (x64 bit) were used for running this code :)

# Seurat is a multimodal single Cell RNA seq analysis algorithm created by
# The Satija Lab. For more information please see: https://satijalab.org/seurat/

# Contributors (by their Github handles):
# @Dragonmasterx87 (Dept. of Cell Biology, UM)
# @msaadsadiq (Dept. of Electrical and Computer Engineering, UM)

# Install plot_ly
#install.packages('plotly')

# Load plot_ly
library(plotly)

# Construct a dataframe using data from your pre-clustered Seurat v3.1.1 object
# Here 'seurat_clusters' is list of numeric cluster identities, you can find it here: yourseuratobject[["seurat_cluster"]], 
# or yourseuratobject$seurat_clusters, where 'yourseuratobject' is a Seurat object created with Seurat v3.1.1 (works for v3.0.0 as well)
yourseuratobject <- combined

# Re-run UMAPs that you have accurate calculations for all UMAP(s)
yourseuratobject <- RunUMAP(yourseuratobject,
                            dims = 1:10,
                            n.components = 3L)

# Extract tSNE information from Seurat Object
umap_1 <- yourseuratobject[["umap"]]@cell.embeddings[,1]
umap_2 <- yourseuratobject[["umap"]]@cell.embeddings[,2]
umap_3 <- yourseuratobject[["umap"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = yourseuratobject, reduction = "umap")

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
fig <- plot_ly(data = plot.data, 
               x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
               color = ~seurat_clusters, 
               colors = custom_color,
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 5, width=2), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names


# Updates stemming from Issue #9 Having a fixed scale on axes while selecting particular clusters
# @rtoddler thanks for the suggestions!
# Before you plot, set the ranges of the axis you desire. This set axis range will be 
# present across all clusters, and plotly will not adjust for axis length anymore
# this axis length will persist even when selecting some clusters

# xaxis
axx <- list(
  nticks = 4,
  range = c(-10,10) #select range of xaxis
)

# yaxis
axy <- list(
  nticks = 4,
  range = c(-10,10) #select range of yaxis
)

#zaxis
axz <- list(
  nticks = 4,
  range = c(-10,10) #select range of zaxis
)

fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
fig_cube <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, aspectmode='cube')) # To maintain cubic aspect
fig
fig_cube

# Say you wanto make a gene-expression 3D plot, where you can plot gene expression against a color scale
# Here using the same seurat object as above, we extract gene expression information for beta-actin 'ACTB'
# Here we concentrate on SCT normalized data, or log normalized RNA NOT raw counts.
# In addition if you want, you may look at normalised-RNA, SCT or integrated slots, to look at gene expression
# Setting your DefaultAssay() will inform R which assay to pick up expression data from.
DefaultAssay(object = yourseuratobject)
DefaultAssay(object = yourseuratobject) <- "RNA"
DefaultAssay(object = yourseuratobject) <- "integrated"
DefaultAssay(object = yourseuratobject) <- "SCT"

# create a dataframe
plot.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "ACTB"), slot = 'data')

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'changed' column of your dataframe
plot.data$changed <- ifelse(test = plot.data$ACTB <1, yes = plot.data$ACTB, no = 1)

# Add the label column, so that now the column has 'cellname-its expression value'
plot.data$label <- paste(rownames(plot.data)," - ", plot.data$ACTB, sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgreen', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5, width=2), 
        text=~label,
        hoverinfo="text"
)

# On running this code the HTML output should appear in RStudio. You can save the output as a
# HTML file. Once you have saved, just open the HTML file in any web browser (double click on the html- file
# and if asked select to open with any web browser like google chrome/safari/mozilla/explorer etc).
# It should be have all of the integrated features you saw in the RStudio output file.

########## #
########## #

# Alternative method as designed by @vertesy (Thanks for the suggestions!)
# create a dataframe
goi <- "TOP2A"
plotting.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Expression"=goi), slot = 'data')

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'Expression' column of your dataframe
# Cutoff <- 2
Cutoff <- quantile(plotting.data[,goi], probs = .95)
plotting.data$"ExprCutoff" <- ifelse(test = plotting.data[,goi] <Cutoff, yes = plotting.data[,goi], no = Cutoff)

# Add the label column, so that now the column has 'cellname-its expression value'
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data[,goi], sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
plot_ly(data = plotting.data,
        # name = goi,
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~ExprCutoff, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('darkgrey', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 1), 
        text=~label,
        hoverinfo="text"
) %>%layout(title=goi)

# Thank you for reading and using this code to further your scRNAseq analysis!
# If you liked it, dont forget to acknowledge, fork and star!
# Citation information is within the Readme, please dont forget to cite!
# Have a wonderful day!!

######### Archive

origin_ident <- as.factor(combined$orig.ident)
tissue_ident <- origin_ident
timepoint_ident <- origin_ident
levels(tissue_ident) <-
  c("CNS",
    "Blood",
    "Blood",
    "CNS",
    "Blood",
    "CNS",
    "Blood",
    "CNS",
    "Blood",
    "CNS")
levels(timepoint_ident) <-
  c(
    "Preclinical",
    "Preclinical",
    "Preclinical",
    "Preclinical",
    "Onset",
    "Onset",
    "Peak",
    "Peak",
    "Late",
    "Late"
  )

combined <- AddMetaData(combined, tissue_ident, col.name = "Tissue")
combined <-
  AddMetaData(combined, timepoint_ident, col.name = "Timepoint")

Idents(combined) <- combined$Timepoint
DimPlot(combined, reduction = "umap", label = TRUE)


Idents(combined) <- combined$Tissue
DimPlot(combined, reduction = "umap", label = TRUE)

Idents(combined) <- combined$seurat_clusters
## VISION analysis
saveRDS(combined, "combined.rds")

