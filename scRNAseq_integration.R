###############################################################
# Title: Integration of scRNA-seq Samples Using CCA-Based Integration
# Author: [Your Name]
# Date: [Workshop Date]
#
# Description:
#   This module demonstrates how to integrate multiple scRNA-seq samples
#   into a single Seurat object using a CCA-based integration approach.
#
#   In this workflow, we:
#     1. Split the RNA assay by sample (using the 'sample' metadata column).
#     2. Process the unintegrated data (normalization, variable feature 
#        selection, scaling, PCA, and preliminary clustering/UMAP).
#     3. Integrate the split layers using the IntegrateLayers function with
#        the CCAIntegration method.
#     4. Re-join the RNA assay layers post-integration.
#     5. Perform clustering and run UMAP on the integrated data.
#     6. Visualize the results to assess sample mixing and cluster consistency.
#
#   Discussion on Integration Methods:
#     - **CCA-Based Integration:** Works well when samples contain similar cell
#       types and batch effects are moderate. It uses canonical correlation 
#       analysis to align shared sources of variation.
#     - **SCTransform Integration:** Uses variance-stabilizing transformation
#       to normalize data and may further reduce technical noise.
#     - **Other Methods:** Tools like Harmony or MNN are alternatives depending
#       on the dataset's heterogeneity and the nature of batch effects.
###############################################################

##########################
# Step 1: Split the RNA Assay by Sample
##########################

View(seurat_object@meta.data)

# The RNA assay is split by the 'sample' metadata, creating separate layers for
# each sample (e.g., 5 samples). This is required for the integration pipeline.
seurat_object[["RNA"]] <- split(seurat_object[["RNA"]], f = seurat_object$sample)

# Confirm the RNA assay is now split into multiple layers.
print(seurat_object)

##########################
# Step 2: Process Unintegrated Data
##########################

# Run the standard preprocessing pipeline on the unintegrated data.
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object))

# Perform preliminary clustering and UMAP visualization based on the PCA reduction.
seurat_object <- FindNeighbors(seurat_object, dims = 1:30, reduction = "pca")
seurat_object <- FindClusters(seurat_object, resolution = 2, cluster.name = "unintegrated_clusters")
seurat_object <- RunUMAP(seurat_object, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# Visualize the unintegrated data. This helps assess differences between samples prior to integration.
DimPlot(seurat_object, reduction = "umap.unintegrated", group.by = c("sample", "seurat_clusters"))

##########################
# Step 3: Integrate the Data Using CCA-Based Integration
##########################

# Integrate the split layers using the CCAIntegration method.
# 'orig.reduction' is set to "pca" from the pre-integration step,
# and the new integrated reduction will be stored as "integrated.cca".
seurat_object <- IntegrateLayers(object = seurat_object, 
                                 method = CCAIntegration, 
                                 orig.reduction = "pca", 
                                 new.reduction = "integrated.cca",
                                 verbose = FALSE)

##########################
# Step 4: Re-join the RNA Assay Layers Post-Integration
##########################

# After integration, we re-join the split RNA assay layers into a single assay.
seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])

##########################
# Step 5: Downstream Analysis on Integrated Data
##########################

# Use the integrated reduction "integrated.cca" for further analysis.
seurat_object <- FindNeighbors(seurat_object, reduction = "integrated.cca", dims = 1:30)
seurat_object <- FindClusters(seurat_object, resolution = 1)

# Run UMAP using the integrated reduction to visualize the integrated data.
seurat_object <- RunUMAP(seurat_object, dims = 1:30, reduction = "integrated.cca")

##########################
# Step 6: Visualization of Integrated Data
##########################

# Visualize the integrated UMAP.
# In this example, we assume the metadata includes columns such as 'sample'
# (e.g., treatment conditions) and 'seurat_annotations' for cluster annotations.
DimPlot(seurat_object, reduction = "umap", group.by = c("sample", "seurat_annotations"))

###############################################################
# End of Integration Module
###############################################################