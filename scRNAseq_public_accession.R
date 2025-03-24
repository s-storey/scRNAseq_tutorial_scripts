###############################################################
# Title: Data Import from GEO and Seurat Object Creation
# Author: [Your Name]
# Date: [Workshop Date]
#
# Description:
#   This module demonstrates how to:
#     1. Download filtered feature-barcode matrix data from GEO using GEOquery.
#     2. Read the data into R with Seurat's Read10X_h5 function.
#     3. Create a Seurat object from the counts matrix.
#     4. Download sample-level metadata from GEO.
#     5. Map barcode suffixes to sample names and add a new metadata column.
#     6. Split the Seurat object into multiple objects (one per sample).
#
#   Each step is commented to explain its purpose and guide parameter adjustments.
###############################################################

##########################
# 0) Environment Setup
##########################

# Install and load required packages.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(GEOquery)

##########################
# 1) Set GEO Accession & Download Supplementary Files
##########################

# Define your GEO accession number (replace with your dataset)
geo_accession <- "GSE226373"  # Example accession

# Download all supplementary files for the GEO series.
# A folder named after the accession (e.g., "GSE226373") will be created.
getGEOSuppFiles(geo_accession)

##########################
# 2) Define File Path & Read the Expression Matrix
##########################

# Define the path to the filtered_feature_bc_matrix.h5 file.
# Ensure the filename exactly matches the supplementary file name.
h5_file_path <- file.path(geo_accession, "GSE226373_filtered_feature_bc_matrix.h5")

# Check if the file exists
if (!file.exists(h5_file_path)) {
  stop("The file filtered_feature_bc_matrix.h5 was not found in the folder ", geo_accession)
}

# Read the HDF5 file using Seurat's Read10X_h5 function.
counts_data <- Read10X_h5(filename = h5_file_path)

# Create a Seurat object from the counts matrix.
seurat_object <- CreateSeuratObject(counts = counts_data, project = geo_accession)

##########################
# 3) Download and Extract Sample-Level Metadata from GEO
##########################

# Download GEO metadata (GSEMatrix format) for the accession.
gse <- getGEO(geo_accession, GSEMatrix = TRUE)

# Extract sample-level metadata (assuming the first element contains the metadata).
sample_metadata <- pData(gse[[1]])

# Optionally, inspect the Seurat object's metadata.
head(seurat_object@meta.data)
tail(seurat_object@meta.data)

##########################
# 4) Map Barcode Suffixes to Sample Names
##########################

# The cell barcodes include a suffix (e.g., "-1", "-2", etc.) that indicates the sample.
# Extract barcodes from the Seurat object's metadata.
barcodes <- rownames(seurat_object@meta.data)

# Use a regular expression to extract the numeric suffix from each barcode.
# For example, "AAACCCAGTATGAGCG-1" will yield "1".
suffixes <- as.numeric(gsub(".*-(\\d+)$", "\\1", barcodes))

# Check the distribution of suffixes (should be between 1 and 5).
table(suffixes)

# Create a new metadata column 'sample' by mapping the suffix index to the corresponding
# sample title in the GEO metadata. (Assumes sample_metadata$title is ordered appropriately.)
seurat_object@meta.data$sample <- sample_metadata$title[suffixes]

# Verify that the new 'sample' column has been added correctly.
head(seurat_object@meta.data)
tail(seurat_object@meta.data)

##########################
# 5) Split the Seurat Object into Individual Sample Objects
##########################

# Split the Seurat object based on the new 'sample' metadata column.
seurat_list <- SplitObject(seurat_object, split.by = "sample")

# Check the names of the resulting list (each corresponds to a unique sample).
names(seurat_list)

# Example: Access individual sample objects -- you can create individual seurat objects this way
s1_uninjured <- seurat_list[[1]]
head(s1_uninjured@meta.data)
tail(s1_uninjured@meta.data)

s2_44hpl <- seurat_list[[2]]
head(s2_44hpl@meta.data)
tail(s2_44hpl@meta.data)




