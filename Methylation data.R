if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("GenomicRanges")
BiocManager::install("SummarizedExperiment")
library(TCGAbiolinks)
library(GenomicRanges)
library(SummarizedExperiment)

# Query for specific genes
query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value"
)

# Download and prepare the data
GDCdownload(query)
data <- GDCprepare(query)

# Load necessary libraries
library(SummarizedExperiment)

# Define your genes of interest
genes_of_interest <- c("SLC8A1", "SLC8A3")

# Get row names (gene symbols) from the data
gene_symbols <- rowData(data)$external_gene_name

# Find indices of the genes of interest
gene_indices <- which(gene_symbols %in% genes_of_interest)

# Extract methylation data for these genes
filtered_methylation_data <- assay(data)[gene_indices, ]

# Extract gene names for filtered data
filtered_gene_names <- gene_symbols[gene_indices]

# Create a DataFrame with the filtered data and gene names
filtered_data_df <- data.frame(Gene = filtered_gene_names, t(filtered_methylation_data))

# Save the filtered data to a CSV file
write.csv(filtered_data_df, "filtered_methylation_data.csv", row.names = FALSE)

