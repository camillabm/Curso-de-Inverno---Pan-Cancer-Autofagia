######################
### Data preparing ###
######################

# Install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("SummarizedExperiment", "biomaRt"))
install.packages("dplyr")

library(SummarizedExperiment)
library(dplyr)
library (biomaRt)

# GBM
  # Metadata
metadata <- as.data.frame(colData(data_gbmR)) #Extract metadata
metadata$project <- "TCGA-GBM"  #Add project name
metadata$sample_type <- substr(metadata$sample, 14, 15) #Extract characters 14 e 15 from sample name
metadata <- metadata[metadata$sample_type %in% c("01", "11"), ] #Filter for samples 01 e 11
metadata$condition <- ifelse(metadata$sample_type == "01", "Tumor", "Normal") #Name condition
metadata$sample_type <- NULL #Remove $sample_type

metadata[] <- lapply(metadata, function(x) { 
  if (is.list(x)) {
    return(sapply(x, function(y) paste(y, collapse = ", ")))
  }
  return(x)
})
write.csv(metadata, "metadata_gbmR.csv", row.names = TRUE) #Save metadata as .csv

  # Expression
matrix <- assay(data_gbmR) #Extract expression matrix

valid_samples <- metadata$barcode #Identify samples in metadata
expression_matrix_filtered <- matrix[, colnames(matrix) %in% valid_samples] 

rownames(expression_matrix_filtered) <- sub("\\..*", "", rownames(expression_matrix_filtered)) #Normalize gene sufixes
rownames(expression_matrix_filtered) <- as.character(rownames(expression_matrix_filtered))
any(duplicated(rownames(expression_matrix_filtered)))  #Verify gene unicity (should return FALSE)
expression_matrix_filtered <- expression_matrix_filtered[!duplicated(rownames(expression_matrix_filtered)), ] #Remove duplicate lines

write.csv(expression_matrix_filtered, "expression_gbm.csv", row.names = TRUE) #Save expression matrix as .csv

# BRCA
  # Metadata
metadata <- as.data.frame(colData(data_brcaR)) #Extract metadata
metadata$project <- "TCGA-BRCA"  #Add project name
metadata$sample_type <- substr(metadata$sample, 14, 15) #Extract characters 14 e 15 from sample name
metadata <- metadata[metadata$sample_type %in% c("01", "11"), ] #Filter for samples 01 e 11
metadata$condition <- ifelse(metadata$sample_type == "01", "Tumor", "Normal") #Name condition
metadata$sample_type <- NULL #Remove $sample_type

metadata[] <- lapply(metadata, function(x) { 
  if (is.list(x)) {
    return(sapply(x, function(y) paste(y, collapse = ", ")))
  }
  return(x)
})
write.csv(metadata, "metadata_gbmR.csv", row.names = TRUE) #Save metadata as .csv

  # Expression
matrix <- assay(data_brcaR) #Extract expression matrix

valid_samples <- metadata$barcode #Identify samples in metadata
expression_matrix_filtered <- matrix[, colnames(matrix) %in% valid_samples] 

rownames(expression_matrix_filtered) <- sub("\\..*", "", rownames(expression_matrix_filtered)) #Normalize gene sufixes
rownames(expression_matrix_filtered) <- as.character(rownames(expression_matrix_filtered))
any(duplicated(rownames(expression_matrix_filtered)))  #Verify gene unicity (should return FALSE)
expression_matrix_filtered <- expression_matrix_filtered[!duplicated(rownames(expression_matrix_filtered)), ] #Remove duplicate lines

write.csv(expression_matrix_filtered, "expression_brca.csv", row.names = TRUE) #Save expression matrix as .csv

# LUAD
  # Metadata
metadata <- as.data.frame(colData(data_luadR)) #Extract metadata
metadata$project <- "TCGA-LUAD"  #Add project name
metadata$sample_type <- substr(metadata$sample, 14, 15) #Extract characters 14 e 15 from sample name
metadata <- metadata[metadata$sample_type %in% c("01", "11"), ] #Filter for samples 01 e 11
metadata$condition <- ifelse(metadata$sample_type == "01", "Tumor", "Normal") #Name condition
metadata$sample_type <- NULL #Remove $sample_type

metadata[] <- lapply(metadata, function(x) { 
  if (is.list(x)) {
    return(sapply(x, function(y) paste(y, collapse = ", ")))
  }
  return(x)
})
write.csv(metadata, "metadata_luadR.csv", row.names = TRUE) #Save metadata as .csv

  # Expression
matrix <- assay(data_luadR) #Extract expression matrix

valid_samples <- metadata$barcode #Identify samples in metadata
expression_matrix_filtered <- matrix[, colnames(matrix) %in% valid_samples] 

rownames(expression_matrix_filtered) <- sub("\\..*", "", rownames(expression_matrix_filtered)) #Normalize gene sufixes
rownames(expression_matrix_filtered) <- as.character(rownames(expression_matrix_filtered))
any(duplicated(rownames(expression_matrix_filtered)))  #Verify gene unicity (should return FALSE)
expression_matrix_filtered <- expression_matrix_filtered[!duplicated(rownames(expression_matrix_filtered)), ] #Remove duplicate lines

write.csv(expression_matrix_filtered, "expression_luad.csv", row.names = TRUE) #Save expression matrix as .csv

# COAD
  # Metadata
metadata <- as.data.frame(colData(data_coadR)) #Extract metadata
metadata$project <- "TCGA-COAD"  #Add project name
metadata$sample_type <- substr(metadata$sample, 14, 15) #Extract characters 14 e 15 from sample name
metadata <- metadata[metadata$sample_type %in% c("01", "11"), ] #Filter for samples 01 e 11
metadata$condition <- ifelse(metadata$sample_type == "01", "Tumor", "Normal") #Name condition
metadata$sample_type <- NULL #Remove $sample_type

metadata[] <- lapply(metadata, function(x) { 
  if (is.list(x)) {
    return(sapply(x, function(y) paste(y, collapse = ", ")))
  }
  return(x)
})
write.csv(metadata, "metadata_luadR.csv", row.names = TRUE) #Save metadata as .csv

  # Expression
matrix <- assay(data_luadR) #Extract expression matrix

valid_samples <- metadata$barcode #Identify samples in metadata
expression_matrix_filtered <- matrix[, colnames(matrix) %in% valid_samples] 

rownames(expression_matrix_filtered) <- sub("\\..*", "", rownames(expression_matrix_filtered)) #Normalize gene sufixes
rownames(expression_matrix_filtered) <- as.character(rownames(expression_matrix_filtered))
any(duplicated(rownames(expression_matrix_filtered)))  #Verify gene unicity (should return FALSE)
expression_matrix_filtered <- expression_matrix_filtered[!duplicated(rownames(expression_matrix_filtered)), ] #Remove duplicate lines

write.csv(expression_matrix_filtered, "expression_luad.csv", row.names = TRUE) #Save expression matrix as .csv
