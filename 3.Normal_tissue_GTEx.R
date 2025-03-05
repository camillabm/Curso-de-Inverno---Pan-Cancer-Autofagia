##########################
### Normal Tissue GTEx ###
##########################

# GTEx data (.gdc) reading
  # GBM
    # Brain cortex (1)
gct_data <- read.table("2.0/GTEx_data/GBM/gene_reads_brain_cortex.gct", skip = 2, header = TRUE, sep = "\t") #Skip the first two lines
rownames(gct_data) <- gct_data$Name

samples_gtex <- colnames(gct_data)[-c(1, 2)]  #Samples extraction
metadata_1 <- data.frame(
  condition = rep("Normal", length(samples_gtex)),  #condition: "Normal"
  project = rep("GTEx", length(samples_gtex)),      #project: "GTEx"
  tissue_or_organ_of_origin = rep("Brain cortex", length(samples_gtex)), #tissue: "Brain cortex"
  barcode = samples_gtex  
)
rownames(metadata_1) <- samples_gtex

gct_data_1 <- gct_data[,-c( 1, 2)]
expression_matrix_1 <- as.matrix(gct_data_1) #expression matrix extraction

    # Brain frontal cortex (2)
gct_data <- read.table("2.0/GTEx_data/GBM/gene_reads_brain_frontal_cortex_ba9.gct", skip = 2, header = TRUE, sep = "\t") #Skip the first two lines
rownames(gct_data) <- gct_data$Name

samples_gtex <- colnames(gct_data)[-c(1, 2)]  #Samples extraction
metadata_2 <- data.frame(
  condition = rep("Normal", length(samples_gtex)),  #condition: "Normal"
  project = rep("GTEx", length(samples_gtex)),      #project: "GTEx"
  tissue_or_organ_of_origin = rep("Brain frontal cortex", length(samples_gtex)), #tissue: "Brain cortex"
  barcode = samples_gtex  
)
rownames(metadata_2) <- samples_gtex

gct_data_2 <- gct_data[,-c( 1, 2)]
expression_matrix_2 <- as.matrix(gct_data_2) #expression matrix extraction

    # Combining 1 and 2
matrix_comb <- cbind(expression_matrix_1, expression_matrix_2)
metadata_comb <- rbind(metadata_1, metadata_2)

write.csv(matrix_comb, "2.0/Expressão Pan/Expression/matrix_gbmR_normal.csv", row.names = TRUE) #Expression matrix
write.csv(metadata_comb, "2.0/Expressão Pan/Metadata/metadata_gbmR_normal.csv", row.names = TRUE) #Metadata

  # BRCA
    # Breast mammary tissue (3)
gct_data <- read.table("2.0/GTEx_data/BRCA/gene_reads_breast_mammary_tissue.gct", skip = 2, header = TRUE, sep = "\t") #Skip the first two lines
rownames(gct_data) <- gct_data$Name

samples_gtex <- colnames(gct_data)[-c(1, 2)]  #Samples extraction
metadata_3 <- data.frame(
  condition = rep("Normal", length(samples_gtex)),  #condition: "Normal"
  project = rep("GTEx", length(samples_gtex)),      #project: "GTEx"
  tissue_or_organ_of_origin = rep("Breast mammary tissue", length(samples_gtex)), #tissue: "Breast mammary tissue"
  barcode = samples_gtex  
)
rownames(metadata_3) <- samples_gtex

gct_data_3 <- gct_data[,-c( 1, 2)]
expression_matrix_3 <- as.matrix(gct_data_3) #expression matrix extraction

write.csv(expression_matrix_3, "2.0/Expressão Pan/Expression/total/matrix_brcaR_normal.csv", row.names = TRUE) #expression matrix
write.csv(metadata_3, "2.0/Expressão Pan/Metadata/total/metadata_brcaR_normal.csv", row.names = TRUE) #metadata

  # LUAD
    # Lung (4)
gct_data <- read.table("2.0/GTEx_data/LUAD/gene_reads_lung.gct", skip = 2, header = TRUE, sep = "\t") #Skip the first two lines
rownames(gct_data) <- gct_data$Name

samples_gtex <- colnames(gct_data)[-c(1, 2)]  #Samples extraction
metadata_4 <- data.frame(
  condition = rep("Normal", length(samples_gtex)),  #condition: "Normal"
  project = rep("GTEx", length(samples_gtex)),      #project: "GTEx"
  tissue_or_organ_of_origin = rep("Lung", length(samples_gtex)), #tissue: "Lung"
  barcode = samples_gtex  
)
rownames(metadata_4) <- samples_gtex

gct_data_4 <- gct_data[,-c( 1, 2)]
expression_matrix_4 <- as.matrix(gct_data_4) #expression matrix extraction

write.csv(expression_matrix_4, "2.0/Expressão Pan/Expression/total/matrix_luadR_normal.csv", row.names = TRUE) #Expression matrix
write.csv(metadata_4, "2.0/Expressão Pan/Metadata/total/metadata_luadR_normal.csv", row.names = TRUE) #Metadata

  # COAD
    # Colon Sigmoid (5)
gct_data <- read.table("2.0/GTEx_data/COAD/gene_reads_colon_sigmoid.gct", skip = 2, header = TRUE, sep = "\t") #Skip the first two lines
rownames(gct_data) <- gct_data$Name

samples_gtex <- colnames(gct_data)[-c(1, 2)]  #Samples extraction
metadata_5 <- data.frame(
  condition = rep("Normal", length(samples_gtex)),  #condition: "Normal"
  project = rep("GTEx", length(samples_gtex)),      #project: "GTEx"
  tissue_or_organ_of_origin = rep("Colon Sigmoid", length(samples_gtex)), #tissue: "Colon Sigmoid"
  barcode = samples_gtex  
)
rownames(metadata_5) <- samples_gtex

gct_data_5 <- gct_data[,-c( 1, 2)]
expression_matrix_5 <- as.matrix(gct_data_5) #expression matrix extraction

    # Colon transverse (6)
gct_data <- read.table("2.0/GTEx_data/COAD/gene_reads_colon_transverse.gct", skip = 2, header = TRUE, sep = "\t") #skip the first two lines
rownames(gct_data) <- gct_data$Name

samples_gtex <- colnames(gct_data)[-c(1, 2)]  #Samples extraction
metadata_6 <- data.frame(
  condition = rep("Normal", length(samples_gtex)),  #condition: "Normal"
  project = rep("GTEx", length(samples_gtex)),      #project: "GTEx"
  tissue_or_organ_of_origin = rep("Colon transverse", length(samples_gtex)), #tissue: "Colon transverse"
  barcode = samples_gtex  
)
rownames(metadata_6) <- samples_gtex

gct_data_6 <- gct_data[,-c( 1, 2)]
expression_matrix_6 <- as.matrix(gct_data_6) #expression matrix extraction

    # Combining 5 and 6
matrix_comb <- cbind(expression_matrix_5, expression_matrix_6)
metadata_comb <- rbind(metadata_5, metadata_6)

write.csv(matrix_comb, "2.0/Expressão Pan/Expression/total/matrix_coadR_normal.csv", row.names = TRUE) #Expression matrix
write.csv(metadata_comb, "2.0/Expressão Pan/Metadata/total/metadata_coadR_normal.csv", row.names = TRUE) #Metadata
