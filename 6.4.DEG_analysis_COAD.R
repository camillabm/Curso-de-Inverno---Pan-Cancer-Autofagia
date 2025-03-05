#####################
### DEGs analysis ###
#####################

######## COAD ########

# Libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "biomaRt", "limma"))
install.packages("dplyr")

library(dplyr)
library(DESeq2)
library (biomaRt)
library (limma)

# Data identification 
expression_analysis_coad <- read.csv("2.0/Expression Pan/Expression/cohort/matrix_coadR_comb.csv", row.names = 1) #Previous expression table (genes x samples)
head(expression_analysis_coad)
metadata_coad <- read.csv("2.0/Expression Pan/Metadata/cohort/metadata_coadR_comb.csv") #Metadata table after manual curation (in .csv) (sample x metadata)
metadata_coad$condition <- as.factor(metadata_coad$condition)
head(metadata_coad)

valid_samples <- metadata_coad$barcode #Identification of samples in metadata after manual curation
colnames(expression_analysis_coad) <- ifelse(    #Sample ID padronization
  grepl("^TCGA", colnames(expression_analysis_coad)),        
  gsub("\\.", "-", colnames(expression_analysis_coad)),      
  colnames(expression_analysis_coad)                         
) 
expression_matrix_filtered <- expression_analysis_coad[, colnames(expression_analysis_coad) %in% valid_samples] #Expression table filtered

# Subset selection
metadata_subset <- metadata_coad %>% filter(project %in% c("TCGA-COAD", "GTEx")) #Project selection (it may not be necessary)
rownames(metadata_subset) <- metadata_subset$barcode
expression_subset <- expression_matrix_filtered[, metadata_subset$barcode] #Expression table filtered by metadata_subset #Expression table filtered by metadata_subset
expression_subset[expression_subset == 0] <- 1

### INITIAL DESEQ2 ANALYSIS ###

dds_coad <- DESeqDataSetFromMatrix(countData = expression_subset, 
                              colData = metadata_subset, 
                              design = ~ condition) #Compare expression data based on condition (normal x tumor)

dds_specific_coad <- DESeq(dds_coad) #DESeq analysis
results_specific <- results(dds_specific_coad, contrast = c("condition", "Tumor", "Normal")) #Results
res_df_coad <- as.data.frame(results_specific) #Results as dataframe
head(res_df_coad)

# Convertion of ensembl ID to HGNC symbol
res_df_coad$ensembl_id <- rownames(res_df_coad) 

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl",
                   host = "https://grch37.ensembl.org") #Connection to ensembl

hgnc_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), #Map matches from ensemnl ID and HGNC symbol
                      filters = "ensembl_gene_id",
                      values = res_df_coad$ensembl_id,
                      mart = ensembl)

res_df_coad <- merge(res_df_coad, hgnc_mapping, by.x = "ensembl_id", by.y = "ensembl_gene_id", all.x = TRUE) #Add HGNC symbol to results
head(res_df_coad)
rownames(res_df_coad) <- res_df_coad$hgnc_symbol #Replace ensembl ID for HGNC symbol
res_df_coad$ensembl_id <- NULL
head(res_df_coad)

# Filtering DEGs (padj < 0.01 & abs(log2FoldChange) > 1.5)
res_sig_coad <- res_df_coad[which(res_df_coad$padj < 0.01 & abs(res_df_coad$log2FoldChange) > 1.5), ]
nrow(res_sig_coad) #Annotate result

# Save differential expression data as .csv
write.csv(res_sig_coad, "2.0/Expression Pan/Results/tables/result_coadR_ARG.csv", row.names = FALSE)
