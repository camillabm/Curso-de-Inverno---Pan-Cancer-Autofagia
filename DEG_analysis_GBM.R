#####################
### DEGs analysis ###
#####################

######## GBM ########

# Libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano", "biomaRt", "limma"))
install.packages("pheatmap")
install.packages("dplyr")
install.packages("RColorBrewer")
install.packages("ggplot2")
install.packages("fmsb")

library(dplyr)
library(DESeq2)
library(pheatmap)
library(readxl)
library (biomaRt)
library(ggplot2)
library (EnhancedVolcano)
library(RColorBrewer)
library(fmsb)

# Data identification 
expression_analysis_gbm <- read.csv("2.0/Expression Pan/Expression/cohort/matrix_gbmR_comb.csv", row.names = 1) #Previous expression table (genes x samples)
head(expression_analysis_gbm)
metadata_gbm <- read.csv("2.0/Expression Pan/Metadata/cohort/metadata_gbmR_comb.csv") #Metadata table after manual curation (in .csv) (sample x metadata)
metadata_gbm$condition <- as.factor(metadata_gbm$condition)
head(metadata_gbm)

valid_samples <- metadata_gbm$barcode #Identification of samples in metadata after manual curation
colnames(expression_analysis_gbm ) <- ifelse(    #Sample ID padronization
  grepl("^TCGA", colnames(expression_analysis_gbm )),        
  gsub("\\.", "-", colnames(expression_analysis_gbm )),      
  colnames(expression_analysis_gbm )                         
) 
expression_matrix_filtered <- expression_analysis_gbm[, colnames(expression_analysis_gbm) %in% valid_samples] #Expression table filtered

# Subset selection
metadata_subset <- metadata_gbm %>% filter(project %in% c("TCGA-GBM", "GTEx")) #Project selection (it may not be necessary)
rownames(metadata_subset) <- metadata_subset$barcode
expression_subset <- expression_matrix_filtered[, metadata_subset$barcode] #Expression table filtered by metadata_subset
expression_subset[expression_subset == 0] <- 1

### INITIAL DESEQ2 ANALYSIS ###

dds_gbm <- DESeqDataSetFromMatrix(countData = expression_subset, 
                              colData = metadata_subset, 
                              design = ~ condition) #Compare expression data based on condition (normal x tumor)

dds_specific_gbm <- DESeq(dds_gbm) #DESeq analysis
results_specific <- results(dds_specific_gbm, contrast = c("condition", "Tumor", "Normal")) #Results
res_df_gbm <- as.data.frame(results_specific) #Results as dataframe
head(res_df_gbm)

# Convertion of ensembl ID to HGNC symbol
res_df_gbm$ensembl_id <- rownames(res_df_gbm) 

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl",
                   host = "https://grch37.ensembl.org") #Connection to ensembl

hgnc_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), #Map matches from ensemnl ID and HGNC symbol
                      filters = "ensembl_gene_id",
                      values = res_df_gbm$ensembl_id,
                      mart = ensembl)

res_df_gbm <- merge(res_df_gbm, hgnc_mapping, by.x = "ensembl_id", by.y = "ensembl_gene_id", all.x = TRUE) #Add HGNC symbol to results
head(res_df_gbm)
rownames(res_df_gbm) <- res_df_gbm$hgnc_symbol #Replace ensembl ID for HGNC symbol
res_df_gbm$ensembl_id <- NULL
head(res_df_gbm)

# Filtering DEGs (padj < 0.01 & abs(log2FoldChange) > 1.5)
res_sig_gbm <- res_df_gbm[which(res_df_gbm$padj < 0.01 & abs(res_df_gbm$log2FoldChange) > 1.5), ]
nrow(res_sig_gbm) #Annotate result

# Save differential expression data as .csv
write.csv(res_sig_gbm, "2.0/Expression Pan/Results/Tables/result_gbmR_ARG.csv", row.names = FALSE)

# Volcano Plot
volcano_plot <- EnhancedVolcano(res_df_gbm, 
                               lab = rownames(res_df_gbm),
                               x = 'log2FoldChange',  
                               y = 'padj',            
                               title = 'Differential Expression Analysis',
                               pCutoff = 0.01,        
                               FCcutoff = 1.5,          
                               pointSize = 3,         
                               labSize = 0,           
                               drawConnectors = FALSE, 
                               widthConnectors = 0.5) 


ggsave("2.0/Expression Pan/Results/Figures/volcano_plot_gbmR_ARG.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)

# UMAP
library(umap)
library(uwot)

vst_gbm <- varianceStabilizingTransformation(dds_specific_gbm, blind = TRUE)
expr_gbm <- assay(vst_gbm)
metadata_gbm <- colData(dds_specific_gbm)
metadata_gbm <- as.data.frame(metadata_gbm)

set.seed(42)  # Para reprodutibilidade

umap_result <- umap(
  t(expr_gbm),   # Matriz de expressão (genes x amostras)
  n_neighbors = 15,  # Ajusta o número de vizinhos (valores menores = mais locais, valores maiores = mais globais)
  min_dist = 0.001,  # Controla a dispersão dos pontos (menor = clusters mais compactos, maior = mais espalhado)
  spread =  2,  # Controla a separação dos grupos
  metric = "correlation",  # Pode testar "euclidean", "correlation" ou "cosine"
  n_components = 2,  # Número de dimensões do UMAP
  learning_rate = 0.1,
  fast_sgd = TRUE,  # Torna o cálculo mais rápido
  verbose = TRUE  # Mostra informações durante o processo
)

umap_df <- data.frame(
  UMAP1 = umap_result[,1],
  UMAP2 = umap_result[,2],
  Condition = metadata_gbm$condition,  # Ou outra variável relevante
  Batch = metadata_gbm$project  # Para verificar efeito do batch
)

colors <- c(
  "Normal" = "gray70",  
  "Tumor" = "#9B59B6")

UMAP <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Condition)) +
  geom_point(size = 2, alpha = 0.7) +  # Pontos maiores e mais suaves
  scale_color_manual(values = colors) +  # Aplicar as cores definidas
  theme_minimal(base_family = "Arial") +  # Fonte Arial para remover serifas
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(
    title = "UMAP - Visualização Melhorada",
    color = "Condição"
  )

ggsave("2.0/Expression Pan/Results/Figures/UMAP_gbmR_ARG.png", plot = UMAP, width = 8, height = 6, dpi = 300, bg = "white")

  # Saving umap groups
umap_df_gbm_groups <- as.data.frame(umap_df)
umap_df_gbm_groups$SampleId <- rownames(umap_df_gbm_groups)
umap_info <-umap_df_gbm_groups[, c("SampleId", "UMAP1", "UMAP2", "Condition")]
head(umap_info)
umap_df_gbm_groups <- umap_df_gbm_groups %>%
  mutate(Group = case_when(
    UMAP1 > 0 & UMAP2 > 0  ~ "1",
    UMAP1 > 0 & UMAP2 <= 0 ~ "2",
    TRUE ~ "Normal"  # Qualquer outro caso será "Normal"
  ))

umap_groups_export <- umap_df_gbm_groups[, c("SampleId", "Group")]
write.csv(umap_groups_export, "2.0/Expression Pan/Results/Tables/groups_gbm_umap.csv", row.names = FALSE)
