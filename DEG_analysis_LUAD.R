#####################
### DEGs analysis ###
#####################

######## LUAD ########

# Libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano"))
install.packages("dplyr")
install.packages("readxl")
install.packages("fpc")
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
expression_analysis_luad <- read.csv("2.0/Expression Pan/Expression/cohort/matrix_luadR_comb.csv", row.names = 1) #Previous expression table (genes x samples)
head(expression_analysis_luad)
metadata_luad <- read.csv("2.0/Expression Pan/Metadata/cohort/metadata_luadR_comb.csv") #Metadata table after manual curation (in .csv) (sample x metadata)
metadata_luad$condition <- as.factor(metadata_luad$condition)
head(metadata_luad)

valid_samples <- metadata_luad$barcode #Identification of samples in metadata after manual curation
colnames(expression_analysis_luad) <- ifelse(    #Sample ID padronization
  grepl("^TCGA", colnames(expression_analysis_luad)),        
  gsub("\\.", "-", colnames(expression_analysis_luad)),      
  colnames(expression_analysis_luad)                         
) 
expression_matrix_filtered <- expression_analysis_luad[, colnames(expression_analysis_luad) %in% valid_samples] #Expression table filtered

# Subset selection
metadata_subset <- metadata_luad %>% filter(project %in% c("TCGA-LUAD", "GTEx")) #Project selection (it may not be necessary)
rownames(metadata_subset) <- metadata_subset$barcode
expression_subset <- expression_matrix_filtered[, metadata_subset$barcode] #Expression table filtered by metadata_subset
expression_subset[expression_subset == 0] <- 1

### INITIAL DESEQ2 ANALYSIS ###

dds_luad <- DESeqDataSetFromMatrix(countData = expression_subset, 
                              colData = metadata_subset, 
                              design = ~ condition) #Compare expression data based on condition (normal x tumor)

dds_specific_luad <- DESeq(dds_luad) #DESeq analysis
results_specific <- results(dds_specific_luad, contrast = c("condition", "Tumor", "Normal")) #Results
res_df_luad <- as.data.frame(results_specific) #Results as dataframe
head(res_df_luad)

# Convertion of ensembl ID to HGNC symbol
res_df_luad$ensembl_id <- rownames(res_df_luad) 

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl",
                   host = "https://grch37.ensembl.org") #Connection to ensembl

hgnc_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), #Map matches from ensemnl ID and HGNC symbol
                      filters = "ensembl_gene_id",
                      values = res_df_luad$ensembl_id,
                      mart = ensembl)

res_df_luad <- merge(res_df_luad, hgnc_mapping, by.x = "ensembl_id", by.y = "ensembl_gene_id", all.x = TRUE) #Add HGNC symbol to results
head(res_df_luad)
rownames(res_df_luad) <- res_df_luad$hgnc_symbol #Replace ensembl ID for HGNC symbol
res_df_luad$ensembl_id <- NULL
head(res_df_luad)

# Filtering DEGs (padj < 0.01 & abs(log2FoldChange) > 1.5)
res_sig_luad <- res_df_luad[which(res_df_luad$padj < 0.01 & abs(res_df_luad$log2FoldChange) > 1.5), ]
nrow(res_sig_luad) #Annotate result

# Save differential expression data as .csv
write.csv(res_sig_luad, "2.0/Expression Pan/Results/tables/result_luadR_ARG.csv", row.names = FALSE)

# Volcano Plot
volcano_plot <- EnhancedVolcano(res_df_luad, 
                                lab = rownames(res_df_luad),
                                x = 'log2FoldChange',  # Eixo X para log2 Fold Change
                                y = 'padj',            # Eixo Y para p-value ajustado
                                title = 'Differential Expression Analysis',
                                pCutoff = 0.01,        # Filtro para p-value ajustado
                                FCcutoff = 1.5,          # Filtro para log2 Fold Change
                                pointSize = 3,         # Tamanho dos pontos
                                labSize = 0,           # Tamanho das labels
                                drawConnectors = FALSE, # Conectar as labels aos pontos
                                widthConnectors = 0.5) # Espessura dos conectores


ggsave("2.0/Expression Pan/Results/figures/volcano_plot_luadR_ARG.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)

# UMAP
library(umap)
library(uwot)

vst_luad <- varianceStabilizingTransformation(dds_specific_luad, blind = TRUE)
expr_luad <- assay(vst_luad)
metadata_luad <- colData(dds_specific_luad)
metadata_luad <- as.data.frame(metadata_luad)

set.seed(42)  # Para reprodutibilidade

umap_result <- umap(
  t(expr_luad),   # Matriz de expressão (genes x amostras)
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
  Condition = metadata_luad$condition,  # Ou outra variável relevante
  Batch = metadata_luad$project  # Para verificar efeito do batch
)

colors <- c(
  "Normal" = "gray70",  
  "Tumor" = "#D8B")

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

ggsave("2.0/Expression Pan/Results/figures/UMAP_luadR_ARG.png", plot = UMAP, width = 8, height = 6, dpi = 300, bg = "white")

