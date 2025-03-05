#######################
### Pan comparisons ###
#######################

#Install libraries
install.packages("dplyr")
library(dplyr)

### UMAP ###
install.packages("umap")   
install.packages("uwot")   
install.packages("ggplot2")  

library(umap)
library(uwot)
library(ggplot2)

vst_gbm <- varianceStabilizingTransformation(dds_gbm, blind = TRUE)
vst_brca <- varianceStabilizingTransformation(dds_brca, blind = TRUE)
vst_luad <- varianceStabilizingTransformation(dds_luad, blind = TRUE)
vst_coad <- varianceStabilizingTransformation(dds_coad, blind = TRUE)

expr_gbm <- assay(vst_gbm)
expr_brca <- assay(vst_brca)
expr_luad <- assay(vst_luad)
expr_coad <- assay(vst_coad)

common_genes <- Reduce(intersect, list(rownames(expr_gbm), rownames(expr_brca), rownames(expr_luad), rownames(expr_coad)))

expr_gbm <- expr_gbm[common_genes, ]
expr_brca <- expr_brca[common_genes, ]
expr_luad <- expr_luad[common_genes, ]
expr_coad <- expr_coad[common_genes, ]

expr_combined <- cbind(expr_gbm, expr_brca, expr_luad, expr_coad)

metadata_gbm <- colData(dds_gbm)
metadata_brca <- colData(dds_brca)
metadata_luad <- colData(dds_luad)
metadata_coad <- colData(dds_coad)

metadata_gbm$Batch <- paste("GBM", metadata_gbm$condition, sep = "_")
metadata_brca$Batch <- paste("BRCA", metadata_brca$condition, sep = "_")
metadata_luad$Batch <- paste("LUAD", metadata_luad$condition, sep = "_")
metadata_coad$Batch <- paste("COAD", metadata_coad$condition, sep = "_")

metadata_gbm_f <- metadata_gbm[, c("barcode", "project", "condition", "Batch")]
metadata_brca_f <- metadata_brca[, c("barcode", "project", "condition", "Batch")]
metadata_luad_f <- metadata_luad[, c("barcode", "project", "condition", "Batch")]
metadata_coad_f <- metadata_coad[, c("barcode", "project", "condition", "Batch")]

metadata_combined <- rbind(metadata_gbm_f, metadata_brca_f, metadata_luad_f, metadata_coad_f)

metadata_combined <- metadata_combined[colnames(expr_combined), ]

all(colnames(expr_combined) == rownames(metadata_combined))  #Should return TRUE

library(sva)
expr_combined_corrected <- ComBat(dat = expr_combined, batch = metadata_combined$Batch)

set.seed(42)  

umap_result <- umap(
  t(expr_combined_corrected), 
  n_neighbors = 10,  
  min_dist = 0.00001,  
  spread =  2,  
  metric = "correlation",  
  n_components = 2,  
  learning_rate = 0.1,
  fast_sgd = TRUE,  
  verbose = TRUE  
)

umap_df <- data.frame(
  UMAP1 = umap_result[,1],
  UMAP2 = umap_result[,2],  
  Condition = metadata_coad$condition,
  Batch = metadata_combined$Batch  
)

  # Basic Plot
custom_colors <- c(
  "GBM_Tumor" = "#9B59B6",
  "GBM_Normal" = "#C78CD7",
  "BRCA_Tumor" = "#4682B4",
  "BRCA_Normal" = "#7CA6D9",
  "LUAD_Tumor" = "#D8B",
  "LUAD_Normal" = "#F0C4D9",
  "COAD_Tumor" = "#B5EAD7",
  "COAD_Normal" = "#E2F7E0")

UMAP <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Batch)) +
  geom_point(size = 2, alpha = 0.7) +  # Pontos maiores e mais suaves
  scale_color_manual(values = custom_colors) +  # Aplicar as cores definidas
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(
    title = "UMAP - Visualização Melhorada",
    color = "Condição"
  )

ggsave("2.0/Expression Pan/Results/figures/UMAP_pan.png", plot = UMAP, width = 8, height = 6, dpi = 300, bg = "white")

  # GBM Plot
umap_df <- as.data.frame(umap_result)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Batch <- metadata_combined$Batch

umap_df$condition_color <- ifelse(grepl("_Normal$", umap_df$Batch), "Normal", 
                                  ifelse(grepl("GBM_Tumor", umap_df$Batch), "GBM", 
                                         "Other Project"))
custom_colors <- c(
  "Other Project" = "gray70",
  "Normal" = "#F8E47B",   
  "GBM" = "#9B59B6"  )

UMAP <- ggplot() +
  geom_point(data = umap_df[umap_df$condition_color != "GBM", ], 
             aes(x = UMAP1, y = UMAP2, color = condition_color), size = 2, alpha = 0.7) +
  geom_point(data = umap_df[umap_df$condition_color == "GBM", ], 
             aes(x = UMAP1, y = UMAP2, color = condition_color), size = 2, alpha = 0.7) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  ggtitle("UMAP com GBM_Tumor Por Cima")

ggsave("2.0/Expression Pan/Results/figures/UMAP_pan_gbm.png", plot = UMAP, width = 8, height = 6, dpi = 300, bg = "white")

### DEGs comparison ###
  # Data loading
df_gbm <- read.csv("2.0/Expression Pan/Results/tables/result_gbmR_ARG.csv")
df_brca <- read.csv("2.0/Expression Pan/Results/tables/result_brcaR_ARG.csv")
df_luad <- read.csv("2.0/Expression Pan/Results/tables/result_luadR_ARG.csv")
df_coad <- read.csv("2.0/Expression Pan/Results/tables/result_coadR_ARG.csv")

  # Extracting column hgnc_symbol
hgnc_gbm <- df_gbm$hgnc_symbol
hgnc_brca <- df_brca$hgnc_symbol
hgnc_luad <- df_luad$hgnc_symbol
hgnc_coad <- df_coad$hgnc_symbol

  # Obtaining DEGs data (no log2FC signal)
unique_values_no_sign <- union(union(union(hgnc_gbm, hgnc_brca), hgnc_luad), hgnc_coad)

intersect_all_no_sign <- Reduce(intersect, list(hgnc_gbm, hgnc_brca, hgnc_luad, hgnc_coad))

intersect_gb_no_sign <- intersect(hgnc_gbm, hgnc_brca)
intersect_gl_no_sign <- intersect(hgnc_gbm, hgnc_luad)
intersect_gc_no_sign <- intersect(hgnc_gbm, hgnc_coad)
intersect_bl_no_sign <- intersect(hgnc_brca, hgnc_luad)
intersect_bc_no_sign <- intersect(hgnc_brca, hgnc_coad)
intersect_lc_no_sign <- intersect(hgnc_luad, hgnc_coad)

intersect_gbl_no_sign <- Reduce(intersect, list(hgnc_gbm, hgnc_brca, hgnc_luad))
intersect_gbc_no_sign <- Reduce(intersect, list(hgnc_gbm, hgnc_brca, hgnc_coad))
intersect_glc_no_sign <- Reduce(intersect, list(hgnc_gbm, hgnc_luad, hgnc_coad))
intersect_blc_no_sign <- Reduce(intersect, list(hgnc_brca, hgnc_luad, hgnc_coad))

list(
  unique_values_no_sign = unique_values_no_sign,
  intersect_all_no_sign = intersect_all_no_sign,
  intersect_gb_no_sign = intersect_gb_no_sign,
  intersect_gl_no_sign = intersect_gl_no_sign,
  intersect_gc_no_sign = intersect_gc_no_sign,
  intersect_bl_no_sign = intersect_bl_no_sign,
  intersect_bc_no_sign = intersect_bc_no_sign,
  intersect_lc_no_sign = intersect_lc_no_sign,
  intersect_gbl_no_sign = intersect_gbl_no_sign,
  intersect_gbc_no_sign = intersect_gbc_no_sign,
  intersect_glc_no_sign = intersect_glc_no_sign,
  intersect_blc_no_sign = intersect_blc_no_sign
)

install.packages("VennDiagram")
library(VennDiagram)

sets <- list(
  GBM = hgnc_gbm, 
  BRCA = hgnc_brca, 
  LUAD = hgnc_luad, 
  COAD = hgnc_coad
)

grid.newpage()

install.packages("extrafont")
library(extrafont)
font_import() 
loadfonts(device = "win")
par(family = "Arial")

png("2.0/Expression Pan/Results/figures/venn_diagram.png", width = 1200, height = 1000, res = 150)
grid.newpage()
venn.plot <- venn.diagram(
  x = sets,
  category.names = c("GBM", "BRCA", "LUAD", "COAD"),
  filename = NULL,  
  output = TRUE,
  fill = c("#9B59B6", "#4682B4", "#D8B", "#B5EAD7"), 
  alpha = 0.6,  
  cat.cex = 1.5,  
  cex = 1.5,  
  cat.pos = c(0, 0, 0, 0),  
  cat.dist = 0.1,  
  main = "Venn Diagram - Interseções de hgnc_symbol",
  main.fontfamily = "Arial",  
  main.fontface = "plain",  
  sub = "Comparação entre os conjuntos de genes",  
  sub.fontfamily = "Arial",  
  sub.fontface = "plain",  
  fontfamily = "Arial",  
  cat.fontfamily = "Arial"  
)
grid.draw(venn.plot)
dev.off() 

  # Obtaining DEGs data (log2FC signal)
check_sign_consistency <- function(df_list) {
  all_genes <- Reduce(union, lapply(df_list, function(df) df$hgnc_symbol)) #All genes present
  consistency_results <- setNames(rep(NA, length(all_genes)), all_genes)
  for (gene in all_genes) { #Sinal consistency verification
    gene_signs <- sapply(df_list, function(df) {
      logfc <- df$log2FoldChange[df$hgnc_symbol == gene]
      if (length(logfc) > 0) sign(logfc) else NA  
    })
    gene_signs <- na.omit(gene_signs)
    consistency_results[gene] <- length(unique(gene_signs)) == 1
  }
  return(data.frame(hgnc_symbol = names(consistency_results), Consistent = consistency_results))
} #DF with results

df_list <- list(GBM = df_gbm, BRCA = df_brca, LUAD = df_luad, COAD = df_coad)

sign_consistency_results <- check_sign_consistency(df_list)

head(sign_consistency_results)

  # Obtaingn individual DEGs from GBM
genes_only_in_gbm <- hgnc_gbm[!hgnc_gbm %in% union(hgnc_brca, union(hgnc_luad, hgnc_coad))]
print(genes_only_in_gbm)

### DEGs GO ###

  # Libraries
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library (clusterProfiler)
library (org.Hs.eg.db)

  # Function to convert HGNC to Entrez
hgnc_to_entrez <- function(hgnc_symbols) {
  gene_ids <- bitr(hgnc_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(gene_ids$ENTREZID)  # Retorna os IDs Entrez
}

  # Function to GO analysis
    #Biological Process
go_analysis_BP <- function(genes) {
  entrez_ids <- hgnc_to_entrez(genes)
  go_results <- enrichGO(entrez_ids, 
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP",   # Você pode escolher BP (Biological Process), MF (Molecular Function), ou CC (Cellular Component)
                         pAdjustMethod = "BH", # Método de ajuste de p-valor
                         qvalueCutoff = 0.05)  # Definir o cutoff para significância
  
  return(go_results)
}

    #Molecular Function
go_analysis_MF <- function(genes) {
  entrez_ids <- hgnc_to_entrez(genes)
  go_results <- enrichGO(entrez_ids, 
                         OrgDb = org.Hs.eg.db, 
                         ont = "MF",   # Você pode escolher BP (Biological Process), MF (Molecular Function), ou CC (Cellular Component)
                         pAdjustMethod = "BH", # Método de ajuste de p-valor
                         qvalueCutoff = 0.05)  # Definir o cutoff para significância
  
  return(go_results)
}

    #Cellular Component
go_analysis_CC <- function(genes) {
  entrez_ids <- hgnc_to_entrez(genes)
  go_results <- enrichGO(entrez_ids, 
                         OrgDb = org.Hs.eg.db, 
                         ont = "CC",   # Você pode escolher BP (Biological Process), MF (Molecular Function), ou CC (Cellular Component)
                         pAdjustMethod = "BH", # Método de ajuste de p-valor
                         qvalueCutoff = 0.05)  # Definir o cutoff para significância
  
  return(go_results)
}

  # Positive DEGs
df_gbm_pos <- df_gbm %>% filter(log2FoldChange > 0)
df_brca_pos <- df_brca %>% filter(log2FoldChange > 0)
df_luad_pos <- df_luad %>% filter(log2FoldChange > 0)
df_coad_pos <- df_coad %>% filter(log2FoldChange > 0)

df_gbm_pos <- df_gbm_pos$hgnc_symbol
df_brca_pos <- df_brca_pos$hgnc_symbol
df_luad_pos <- df_luad_pos$hgnc_symbol
df_coad_pos <- df_coad_pos$hgnc_symbol

intersect_all_no_sign_pos <- Reduce(intersect, list(df_gbm_pos, df_brca_pos, df_luad_pos, df_coad_pos))
all_other_genes_pos <- unique(c(df_brca_pos, df_luad_pos, df_coad_pos))
genes_only_in_gbm_pos <- df_gbm_pos[!df_gbm_pos %in% all_other_genes_pos]

intersect_all_no_sign_pos<- go_analysis_CC(intersect_all_no_sign_pos)
df_gbm_pos <- go_analysis_CC(df_gbm_pos)
genes_only_in_gbm_pos <- go_analysis_CC(genes_only_in_gbm_pos)

barplot(intersect_all_no_sign_pos, showCategory = 5)
barplot(df_gbm_pos, showCategory = 5)
barplot(genes_only_in_gbm_pos, showCategory = 5)

  # Negative DEGs
df_gbm_neg <- df_gbm %>% filter(log2FoldChange < 0)
df_brca_neg <- df_brca %>% filter(log2FoldChange < 0)
df_luad_neg <- df_luad %>% filter(log2FoldChange < 0)
df_coad_neg <- df_coad %>% filter(log2FoldChange < 0)

df_gbm_neg <- df_gbm_neg$hgnc_symbol
df_brca_neg <- df_brca_neg$hgnc_symbol
df_luad_neg <- df_luad_neg$hgnc_symbol
df_coad_neg <- df_coad_neg$hgnc_symbol

intersect_all_no_sign_neg <- Reduce(intersect, list(df_gbm_neg, df_brca_neg, df_luad_neg, df_coad_neg))
all_other_genes_neg <- unique(c(df_brca_neg, df_luad_neg, df_coad_neg))
genes_only_in_gbm_neg <- df_gbm_neg[!df_gbm_neg %in% all_other_genes_neg]

intersect_all_no_sign_neg <- go_analysis_BP(intersect_all_no_sign_neg)
df_gbm_neg <- go_analysis_BP(df_gbm_neg)
genes_only_in_gbm_neg <- go_analysis_BP(genes_only_in_gbm_neg)

barplot(intersect_all_no_sign_neg, showCategory = 5)
barplot(df_gbm_neg, showCategory = 5)
barplot(genes_only_in_gbm_neg, showCategory = 5)
