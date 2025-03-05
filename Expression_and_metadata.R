### PREPARAÇÃO DOS DADOS PARA AS ANÁLISES ###

# Instalar pacotes e chamar bibliotecas
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("SummarizedExperiment", "biomaRt"))
install.packages("dplyr")

library(SummarizedExperiment)
library(dplyr)
library (biomaRt)

# Salvar matriz de expressão e metadados

  #Metadados
metadata <- as.data.frame(colData(data_coadR)) # extrai metadados
metadata$project <- "TCGA-COAD"  # adiciona o nome do projeto
metadata$sample_type <- substr(metadata$sample, 14, 15) #extrai os caracteres 14 e 15
metadata <- metadata[metadata$sample_type %in% c("01", "11"), ] #filtra apenas amostras 01 e 11
metadata$condition <- ifelse(metadata$sample_type == "01", "Tumor", "Normal") #atribui condição
metadata$sample_type <- NULL #remove a coluna sample_type

metadata[] <- lapply(metadata, function(x) { # Converte listas para caracteres (ou outro tipo, conforme necessário)
  if (is.list(x)) {
    return(sapply(x, function(y) paste(y, collapse = ", ")))
  }
  return(x)
})
write.csv(metadata, "metadata_coadR.csv", row.names = TRUE) #salva metadados em csv

  #Valores de expressão
matrix <- assay(data_coadR) #extrai a matriz de expressão

valid_samples <- metadata$barcode #identifica amostras que estão nos metadados
expression_matrix_filtered <- matrix[, colnames(matrix) %in% valid_samples] #filtra amostras que estão nos metadaos na matriz de expressão

rownames(expression_matrix_filtered) <- sub("\\..*", "", rownames(expression_matrix_filtered)) #normaliza sufixos dos genes
rownames(expression_matrix_filtered) <- as.character(rownames(expression_matrix_filtered)) #converte os genes como character
any(duplicated(rownames(expression_matrix_filtered)))  #verifica a unicidade dos genes (deve retornar FALSE)
expression_matrix_filtered <- expression_matrix_filtered[!duplicated(rownames(expression_matrix_filtered)), ] #remove linhas duplicadas

write.csv(expression_matrix_filtered, "expression_coad.csv", row.names = TRUE) #salva martiz de expressão em csv

# Salvar matriz de metilação e metadados

  # Metadados
metadata <- as.data.frame(colData(data_brcaM)) # extrai metadados
metadata$project <- "TCGA-COAD"  # adiciona o nome do projeto
metadata$sample_type <- substr(metadata$sample, 14, 15) #extrai os caracteres 14 e 15
metadata <- metadata[metadata$sample_type %in% c("01", "11"), ] #filtra apenas amostras 01 e 11
metadata$condition <- ifelse(metadata$sample_type == "01", "Tumor", "Normal") #atribui condição
metadata$sample_type <- NULL #remove a coluna sample_type

metadata[] <- lapply(metadata, function(x) { # Converte listas para caracteres (ou outro tipo, conforme necessário)
  if (is.list(x)) {
    return(sapply(x, function(y) paste(y, collapse = ", ")))
  }
  return(x)
})
write.csv(metadata, "metadata_coadM.csv", row.names = TRUE) #salva metadados em csv

  # Metilação
beta_values <- assay(data_brcaM) #extrai os beta values
beta_values <- beta_values[, metadata$barcode] #filtra só as amostras que tem nos metadados
rownames(beta_values) <- sub("\\..*", "", rownames(beta_values)) #normaliza sufixos dos genes

write.csv(beta_values, "methylation_coad.csv")

  # Caso BRCA
methylation_data <- load("BRCA_methylation_data.rda")
head(methylation_data)
str(methylation_data)
class(methylation_data)  # Check the object class
typeof(methylation_data)  # Check the object type


# Identificar as amostras comuns e filtrar
  # Metadados
meta1 <- read.csv("D:/TCC/Documentos finais/exp/metadata_coadR_f.csv", sep = ";", stringsAsFactors = FALSE) #converte ; para ,
meta2 <- read.csv("D:/TCC/Documentos finais/met/metadata_coadM_f.csv", sep = ";", stringsAsFactors = FALSE) #converte ; para ,

linhas_comuns <- inner_join(meta1, meta2, by = "sample") #encontra as linhas comuns
cat("Número de linhas em comum:", nrow(linhas_comuns), "\n")

meta1_filtrado <- meta1 %>% filter(sample %in% linhas_comuns$sample)
meta2_filtrado <- meta2 %>% filter(sample %in% linhas_comuns$sample)

table(meta2_filtrado$condition)

write.csv(meta1_filtrado, "metadata_coadR_comb.csv", row.names = FALSE)
write.csv(meta2_filtrado, "metadata_coadM.comb.csv", row.names = FALSE)

  # Matrizes
arquivo_filtros1 <- "D:/TCC/Documentos finais/exp/expression_coad.csv"
arquivo_filtros2 <- "D:/TCC/Documentos finais/met/methylation_coad.csv"

dados1 <- read.csv(arquivo_filtros1, sep = ",", stringsAsFactors = FALSE)
dados2 <- read.csv(arquivo_filtros2, sep = ",", stringsAsFactors = FALSE)

dados1 <- dados1[, meta1$barcode]
dados2 <- dados2[, meta2$barcode]

write.csv(dados1, "expression_coad_comb.csv", row.names = FALSE)
write.csv(dados2, "methylation_coad_comb.csv", row.names = FALSE)

