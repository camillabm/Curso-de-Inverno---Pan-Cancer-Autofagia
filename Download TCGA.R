### DOWNLOAD DADOS TCGA ###

# Instalar pacotes e chamar bibliotecas
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

# Obter dados de RNA-seq de vários projetos (pancancer)

  #GBM
GBM <- c("TCGA-GBM")
query_rna_gbm <- GDCquery(
  project = GBM,
  experimental.strategy = "RNA-Seq",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query = query_rna_gbm, method = "api", files.per.chunk = 1) #se já estiver baixado ele não baixa de novo
data_gbmR <- GDCprepare(query_rna_gbm) 

  #BRCA
BRCA <- c("TCGA-BRCA") #coloque aqui os TCGA escolhidos
query_rna_brca <- GDCquery(
  project = BRCA,
  experimental.strategy = "RNA-Seq",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query = query_rna_brca, method = "api", files.per.chunk = 1) #se já estiver baixado ele não baixa de novo
data_brcaR <- GDCprepare(query_rna_brca) 

  #LUAD
LUAD <- c("TCGA-LUAD") #coloque aqui os TCGA escolhidos
query_rna_luad <- GDCquery(
  project = LUAD,
  experimental.strategy = "RNA-Seq",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query = query_rna_luad, method = "api", files.per.chunk = 1) #se já estiver baixado ele não baixa de novo
data_luadR <- GDCprepare(query_rna_luad) 

#COAD
COAD <- c("TCGA-COAD") #coloque aqui os TCGA escolhidos
query_rna_coad <- GDCquery(
  project = COAD,
  experimental.strategy = "RNA-Seq",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query = query_rna_coad, method = "api", files.per.chunk = 1) #se já estiver baixado ele não baixa de novo
data_coadR <- GDCprepare(query_rna_coad) 

# Obter dados de metilação de DNA de vários projetos (pancancer)
BiocManager::install("sesameData")
BiocManager::install("sesame")
BiocManager::install("GenomeInfoDb")
BiocManager::install("GenomicRanges")


library (sesameData)
library (sesame)
sesameDataCache() #rodar apenas na primeira vez que instala o pacote
library (GenomeInfoDb)
library (GenomicRanges)

  #GBM
GBM <- c("TCGA-GBM")
query_met_gbm <- GDCquery(
  project = GBM,
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)

GDCdownload(query = query_met_gbm, method = "api", files.per.chunk = 1) #deixou mais leve pra baixar
data_gbmM <- GDCprepare(query_met_gbm) 

  #BRCA
BRCA <- c("TCGA-BRCA") #coloque aqui os TCGA escolhidos
query_met_brca <- GDCquery(
  project = BRCA,
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)

GDCdownload(query = query_met_brca, method = "api", files.per.chunk = 1) #deixou mais leve pra baixar
data_brcaM <- GDCprepare(query_met_brca, 
                         summarizedExperiment = FALSE, #tentativa de otimizar o código
                         remove.files.prepared = FALSE,
                         save = FALSE, 
                         save.filename = "BRCA_methylation_data.rda") 

  #LUAD
LUAD <- c("TCGA-LUAD") #coloque aqui os TCGA escolhidos
query_met_luad <- GDCquery(
  project = LUAD,
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450" 
)

GDCdownload(query = query_met_luad, method = "api", files.per.chunk = 1) #deixou mais leve pra baixar
data_luadM <- GDCprepare(query_met_luad) 

  #COAD
COAD <- c("TCGA-COAD") #coloque aqui os TCGA escolhidos
query_met_coad <- GDCquery(
  project = COAD,
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)

GDCdownload(query = query_met_coad, method = "api", files.per.chunk = 1) #deixou mais leve pra baixar
data_coadM <- GDCprepare(query_met_coad) 
