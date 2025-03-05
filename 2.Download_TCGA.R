###########################
### DOWNLOAD TCGA DATA  ###
###########################

# Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)

# RNA-seq Data from TCGA

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
