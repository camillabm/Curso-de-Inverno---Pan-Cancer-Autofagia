#######################
### ARGs definition ###
#######################

# ARGs table as csv
df <- read_excel("2.0/ARG.xlsx")
write.csv(df, "2.0/ARG.csv", row.names = FALSE)
ARG <- read.csv("2.0/ARG.csv")

# Ensemble ID annotation
library (biomaRt)

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", 
                   dataset = "hsapiens_gene_ensembl",
                   host = "https://grch37.ensembl.org") #Connection to ensembl

hgnc_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"), 
  filters = "hgnc_symbol",  # Agora filtramos pelo HGNC.code
  values = ARG$HGNC.code,  # Passamos os HGNC codes
  mart = ensembl
)

ARG <- merge(ARG, hgnc_mapping, by.x = "HGNC.code", by.y = "hgnc_symbol", all.x = TRUE) #Add HGNC symbol to results
ARG$Ensembl_ID <- NULL

ARG <- ARG[!is.na(ARG$ensembl_gene_id), ]

# File saving as .csv
write.csv(ARG, "2.0/ARG_v1.csv", row.names = FALSE)

