df<-read.csv("data/UTR_overlapping_intervals.csv", header = TRUE)

#initially need to remove any numbers after the dots, and re-name the GENCODE_ID column to ensembl_transcript_id

df$GENCODE_ID <- gsub("\\.\\d+$", "", df$GENCODE_ID)
names(df)[names(df) == 'GENCODE_ID'] <- 'ensembl_transcript_id'

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library("BiocManager")

#BiocManager::install("biomaRt")

library("biomaRt")

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
res <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name", "description"),
  filters = "ensembl_transcript_id",
  values = df$ensembl_transcript_id,
  mart = ensembl)
res

#now need to join the data from res (external_gene_id, external_gene_name, and description) to the data in df- match according to transcript_id
#When merging the dataframes, make sure they have matching column names

#res2 <- res$external_gene_name
#write.csv(res2, "data/gene_list.csv")

df2 <- merge(df, res, by = "ensembl_transcript_id")

#order by start coordinate, re-order columns, and drop 'X'

library(dplyr)
df2 <- arrange(df2, start)
df2$ensembl_gene_id <- as.factor(df2$ensembl_gene_id)
df2 <- relocate(df2, ensembl_gene_id, .before = X)

#export df2 as .csv

write.csv(df2,"output/overlapping_intervals_gene_names.csv", row.names = FALSE)
