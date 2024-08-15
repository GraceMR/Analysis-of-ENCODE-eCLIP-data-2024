###Gene ontology analysis###
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library("AnnotationDbi")
library(org.Hs.eg.db)
library(ggplot2)

df <- read.csv("data/DDX6_3UTR_binding_sites_and_transcript_info.csv", header = TRUE)
DDX6_genes <- df$external_gene_name
DDX6_genes <- unique(DDX6_genes)

df2 = bitr(DDX6_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

#'BIOLOGICAL PROCESS' GO TERMS

BPEGO = enrichGO(gene = df2$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.05, 
                  pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
BPEGO

BPSIMGO = simplify(BPEGO, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", 
                    semData = NULL)
nrow(BPSIMGO)

png(paste0("output/clusterProfiler_GO-BP_ORA_simplify", ".png"), width = 9, height = 6, units = "in", res = 300)

barplot(BPSIMGO, showCategory = 20) + 
  ggtitle(paste0("GO-BP ORA of DDX6-bound transcripts")) + 
  xlab("Enriched terms") + ylab("Count") +
  theme(axis.text.y = element_text(size = 9))

invisible(dev.off())

GO_BP_summary <- data.frame(BPEGO)
write.csv(GO_BP_summary, file = "output/GO_BP_Summary.csv", row.names=FALSE)

#'CELLULAR COMPONENT' GO TERMS

CCEGO = enrichGO(gene = df2$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "CC", pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
CCEGO

CCSIMGO = simplify(CCEGO, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", 
                   semData = NULL)
nrow(CCSIMGO)

png(paste0("output/clusterProfiler_GO-CC_ORA_simplify", ".png"), width = 9, height = 6, units = "in", res = 300)

barplot(CCSIMGO, showCategory = 20) + 
  ggtitle(paste0("GO-CC ORA of DDX6-bound transcripts")) + 
  xlab("Enriched terms") + ylab("Count") +
  theme(axis.text.y = element_text(size = 9))

invisible(dev.off())

GO_CC_summary <- data.frame(CCEGO)
write.csv(GO_CC_summary, file = "output/GO_CC_Summary.csv", row.names=FALSE)

#'MOLECULAR FUNCTION' GO TERMS

MFEGO = enrichGO(gene = df2$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "MF", pvalueCutoff = 0.05, 
                 pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
MFEGO

MFSIMGO = simplify(MFEGO, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", 
                   semData = NULL)
nrow(MFSIMGO)

png(paste0("output/clusterProfiler_GO-MF_ORA_simplify", ".png"), width = 9, height = 6, units = "in", res = 300)

barplot(MFSIMGO, showCategory = 20) + 
  ggtitle(paste0("GO-MF ORA of DDX6-bound transcripts")) + 
  xlab("Enriched terms") + ylab("Count") +
  theme(axis.text.y = element_text(size = 9))

invisible(dev.off())

GO_MF_summary <- data.frame(MFEGO)
write.csv(GO_MF_summary, file = "output/GO_MF_Summary.csv", row.names=FALSE)
