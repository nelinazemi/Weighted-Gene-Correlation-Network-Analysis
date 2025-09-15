library(DESeq2)
library(tibble)
library(tidyverse)
library(dplyr)
library(EnhancedVolcano)
library(apeglm)
library(biomartr)
library(clusterProfiler)
library(tidyverse)
library(enrichplot)
suppressPackageStartupMessages(library("org.At.tair.db"))
library(biomaRt)
library(pheatmap)

xp_design <- read.delim2("experimental_design_modified.txt",
                         header = T,
                         stringsAsFactors = F,
                         colClasses = rep("character",4))

xp_design_filtered <- xp_design %>%
  filter(seed == "MgCl2", dpi == "7")

raw_counts <- read.delim2("counts.txt",
                          header = T,
                          stringsAsFactors = F) %>%
  column_to_rownames("gene")

raw_counts <- raw_counts[, xp_design$sample]

raw_counts_filtered <- raw_counts[, colnames(raw_counts) %in% xp_design_filtered$sample]

dds <- DESeqDataSetFromMatrix(countData = raw_counts_filtered,
                              colData = xp_design_filtered,
                              design = ~ infected)
dds$infected

dds <- DESeq(dds)

results <- results(dds)

all_genes_results <- results(dds, contrast = c("infected",
                                               "Pseudomonas_syringae_DC3000",
                                               "mock"))
head(all_genes_results)

all_genes_results %>%
  as.data.frame() %>% 
  filter(padj < 0.01) %>% 
  dim()

all_genes_results %>%
  as.data.frame() %>% 
  filter(padj < 0.001) %>% 
  dim()

hist(all_genes_results$padj, col="lightblue", main = "Adjusted p-value distribution")
hist(all_genes_results$pvalue, col="lightblue", main = "Undjusted p-value distribution")

diff_expressed_genes = all_genes_results %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(padj < 0.01) %>% 
  arrange(desc(log2FoldChange), 
          desc(padj))

head(diff_expressed_genes)

resLFC <- lfcShrink(dds = dds, 
                    res = all_genes_results,
                    type = "normal",
                    coef = "infected_Pseudomonas_syringae_DC3000_vs_mock")

EnhancedVolcano(toptable = resLFC,
                x = "log2FoldChange",
                y = "padj",
                lab = rownames(resLFC)
)

nrow(diff_expressed_genes)

diff_expressed_genes %>% 
  filter(log2FoldChange > quantile(log2FoldChange, c(0.75))) %>% # keeping fold changes above the 75th percentile
  dplyr::select(gene) %>% 
  write.table(., file = "diff_genes_for_agrigo.tsv", row.names = FALSE, quote = FALSE)



all_genes <- raw_counts <- read.delim2("counts.txt",
                                       header = T,
                                       stringsAsFactors = F)[,1]
head(all_genes)






biomartr::organismBM(organism = "Arabidopsis thaliana")

arabido_attributes = 
  biomartr::organismAttributes("Arabidopsis thaliana") %>% 
  filter(dataset == "athaliana_eg_gene")

arabido_attributes

attributes_to_retrieve = c("tair_symbol","entrezgene_id", "uniprotswissprot")

options(timeout = 300)

all_genes_annotated <- biomartr::biomart(
  genes = all_genes,
  mart = "plants_mart",
  dataset = "athaliana_eg_gene",
  attributes = attributes_to_retrieve,
  filters = "ensembl_gene_id"
)

options(timeout = 3000)

diff_expressed_genes$gene <- as.character(diff_expressed_genes$gene)

diff_expressed_genes_annotatedd <- biomartr::biomart(genes = diff_expressed_genes$gene,
                                                     mart       = "plants_mart",                 
                                                     dataset    = "athaliana_eg_gene",           
                                                     attributes = attributes_to_retrieve,        
                                                     filters =  "ensembl_gene_id" )  

diff_expressed_genes_annotatedd

ora_analysis_bp <- enrichGO(gene = diff_expressed_genes$entrezgene_id, 
                            universe = all_genes$entrezgene_id, 
                            OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                            keyType = "ENTREZID",
                            ont = "BP",              # either "BP", "CC" or "MF",
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            readable = TRUE, 
                            pool = FALSE)

ora_analysis_bp_simplified <- clusterProfiler::simplify(ora_analysis_bp)

write_delim(x = as.data.frame(ora_analysis_bp_simplified@result), 
            path = "go_results.tsv", 
            delim = "\t")

ora_analysis_bp_simplified@result[1:5,1:8]

