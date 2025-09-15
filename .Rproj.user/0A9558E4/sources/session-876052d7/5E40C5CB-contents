library(WGCNA)
library(DESeq2)
library(dplyr)
library(tibble)
library(GEOquery)
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(devtools)
library(gridExtra)
remotes::install_github("kevinblighe/CorLevelPlot")
library(CorLevelPlot)

data <- read.delim('GSE152418_p20047_Study1_RawCounts.txt', header = T)

geo_id <- "GSE152418"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))
head(phenoData)

data_modified <- data %>% 
  gather(key = 'samples', value = 'counts', -ENSEMBLID) %>% 
  mutate(samples = gsub('\\.', '-', samples)) %>% 
  inner_join(., phenoData, by = c('samples' = 'title')) %>% 
  select(1,3,4) %>% 
  spread(key = 'geo_accession', value = 'counts') %>% 
  column_to_rownames(var = 'ENSEMBLID')

data_modified %>%
  as.data.frame() %>% 
  rownames_to_column("ENSEMBLID") %>% 
  pivot_longer(cols = -ENSEMBLID, names_to = "sample", values_to = "counts") %>% 
  group_by(ENSEMBLID) %>%
  summarise(gene_average = mean(counts), gene_stdev = sd(counts)) %>% 
  ungroup() %>% 
  ggplot(aes(x = log10(gene_average), y = log10(gene_stdev))) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") +
  labs(x = "Gene count average (log10 scale)",
       y = "Gene count standard deviation (log10 scale)") +
  ggtitle("Mean - Standard deviation relationship\n(no variance stabilisation)")

gsg <- goodSamplesGenes(t(data_modified))
summary(gsg)
data_modified <- data_modified[gsg$goodGenes == TRUE,]

colData <- phenoData
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))

all(rownames(colData) %in% colnames(data_modified))
all(rownames(colData) == colnames(data_modified))

dds <- DESeqDataSetFromMatrix(countData = data_modified,
                              colData = colData,
                              design = ~ 1)

dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]

vsd = varianceStabilizingTransformation(object = dds75, 
                                        blind = TRUE,
                                        fitType = "parametric")

variance_stabilised_counts <- assay(vsd)

norm.counts <- assay(vsd) %>% 
  t()

variance_stabilised_counts %>% 
  as.data.frame() %>% 
  rownames_to_column("ENSEMBLID") %>% 
  pivot_longer(cols = -ENSEMBLID, names_to = "sample", values_to = "counts") %>% 
  group_by(ENSEMBLID) %>% 
  summarise(gene_average = mean(counts), gene_stdev = sd(counts)) %>% 
  ungroup() %>% 
  ggplot(aes(x = gene_average, y = gene_stdev)) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") +
  labs(x = "Gene count average (variance stabilised)", 
       y = "Gene count standard deviation (variance stabilised)") +
  ggtitle("Mean - Standard deviation relationship\n(after variance stabilisation)")

power <- c(c(1:10), seq(from = 12, to = 50, by = 2))
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 18
cor <- WGCNA::cor

bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

module_eigengenes <- bwnet$MEs
head(module_eigengenes)
table(bwnet$colors)

plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

traits <- colData %>% 
  select(47) %>%
  mutate(disease_state = ifelse(grepl('COVID-19', disease_state), 1, 0))
  
colData$severity <- factor(colData$severity, levels = c("Healthy", "Convalescent", "ICU", "Moderate", "Severe"))
severity.out <- binarizeCategoricalColumns(colData$severity,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)
traits <- cbind(traits, severity.out)
traits

nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[17:21],
             y = names(heatmap.data)[1:17],
             col = c("blue1", "skyblue", "white", "pink", "red"))

module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'turquoise') %>% 
  rownames()

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:10,1:10]

gene.signf.corr <- cor(norm.counts, traits$data.Severe.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)
