library(limma)
library(edgeR)
library(AnnotationHub)
library(tidyverse)
library(magrittr)
library(scales)
library(pander)
library(ggrepel)
library(ensembldb)
library(rtracklayer)

library(reshape2)
library(readxl)

library(corrplot)
library(RColorBrewer)
library(grid)
library(scales)
library(knitr)
library(kableExtra)
library(pheatmap)
library(ggrepel)
library(GGally) #for ggnet

# Viewports for grid
vp_left <- viewport(x = 0, y = 0, width = 0.5, height = 1, just = c(0, 0))
vp_right <- viewport(x = 0.5, y = 0, width = 0.5, height = 1, just = c(0,0))


# Set working directory ---------------------------------------------------
setwd("~/Desktop/20190717_Lardelli_RNASeq_Larvae-master/R")

# Load file --------------------------------------------------------------
counts <- read_tsv("../2_alignedData/featureCounts/counts.out", skip = 1) 


colnames(counts) <- colnames(counts) %>%
  basename() %>%
  str_remove("LardelliAligned.sortedByCoord.out.bam")

counts
  
dgeList <- counts %>%
  dplyr::select(-(c("Chr", "Start", "End", "Strand", "Length"))) %>%
  as.data.frame() %>%
  column_to_rownames("Geneid") %>%
  DGEList() %>%
  calcNormFactors()

dgeList$samples

dgeList$samples$group <- colnames(dgeList) %>%
  str_extract("_(1|2)") %>%
  factor(levels = c("1", "2"))

ah <- AnnotationHub()
ah %>%
  subset(species == "Danio rerio") %>%
  subset(dataprovider == "Ensembl") %>%
  subset(rdataclass == "EnsDb")
ensDb <- ah[["AH64906"]]
ensDb
genesGR <- genes(ensDb)
granges(genesGR)
seqinfo(genesGR)

transcriptsGR <- transcripts(ensDb)
granges(transcriptsGR)
seqinfo(transcriptsGR)
mcols(transcriptsGR)
mcols(genesGR)
mcols(genesGR) <- mcols(genesGR)[c("gene_id", "gene_name", "gene_biotype", "entrezid")]
dgeList$genes <- genesGR[rownames(dgeList),]
dgeList[1:4,]
#
dgeList$counts %>% 
  rowSums() %>%
  is_greater_than(0) %>%
  table()

dgeList %>%
  cpm() %>%
  is_greater_than(1) %>%
  rowSums() %>%
  is_greater_than(3) %>%
  table()

genes2keep <- dgeList %>%
  cpm() %>%
  is_greater_than(1) %>%
  rowSums() %>%
  is_greater_than(3)
dgeFilt <- dgeList[genes2keep,] %>% calcNormFactors()

par(mfrow = c(1,2))
dgeList %>%
  cpm(log = TRUE) %>%
  plotDensities(legend = FALSE, main = "A. Before Filtering")
dgeFilt %>%
  cpm(log = TRUE) %>%
  plotDensities(legend = FALSE, main = "B. After Filtering")

# Library size ------------------------------------------------------------

dgeFilt$samples %>%
  ggplot(aes(group, lib.size, fill = group)) +
  geom_boxplot() +
  scale_y_continuous(labels = comma) +
  labs(x = "Timepoint", y = "Library Size") +
  theme_bw() 


# PCA ---------------------------------------------------------------------

pca <- dgeFilt %>%
  cpm(log = TRUE) %>%
  t() %>%
  prcomp() 
summary(pca)$importance %>% pander(split.tables = Inf)

plotly::ggplotly(
  pca$x %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    as_tibble() %>%
    dplyr::select(sample, PC1, PC2) %>%
    left_join(rownames_to_column(dgeFilt$samples, "sample")) %>%
    ggplot(aes(PC1, PC2, colour = group, label = sample)) +
    geom_point(size = 3) +
    theme_bw()
)

# Differential Expression -------------------------------------------------
voomData <- voom(dgeFilt)

topTable <- voomData %>% 
  lmFit() %>%
  eBayes %>%
  topTable(coef = "group24mth", n = Inf) %>%
  as_tibble()