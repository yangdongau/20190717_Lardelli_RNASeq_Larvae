if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GO.db", "org.Hs.eg.db"))
BiocManager::install("AnnotationHub")
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("ensembldb")
BiocManager::install("rtracklayer")
BiocManager::install("WGCNA")
BiocManager::install("RUVSeq")

install.packages("stringr")
install.packages("tidyverse")
install.packages("magrittr")
install.packages("scales")
install.packages("pander")
install.packages("ggrepel")


