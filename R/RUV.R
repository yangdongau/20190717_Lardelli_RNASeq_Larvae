library(RUVSeq)

# Find some 'unchanged' genes. Let's just grab the lowest ranked 5000
genes4Control <- deGenes %>%
  arrange(desc(P.Value)) %>%
  dplyr::slice(1:5000) %>%
  .[["ensembl_gene_id"]]
length(genes4Control)

## Let's try to remove 1 unwanted source of variation
k <- 1
ruv <- RUVg(
  dgeList$counts, 
  cIdx = genes4Control,
  k = k, 
  isLog = FALSE,
  round = TRUE
)
dgeList$samples <- cbind(dgeList$samples, ruv$W)

## Check the PCA
pcaRUV <- ruv$normalizedCounts %>%
  cpm(log = TRUE) %>%
  t() %>%
  prcomp() 
pcaRUV$x %>%
  cbind(dgeList$samples[rownames(.),]) %>%
  ggplot(aes(PC1, PC2, colour = Genotype)) +
  geom_point() +
  stat_ellipse(aes(fill = Genotype), geom = "polygon", alpha = 0.2) +
  xlab(
    paste0(
      "PC1 (", percent(summary(pcaRUV)$importance[2, "PC1"]), ")"
    )
  ) +
  ylab(
    paste0(
      "PC2 (", percent(summary(pcaRUV)$importance[2, "PC2"]), ")"
    )
  )

## Setup a design matrix using each pool as it's own intercept, with the genotype as the common difference
expDesign <- model.matrix(~0 + pair + Genotype + W_1, dgeList$samples) %>%
  set_colnames(str_replace(colnames(.), pattern = "Genotyp.+", "Mutant"))

## Recalculate the dispersions
dgeList %<>%
  estimateGLMCommonDisp(expDesign) %>%
  estimateGLMTagwiseDisp(expDesign)

## Now run DE Analysos after RUV
minLfc <- 0.75
topTable <- dgeList %>%
  glmFit(expDesign) %>%
  glmLRT(coef = "Mutant") %>%
  topTags(n = Inf) %>%
  as.data.frame() %>%
  set_colnames(gsub("ID.", "", colnames(.))) %>%
  as_tibble() %>%
  mutate(
    DE = FDR < 0.05 & abs(logFC) > minLfc
  ) %>%
  dplyr::select(
    ensembl_gene_id, external_gene_name, logFC, logCPM, LR,
  PValue, FDR, DE
  ) %>%
  arrange(PValue)

## Check the MA plot
tested <- c("si:ch211-213a13.2", "si:ch211-11p18.6", "mdh1ab", "CABZ01034698.2")
topTable %>%
  ggplot(aes(logCPM, logFC)) +
  geom_point(aes(colour = DE), alpha = 0.5) +
  geom_smooth(se = FALSE) +
  geom_hline(yintercept = c(-1, 1)*minLfc, linetype = 2) +
  geom_text_repel(
    aes(label = external_gene_name, colour = DE),
    data = . %>%
      dplyr::filter(abs(logFC) > 1.2 | external_gene_name %in% tested),
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c("grey50", "red"))

## Check the volcano plot
topTable %>%
  ggplot(aes(logFC, -log10(PValue))) +
  geom_point(aes(colour = DE), alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1)*minLfc, linetype = 2) +
  geom_text_repel(
    aes(label = external_gene_name, colour = DE),
    data = . %>%
      dplyr::filter(abs(logFC) > 1.2 & DE | external_gene_name %in% tested),
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c("grey50", "red"))

topTable %>%
  ggplot(aes(PValue)) +
  geom_histogram(fill = "grey50", colour = "black", bins = 100)

deGenes %>%
  ggplot(aes(P.Value)) +
  geom_histogram(fill = "grey50", colour = "black", bins = 100)

## Compare the two lists
topTable %>%
  dplyr::select(ensembl_gene_id, RUV = logFC, DE) %>%
  left_join(deGenes) %>%
  dplyr::rename(Voom = logFC) %>%
  mutate(
    Status = case_when(
      Sig & DE ~ "Both",
      Sig & !DE ~ "Voom Only",
      !Sig & DE ~ "RUV Only",
      !Sig & !DE ~ "Not DE"
    )
  ) %>%
  ggplot(aes(Voom, RUV,)) +
  geom_point(aes(colour = Status), alpha = 0.5) +
  geom_abline(slope= 1, linetype = 2) +
  geom_vline(xintercept = c(-1, 1)*minLfc) +
  geom_hline(yintercept = c(-1, 1)*minLfc) +
  scale_colour_manual(values = c("green", "grey50", "red", "blue"))

## Try running GSEA
ranks <- topTable %>%
  with(structure(sign(logFC)*LR, names = ensembl_gene_id)) %>%
  sort()
set.seed(22)
# Run GSEA for hallmark
fgseaHallmark <- fgsea(hallmark, ranks, nperm=1e5) %>%
  as_tibble() %>%
  dplyr::rename(FDR = padj) %>%
  mutate(padj = p.adjust(pval, "bonferroni")) %>%
  dplyr::arrange(pval)

if (interactive()) grid::grid.newpage()
plotGseaTable(
  hallmark[dplyr::filter(fgseaHallmark, padj < 0.05)$pathway], ranks, fgseaHallmark, gseaParam = 0.5
)

## Try running without direction
ranks <- topTable %>%
  with(structure(sign(logFC)*LR, names = ensembl_gene_id)) %>%
  sort()
set.seed(22)
# Run GSEA for hallmark
fgseaHallmark <- fgsea(hallmark, ranks, nperm=1e5) %>%
  as_tibble() %>%
  dplyr::rename(FDR = padj) %>%
  mutate(padj = p.adjust(pval, "bonferroni")) %>%
  dplyr::arrange(pval)

if (interactive()) grid::grid.newpage()
plotGseaTable(
  hallmark[dplyr::filter(fgseaHallmark, padj < 0.5)$pathway], ranks, fgseaHallmark, gseaParam = 0.5
)


fgseaKEGG <- fgsea(kegg, ranks, nperm=1e5) %>%
  as_tibble() %>%
  dplyr::rename(FDR = padj) %>%
  mutate(padj = p.adjust(pval, "bonferroni")) %>%
  dplyr::arrange(pval)
