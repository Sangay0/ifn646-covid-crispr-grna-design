##Limma zoom checking 

library(edgeR)
library(limma)

ge <- DGEList(counts = as.matrix(data))
design <- model.matrix(~ 0 + sampletype, data = meta)
colnames(design) <- sub("sampletype", "", colnames(design))
keep <- filterByExpr(dge, design)
table(keep)

dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
dge_v <- voom(dge, design, plot = TRUE)

contrast_matrix <- makeContrasts(SARS - mock, levels = design)
contrast_matrix

fit <- lmFit(dge_v, design)
fit_2 <- contrasts.fit(fit, contrast_matrix)
fit_2 <- eBayes(fit_2)
results_limma <- topTable(fit_2, coef = 1, number = Inf, sort.by = "p")
head(results_limma)

sig_limma <- results_limma %>% filter(adj.P.Val < padj.cutoff & abs(logFC) > lfc.cutoff)
nrow(sig_limma)

sum(results_limma$P.Value < 0.05 & abs(results_limma$logFC) > lfc.cutoff)

comb <- sigSARS %>%
  +     select(gene, deseq2_log2FC = log2FoldChange, deseq2_padj = padj) %>%
  +     left_join(limma_tb, by = "gene") %>%
  +     mutate(direction = ifelse(deseq2_log2FC > 0, "up", "down"),
               +            consensus = !is.na(limma_p) & limma_p < 0.05 & abs(limma_logFC) > lfc.cutoff &
                 +                sign(limma_logFC) == sign(deseq2_log2FC))
consensus_genes <- comb %>% filter(consensus) %>% arrange(deseq2_padj)
consensus_genes %>% select(gene, deseq2_log2FC, deseq2_padj, limma_p) %>% print(n = Inf)

gene_info <- count_new_data %>%
  select(gene_id, gene_name, location) %>%
  distinct()
consensus_genes <- consensus_genes %>%
  left_join(gene_info, by = c("gene" = "gene_id"))

consensus_genes %>%
  filter(grepl("^chr[12]:", location)) %>%
  select(gene_name, location, deseq2_padj, limma_p)



