library("DESeq2")
library("ggplot2")
library("pheatmap")
count_data <- read.csv("GSE198063_raw_count_table.csv")
col_data <- read.csv("GSE198063_col_table.csv")

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~bin)
dds

### Removing genes which have low count (<100)
keep <- rowSums(counts(dds)) >= 100
keep

### Modifying dds
dds <- dds[keep,]
dds

prdds <- DESeq(dds)
prdds

res <- results(prdds, alpha = 0.1)
res

### Removing NA values in res
res <- res[!is.na(res$padj),]
res

### Seperating genes which are differentially expressed (p-adjusted < 0.1)
DE_genes <- res[res$padj < 0.1,]
DE_genes

DE_genes_list <- DE_genes[, c("log2FoldChange", "padj")]
DE_genes_list

### MA plot
plotMA(res, ylim=c(-5,5))
dev.copy2pdf(file='MA_plot.pdf', out.type = 'pdf')

### Volcano plot
vol_data <- data.frame(logFC=res$log2FoldChange, adjP = res$padj)

vol_data$sig <- "not-significant"
vol_data$sig[which((vol_data$adjP < 0.05) & (vol_data$logFC > 1))] <- "up"
vol_data$sig[which((vol_data$adjP < 0.05) & (vol_data$logFC < -1))] <-  "down"

ggplot(data = vol_data,aes(x=logFC, y=-1*log10(adjP), color = sig)) + geom_point()

dev.copy2pdf(file='vol_plot.pdf', out.type = 'pdf')

### Hierarchical Heatmap
vst_res <- vst(dds, blind = TRUE)
vsd_mat_res2 <- assay(vst_res)
vsd_cor_res2 <- cor(vsd_mat_res2)

### Plot the heatmap
pheatmap(vsd_cor_res2)

dev.copy2pdf(file='Heatmap.pdf', out.type = 'pdf')

### Dispersion plot
res_disp <- DESeq(dds)
plotDispEsts(res_disp)

dev.copy2pdf(file='disp_plot.pdf', out.type = 'pdf')




