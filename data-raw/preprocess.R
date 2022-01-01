## load in packages (commented out for markdown)
library(GAclustenhancer)
library(DESeq2)
library(SummarizedExperiment)

## set paths (change paths before using)
setwd("C:/Users/frank/OneDrive/Documents/GAclustenhancer")

in_counts <- "./data-raw/BeatAML_count.rds" #raw gene expression counts
in_l2fc <- "./data-raw/deseq_results_tcga_flt3.rds" #l2fc results calculated from DESeq2

## read in L2fc and subset to only significant genes
data_l2fc <- readRDS(in_l2fc)
data_l2fc <- data_l2fc[data_l2fc$padj < 0.05, ]
input_lfc <- data_l2fc$log2FoldChange
names(input_lfc) <- rownames(data_l2fc)

## read in count matrix and convert from ExpressionSet to Summarized Experiment
count_data <- readRDS(in_counts)
feature_data <- fData(count_data)

## convert to DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(exprs(count_data)), colData = pData(count_data), design = ~ 1)

## apply vst transform (more information: https://bioconductor.org/packages/devel/bioc/manuals/DESeq2/man/DESeq2.pdf)
gene_counts_vst_norm <- varianceStabilizingTransformation(dds, blind = TRUE)

## drop genes with duplicated gene names
gene_counts_vst_norm <- gene_counts_vst_norm[, !is.na(gene_counts_vst_norm$FLT3.ITD)]
gene_counts_vst_norm <- gene_counts_vst_norm[!duplicated(feature_data$gene_name), ]
feature_data <- feature_data[!duplicated(feature_data$gene_name), ]

## reassign rownames from ENSEMBL gene IDs to fData mapping
rownames(gene_counts_vst_norm) <- feature_data$gene_name

## subset only commongenes in each dataset
commongene <- intersect(rownames(gene_counts_vst_norm), names(input_lfc))

## subset gene counts to only include commongenes
gene_counts_vst_norm <- gene_counts_vst_norm[commongene, ]
input_lfc <- input_lfc[commongene]

## finalize data and transposed data
count_data <- as.data.frame(assay(gene_counts_vst_norm))
t_count_data <- as.data.frame(t(assay(gene_counts_vst_norm)))
