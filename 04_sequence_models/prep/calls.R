# DESeq2 calls for every consensus peak, exactly as in Step 3 (03_differential/differential_atac.ipynb, section 4):
# size factors x copy number as normalization factors, no GC correction, tested for a change bigger than 2-fold.
# Writes work/calls.tsv, used by prepare_bundle.py to pick the featured peaks.
suppressPackageStartupMessages(library(DESeq2))
peaks <- read.delim("../../03_differential/data/peaks.tsv.gz")
samples <- c("COLO320DM_REP1", "COLO320DM_REP2", "COLO320HSR_REP1", "COLO320HSR_REP2")
counts <- as.matrix(peaks[, samples]); rownames(counts) <- peaks$peak
sample_info <- data.frame(line = factor(c("DM", "DM", "HSR", "HSR"), levels = c("HSR", "DM")), row.names = samples)
dds <- DESeq(DESeqDataSetFromMatrix(counts, sample_info, design = ~ line), fitType = "local", quiet = TRUE)

copy_number <- cbind(peaks$copy_ratio_DM, peaks$copy_ratio_DM, peaks$copy_ratio_HSR, peaks$copy_ratio_HSR)
has_cn <- complete.cases(copy_number)
dds_cn <- DESeqDataSetFromMatrix(counts[has_cn, ], sample_info, design = ~ line)
factors <- matrix(sizeFactors(dds), nrow(counts), 4, byrow = TRUE)[has_cn, ] * copy_number[has_cn, ]
normalizationFactors(dds_cn) <- factors / exp(rowMeans(log(factors)))
dds_cn <- DESeq(dds_cn, fitType = "local", quiet = TRUE)

res <- results(dds_cn, alpha = 0.05, lfcThreshold = 1)
call <- ifelse(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1, "up in DM",
        ifelse(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1, "up in HSR", "not significant"))
x <- peaks[has_cn, ]
dir.create("work", showWarnings = FALSE)
write.table(data.frame(peak = x$peak, chr = x$chr, start = x$start, end = x$end, region_class = x$region_class,
                       log2FC = res$log2FoldChange, padj = res$padj, call = call),
            "work/calls.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(table(call))
