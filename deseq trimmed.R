setwd("D:\\Bioinfo Project\\bioconductor\\deseq trimmed")

# BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)


#input of the count matrix as csv file
Matrix = read.csv("BLCA_HTseq count.csv")  
attach(Matrix)

#setting the row names (genes) for each samples

BRCA_matrix = Matrix[,-1]
rownames(BRCA_matrix) = Matrix[,1]

#input of the metadata
Metadata = read.csv("Metadata_tumor.csv", row.names = 1)

#if the row names of metadata is in same seq as column of matrix
all(rownames(Metadata) == colnames(BRCA_matrix))

# Create DESeq object
dds_wt <- DESeqDataSetFromMatrix(countData = round(BRCA_matrix),
                                 colData = Metadata,
                                 design = ~ Group)
View(counts(dds_wt))
#Count normalization---> calculation
dds_wt <- estimateSizeFactors(dds_wt)
sizeFactors(dds_wt)

#Count normalization---> extraction
normalized_wt_counts <- counts(dds_wt, normalized=TRUE)
View(normalized_wt_counts)
write.csv(normalized_wt_counts,'bladder_tumor_normalized.csv')

#------------------------------------------------------------------
# # Calculating mean for each gene (each row)
# mean_counts <- apply(PRAD_matrix[, 496:547], 1, mean)
# View(mean_counts)
# 
# # Calculating variance for each gene (each row)
# variance_counts <- apply(PRAD_matrix[, 496:547], 1, var)
# 
# # Creating data frame with mean and variance for every gene
# df <- data.frame(mean_counts, variance_counts)
# ggplot(df) +
#   geom_point(aes(x=mean_counts, y=variance_counts)) +
#   scale_y_log10() +
#   scale_x_log10() +
#   xlab("Mean counts per gene") +
#   ylab("Variance per gene")
#-------------------------------------------------------------------

# Plot dispersion estimates
# plotDispEsts(dds_wt)

#===================================================================

# Run analysis
dds_wt <- DESeq(dds_wt)

results(dds_wt, alpha = 0.05) #sample result

wt_res <- results(dds_wt,
                  contrast = c("Group", "C1",
                               "C2"),
                  alpha = 0.05,
                  lfcThreshold = 1)
wt_res

#DESeq2 LFC shrinkage
# plotMA(wt_res, ylim=c(-8,8))
# 
# wt_res <- lfcShrink(dds_wt,
#                     contrast=c("Group", "TP", "NT"),
#                     res=wt_res)
# plotMA(wt_res, ylim=c(-8,8))

#DESeq2 results table
mcols(wt_res)
head(wt_res, n=10)
summary(wt_res)

result = data.frame(wt_res)

write.csv(result, file = 'deg_C1_vs_C2_lfc1.csv')

q(save = "no")
