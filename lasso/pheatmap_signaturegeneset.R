setwd("D:\\Bioinfo Project\\bioconductor\\deseq trimmed\\lasso")
blcafpkm <- read.csv("heatmap 6 genes.csv",row.names = 1)


ZScore <- scale((blcafpkm))
write.csv(ZScore,"ZScore_blca.csv")
library(pheatmap)
ZScore <- read.csv("ZScore_blca.csv", row.names = 1)



annotationall <- read.csv("heatmap annotation.csv", row.names = 1)


####tonmoy code pheatmap#####

breaksList = seq(-2, 2, by = .1)

pheatmap(t(ZScore),  show_colnames = F, 
         show_rownames = T,
         fontsize_row = 4, fontsize_col = 4,
         angle_col = 90,
         # annotation_row = annotation,
         annotation_col = annotationall,
         gaps_col= 203,
         annotation_colors = NA,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         breaks = breaksList,
         colorRampPalette(c("#292C6D", "white", "#EC255A"))(length(breaksList)),
         cluster_cols = F,
         cluster_rows = F,
         border_color = NA)
