setwd("D:\\Bioinfo Project\\bioconductor\\deseq trimmed")
deseq <- read.csv("deg_C1_vs_C2_lfc1.csv")
library(ggplot2)
library(viridis)
# library(ggplot2)

 # install.packages("viridis")                                                                       

 library(ggrepel)
 
 
 
 ggplot2::ggplot(deseq,aes(x=log2FoldChange, y=-log10(padj)))+
   geom_point(aes(color=Category),alpha=0.4, size=2.75, 
show.legend = F)+scale_color_viridis(discrete = TRUE, option = "E")+theme_bw()+
   geom_vline(xintercept=c(-1, 1), col="black",linetype = "dashed") +
   geom_hline(yintercept=-log10(0.05), col="black",linetype = "dashed")
 
 ?geom_vline
 
 
 