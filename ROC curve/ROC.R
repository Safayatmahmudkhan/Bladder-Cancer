setwd("D:\\Bioinfo Project\\bioconductor\\deseq trimmed\\ROC curve")
library(survival)
install.packages("timeROC")
library(timeROC)
data <- read.csv("signature gene subtype.csv")
ROC.DSST <- timeROC(T=data$time,delta=data$OS,
                  marker=data$Risk.Score,cause=1,
                  weighting="marginal",
                  times=c(1,3,5),ROC=T,iid= TRUE)
ROC.DSST 

three <- plot(ROC.DSST,time = 3, title = F,lwd=2)      
one <- plot(ROC.DSST,time = 1,add = T, col = "blue",lwd=2)
five <- plot(ROC.DSST,time = 5, col="green",add=T,lwd=2)
legend("right",c("1 year(0.6497)","3 years(.7394)","(5 years(.7365)"),col=c("blue","red", "green"),lwd=2)


plotAUCcurve(ROC.DSST,FP=2,conf.int = TRUE,
             conf.band = TRUE) 
?legend
