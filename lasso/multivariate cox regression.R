setwd("D:\\Bioinfo Project\\bioconductor\\deseq trimmed")

library(readxl)
library(openxlsx)

sur = read.csv("survival lfc1 xena dseqcsv.csv")

dim(sur)
covariates = colnames(sur[,4:17])
# install.packages(c("survival", "survminer"))
library("survival")
library("survminer")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS.time, OS)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = sur)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res.cox <- coxph(Surv(OS.time, OS) ~ DSC1+IL9RP3+TCHH+MDGA2+DLX1+GGT2+GNLY+UGT2B4+GPC2+CLIC6+ZPBP2+AC016730.1+AC079949.1   
                 +CYP4F24P, data = sur )
summary(res.cox)










library(ggplot2)
library(RColorBrewer)
setwd("D:\\Bioinfo Project\\bioconductor\\deseq trimmed\\lasso")
risk <- read.csv("signature gene subtype.csv")
attach(risk)
as.matrix(risk)
#####rank plot####
ggplot2::ggplot(risk,aes(x=Rank, y=OS.time/365))+
  geom_point(aes(color=OS), size=2)+scale_color_brewer(palette = 'Set2')+geom_point(alpha=.1)+
  theme_bw()

?geom_point


shapiro.test(risk$OS.1..dead)
wilcox.test(risk$OS.1..dead)

chisq.test(risk,Category~OS.time)
