setwd("D:\\Bioinfo Project\\bioconductor\\Independent dataset\\Final 1")

library(readxl)
library(openxlsx)

sur = read.csv("chungbuk signature gene tumor tissue.csv")

dim(sur)
covariates = colnames(sur[,3:6])
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
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

write.csv(res,"univariate lfc1 subtype.csv", row.names = TRUE)


res.cox <- coxph(Surv(OS.time, OS) ~ MDGA2+DLX1+GNLY+UGT2B4
                 , data = sur )
summary(res.cox)




library(ggfortify)



fit <- survfit(Surv(OS.time, OS) ~ Category, data = sur)
ggsurvplot(fit, pval=TRUE, data = sur)
