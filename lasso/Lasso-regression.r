setwd("D:\\Bioinfo Project\\bioconductor\\deseq trimmed\\lasso")

library(readxl)
library(glmnet) 
library(faraway)

library(dplyr)

# uni <- read.csv("D:\\Bioinfo Project\\bioconductor\\deseq2\\univariate lfc1.5 subtype.csv")
# uniinput <- read.csv("D:\\Bioinfo Project\\bioconductor\\deseq2\\survival lfc1.5 xena desq2.csv")
# uni <- uni %>% filter(p.value < 0.05)
# uniinput_2 <- cbind(uniinput[,2:3], uni %>% select(uniinput %in% uni$X))



Data = read.csv("genes for lasso lfc1.csv")
data.matrix(Data)
dim (Data)
Data = Data[,4:363]
set.seed(1)

X <- model.matrix(time ~ ., data=Data)[,-1]

Y <- Data[,"time"] 
Y=as.numeric(unlist(Y))



cv.lambda.lasso = cv.glmnet(x=X, y=Y, alpha = 1) 

plot(cv.lambda.lasso) 


plot(cv.lambda.lasso$glmnet.fit, xvar = "lambda", label=FALSE)

l.lasso.min <- cv.lambda.lasso$lambda.min
l.lasso.min 

lasso.model <- glmnet(x=X, y=Y,
                      alpha  = 1,  lambda = l.lasso.min)
                    
lasso.model$beta  
tfit <- glmnet(x=X, y=Y, nfold=10,lower.limits = -4, upper.limits = 4)
plot(tfit)

######
ols.model <- glm(time ~DSC1+IL9RP3+TCHH+MDGA2+DLX1+GGT2+GNLY+UGT2B4+GPC2+CLIC6+ZPBP2+AC016730.1+AC079949.1   
                 +CYP4F24P, data=Data)
summary(ols.model) 




