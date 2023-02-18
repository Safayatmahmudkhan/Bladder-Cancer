setwd("D:\\Bioinfo Project\\bioconductor\\deseq trimmed\\lasso")
library(dplyr)
library(ggplot2)
library(tidyverse)
risk <- read.csv("signature gene subtype.csv")

risk %>%
  ggplot( aes(x=Category, y=MDGA2)) +
  geom_boxplot(aes(fill=Category), outlier.shape = 15, outlier.colour = "black",width=0.5, alpha=.5) +
  
  geom_jitter(aes(color=Allele), size=0.6, alpha=0.9)+
  xlab("Allele")+theme_classic()
t.test(UGT2B4~Category,select)
select <- select(risk, c(2,5:10))

select %>%
  pivot_longer(-Category) %>%
  mutate(real = factor(Category)) %>%
  ggplot(aes(Category, value)) + 
  geom_boxplot(aes(fill=Category)) + 
  
  facet_wrap(~name)+ theme_bw()+ ylim(0,4)


select %>%
  pivot_longer(-Category) %>%
  mutate(real = factor(Category)) %>%
  ggplot(aes(Category, value)) + 
  geom_boxplot(aes(fill=Category)) + 
facet_grid(~name,scales = "free", space = "free_x")






######shapiro test#####

high <- read.csv("high genes.csv")
low <- read.csv("low genes.csv")

shapiro.test(low$UGT2B4)
wilcox.test(UGT2B4~Category,select)
wilcox.test(select$DSC1~Category,select)
wilcox.test(select$MDGA2~Category,select)
wilcox.test(select$DLX1~Category,select)
wilcox.test(select$GGT2~Category,select)
wilcox.test(select$GNLY~Category,select)
