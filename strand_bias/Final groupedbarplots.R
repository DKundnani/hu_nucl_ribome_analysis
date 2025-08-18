
setwd('C:/Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/Strand_bias/template_non-template')

library(ggthemes)
library(lemon)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(tidyr)
library(stringr)
library(ggpubr)
library(matrixStats)
library(rstatix)
library(ggplot2)
#library(gg.gap)
library(scales)
library(stringr) #wrapping Labels
library(tidyr)
library(ggpubr)
#if(!require(weatherData)) {install.packages("ggbreak", repos = "http://cran.us.r-project.org")}
library(ggbreak) 
library(cowplot)
library(grid)
library(gridExtra) 
library("extrafont")


#################################################
#file='Template Non-template for high low mod expression - 1kbdown.csv'
file='Template Non-template for high low mod expression - 1kbdown avg.csv'
#file='Template Non-template for high low mod expression - 1kbdown med.csv'
#file='Template Non-template for high low mod expression - 1kbdown min.csv'
#file='Template Non-template for high low mod expression - 1kbdown max.csv'
#file='Template Non-template for high low mod expression - 500down.csv'
file='Template Non-template for high low mod expression - TSS - TTS.csv'

df<-as.data.frame(read.table(file, sep=",", header = T, quote="", fill=TRUE, check.names = F))
colnames(df)<-c("strand","Celltype","High", "Low", "Mod")
mat<-df %>%
  gather(key = "Exp", value = "EF", -c("strand","Celltype"))
mat$Exp <- factor(mat$Exp, levels = c("Low", "Mod", "High"))
mat$Celltype <- factor(mat$Celltype, levels = c( "CD4T","hESC","HEK293T","RNH2AKO-T3-8","RNH2AKO-T3-17"))

stat.test <- mat %>%
  group_by(Celltype, Exp) %>%
  wilcox_test(EF ~ strand) %>% 
  add_significance("p") %>%
  adjust_pvalue(method = "bonferroni")   %>%
  add_xy_position(x = "Celltype", dodge = 0.8)

diff=mat[mat$strand == "template","EF"]- mat[mat$strand == "non-template","EF"]

finalmat = cbind(mat[mat$strand == "non-template",c("Celltype","Exp")],EF=diff)
finalmat = finalmat[-c(24,74,128),] #Removing FS310 only for 1kb down of TSS
finalmat = finalmat[-which(rownames(finalmat) == "74"),]

stat.test <- finalmat %>%
  group_by(Celltype) %>%
  wilcox_test(EF ~ Exp) %>% 
  add_significance("p") %>%
  adjust_pvalue(method = "bonferroni")   %>%
  add_xy_position(x = "Celltype", dodge = 0.8)

#stat.test$p.signif= mystarformat(stat.test$p)
stat.test$p.adj.signif=stat.test$p.signif

print(stat.test)
svg(paste(file,"diff.svg",sep=""),width=length(unique(mat$Celltype))*1.4, height=6)

png(paste(file,"diff.png",sep=""),width=length(unique(mat$Celltype))*1.4, height=6,units= "in",  res=600)
obj=ggboxplot(finalmat, x = "Celltype", y = "EF", fill="Exp",
          palette = c("#f3a712","#ef5b5b","#a31621"),
          alpha=0.55, xlab="", ylab="Difference", outlier.shape = NA, color = 'black')+
  geom_point(finalmat,mapping=aes(Celltype,EF, color=Exp,fill=Exp),alpha =0.55, position=position_dodge(0.8),size=1.25,pch=21,stroke = 1)+
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-1,0.6) , breaks = round(seq(-1,0.5,0.5), digits = 2)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12))+
  stat_pvalue_manual( stat.test[stat.test$group1 == "Low" & stat.test$group2=="Mod",],label='{p.signif}', y.position = 0.5,hide.ns = TRUE,  color = "black")+ #,   label = "p.adj", tip.length = 0.01)+
  stat_pvalue_manual( stat.test[stat.test$group1 == "Low" & stat.test$group2=="High",],label='{p.signif}', y.position = 0.6, hide.ns = TRUE, color = "black")+ #,   label = "p.adj", tip.length = 0.01)+
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black"),
        axis.text.y = element_text(color = "black"))

p=ggpar(obj, legend=c("right"), legend.title = "")
print(p)
dev.off()

################################
bias=diff/(mat[mat$strand == "non-template","EF"]+ mat[mat$strand == "non-template","EF"])

finalmat = cbind(mat[mat$strand == "non-template",c("Celltype","Exp")],EF=bias)
finalmat=finalmat[-c(24,74,128),]

stat.test <- finalmat %>%
  group_by(Celltype) %>%
  #wilcox_test(EF ~ Exp) %>% 
  wilcox_test(EF ~ Exp, alternative = "less") %>% 
  add_significance("p") %>%
  adjust_pvalue(method = "bonferroni")   %>%
  add_xy_position(x = "Celltype", dodge = 0.8)

#stat.test$p.signif= mystarformat(stat.test$p)
stat.test$p.adj.signif=stat.test$p.signif

png(paste(file,"bias.png",sep=""),width=length(unique(mat$Celltype))*1.4, height=6,units= "in",  res=600)
obj=ggboxplot(finalmat, x = "Celltype", y = "EF", fill="Exp",
              palette = c("#f3a712","#ef5b5b","#a31621"),
              alpha=0.55, xlab="", ylab="Bias", outlier.shape = NA, color = 'black')+
  geom_point(finalmat,mapping=aes(Celltype,EF, color=Exp,fill=Exp),alpha =0.55, position=position_dodge(0.8),size=1.25,pch=21,stroke = 1)+
  geom_hline(yintercept = 0, linetype = "dashed") +
  #scale_y_continuous(limits = c(-0.25,0.5) , breaks = round(seq(-0.24,0.24,0.08), digits = 2)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12))+
  stat_pvalue_manual( stat.test[stat.test$group1 == "Low" & stat.test$group2=="Mod",],label='{p.signif}', y.position = 0.26,hide.ns = TRUE,  color = "black")+ #,   label = "p.adj", tip.length = 0.01)+
  stat_pvalue_manual( stat.test[stat.test$group1 == "Low" & stat.test$group2=="High",],label='{p.signif}', y.position = 0.29, hide.ns = TRUE, color = "black")+ #,   label = "p.adj", tip.length = 0.01)+
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black"),
        axis.text.y = element_text(color = "black"))

p=ggpar(obj, legend=c("right"), legend.title = "")
print(p)
dev.off()

################################

ratio=mat[mat$strand == "template","EF"]/ mat[mat$strand == "non-template","EF"]
finalmat = cbind(mat[mat$strand == "non-template",c("Celltype","Exp")],EF=ratio)
finalmat=finalmat[-c(24,76,128),]

stat.test <- finalmat %>%
  group_by(Celltype) %>%
  #wilcox_test(EF ~ Exp) %>% 
  wilcox_test(EF ~ Exp, alternative = "less") %>% 
  add_significance("p") %>%
  adjust_pvalue(method = "bonferroni")   %>%
  add_xy_position(x = "Celltype", dodge = 0.8)

#stat.test$p.signif= mystarformat(stat.test$p)
stat.test$p.adj.signif=stat.test$p.signif

png(paste(file,"ratio.png",sep=""),width=length(unique(mat$Celltype))*1.4, height=6,units= "in",  res=600)
obj=ggboxplot(finalmat, x = "Celltype", y = "EF", fill="Exp",
              palette = c("#f3a712","#ef5b5b","#a31621"),
              alpha=0.55, xlab="", ylab="Ratio", outlier.shape = NA, color = 'black')+
  geom_point(finalmat,mapping=aes(Celltype,EF, color=Exp,fill=Exp),alpha =0.55, position=position_dodge(0.8),size=1.25,pch=21,stroke = 1)+
  geom_hline(yintercept = 1, linetype = "dashed") +
  #scale_y_continuous(limits = c(-0.25,0.5) , breaks = round(seq(-0.24,0.24,0.08), digits = 2)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12))+
  stat_pvalue_manual( stat.test[stat.test$group1 == "Low" & stat.test$group2=="Mod",],label='{p.signif}', y.position = 1.5,hide.ns = TRUE,  color = "black")+ #,   label = "p.adj", tip.length = 0.01)+
  stat_pvalue_manual( stat.test[stat.test$group1 == "Low" & stat.test$group2=="High",],label='{p.signif}', y.position = 1.6, hide.ns = TRUE, color = "black")+ #,   label = "p.adj", tip.length = 0.01)+
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black"),
        axis.text.y = element_text(color = "black"))

p=ggpar(obj, legend=c("right"), legend.title = "")
print(p)
dev.off()

