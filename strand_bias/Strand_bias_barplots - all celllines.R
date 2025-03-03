setwd('C:/Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/Strand_bias/template_non-template')

top='same_genes_pc_noXY_subtype_percent_normalized.tsv'
bottom='opp_genes_pc_noXY_subtype_percent_normalized.tsv'

order='files_short'
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


same=as.data.frame(read.table(top, sep="\t", fill = TRUE,header =T,stringsAsFactors=FALSE, quote="", check.names = FALSE)); same$Strand='Non-template'
opp=as.data.frame(read.table(bottom, sep="\t", fill = TRUE,header =T,stringsAsFactors=FALSE, quote="" , check.names = FALSE)); opp$Strand='Template'
merged=rbind(same,opp)

order=as.data.frame(read.table(order, sep="\t",header = F, quote=""))
colnames(order)=c("num", "Filename","Library","Cellline", "Rep")
order=order[order(order$num),]
order$lab=paste(order$Library,order$Cellline,order$Rep,sep="-")

group=gather(merged[,c(-2)], lib, Enrich,-c(Strand,Group.1))
for (lib in unique(group$lib)) {group[group$lib==lib,'Cellline']=order[order$lab==lib,'Cellline']}
colnames(group)[1]="Subtype"
group$Enrich=group$Enrich*2
ym=ceiling(ceiling(max(group$Enrich))/2)*2
ym=3
group=group %>% mutate(Cellline = fct_relevel(Cellline,unique(group$Cellline)))



png("All_cellline_strand_bias.png",width=6, height=6,units= "in", family = "Arial", res=600)
stat.test <- group %>%
  group_by(Cellline) %>%
  wilcox_test(Enrich ~ Strand, ref.group = "Non-template") %>% 
  add_significance("p") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_xy_position(x = "Cellline", dodge = 0.8)
  
obj=ggbarplot(group, x = "Cellline", y = "Enrich", fill="Strand",
              add = c("mean_se"),
              add.params=list(width=0.4,size=0.5),
              palette = c("#60C2B1","#A77AE6"),
              position = position_dodge(0.75), 
              xlab="", ylab="rNMP EF", 
              width = 0.75,
              size = 0.5) +
  geom_point(group,mapping=aes(Cellline,Enrich,color=Strand), position=position_dodge(0.75),size=1)+guides(color = FALSE) +
  scale_color_manual(values=c("black","black"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,ym), breaks = seq(0,ym,1)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 6))+
  #stat_pvalue_manual(stat.test,label='{p} {p.signif}', y.position = 2.2, hide.ns = TRUE) #,   label = "p.adj", tip.length = 0.01)
  stat_pvalue_manual(stat.test,label='{p.signif}', y.position = 1.75, hide.ns = TRUE,label.size = 8,bracket.size = 0.8,tip.length = 0.05)+ #,   label = "p.adj", tip.length = 0.01)+
  theme_classic(base_size = 25)+
  theme(legend.position = "none",
        axis.text.x=element_text(color="black",size=10,angle = 45, vjust = 1, hjust=1),axis.text.y=element_text(color="black",size=25),
        axis.ticks = element_line(size=0.8,color = "black"),
        axis.line = element_line(size=0.8,color = "black"),
        plot.margin = unit(c(0.5, 0.0, 0.0, 0.0), "cm") #t, r, b, l
    ) 
print(obj)
dev.off()
  


