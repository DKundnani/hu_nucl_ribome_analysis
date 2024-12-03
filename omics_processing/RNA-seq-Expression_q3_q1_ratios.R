#PCA.R
library(dplyr)
library(factoextra)
library("ggplot2")
library("ggfortify")

colors <- c(
  #"CD4T"="#FF7F0E",
  #"hESC"="#2CA02C",
  "HEK_TPM" = "#D62728",
  
  "HEK293T.KO.T3.8."="#E377C2",
  "HEK293T.KO.T3.17."="#9467BD",
  "HEK293T.WT."="#D62728"
)

colors <- c(
  "HEK293T.WT.1"="#D62728",
  "HEK293T.KO.T3.8.1"="#E377C2",
  "HEK293T.KO.T3.17.1"="#9467BD",
  "HEK293T.WT.2"="#D62728",
  "HEK293T.KO.T3.8.2"="#E377C2",
  "HEK293T.KO.T3.17.2"="#9467BD"
)

setwd('C:/Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/RNA-seq') #nolint
HEKExp_open<-read.table("hek_TPM.geneid.TSV", sep="\t", header=TRUE)

HEKExp<-read.table("C:/Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/RNA-seq/HEK_tpm.csv", sep=",", header=TRUE); out='exp'
colnames(HEKExp_open)<-c("Geneid","HEK_TPM1","HEK_TPM2","HEK_TPM3")
HEKExp_merged=merge(HEKExp_open, HEKExp, by="Geneid", all.x=TRUE) %>% distinct()


pca_mat <- na.omit(HEKExp_merged %>% select(matches("HEK")))
pca_mat <- log2(pca_mat)
pca_mat <- pca_mat[!apply(pca_mat, 1, function(row) any(is.infinite(row))), ]
pca_mat <- as.data.frame(lapply(pca_mat, function(x) x - median(x)+mean(apply(pca_mat, 2, median))))

order <- data.frame(
  colnames = colnames(pca_mat),
  Group = unlist(strsplit(colnames(pca_mat), "1$|2$|3$")),
  col = c("#D62728","#D62728","#D62728","#9467BD","#9467BD","#E377C2","#E377C2","#D62728","#D62728")
)

res.pca <- prcomp(t(pca_mat))
pc=1

png(paste(out,"_pca",pc,".png", sep=""), width =7, height = 5, unit ='in', res=300)

autoplot(res.pca, data = order,color='Group', size = 3, x=pc, y=(pc+1))
  scale_fill_manual(values = colors)+
  scale_color_manual(values=colors)+
  theme_classic()+
  theme(#legend.position = "none",
    axis.ticks = element_line(size=0.8,color = "black"),
    axis.ticks.length=unit(.25, "cm"),
    axis.line = element_line(size=0.8,color = "black"),
    axis.text = element_text(size=25,color="black"),
    plot.margin = unit(c(0.0, 0.0, 0.0, 0.0), "cm"), #t, r, b, l
    axis.title=element_text(size=25))
dev.off()


#######################################
library(tidyr)
library(stringr)
library(ggplot2)

mat_centered=pca_mat

mat <- mat_centered %>%
  gather(key = "Library", value = "TPM")

mat$Library <- factor(mat$Library, levels = c("HEK_TPM1","HEK_TPM2","HEK_TPM3","HEK293T.WT.1", "HEK293T.WT.2", "HEK293T.KO.T3.8.1","HEK293T.KO.T3.8.2","HEK293T.KO.T3.17.1","HEK293T.KO.T3.17.2"))

mat$group <- str_replace(mat$Library,"1$|2$|3$", "" )
#mat$group <- mat$Library

png("boxplot.png", width = 6, height = 6, units = "in", res=300, type="cairo", bg="transparent")
dodge=position_dodge(width = 0.8)

ggplot(mat, aes(x=Library,y=TPM,fill = group))+
  geom_violin(width=1,color=NA,lwd = 0.3,alpha=0.25,position=dodge, trim = FALSE)+
  geom_boxplot(width=0.3, lwd = 0.5,alpha=0.5,outlier.size= 0.8)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 12))+
  scale_y_continuous(n.breaks=6)+
  #stat_summary(fun = "mean", geom="text", size = 2, vjust=-50, label=round(gene_stats[[1]]$mean,2))+
  scale_fill_manual(values = colors)+   # Custom colors
  theme_classic(base_size=20)+
  theme(axis.title=element_blank(),
      axis.ticks = element_line(size=0.6,color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.line = element_line(size=0.6,color = "black"),
      axis.text = element_text(color="black"),
      panel.background = element_rect(fill = "transparent",colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "transparent",colour = NA),
      strip.background = element_blank(),
      axis.text.x=element_blank(),legend.position = "none",panel.border = element_blank(),
      #axis.text.y=element_blank(),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
dev.off()
'''
library(ggridges)

png("boxplot.png", width = 3.5, height = 3, units = "in", res=600, type="cairo", bg="transparent")
dodge=position_dodge(width = 0.8)

mat$Library <- factor(mat$Library, levels = c("HEK293T.KO.T3.17.1","HEK293T.KO.T3.17.2", "HEK293T.KO.T3.8.1","HEK293T.KO.T3.8.2","HEK293T.WT.1", "HEK293T.WT.2"))

ggplot(mat, aes(x=TPM,group = Library,color=group))+
  #geom_density_ridges(scale = 0.9, rel_min_height = 0.05, alpha = 0.5,stat="binline",bins=25)+
  #geom_density_ridges(scale = 0.9, rel_min_height = 0.05, alpha = 0.5)+
  geom_density(alpha=0.2, size = 0.9,adjust = 0.9)+ 
  scale_x_continuous(limits = c(-10, 6), breaks = seq(-10, 6, 2))+
  scale_colour_manual(values = colors)+
  theme_classic(base_size=10)# Custom colors
theme(legend.position = "none",panel.border = element_blank(),
      axis.title=element_blank(),
      axis.ticks = element_line(size=0.6,color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.line = element_line(size=0.6,color = "black"),
      axis.text = element_text(color="black"),
      panel.background = element_rect(fill = "transparent",colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "transparent",colour = NA),
      strip.background = element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

out=data.frame()

for (e in seq(1,ncol(mat_centered))) { 
  df=mat_centered[!is.na(mat_centered[,e]>0),]
  df[,paste(colnames(df)[e], "_exp", sep="")]<-cut(df[,e], breaks=quantile(df[,e], probs = c(0,0.25,0.75,1), na.rm = T), include.lowest = T, labels=c("Low", "Mod", "High"))
  for (exp in c("Low", "Mod", "High")) {
    m=m[m[paste(colnames(m)[e], "_exp", sep="")]== exp,]
    out[e,exp]=exp;out[e,exp]=exp;
    write.table(m[,c(2:4,1,5:ncol(m))],paste(colnames(df)[e], exp,".bed", sep=""),sep='\t', col.names=F,row.names=F,quote=F)
  }
}
'''

std_dev <- apply(mat_centered, 2, sd)
var <- apply(mat_centered, 2, var)
stats=apply(mat_centered, 2, quantile)
ratio=data.frame('values'=stats["75%",]/stats["25%",], 'cell_type' =str_replace(colnames(stats),"1$|2$|3$", "" ) )


ratio.summary <- ratio %>%
  group_by(cell_type) %>%
  summarise(
    mean_values = mean(values, na.rm = TRUE),
    sd = sd(values, na.rm = TRUE), # / sqrt(n()),
  )
colnames(ratio.summary)[2]<-"values"

png("ratios.png", width = 4, height = 3, units = "in", res=300, type="cairo", bg="transparent")
dodge=position_dodge(width = 0.8)
ratio$cell_type <- factor(ratio$cell_type, levels = c("HEK_TPM","HEK293T.WT.","HEK293T.KO.T3.8.","HEK293T.KO.T3.17."))
ratio.summary$cell_type <- factor(ratio.summary$cell_type, levels = c("HEK_TPM","HEK293T.WT.","HEK293T.KO.T3.8.","HEK293T.KO.T3.17."))

ggplot(ratio, aes(x=cell_type,y=values,group=cell_type,color=cell_type, fill=cell_type))+ 
  #geom_bar(stat = "identity",position= dodge,color = "black",  aes(color = cell_type), width = 0.8)+
  geom_point(size = 4,alpha = 0.6, position =  position_dodge2(w = 0.7), 
             aes(fill=cell_type, color=cell_type),pch=21, stroke = 1)+
  #geom_boxplot(width=0.5, lwd = 0.2,alpha=0.5,outlier.size= 0.4)+
  #geom_pointrange(data = ratio.summary, aes(ymin = values-sd, ymax = values+sd), size = 0)+
  geom_errorbar(data = ratio.summary, aes(ymin = values-sd, ymax = values+sd), color = "black", width = 0.6)+
  scale_y_continuous( limits = c(2.0,2.5),n.breaks=5)+
  #stat_summary(fun = "mean", geom="text", size = 2, vjust=-50, label=round(gene_stats[[1]]$mean,2))+
    scale_fill_manual(values =  c("#D62728","#D62728", "#E377C2","#9467BD"))+
  scale_color_manual(values =  c("#D62728","#D62728", "#E377C2","#9467BD"))+   # Custom colors
  theme_classic(base_size=20)+
  theme(axis.title=element_blank(),
        axis.ticks = element_line(size=0.6,color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(size=0.6,color = "black"),
        axis.text = element_text(color="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_blank(),
        axis.text.x=element_blank(),legend.position = "none",panel.border = element_blank(),
        #axis.text.y=element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
dev.off()


for (e in seq(1,ncol(mat_centered))) { 
  df=mat_centered[!is.na(mat_centered[,e]>0),]
  quantile(df[,e], probs = c(0,0.25,0.75,1), na.rm = T)
}

####################### AGS Exp
#AGSexp1<-read.table("C:/Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/RNA-seq/AGS_tpm.csv", sep=",", header=TRUE); rownames(AGSexp) = AGSexp$GeneID;AGSexp = AGSexp1

AGSexp<-read.table("C:/Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/RNA-seq/GSE57353_norm_counts_TPM_GRCh38.p13_NCBI.tsv", sep="\t", header=TRUE) ; rownames(AGSexp) = AGSexp$GeneID; AGSexp=AGSexp[,-1]
colnames(AGSexp)<-c("Control1-1", "Control1-2", "AGS1-1", "AGS1-2", "AGS2P1-1", "AGS2P1-2","AGS2P2-1", "AGS2P2-2", "AGS4P1-1", "AGS4P1-2", "AGS4P2-1", "AGS4P2-2","AGS5-1", "AGS5-2" )

library(DESeq2)
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID
merged=merge(AGSexp, annot[rownames(AGSexp),c('Symbol','EnsemblGeneID','GeneType','Length')], by=0, all=TRUE)
merged=merged[merged$GeneType == "protein-coding",-1]
rownames(merged) <- merged$Symbol
pca_mat <- na.omit(merged %>% select(matches("Control|AGS")))

#pca_mat <- na.omit(merged %>% select(matches("Control|AGS1|AGS5|AGS4P2")))
pca_mat <- na.omit(AGSexp %>% select(matches("GSM")))
pca_mat <- log2(pca_mat)
pca_mat <- pca_mat[!apply(pca_mat, 1, function(row) any(is.infinite(row))), ]
pca_mat <- as.data.frame(lapply(pca_mat, function(x) x - median(x)+mean(apply(pca_mat, 2, median))))
pca_mat <- pca_mat %>% select(matches("Control|AGS1|AGS5|AGS4P2"))

order <- data.frame(
  colnames = colnames(pca_mat),
  Group = colnames(pca_mat))

res.pca <- prcomp(t(pca_mat))
pc=1
autoplot(res.pca, data = order, size = 3, x=pc, y=(pc+1),   label = TRUE)

mat_centered=pca_mat
stats=apply(mat_centered, 2, quantile)
ratio=data.frame('values'=stats["75%",]/stats["25%",], 'cell_type' =str_replace(colnames(stats),".1$|.2$|.3$", "" ) )

ration_sum <- ratio %>%
  group_by(cell_type) %>%
  summarise(
    mean = mean(values, na.rm = TRUE),
    se = sd(values, na.rm = TRUE) / sqrt(n()),
  )


png("AGS4_ratios.png", width = 4, height = 2.5, units = "in", res=300, type="cairo", bg="transparent")
dodge=position_dodge(width = 0.8)
ratio$cell_type <- factor(ratio$cell_type, levels = c("Control1","AGS1","AGS5","AGS4P2"))

ggplot(ration_sum, aes(x=cell_type,y=mean,color=cell_type))+ 
  geom_point(size = 2.5)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, color = "black")+
  scale_y_continuous(limits = c(2.5,3.5),n.breaks=3)+
  #stat_summary(fun = "mean", geom="text", size = 2, vjust=-50, label=round(gene_stats[[1]]$mean,2))+
  scale_color_manual(values =  c("#D62728","#E377C2","#D62728","#9467BD"))+   # Custom colors
  theme_classic(base_size=20)+
  theme(axis.title=element_blank(),
        axis.ticks = element_line(size=0.6,color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(size=0.6,color = "black"),
        axis.text = element_text(color="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_blank(),
        axis.text.x=element_blank(),legend.position = "none",panel.border = element_blank(),
        #axis.text.y=element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
dev.off()




####################### Mouse AGS Exp
AGSexp1<-read.table("C:/Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/RNA-seq/Mice_AGS_G37S_RPKM_counts.csv", sep=",", header=TRUE,check.names=FALSE); rownames(AGSexp1) = AGSexp1$'Feature ID'; AGSexp1 = AGSexp1[,-1]

pca_mat <- na.omit(AGSexp1)
pca_mat <- log2(pca_mat)
pca_mat <- pca_mat[!apply(pca_mat, 1, function(row) any(is.infinite(row))), ]
pca_mat <- as.data.frame(lapply(pca_mat, function(x) x - median(x)+mean(apply(pca_mat, 2, median))))

order <- data.frame(
  colnames = colnames(pca_mat),
  Group = colnames(pca_mat))

mat_centered=pca_mat
stats=apply(mat_centered, 2, quantile)
ratio=data.frame('values'=stats["75%",]/stats["25%",], 'cell_type' = colnames(stats) )
ratio$cell_type <- factor(ratio$cell_type, levels = c("WT","G37S..","G37S.G37S"))

ration_sum <- ratio %>%
  group_by(cell_type) %>%
  summarise(
    mean = mean(values, na.rm = TRUE),
    se = sd(values, na.rm = TRUE) / sqrt(n()),
  )


png("G37S_ratios.png", width = 4, height = 2.5, units = "in", res=300, type="cairo", bg="transparent")
dodge=position_dodge(width = 0.8)

ggplot(ration_sum, aes(x=cell_type,y=mean,color=cell_type))+ 
  geom_point(size = 2.5)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, color = "black")+
  scale_y_continuous(limits = c(3.5,4.5),n.breaks=3)+
  #stat_summary(fun = "mean", geom="text", size = 2, vjust=-50, label=round(gene_stats[[1]]$mean,2))+
  scale_color_manual(values =  c("#D62728","#E377C2","#9467BD","#D62728"))+   # Custom colors
  theme_classic(base_size=20)+
  theme(axis.title=element_blank(),
        axis.ticks = element_line(size=0.6,color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(size=0.6,color = "black"),
        axis.text = element_text(color="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_blank(),
        axis.text.x=element_blank(),legend.position = "none",panel.border = element_blank(),
        #axis.text.y=element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
dev.off()

mat_centered=pca_mat
mat_centered <- mat_centered[is.finite(rowSums(mat_centered)),]
#mat_centered$Group.1=rownames(mat_centered); mat_agg = mat_centered
library(Hmisc); mat_centered$groups<-as.numeric(cut2(mat_centered$WT, g=25))
library(stats); mat_agg = aggregate(mat_centered[1:3], list(mat_centered$groups), mean)
mat=mat_agg %>% gather(key = "Library", value = "TPM",-Group.1 )
mat$Library <- factor(mat$Library, levels = c("WT","G37S..","G37S.G37S"))

ggplot(mat, aes(x = Library, y = TPM)) + 
  geom_point(aes(color = Library), alpha=0.4) +  
  #geom_jitter(aes(color = Library), width = 0.05, height = 0.01) +
  geom_line(aes(group = Group.1),alpha=0.4) +
  scale_color_manual(values = c("#D62728","#E377C2","#9467BD")) +
  theme_classic(base_size=20)+
  theme(axis.title=element_blank(),
        axis.ticks = element_line(size=0.6,color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(size=0.6,color = "black"),
        axis.text = element_text(color="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_blank(),
        axis.text.x=element_blank(),legend.position = "none",panel.border = element_blank(),
        #axis.text.y=element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))


