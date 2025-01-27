#!/usr/bin/env Rscript

#Taking in arguments
library("optparse")
 
option_list = list(
  make_option(c("-m", "--mat"), type="character", default=NULL, 
              help="matrix or data with headers [default %default]", metavar="filetype"),
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="name of analysis [default %default]", metavar="filetype"),
  make_option(c("-c", "--color"), type="character", default=NULL,
              help="colors used hex codes with # [default %default]", metavar="filetype"),
  make_option(c("-e", "--exp"), type="integer", default=NULL,
              help="column for the X axis [default %default]", metavar="filetype"),
  make_option(c("-r", "--rNMP"), type="integer", default=NULL,
              help="column for Y axis [default %default]", metavar="filetype"),
  make_option(c("-s", "--switch"), action = "store_true", default = FALSE,
              help="to switch the direction of x axis (mulitple by -1 for the Xaxis values) [default %default]", metavar="filetype"),
  make_option(c("-a", "--abs"), action = "store_true", default = FALSE,
              help="take absolute (mulitple by -1 for the Xaxis values) [default %default]", metavar="filetype"),
  make_option(c("-p", "--perc"), action = "store_true", default = FALSE,
              help="group based on percentiles with equal number of datapoints in each group [default %default]", metavar="filetype"),
  make_option(c("-y", "--y_max"), type="integer", default=4, 
              help="max y limit of plot [default %default]", metavar="integer"),
  make_option(c("-o", "--output_prefix"), type="character", default="out", 
              help="output file name [default %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Checking for required input
writeLines("\n...Checking input...\n")
if (is.null(opt$mat)){
  print_help(opt_parser)
  stop("Please specify input matrix.n", call.=FALSE)
}
m=opt$mat #mat="sorted_nucl_mono_0"

if (is.null(opt$name)){
  print_help(opt_parser)
  stop("Please specify name.n", call.=FALSE)
}
name=opt$name #mat="sorted_nucl_mono_0"


if (is.null(opt$color)){
  print_help(opt_parser)
  stop("Please specify color.n", call.=FALSE)
}
col=opt$color #mat="sorted_nucl_mono_0"
out=opt$output_prefix
exp=opt$exp
rNMP=opt$rNMP
ymax=opt$y_max
#writeLines(paste("Bin width used:", opt$bin_width))
#writeLines(paste("Maximum limit of yaxis:", opt$y_max))

#Usage:source(paste(getwd(),"/hotspot_analysis.R", sep=""))
suppressMessages(library(Hmisc, quietly = T) )
suppressMessages(library(gsubfn, quietly = T))
suppressMessages(library(ggplot2, quietly = T))
suppressMessages(library(ggpubr, quietly = T))
suppressMessages(library(tidyr, quietly = T))
suppressMessages(library(plyr, quietly = T))
suppressMessages(library(tools,quietly = T)) #for getting file basename

getstats <- function(df,name,featurex='featurex', featurey='featurey', bins=10, percentile=TRUE){
  data=df[c(featurex,featurey)]
  colnames(data)<-c('x', 'y')
  #Creating groups
  data=na.omit(data) ; if (max(data['x']) > 50) {data['x']=log2(data['x']+1); data=data[data[1]>0,]}   
  if (percentile) {data$groups<-as.numeric(cut2(data$x, g=bins)); data$subgroups=as.numeric(ecdf(data$x)(data$x))} else {data$groups=ceiling(data$x/max(data$x, na.rm=TRUE)*bins);data$subgroups=data$x}
  #Creating labels for groups
  range<-as.data.frame(gsub(')',']', levels(cut2(seq(0,100), g=bins))))
  range<-separate(range,1,c("min","max"),sep=',')
  range[,1]<-gsub('[^[:alnum:] ]', " ",range[,1])
  range[,2]<-gsub('[^[:alnum:] ]', " ",range[,2])
  range[1:(bins-1),2]<-as.numeric(range[1:(bins-1),2])-1
  range[,3]=paste(range[,1],range[,2], sep="-")
  label=as.list(range[,3])
  #Getting Stats
  
  stats<-ddply(data, c("groups"), summarise,
               N    = length(y),
               Featurex = mean(x),
               ymean = mean(y),
               ymedian = median(y),
               sd   = sd(y),
               se   = sd / sqrt(N),
               min   = min(y),
               max   = max(y),
               CI25 =quantile(y, probs = 0.25),
               CI75 =quantile(y, probs = 0.75))
  
  stats$low<-stats$ymean-stats$se
  #stats$low<-stats$CI25
  stats$high<-stats$ymean+stats$se
  #stats$high<-stats$CI75
  stats$name<-name
  #if (percentile) {stats$label<-range[,3]} else {stats$label<-stats$groups;stats$Featurex<-stats$groups}
  stats$label <- stats$groups
  stats<-stats[order(stats$Featurex),]
  coeff<-lm(ymedian ~ Featurex , data = stats)
  #print(summary(coeff))
  #print(summary(coeff)$adj.r.squared)
  p<-cor.test(as.numeric(stats$ymean), as.numeric(stats$Featurex), method='pearson')
  s<-cor.test(as.numeric(stats$ymean), as.numeric(stats$Featurex), method='spearman')
  #print(p$estimate)
  #print(p$p.value)
  #print(s$estimate)
  #print(s$p.value)
  return(list(data,label, stats, summary(coeff),summary(coeff)$adj.r.squared,p$estimate, p$p.value, s$estimate, s$p.value))
}

mat=read.table(m, sep="\t", header=TRUE)
if (opt$switch) {mat[exp] = (mat[exp]*-1)}
if (opt$abs) {mat[exp] = abs(mat[exp])}
if (opt$perc) {
  stats=getstats(mat, name, colnames(mat)[exp],colnames(mat)[rNMP],bins=10,percentile=TRUE)
} else {
  stats=getstats(mat, name, colnames(mat)[exp],colnames(mat)[rNMP],bins=10,percentile=FALSE)
}
#stats[3][[1]] = stats[3][[1]][-1,]; rownames(stats[3][[1]])<- stats[3][[1]]$label

plot_stats=stats[3][[1]]; write.table(plot_stats,paste(out,'/',name,"_stats.tsv",sep=""), append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)
data=stats[1][[1]]; write.table(data,paste(out,'/',name,"_data.tsv",sep=""), append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)

png(paste(out,'/',name,"bins.png",sep=""), width = 2, height = 2.1, units = "in", res=600, type="cairo", bg="transparent")
data$groups<-as.factor(data$groups)
dodge=position_dodge(width = 0.8)
ymin=0; if(min(data$y)<0) {ymin=-ymax}
p<-ggplot(data, aes(x=groups,y=y, fill=groups))+
  geom_violin(width=1,color="grey30",lwd = 0.1,position=dodge)+
  geom_boxplot(width=0.2, color="grey30", lwd = 0.1,alpha=0.5, outlier.shape = NA, position=dodge)+
  scale_x_discrete( labels=stats[3][[1]]$FeatureX) + #Alternate to have FeatureX: stats[3][[1]]$Featurex or #stats[2][[1]]
  scale_y_continuous(limits = c(ymin,ymax),n.breaks=4, name = "rNMP Enrichment")+
  theme_classic(base_size=20)+
  geom_abline(intercept = stats[4][[1]]$coefficients[1], slope = stats[4][[1]]$coefficients[2], color = "black", size = 0.2)+
  stat_cor(label.x = 3, label.y = 20)+
  #stat_summary(fun = "mean", geom="text", size = 2, hjust=0,vjust=20, label=round(stats[[3]]$ymean,2))+
  stat_summary(fun = "median", geom="text", size = 1.2, vjust=5,label=round(stats[[3]]$ymean,2))+
  #scale_fill_manual(values=get_brewer_pal("Spectral", n= bins, contrast = c(0.3, 0.8), stretch = F, plot = F))+
  #scale_fill_manual(values=colorRampPalette(c("#F8C6CB","#CB4A42"),  bias=1)(10))+
  scale_fill_manual(values=colorRampPalette(c(colorRampPalette(c("#FFFFFF",col),  bias=1)(11)[4],col),  bias=1)(nrow(data)))+
  #scale_fill_manual(values=get_brewer_pal("BuPu", n= bins, contrast = c(0.3, 0.6), stretch = F, plot = F))+
  #guides(colour = "colorbar", size = "legend", shape = "none")+
  #guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10))+
  ggtitle(bquote(italic(adj.R)^2 == .(format(stats[[5]], digits = 3))~";"~italic(p) == .(format(stats[[7]], digits = 3))~";"~italic(slope) == .(format(stats[[4]]$coefficients[2], digits = 3))))+
  theme(legend.position = "none",panel.border = element_blank(),
        plot.title = element_text(color="black",size=7,hjust=0.5),
        axis.title= element_blank(),
        axis.ticks = element_line(linewidth=0.3,color = "black"),
        axis.ticks.length = unit(.05, "cm"),
        axis.line = element_line(linewidth=0.3,color = "black"),
        axis.text = element_text(color="black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        strip.background = element_blank(),
        axis.text.x=element_text(size=8,  angle = 90, hjust=0.75),
        #axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

print(p)
dev.off()

mat=na.omit(mat); mat=mat[mat[exp]>0,]
if (max(mat[exp])>50) {mat[exp]=log2(mat[exp]+1)}
scatter_stats<-lm(rNMP ~ exp , data = mat)
colnames(mat)[exp]='exp'
png(paste(out,'/',name,"_scatter.png",sep=""), width = 4, height = 4, units = "in", res=600, type="cairo", bg="transparent")
ggscatter(mat, x = 'exp' , y = colnames(mat)[rNMP], size = 0.65, color = col, margin.params = list(fill = "grey30"), add = "reg.line", add.params = list(color='grey30', fill='grey30'),conf.int = TRUE, cor.coef = T, cor.coeff.args = list(method='spearman',label.x.npc = 0,label.y.npc = 1, size=5)) + theme_classic(base_size=20) + theme(axis.ticks = element_line(size=0.4,color = "black"),  axis.line = element_line(size=0.4,color = "black"), axis.title=element_blank(),)
dev.off()


writeLines("............... Stats .............")
writeLines(paste("lmcoeff","Pearson'sR", "Pearson'sPvalue", "Spearman'srho", "SpearmansPvalue", sep="\t"))
print(unlist(stats[5:9]))
corr=data.frame(stats=c("lmcoeff","Pearson'sR", "Pearson'sPvalue", "Spearman'srho", "SpearmansPvalue", "intercept", "slope", "scatter_intercept","scatter_slope"), value=c(unlist(stats[5:9]),stats[[4]]$coefficients[1:2],scatter_stats$coefficients[1:2]))
write.table(corr,paste(out,'/',name,"_corr.tsv",sep=""), append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)
writeLines("...................................")

