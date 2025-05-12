#plot_REZ
library(karyoploteR)
library(bedr)
library(dplyr)
library(data.table)
library(regioneR)
library(stringr)
library(zoo)
library(httpgd)
library(GenomicRanges)
#detach('cli', unload=TRUE)
#remove.packages("xfun")
#install.packages("xfun")
#install.packages("cli", dependencies = TRUE, update = TRUE)

colors <- c(
  "CD4T"="#FF7F0E",
  "hESC"="#2CA02C",
  "HEK293T"="#D62728",
  "RNH2A-KO T3-8"="#E377C2",
  "RNH2A-KO T3-17"="#9467BD"
)

plot <- function(kp) {
  kpAxis(kp,cex = 0.8,r0=0.2, r1=0.85, label.margin=-20e5,tick.len=10e5,ymin=0, ymax=4, numticks=2, data.panel=1)
  kpAxis(kp,cex = 0.8,r0=0.2, r1=0.85, label.margin=-20e5,tick.len=10e5,ymin=0, ymax=4, numticks=2, data.panel=2)
  
  #kpAbline(kp, h=1, lty=2,col="#666666",ymin=ymin, ymax=ymax)
  
  kpPlotRegions(kp, data = REZ.split$`+`, r0 = 0.2, r1=0.85, col=transparent('#86ECFA', amount=0.4),border=NA, data.panel=1)
  kpPlotRegions(kp, data = REZ.split$`-`, r0 = 0.2, r1=0.85, col=transparent('#86ECFA', amount=0.4),border=NA, data.panel=2)
  
  kpPlotRegions(kp,annot.range,r0 = -0.17, r1=0.17,col=transparent('#9E1717', amount=0.90),avoid.overlapping=FALSE,data.panel=2)
  #kpPlotRegions(kp1,annot.range,r0 = 0, r1=0.17,col=transparent('#9E1717', amount=0.90),data.panel=1)
  #kpPlotRegions(kp1,annot.range,r0 = 0, r1=0.17,col=transparent('#9E1717', amount=0.90),data.panel=2)
  
  for (n in names(colors)){
    EF.split=EF.range[which(elementMetadata(EF.range)[,'Celltype'] == n)]
    gr_pos=sort(split(EF.split, strand(EF.split))$`+`); gr_neg=sort(split(EF.split, strand(EF.split))$`-`)
    kpLines(kp, data=gr_pos, y=gr_pos$EF, data.panel = 1, col=unname(colors[n]), ymin=0, ymax=4, r0=0.2, r1=0.85,lwd=1)
    kpLines(kp, data=gr_neg, y=gr_neg$EF, data.panel = 2, col=unname(colors[n]), ymin=0, ymax=4, r0=0.2, r1=0.85,lwd=1)
  }
}

anno='hg38_cpg_islands.bed'

fract=0.9; bin=10;slide=bin-1; en_thresh=1.2

n=length(unique(EF$Sample)); sample_thresh=fract*n

setwd('C:/Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/REZ') #nolint
#Genome
genome<-read.table('filtered_hg38-nucleus-noXY.fa.fai',sep="\t",header=F)
custom.genome <- toGRanges(data.frame(chr=genome[1], start=rep(1, nrow(genome)), end=genome[2]))
#rNMP Data processing

EF=read.table('REZ_500000_1.2_0.8_5_ef.tsv',sep="\t",header=T); REZ=read.table('REZ_500000_1.2_0.8_5_common_rezs.tsv',sep="\t",header=T)

head(EF)
#'''

EF_mean=aggregate(EF$EF, list(EF$Chromosome,EF$Start,EF$End,EF$Strand), FUN=mean); colnames(EF_mean)<-c("chr","start","stop","strand","EF_mean")
EF$count=1;EF$count[EF$EF<en_thresh]=0; EF_count=aggregate(EF$count, list(EF$Chromosome,EF$Start,EF$End,EF$Strand), FUN=sum); colnames(EF_count)<-c("chr","start","stop","strand","count")
EF_merge=merge(EF_mean,EF_count); EF_merge[order(EF_merge$strand, EF_merge$chr, EF_merge$start),]

for (chr in unique(EF_merge$chr)) {
  for (strand in unique(EF_merge$strand)) {
    print(chr); print(strand)
    df=EF_merge[EF_merge$chr == chr & EF_merge$strand == strand,]
    df=df[order(df$start),]
    for (row in seq(1,nrow(df)-9)) {
      df[seq(1,1+9),df$EF>1.2 & df$count>sample_thresh]
    }
  }
}

common_REZ

#'''
EF=aggregate(EF$EF, list(EF$Chromosome,EF$Start,EF$End,EF$Strand,EF$Celltype), FUN=mean); colnames(EF)<-c("chr","start","stop","strand","Celltype","EF")
EF.range=makeGRangesFromDataFrame(EF, keep.extra.columns=T, starts.in.df.are.0based=T)
REZ.range=makeGRangesFromDataFrame(REZ, keep.extra.columns=T, starts.in.df.are.0based=T); REZ.split=sort(split(REZ.range, strand(REZ.range)))
#Annotation Data processing

annot=read.table(anno,sep="\t",header=F); colnames(annot)[1:7]<-c("chr","start","stop","id","length","strand", "label");annot.range=makeGRangesFromDataFrame(annot[1:7], keep.extra.columns=T, starts.in.df.are.0based=T); strand(annot.range)='*'

pp <- getDefaultPlotParams(plot.type=2)
pp$ideogramheight=0
pp$data1height=500
pp$data2height=500
pp$leftmargin=0.075
pp$data1inmargin=0
pp$data2inmargin=0

png(paste(anno,'1.2_0.8_5','_22_chr_largel.png', sep=""), height = 15, width = 11,units='in',res=1000)
png(paste(anno,'1.2_0.8_5','_22_chr_large_horizontal.png', sep=""), height = 10, width = 20,units='in',res=1000)
par(mfrow=c(1,1), mar=c(0,0,0,0))
kp1 <- plotKaryotype(genome = custom.genome, plot.type=2, plot.params = pp, labels.plotter = NULL)
kpAddChromosomeNames(kp1, cex = 1.2, xoffset=-0.01,yoffset=0,srt=0)
#kp<-plotKaryotype(genome = "hg38", chromosomes=paste('chr',seq(1,22), sep=""))
plot(kp1)
dev.off()

png(paste(anno,'1.2_0.8_5','_chr19.png', sep=""), height = 1.5, width = 3,units='in',res=1000)
par(mfrow=c(1,1), mar=c(0,0,0,0))
kp2 <- plotKaryotype(genome = custom.genome, chromosomes="chr19", plot.type=2, plot.params = pp, labels.plotter = NULL)
#kpAddChromosomeNames(kp2, cex = 1, xoffset=-0.03,yoffset=0,srt=0)
#kp<-plotKaryotype(genome = "hg38", chromosomes=paste('chr',seq(1,22), sep=""))
plot(kp2)
dev.off()

png('legend.png', height = 1, width = 11,units='in',res=1000)
png('legend_horizontal.png', height = 1, width = 20,units='in',res=1000)

par(mfrow=c(1,1), mar=c(0,0,0,0))
kp3 <- plotKaryotype(genome = custom.genome, chromosomes="chr1", plot.type=2, plot.params = pp, labels.plotter = NULL)
kpAddBaseNumbers(kp3, tick.dist = 5e7, tick.len = 50, tick.col="black", cex=1,minor.tick.dist = 1e7, minor.tick.len = 25, minor.tick.col = "black", units='Mb')
kpPlotRegions(kp3, data = REZ.split$`+`[1], r0 = 0.2, r1=0.85, col=transparent('#86ECFA', amount=0.4),border=NA, data.panel=1)
kpPlotRegions(kp3, data = REZ.split$`+`[1], r0 = 0.2, r1=0.85, col=transparent('#86ECFA', amount=0.4),border=NA, data.panel=1)
kpPlotRegions(kp3,annot.range,r0 = 0, r1=0.34,col=transparent('#9E1717', amount=0.90),avoid.overlapping=FALSE,data.panel=2)
dev.off()

png('legend_chr19.png', height = 0.25, width = 3,units='in',res=1000)
par(mfrow=c(1,1), mar=c(0,0,0,0))
kp4 <- plotKaryotype(genome = custom.genome, chromosomes="chr19", plot.type=2, plot.params = pp, labels.plotter = NULL)
kpAddBaseNumbers(kp4, tick.dist = 1e7, tick.len = 90, tick.col="black", cex=0.6,units='Mb')
dev.off()
