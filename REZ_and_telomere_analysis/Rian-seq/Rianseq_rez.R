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

colors <- c(
  "CD4T"="#FF7F0E",
  "hESC"="#2CA02C",
  "HEK293T"="#D62728",
  "RNH2A-KO T3-8"="#E377C2",
  "RNH2A-KO T3-17"="#9467BD"
)

setwd('C:Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/REZ') #nolint
#anno='GSE254470_RAW/peaks/bgcall_overlaps.bed'
anno='C:/Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/REZ/RIANseq_redo/peaks/bgcall_overlaps.bed'


#Genome
genome<-read.table('C:/Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/REZ/filtered_hg38-nucleus-noXY.fa.fai',sep="\t",header=F)
custom.genome <- toGRanges(data.frame(chr=genome[1], start=rep(1, nrow(genome)), end=genome[2]))
#rNMP Data processing
EF=read.table('C:/Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/REZ/REZ_500000_1.2_0.8_5_ef.tsv',sep="\t",header=T); REZ=read.table('C:/Users/deepa/GaTech Dropbox/Deepali Kundnani/0.IMPBACKUP/data/REZ/REZ_500000_1.2_0.8_5_common_rezs.tsv',sep="\t",header=T)

EF=aggregate(EF$EF, list(EF$Chromosome,EF$Start,EF$End,EF$Strand,EF$Celltype), FUN=mean); colnames(EF)<-c("chr","start","stop","strand","Celltype","EF")
EF.range=makeGRangesFromDataFrame(EF, keep.extra.columns=T, starts.in.df.are.0based=T)
REZ.range=makeGRangesFromDataFrame(REZ, keep.extra.columns=T, starts.in.df.are.0based=T); REZ.split=sort(split(REZ.range, strand(REZ.range)))
#Annotation Data processing

annot=read.table(anno,sep="\t",header=F); colnames(annot)<-c("chr","start","stop","id","length","strand"); annot.range=makeGRangesFromDataFrame(annot, keep.extra.columns=T, starts.in.df.are.0based=T)#; strand(annot.range)='*'

pp <- getDefaultPlotParams(plot.type=2)
pp$ideogramheight=0
pp$data1height=500
pp$data2height=500
pp$leftmargin=0.075
pp$data1inmargin=0
pp$data2inmargin=0

svg(paste(anno,'1.2_0.8_5','_22_chr_large.svg', sep=""), height = 15, width = 11)
#library(svglite)
#svglite(anno,'1.2_0.8_5','_22_chr_large.svg', sep=""),  height = 15, width = 11)

#png('rianseq_22_chr_large.png', height = 15, width = 11,units='in',res=1600)
par(mfrow=c(1,1), mar=c(0,0,0,0))
kp1 <- plotKaryotype(genome = custom.genome, plot.type=2, plot.params = pp, labels.plotter = NULL)
kpAddChromosomeNames(kp1, cex = 1.2, xoffset=-0.01,yoffset=0,srt=0)
kpAxis(kp1,cex = 0.8,r0=0.2, r1=0.85, label.margin=-20e5,tick.len=10e5,ymin=0, ymax=4, numticks=2, data.panel=1)
kpAxis(kp1,cex = 0.8,r0=0.2, r1=0.85, label.margin=-20e5,tick.len=10e5,ymin=0, ymax=4, numticks=2, data.panel=2)
kpPlotRegions(kp1, data = REZ.split$`+`, r0 = 0.2, r1=0.85, col=transparent('#86ECFA', amount=0.4),border=NA, data.panel=1)
kpPlotRegions(kp1, data = REZ.split$`-`, r0 = 0.2, r1=0.85, col=transparent('#86ECFA', amount=0.4),border=NA, data.panel=2)
for (n in names(colors)){
  EF.split=EF.range[which(elementMetadata(EF.range)[,'Celltype'] == n)]
  gr_pos=sort(split(EF.split, strand(EF.split))$`+`); gr_neg=sort(split(EF.split, strand(EF.split))$`-`)
  kpLines(kp1, data=gr_pos, y=gr_pos$EF, data.panel = 1, col=unname(colors[n]), ymin=0, ymax=4, r0=0.2, r1=0.85,lwd=1)
  kpLines(kp1, data=gr_neg, y=gr_neg$EF, data.panel = 2, col=unname(colors[n]), ymin=0, ymax=4, r0=0.2, r1=0.85,lwd=1)
}
kpPlotRegions(kp1,annot.range,r0 = -0.17, r1=0.17,col=transparent('#9519d8', amount=0.95),avoid.overlapping=FALSE,data.panel=2)
dev.off()

svg(paste(anno,'1.2_0.8_5','_chr19.svg', sep=""), height = 1.5, width = 3)

#png(paste(anno,'1.2_0.8_5','_chr19.png', sep=""), height = 1.5, width = 3,units='in',res=1000)
par(mfrow=c(1,1), mar=c(0,0,0,0))
kp2 <- plotKaryotype(genome = custom.genome, chromosomes="chr19", plot.type=2, plot.params = pp, labels.plotter = NULL)
kpAxis(kp2,cex = 0.8,r0=0.2, r1=0.85, label.margin=-20e5,tick.len=10e5,ymin=0, ymax=4, numticks=2, data.panel=1)
kpAxis(kp2,cex = 0.8,r0=0.2, r1=0.85, label.margin=-20e5,tick.len=10e5,ymin=0, ymax=4, numticks=2, data.panel=2)
kpPlotRegions(kp2, data = REZ.split$`+`, r0 = 0.2, r1=0.85, col=transparent('#86ECFA', amount=0.4),border=NA, data.panel=1)
kpPlotRegions(kp2, data = REZ.split$`-`, r0 = 0.2, r1=0.85, col=transparent('#86ECFA', amount=0.4),border=NA, data.panel=2)
for (n in names(colors)){
  EF.split=EF.range[which(elementMetadata(EF.range)[,'Celltype'] == n)]
  gr_pos=sort(split(EF.split, strand(EF.split))$`+`); gr_neg=sort(split(EF.split, strand(EF.split))$`-`)
  kpLines(kp2, data=gr_pos, y=gr_pos$EF, data.panel = 1, col=unname(colors[n]), ymin=0, ymax=4, r0=0.2, r1=0.85,lwd=1)
  kpLines(kp2, data=gr_neg, y=gr_neg$EF, data.panel = 2, col=unname(colors[n]), ymin=0, ymax=4, r0=0.2, r1=0.85,lwd=1)
}
kpPlotRegions(kp2,annot.range,r0 = -0.17, r1=0.17,col=transparent('#9519d8', amount=0.95),avoid.overlapping=FALSE,data.panel=2)

dev.off()
