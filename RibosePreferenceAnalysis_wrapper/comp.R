#!/usr/bin/env Rscript

#Taking in arguments
library("optparse")
 
option_list = list(
  make_option(c("-m", "--mat"), type="character", default=NULL, 
              help="sorted matrix with labels on one site and things to be plotted on another", metavar="filetype"),
  make_option(c("-t", "--trans"), action = "store_true", default = FALSE,
              help="columns are x axis and rows are y axis", metavar="filetype"),
  make_option(c("-c", "--hex"), type="character", default="#D55E00,#0072B2,#F0E442,#009E73", 
              help="hexcodes for colors to be used for the legend [default %default]", metavar="integer")
  #make_option(c("-y", "--y_max"), type="integer", default=300, 
              #help="max limit of y axis for all chromosomes [default %default]", metavar="integer"),
  #make_option(c("-o", "--output_prefix"), type="character", default="out", 
              #help="output file name [default %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Checking for required input
writeLines("\n...Checking input...\n")
if (is.null(opt$mat)){
  print_help(opt_parser)
  stop("Please specify input matrix.n", call.=FALSE)
}
mat=opt$mat #mat="sorted_nucl_mono_0"

transp=opt$t #transp=NULL
writeLines(paste("Input matrix name:", opt$mat))
if (is.null(opt$mat)){
  print_help(opt_parser)
  stop("Please specify hex codes for legend with a #.n", call.=FALSE)
}

#writeLines(paste("Bin width used:", opt$bin_width))
#writeLines(paste("Maximum limit of yaxis:", opt$y_max))


#####Main Script
#1.Calling required Packages
writeLines("\n...Calling required Package(s)...\n")
library(ggplot2)
#library(gg.gap)
library(scales)
library(stringr) #wrapping Labels
library(tidyr)
library(ggpubr)
if(!require(weatherData)) {install.packages("ggbreak", repos = "http://cran.us.r-project.org")}
library(ggbreak) 
library(cowplot)
library(grid)
library(gridExtra) 

#2.Defining functions if any
writeLines("\n...Defining required functions...\n")

#3. Preprocessing input 
writeLines("...Preprocessing input(s)...")
QC<-read.table(mat,row.names=1,sep="\t",header=T)
if (transp) {QC=data.frame(t(QC))}

col=opt$hex #col="#D55E00,#0072B2,#F0E442,#009E73"
col=unlist(str_split(col,","))


#QC=QC/rowSums(QC)*100
QC$Lab=rownames(QC)
#4.Main code
writeLines("\n... Executing Main Code...\n")
perc=gather(QC,"key","value",-Lab)

#5.Plotting
#writeLines("...Plotting...")


png("Horizontal_stacked.png", width=12, height=8,units= "in",  res=1000,bg="transparent" )
ggplot(perc, aes(x=Lab, y=value, fill=forcats::fct_rev(key)))+
  geom_bar(stat="identity", position="fill", width = 0.9)+
  #geom_text(aes(label=paste(round(perc$value, digits=2),"%", sep="")), position=position_stack(vjust=0.5), size=3, fontface="bold")+
  geom_text(aes(label=scales::percent(value, accuracy=0.1)), position=position_stack(vjust=0.5), size=4.5, fontface="bold")+
  #geom_text(aes(label=value),nudge_y= -.01, color="black",size = 5,fontface="bold", )+
  coord_flip(clip = "off")+
  scale_fill_manual(labels = rev(levels(factor(perc$key))),values=rev(col))+#c("#517AC9","#C05D5D"))+
  scale_x_discrete(limits=rev(perc$Lab),labels= rev(perc$Lab))+
  scale_y_continuous(expand = c(0,0),labels = scales::percent_format())+
  ggtitle(" ")+
  ylab(" ")+
  xlab("")+
  theme_classic()+
  theme(legend.background = element_rect(fill="transparent"), #element_rect(fill = "transparent",colour = NA),
        legend.key=element_rect(colour="transparent"), #legend.key = element_rect(fill = "transparent"), 
        legend.key.size = unit(0.5, "cm"),
        panel.border = element_blank(), 
        axis.text = element_text(size=16, colour = "black"),
        axis.ticks = element_line(size=0.8,color = "black"),
        axis.line = element_line(size=0.6,color = "black"),
        plot.margin = unit(c(1, 1, 0.0, 0.0), "cm")) #t, r, b, l
dev.off()

writeLines("\nJob Finished.\n")
writeLines(paste("Output files in current folder"))