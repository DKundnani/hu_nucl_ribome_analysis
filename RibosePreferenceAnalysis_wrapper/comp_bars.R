#!/usr/bin/env Rscript

#Taking in arguments
library("optparse")
 
option_list = list(
  make_option(c("-m", "--mat"), type="character", default=NULL, 
              help="sorted matrix with labels on one site and things to be plotted on another", metavar="filetype"),
   make_option(c("-l", "--libinfo"), type="character", default=NULL, 
              help="lib \t cell \t RE", metavar="filetype"),
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
if (is.null(opt$libinfo)){
  print_help(opt_parser)
  stop("Please specify hex codes for legend with a #.n", call.=FALSE)
}
libinfo=opt$libinfo
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
comp<-read.table(mat,row.names=1,sep="\t",header=T)
comp=comp*100
lib<-read.table(libinfo,row.names=1,sep="\t",header=F)[1]
colnames(lib)="Group"
df=merge(comp,lib,by.x = 0, by.y=0)
df$Group <- factor(df$Group, levels = unique(lib$Group))
col=opt$hex #col="#D55E00,#0072B2,#F0E442,#009E73"
col=unlist(str_split(col,","))

#4.Main code
writeLines("\n... Executing Main Code...\n")
data=gather(df,key="rNTP", value = "Percentage", -c(Group,Row.names))
data[data$rNTP == 'T','rNTP'] ='U'
data$rNTP=paste('r',data$rNTP, sep="")
#5.Plotting
#writeLines("...Plotting...")

png(paste(mat,"_barplts.png",sep=""),width=length(unique(lib$Group))*2.2, height=7,units= "in",  res=600,bg="transparent")
obj<-ggbarplot(data, x = "Group", y = "Percentage", fill="rNTP",
              add = c("mean_se"),
              add.params=list(width=0.4,size=0.5),
              palette = c("#D55E00","#0072B2", "#F0E442","#009E73"),
              position = position_dodge(0.75), 
              xlab="", ylab="Normalized Percentage", 
              width = 0.75,
              size = 0.5) +
  geom_point(data,mapping=aes(Group,Percentage,color=rNTP), position=position_dodge(0.75),size=1)+guides(color = "none")+
  scale_color_manual(values=c("black","black","black","black"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,100), breaks = seq(0,100,10)) +
  theme_classic(base_size = 20, base_family = "Arial")+ 
  theme(#text=element_text(colour ="black", family = 'Arial', size=25), axis.text = element_text(colour = "black"),
  axis.text.x=element_text(color="black",size=20,angle = 0, vjust = 1, hjust=0.5),axis.text.y=element_text(color="black",size=25),
        plot.margin = unit(c(0.5, 0.0, 0.0, 0.0), "cm") #t, r, b, l 
        ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
p=ggpar(obj, legend=c("right"), legend.title = "")
print(p)
dev.off()

writeLines("\nJob Finished.\n")
writeLines(paste("Output files in same folder as the input matrix file"))