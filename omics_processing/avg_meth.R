#!/usr/bin/env Rscript

#Taking in arguments
library("optparse")
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="order in which user prefers to align the matrix",metavar="character"),
  make_option(c("-c", "--col"), type="integer", default=7, 
              help="order in which user prefers to align the matrix",metavar="integer"),
  make_option(c("-n", "--num"), type="integer", default=6, 
              help="number of columns to retain",metavar="integer"),
  make_option(c("-o", "--out"), type="character", default="result.bed", 
              help="output file with path [default %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser)

file=args$file
out=args$out
col=args$col
num=args$num
df=read.table(file, sep='\t', header=FALSE)
mat=aggregate(as.numeric(as.character(unlist(df[col]))),by=as.list(df[1:num]),FUN=mean)
#mat=aggregate(as.numeric(as.character(unlist(df[col]))),by=list(df$V1,df$V2,df$V3),FUN=mean)
mat=mat[order(mat[,1],mat[,2]),]

write.table(mat, out , sep = "\t",quote=F,row.names=F, col.names=F)
#print(mat , sep = "\t",quote=F,row.names=F, col.names=T)
