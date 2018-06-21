#!/usr/bin/env Rscript
#Date:June 8 2018
#Author:Kiranmayee Bakshy

# A program to plot densities of fraction missing per variant of 2 datasets
# The input and output files have to be passed as arguments to this program
# Input = 2 stat files and output = fileplot.pdf

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least two arguments must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "out.pdf"
}

#install.packages("vcfR")
library(ggplot2)
library(data.table)
library(reshape2)

filtered<-read.table(args[1], header=T, stringsAsFactors = F)
f<-filtered[,c(2,6)]

combined<-read.table(args[2], header=T, stringsAsFactors = F)
c<-combined[,c(2,6)]

df<-merge(f,c, by="IID", all.y=T)
df_<-melt(df, id.vars = "IID")
p<-ggplot(df_, aes(value, fill=variable)) + geom_density(alpha=0.2)+ggtitle("Fraction_missing_per_sample")+ylab("Density")+labs(fill="")+scale_fill_discrete(labels=c("Filtered", "Combined"))

pdf(file=args[3], onefile=T) 
p
dev.off()
