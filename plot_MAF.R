#!/usr/bin/env Rscript
#Date:June 8 2018
#Author:Kiranmayee Bakshy

# A program to plot densities of MAF of 2 datasets
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

filtered<-fread(args[1], header=T, stringsAsFactors = F, verbose=T, select=c("SNP", "MAF"), key = c("SNP"))
combined<-fread(args[2], header=T, stringsAsFactors = F, verbose=T, select=c("SNP", "MAF"), key = c("SNP"))

df<-merge(filtered,combined, by="SNP", all.y=T)
df_<-melt(df, id.vars = "SNP")
p<-ggplot(df_, aes(value, fill=variable)) + geom_density(alpha=0.2)+ggtitle("Minor Allele Frequency")+ylab("Density")+labs(fill="")+scale_fill_discrete(labels=c("Filtered", "Combined"))

pdf(file=args[3], onefile=T) 
p
dev.off()
