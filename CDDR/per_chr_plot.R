#!/usr/bin/env Rscript
#Date:June 21 2018
#Author:Kiranmayee Bakshy

# A program to plot variants per chr
# The output files have to be passed as arguments to this program
# Input = RData output is a plot.pdf

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "out.pdf"
}

library(ggplot2)
library(reshape2)

cpc<-read.table(args[1], header=F, stringsAsFactors = F)
names(cpc)<-c("Chromosome", "Combined")
fpc<-read.table(args[2], header=F, stringsAsFactors = F)
names(fpc)<-c("Chromosome", "Filtered")
df<-merge(fpc,cpc, by="Chromosome", sort=F)
df_=melt(df, id.vars=c("Chromosome")) 

pdf(file=args[3], onefile=T, paper='A4r')  

ggplot(df_, aes(Chromosome, value, fill=variable)) + geom_bar(stat='Identity', position=position_dodge())+ggtitle("No. of SNPs per Chromosome")+ylab("SNPs count")+theme(legend.position="top")+labs(fill="") +theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_x_discrete(limits=df$Chromosome)

dev.off()