#!/usr/bin/env Rscript
#Date:June 13 2018
#Author:Kiranmayee Bakshy

# A program to check normality certain info fields in a VCF file
# The input files have to be passed as arguments to this program
# Input = RData output is a plot.pdf

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.pdf"
}

library(ggpubr)
library(gridExtra)

setwd("/mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/filtration")

outname=grep("(\\w+)\\.", args[1], perl=T, value=T)

load(args[1])

print("loaded dataframes")

cols<-c("MQB", "MQSB" ,"RPB")

print("now plotting correlations...")


pdf(file=args[2], onefile=T)  

bqb<-ggqqplot(plot_info$BQB, ylab = "BQB")
mqb<-ggqqplot(plot_info$MQB, ylab = "MQB")
mqsb<-ggqqplot(plot_info$MQSB, ylab = "MQSB")
rpb<-ggqqplot(plot_info$RPB, ylab = "RPB")
  
grid.arrange(arrangeGrob(bqb,mqb,mqsb,rpb, ncol = 2), nrow = 2)

dev.off()

