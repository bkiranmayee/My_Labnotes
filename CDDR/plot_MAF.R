#!/usr/bin/env Rscript
#Date:June 8 2018
#Author:Kiranmayee Bakshy

# A program to plot densities of MAF/fraction missing per variant of 2 datasets (designed for plink output files like file.frq and file.lmiss)
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

library(ggplot2)
library(reshape2)

combined<-readRDS(args[1])
filtered<-readRDS(args[2])

pdf(file=args[3], onefile=T, paper='A4r') 

plot(density(combined[[2]], na.rm=T), col="blue", border="black", main="Density plot", xlab=names(filtered[2]))
  lines(density(filtered[[2]], na.rm=T), col="red", border="black")
  legend("top", c("Combined","Filtered"), lty = c(1,1), col = c("blue","red"))

dev.off()

