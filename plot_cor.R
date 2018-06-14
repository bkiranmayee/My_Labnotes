#!/usr/bin/env Rscript
#Date:June 13 2018
#Author:Kiranmayee Bakshy

# A program to calculate correlation between certain info fields in a VCF file
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

library(ggpubr)
library(gridExtra)

setwd("/mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/filtration")

outname=grep("(\\w+)\\.", args[2], perl=T, value=T)

load(args[1])
combined<-plot_info

rm(plot_info)

load(args[2])
filtered<-plot_info

rm(plot_info)


print("loaded dataframes")

cols<-c("MQB", "MQSB" ,"RPB")

print("now plotting correlations...")


pdf(file=args[3], onefile=T)  

for (i in 1:length(cols)){
  c[i]<-ggscatter(combined, x = "BQB", y = cols[i],
               add = "reg.line", conf.int = TRUE,
               cor.coef = TRUE, cor.method = "pearson", main="Combined",
               xlab = "BQB", ylab = cols[i])
  
  f[i]<-ggscatter(filtered, x = "BQB", y = cols[i],
               add = "reg.line", conf.int = TRUE,
               cor.coef = TRUE, cor.method = "pearson", main="Filtered",
               xlab = "BQB", ylab = cols[i])
  
  grid.arrange(arrangeGrob(c[i], f[i], ncol = 2), nrow =1)
}

dev.off()

