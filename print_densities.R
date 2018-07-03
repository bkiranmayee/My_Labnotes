#!/usr/bin/env Rscript
#Date:June 11 2018
#Author:Kiranmayee Bakshy

# A program to plot densities of info fields in a VCF file
# The output files have to be passed as arguments to this program
# Input = RData output of print_plots.R program

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "out.pdf"
}

setwd("/mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/filtration")

outname=grep("(\\w+)\\.", args[2], perl=T, value=T)

load(args[1])
combined<-plot_info

rm(plot_info)

filtered<-readRDS(args[2])

print("loaded dataframes")

infocols<-c("AC","AN","BQB", "DP","HOB","ICB","MQ", "MQ0F", "MQB", "MQSB" ,"RPB", "SGB", "VDB")

print("now plotting...")

pdf(file=args[3], onefile=T)  

for (i in 2:ncol(combined)){
  plot(density(combined[[i]], na.rm=T), col="blue", border="black", main=paste("Density plot", outname, sep=' '), xlab=names(combined[i]))
  lines(density(filtered[[i]], na.rm=T), col="red", border="black")
  legend("top", c("Combined","Filtered"), lty = c(1,1), col = c("blue","red"))
}

dev.off()

