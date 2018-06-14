#!/usr/bin/env Rscript
#Date:June 11 2018
#Author:Kiranmayee Bakshy

# A program to plot densities of info fields in a VCF file
# The output files have to be passed as arguments to this program
# Input = RData output of print_plots.R program

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  args[1] = "out.pdf"
} 

setwd("/mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated")

#outname=grep("(\\w+)\\.", args[1], perl=T, value=T)

load("combined.ann.bgzip.vcf.gz.RData")
combined<-plot_info

rm(plot_info)

load("filtration/filtered_no_singletons.vcf.gz.RData")
filtered<-plot_info

rm(plot_info)

print("loaded dataframes")

infocols<-c("AC","AN","BQB", "DP","HOB","ICB","MQ", "MQ0F", "MQB", "MQSB" ,"RPB", "SGB", "VDB")

print("now plotting...")

pdf(file=args[1], onefile=T)  

for (i in 2:ncol(combined)){
  plot(density(combined[[i]], border="blue",xlab=names(combined[i])))
  lines(density(filtered[[i]], border="red"))
  legend("topright", c("Combined","Filtered"), lty = c(1,1), col = c("blue","red"))
}

dev.off()

