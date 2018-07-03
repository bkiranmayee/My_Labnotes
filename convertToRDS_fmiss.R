#!/usr/bin/env Rscript
#Date:June 1 2018
#Author:Kiranmayee Bakshy

# A program to extract INFO feilds of a VCF file and save it as a dataframe for plotting
# The input and output files have to be passed as arguments to this program
# Input = file and output = file.rds

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} 

#setwd("/mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/filtration/filtered_run1")

wd = "/mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/filtration/filtered_run1/"


library(data.table)

for (i in length(args)) {
	df[i]<-fread(args[i], sep="auto", header="auto", select=c(2, ncols), data.table=F)
	outname=grep("(\\w+)\\.", args[i], perl=T, value=T)
	file=paste(wd, outname)
	saveRDS(df[i], file=paste(file, "rds", sep="."))
	}


	

