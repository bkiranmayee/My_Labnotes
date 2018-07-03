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

library(gridExtra)
library(Hmisc)
library(corrplot)

setwd("/mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/filtration")

outname=grep("(\\w+)\\.", args[2], perl=T, value=T)

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

load(args[1])
combined<-plot_info
rm(plot_info)
com<-rcorr(as.matrix(combined[,c(2:14)]),type="pearson")
com_flat<-flattenCorrMatrix(com$r, com$P)

print("processed file 1")

filtered<-readRDS(args[2])
fil<-rcorr(as.matrix(filtered[,c(2:14)]), type="pearson")
fil_flat<-flattenCorrMatrix(fil$r, fil$P)

print("now writing results")

write.table(com_flat, file="combined_cor.txt", row.names = F, quote = F)
write.table(fil_flat, file="filtered_cor.txt", row.names = F, quote = F)

print("started plotting")


pdf(file=args[3], onefile=T, paper='A4r')  

par(mfrow = c(1, 2))
corrplot(com$r, type="upper", title="Combined dataset",
            p.mat = com$P, sig.level = 0.01, insig = "blank",  mar=c(0,0,1,0))

corrplot(fil$r, type="upper", title="Filtered_dataset",
            p.mat = fil$P, sig.level = 0.01, insig = "blank",  mar=c(0,0,1,0))

dev.off()

