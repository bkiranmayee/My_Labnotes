#!/usr/bin/env Rscript
#Date:June 1 2018
#Author:Kiranmayee Bakshy

# A program to plot histograms of info fields in a VCF file
# The input and output files have to be passed as arguments to this program
# Input = file.vcf.gz and output = fileplot.pdf

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.pdf"
}

#install.packages("vcfR")
library(vcfR)

dat<-read.vcfR(file= args[1])
chrom<-create.chromR(name='chr', vcf=dat)

print("now plotting graphs")

pdf(file=args[2])  
plot(chrom)
chromoqc(chrom, dp.alpha=20)

dev.off()
