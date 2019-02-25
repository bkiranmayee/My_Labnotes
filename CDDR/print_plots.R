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

packages("vcfR", "ggplot2", "tidyr", "plyr")


dat<-read.vcfR(file= args[1], cols = c(1:8), convertNA = T, check_keys = T, verbose = T)
#m<-metaINFO2df(dat,field="INFO")
infocols<-c("AC","AN","BQB", "DP","HOB","ICB","MQ", "MQ0F", "MQB", "MQSB" ,"RPB", "SGB", "VDB")
plot_info<-extract_info_tidy(x = test, info_fields = infocols, info_types = T)
chr<-c("AC","AN","DP")
plot_info$Key=as.numeric(plot_info$Key)
plot_info$AC=as.numeric(plot_info$AC)
plot_info$AN=as.numeric(plot_info$AN)
plot_info$DP=as.numeric(plot_info$DP)

myPath = "U:\\Stats_plots\\test.pdf"
pdf(file=myPath, onefile=T)  

for (col in 2:ncol(plot_info)){
    hist(plot_info[[col]], main = paste("Histogram of", "test"),col="blue", border="black",xlab=names(plot_info[col]))
  }

dev.off()

