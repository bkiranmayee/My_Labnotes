#!/usr/bin/env Rscript
#Date:June 4 2018
#Author:Kiranmayee Bakshy

# A program to extract info fields in a VCF file and save as a dataframe
# The input file has to be passed as argument to this program
# Input = file.vcf.gz

args = commandArgs(trailingOnly=TRUE)

output_name=grep("(\\w+)\\.", args[1], perl=T, value=T)

library("vcfR", "tidyr")


dat<-read.vcfR(file= args[1], cols = c(1:8), convertNA = T, check_keys = T, verbose = T)
infocols<-c("AC","AN","BQB", "DP","HOB","ICB","MQ", "MQ0F", "MQB", "MQSB" ,"RPB", "SGB", "VDB")
plot_info<-extract_info_tidy(x = test, info_fields = infocols, info_types = T)
plot_info$Key=as.numeric(plot_info$Key)
plot_info$AC=as.numeric(plot_info$AC)
plot_info$AN=as.numeric(plot_info$AN)
plot_info$DP=as.numeric(plot_info$DP)

save(plot_info, file=paste(output_name, "RData", sep="."))