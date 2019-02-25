#!/usr/bin/env Rscript
#Date:Jul 13 2018
#Author:Kiranmayee Bakshy

# A program to make QQplot and manhattan plot of GWAS results obtained using GEMMA package
library(qqman)
library(dplyr)

gwasresults<-read.delim("output/merge1_gwas.assoc.txt", sep="\t", header=T, stringsAsFactors=F)
gwasresults1<-read.delim("output/pheno1_gwas.assoc.txt", sep="\t", header=T, stringsAsFactors=F)
gwasresults2<-read.delim("output/pheno2_gwas.assoc.txt", sep="\t", header=T, stringsAsFactors=F)
gwasresults3<-read.delim("output/pheno3_gwas.assoc.txt", sep="\t", header=T, stringsAsFactors=F)

#SnpsOfInterest
soi <- c("ARS_PIRBRIGHT_5_98985645","ARS_PIRBRIGHT_5_99036250","ARS_PIRBRIGHT_5_99333595","ARS_PIRBRIGHT_5_99445508","ARS_PIRBRIGHT_5_99780983","ARS_PIRBRIGHT_18_62460291","ARS_PIRBRIGHT_18_62559417","ARS_PIRBRIGHT_18_62644240","ARS_PIRBRIGHT_18_62670367","ARS_PIRBRIGHT_18_62774779","ARS_PIRBRIGHT_18_62812297","ARS_PIRBRIGHT_18_63089154","ARS_PIRBRIGHT_23_28535354","ARS_PIRBRIGHT_23_28651088","ARS_PIRBRIGHT_CH240_391K10_KIR_41480","ARS_PIRBRIGHT_CH240_370M3_LILR_LRC_61352","ARS_PIRBRIGHT_LIB14413_LRC_65014","ARS_PIRBRIGHT_CH240_370M3_LILR_LRC_72198","ARS_PIRBRIGHT_LIB14413_LRC_81028","ARS_PIRBRIGHT_CH240_370M3_LILR_LRC_98785","ARS_PIRBRIGHT_LIB14413_LRC_106729","ARS_PIRBRIGHT_CH240_370M3_LILR_LRC_174904","ARS_PIRBRIGHT_LIB14427_MHC_9213","ARS_PIRBRIGHT_LIB14427_MHC_43656","ARS_PIRBRIGHT_LIB14427_MHC_59013","ARS_PIRBRIGHT_LIB14427_MHC_86084","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_115082","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_143922","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_154399","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_208321","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_260913","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_286137","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_317666","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_324231","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_380177","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_400205")
   
pdf("merge1_manhattan.pdf", paper='A4r')
manhattan(gwasresults, chr="chr", bp="ps", snp="rs", p="p_wald", highlight=soi)
dev.off()

pdf("merge1_QQplot.pdf", paper="A4r")
qq(gwasresults$p_wald)
dev.off()

pdf("pheno1_manhattan.pdf", paper='A4r')
manhattan(gwasresults1, chr="chr", bp="ps", snp="rs", p="p_wald", highlight=soi)
dev.off()

pdf("pheno1_QQplot.pdf", paper="A4r")
qq(gwasresults1$p_wald)
dev.off()

pdf("pheno2_manhattan.pdf", paper='A4r')
manhattan(gwasresults2, chr="chr", bp="ps", snp="rs", p="p_wald", highlight=soi)
dev.off()

pdf("pheno2_QQplot.pdf", paper="A4r")
qq(gwasresults2$p_wald)
dev.off()

pdf("pheno3_manhattan.pdf", paper='A4r')
manhattan(gwasresults3, chr="chr", bp="ps", snp="rs", p="p_wald", highlight=soi)
dev.off()

pdf("pheno3_QQplot.pdf", paper="A4r")
qq(gwasresults3$p_wald)
dev.off()

