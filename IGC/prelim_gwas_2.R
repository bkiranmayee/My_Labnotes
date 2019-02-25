#!/usr/bin/env Rscript
#Date:Jul 05 2018
#Author:Kiranmayee Bakshy

# A program to run GWAS using rrBLUP package and plot QQ plots and Manhattan plots

library(rrBLUP)
library(corrplot)

load("set1.RData")

gwasresults1<-GWAS(pheno,geno.gwas,fixed=colnames(pheno)[2:4], K=NULL, n.PC=0, n.core>1, plot=TRUE)
gwasresults2<-GWAS(pheno,geno.gwas,fixed=colnames(pheno)[2:4], K=NULL, n.PC=2, n.core>1, plot=TRUE)
gwasresults3<-GWAS(pheno,geno.gwas,fixed=colnames(pheno)[2:4], K=K.mat, n.PC=0, n.core>1, plot=TRUE)
gwasresults4<-GWAS(pheno,geno.gwas,fixed=colnames(pheno)[2:4], K=K.mat, n.PC=2, n.core>1, plot=TRUE)


###################################################################################
#################################*** QQ PLOT*****##################################
###################################################################################

pdf("set1_QQplot.pdf",width = 7)
par(mfrow=c(2,2))
N <- length(gwasresults1$STATUS)
expected.logvalues <- sort( -log10( c(1:N) * (1/N) ) )
observed.logvalues <- sort( gwasresults1$STATUS)

plot(expected.logvalues , observed.logvalues, main="Naive model(K=NULL,n.PC=0)", 
     xlab="expected -log pvalue ", 
     ylab="observed -log p-values",col.main="blue",col="coral1",pch=20)
abline(0,1,lwd=3,col="black")


N1 <- length(gwasresults2$STATUS)
expected.logvalues1 <- sort( -log10( c(1:N1) * (1/N1) ) )
observed.logvalues1 <- sort( gwasresults2$STATUS)

plot(expected.logvalues1 , observed.logvalues1, main="Q model (K=NULL,n.PC=2)", 
     xlab="expected -log pvalue ", 
     ylab="observed -log p-values",col.main="blue",col="coral1",pch=20)
abline(0,1,lwd=2,col="black")


N2 <- length(gwasresults3$STATUS)
expected.logvalues2 <- sort( -log10( c(1:N2) * (1/N2) ) )
observed.logvalues2 <- sort( gwasresults3$STATUS)

plot(expected.logvalues2 , observed.logvalues2, main="K model (n.PC=0)", 
     xlab="expected -log pvalue ", 
     ylab="observed -log p-values",col.main="blue",col="coral1",pch=20)
abline(0,1,lwd=2,col="black")

N3 <- length(gwasresults4$STATUS)
expected.logvalues3 <- sort( -log10( c(1:N3) * (1/N3) ) )
observed.logvalues3 <- sort( gwasresults4$STATUS)

plot(expected.logvalues3 , observed.logvalues3, main="Q+K model (n.PC=2)", 
     xlab="expected -log pvalue ", 
     ylab="observed -log p-values",col.main="blue",col="coral1",pch=20)
abline(0,1,lwd=2,col="black")

dev.off()

###################################################################################
#################################*** MANHATTAN PLOT*****###########################
###################################################################################
#False Discovery Rate Function

FDR<-function(pvals, FDR){
  pvalss<-sort(pvals, decreasing=F)
  m=length(pvalss)
  cutoffs<-((1:m)/m)*FDR
  logicvec<-pvalss<=cutoffs
  postrue<-which(logicvec)
  print(postrue)
  k<-max(c(postrue,0))
  cutoff<-(((0:m)/m)*FDR)[k+1]
  return(cutoff)
}

alpha_bonferroni=-log10(0.05/length(gwasresults1$STATUS)) ###This is Bonferroni correcton
alpha_FDR_STATUS <- -log10(FDR(10^(-gwasresults1$STATUS),0.05))## This is FDR cut off



#################################*** MANHATTAN PLOT*****###########################

pdf("set1_Manhattan.pdf",width=8,height=8)
par(mfrow=c(2,2))
plot(gwasresults1$STATUS,col=gwasresults1$Chr,ylab="-log10.pvalue",
     main="Naïve model (K=NULL,n.PC=0)",xaxt="n",xlab="Position",ylim=c(0,14))
axis(1,at=c(1:length(unique(gwasresults1$Chr))),labels=unique(gwasresults1$Chr))
#axis(1,at=c(0,440,880,1320,1760))
abline(a=NULL,b=NULL,h=alpha_bonferroni,col="blue",lwd=2)
abline(a=NULL,b=NULL,h=alpha_FDR_Yield,col="red",lwd=2,lty=2)
legend(1,13.5, c("Bonferroni","FDR") , 
       lty=1, col=c('red', 'blue'), bty='n', cex=1,lwd=2)

plot(gwasresults2$STATUS,col=gwasresults2$Chr,ylim=c(0,14),ylab="-log10.pvalue",
     main="Q model (K=NULL,n.PC=2)",xaxt="n",xlab="Position")
axis(1,at=c(0,440,880,1320,1760))
abline(a=NULL,b=NULL,h=alpha_bonferroni,col="blue",lwd=2)
abline(a=NULL,b=NULL,h=alpha_FDR_Yield,col="red",lwd=2,lty=2)
legend(1.5,13.5, c("Bonferroni","FDR") , 
       lty=1, col=c('red', 'blue'), bty='n', cex=1,lwd=2)

plot(gwasresults3$STATUS,col=gwasresults3$Chr,ylim=c(0,14),ylab="-log10.pvalue",
     main="K model (K=K.mat,n.PC=0)",xaxt="n",xlab="Position")
axis(1,at=c(0,440,880,1320,1760))
abline(a=NULL,b=NULL,h=alpha_bonferroni,col="blue",lwd=2)
abline(a=NULL,b=NULL,h=alpha_FDR_Yield,col="red",lwd=2,lty=2)
legend(1.5,13.5, c("Bonferroni","FDR") , 
       lty=1, col=c('red', 'blue'), bty='n', cex=1,lwd=2)

plot(gwasresults4$STATUS,col=gwasresults4$chr,ylim=c(0,14),ylab="-log10.pvalue",
     main="Q+K model (K=K.mat,n.PC=2)",xaxt="n",xlab="Position")
axis(1,at=c(0,440,880,1320,1760))
abline(a=NULL,b=NULL,h=alpha_bonferroni,col="blue",lwd=2)
abline(a=NULL,b=NULL,h=alpha_FDR_Yield,col="red",lwd=2,lty=2)##FDR gives inf for Yield
legend(1,13.5, c("Bonferroni","FDR") , 
       lty=1, col=c('red', 'blue'), bty='n', cex=1,lwd=2)
dev.off()

