#!/usr/bin/env Rscript
#Date:Jul 05 2018
#Author:Kiranmayee Bakshy

# A program to run GWAS using rrBLUP package and plot QQ plots and Manhattan plots

library(rrBLUP)
library(corrplot)

geno<-read.delim("set1_ld_filtered_recoded_geno.tped", sep="\t", header=F, stringsAsFactors=F)
geno<-geno[,c(2,1,4,5:856)]
pheno<-read.delim("set1_ld_filtered_recoded.tfam", sep="\t", header=F, stringsAsFactors=F, row.names=1)
colnames(pheno)<-c("GID", "MID", "PID", "SEX", "STATUS")
colnames(geno)<-c("Marker","Chr","Pos",rownames(pheno))


##########################***IMPUTATION***###############################

### rrBLP program will make imputation. For the simplicity , we impute using 
# the mean but EM algorithm can be also used. 
## rrBLUP also allows to remove markers depending on the Minor allele frequency (MAF),
## in our example we remove those markers with MAF less than 0.05.

Imputation <- A.mat(geno[,c(4:855)],impute.method="EM",return.imputed=T,min.MAF=0.05)
K.mat <- Imputation$A ### KINSHIP matrix
geno.gwas <- Imputation$imputed #NEW geno data.
geno.scale <- scale(geno.gwas,center=T,scale=F) # Data needs to be center.
svdgeno <- svd(geno.scale) 
PCA <- geno.scale%*%svdgeno$v #Principal components

## Screeplot to visualize the proportion of variance explained by PCA
pdf("set1_Screeplot.pdf")
plot(round((svdgeno$d)^2/sum((svdgeno$d)^2),d=7)[1:10],type="o",main="Screeplot",xlab="PCAs",ylab="% variance")
dev.off()

##Proportion of variance explained by PCA1 and PCA2
PCA1 <- 100*round((svdgeno$d[1])^2/sum((svdgeno$d)^2),d=3); PCA1
PCA2 <- 100*round((svdgeno$d[2])^2/sum((svdgeno$d)^2),d=3); PCA2

### Plotting Principal components.
pdf("set1_PCAs.pdf")
plot(PCA[,1],PCA[,2],xlab=paste("Pcomp:",PCA1,"%",sep=""),ylab=paste("Pcomp:",PCA2,"%",sep=""),pch=20,cex=0.7)
dev.off()

### Plotting depending on clustering. 
Eucl <- dist(geno.gwas) ###Euclinean distance
fit <- hclust(Eucl,method="ward.D2") ### Ward criterion makes clusters with same size.
groups3 <- cutree(fit,k=3) ### Selecting three clusters.
table(groups3)# Number of individuals per cluster.
pdf("set1_clustering.pdf")
plot(PCA[,1],PCA[,2],xlab=paste("Pcomp:",PCA1,"%",sep=""),ylab=paste("Pcomp:",PCA2,"%",sep=""),pch=20,cex=0.7,col=groups2)
legend("bottomright",c("Subpop1: ?","Subpop2: ??"),pch=20,col=(c("black","red")),lty=0,bty="n",cex=1)
dev.off()

save.image(file="set1.RData")
