    ## Merging neogen results markers and HD bovine markers from Roslin institute ##


**Subset the main ped and map files to include only the animals that were genotyped by neogen** 

    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile datConCas/datConCasBed --cow --keep /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results/LD/neogen_cleaned_imputed.fam --make-bed --out trial2/set1

**Apply filters similar to the ones mentioned in the paper with 1 extra step: LD**

**Apply HWE filter**

**Apply LD filter**


    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile trial2/set1 --cow --allow-extra-chr --mind 0.1 --maf 0.05 --geno 0.1 --hwe 8.1e-6 --indep 10 5 4 --make-bed --out trial2/set1_clean

**Extract the SNPs that passed the LD filter**

    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile trial2/set1_clean --cow --extract trial2/set1_clean.prune.in --make-bed --out trial2/set1_ld_filtered

**Merging both the datsets**
         
    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile trial2/set1_ld_filtered --bmerge /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results/LD/neogen_cleaned_imputed --cow --out trial2/merge1
    
**Now ready to prepare the 3 pheno groups mentioned in the thesis by Raphaka from University of Edinburg**

The 3 pheno groups are manually made in excel sheet NI_samples_forJHDB recoded_formatted.xlsx
pheno1: VL+P
pheno2: VL+P+NVL
pheno3: VL+P+NVL+I+N (2.6% I and 0.1% N)



Its easy to prepare all three groups using plink then the phenotype coding can be changed in R

    
    [kiranmayee.bakshy@assembler2 trial2]$ plink --bfile merge1 --cow --keep-fam Pheno1.txt --make-bed --out pheno1
    [kiranmayee.bakshy@assembler2 trial2]$ plink --bfile merge1 --cow --keep-fam Pheno2.txt --make-bed --out pheno2
    [kiranmayee.bakshy@assembler2 trial2]$ plink --bfile merge1 --cow --keep-fam Pheno3.txt --make-bed --out pheno3


**Phenotype recoding to controls=0, cases=1**

    df<-read.delim("merge1.fam", sep="\t",header=F,stringsAsFactors=F)
    p<-read.delim("pheno_All.txt", sep="\t",header=T,stringsAsFactors=F)
    dim(p)
    colnames(df)<-c("id","V2","V3","V4","V5","pheno")
    head(df)
    merge<-merge(df,p,by="id")
    head(merge)
    dim(merge)
    merge$pheno<-ifelse(merge$TB_pheno=="VL",1,ifelse(merge$TB_pheno=="NVL",1,0))
    head(merge)
    tail(merge)
    merge$TB_pheno<-NULL
    write.table(merge,"merge1.fam",sep=" ", col.names=F, row.names=F, quote=F)

Similar steps were followed for all the 3 phenotype groups.

**Preparing Covariate files for each group using R**

First load covariates and Q covariates files of all the samples (1966):

    c<-read.delim("/mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/datConCas/datConCasCovarFile.txt",sep=" ", header=F, stringsAsFactors=F)
    
    q<-read.delim("/mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/datConCas/datConCasQcovarFile.txt",sep=" ", header=F, stringsAsFactors=F)

    colnames(c)<-c("id","V2","COV1","COV2","COV3","COV4","COV5")
    colnames(q)<-c("id","V2","QCOV")
    phe<-read.delim("pheno3.fam",sep=" ",header=F,stringsAsFactors=F)
    colnames(phe)<-c("id","V2","V3","V4","V5","pheno")
    head(phe)
    dim(phe)
    m<-merge(phe,c,by="id")
    m<-merge(m,q,by="id")
    m$COV<-rep(1,nrow(m))
    m<-m[,c("COV","COV1","COV2","COV3","COV4","COV5","QCOV")]
    head(m)
    dim(m)
    write.table(m,"pheno3_cov.txt",sep="\t",col.names=F,row.names=F,quote=F)


**Ready to run GWAS using GEMMA**

    #!/bin/bash
    
    #SBATCH --nodes=1
    #SBATCH --mem=20000
    #SBATCH --partition=assemble2
    #SBATCH --job-name=gemma
    #SBATCH -o gemma
    #SBATCH -e gemma
    #SBATCH --ntasks-per-node=1
    #SBATCH --cpus-per-task=1
    
    cd $SLURM_SUBMIT_DIR
    echo $PWD
    
    srun /mnt/nfs/nfs2/bickhart-users/binaries/gemma/bin/gemma -bfile /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/trial2/merge1 -gk 1 -o merge1
    
    srun /mnt/nfs/nfs2/bickhart-users/binaries/gemma/bin/gemma -bfile /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/trial2/pheno1 -gk 1 -o pheno1 &
    
    srun /mnt/nfs/nfs2/bickhart-users/binaries/gemma/bin/gemma -bfile /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/trial2/pheno2 -gk 1 -o pheno2 &
    
    srun /mnt/nfs/nfs2/bickhart-users/binaries/gemma/bin/gemma -bfile /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/trial2/pheno3 -n 2 -gk 1 -o pheno3 &
    
    
    
    srun /mnt/nfs/nfs2/bickhart-users/binaries/gemma/bin/gemma -bfile /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/trial2/merge1 -k output/merge1_ca0_co1.cXX.txt -c merge1_cov.txt -lmm 1 -o merge1_gwas
    
    srun /mnt/nfs/nfs2/bickhart-users/binaries/gemma/bin/gemma -bfile /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/trial2/pheno1 -k output/pheno1.cXX.txt -c pheno1_cov.txt -lmm 1 -o pheno1_gwas &
    
    srun /mnt/nfs/nfs2/bickhart-users/binaries/gemma/bin/gemma -bfile /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/trial2/pheno2 -k output/pheno2.cXX.txt -c pheno2_cov.txt -lmm 1 -o pheno2_gwas &
    
    srun /mnt/nfs/nfs2/bickhart-users/binaries/gemma/bin/gemma -bfile /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/trial2/pheno3 -n 2 -k output/pheno3.cXX.txt -c pheno3_cov.txt -lmm 1 -o pheno3_gwas &
    
    
    wait
    
  **Now plot the GWAS results**

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
    
