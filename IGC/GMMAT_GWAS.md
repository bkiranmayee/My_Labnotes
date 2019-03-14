Labnotes on running GWAS for second round IGC markers using GMMAT
=================================================================

Since GWAS using GEMMA package gave promising results, we are now skeptical and want to re-confirm because the GEMMA package uses linear mixed model (LMM) which is not suitable for association studies with binary phenotype. This is because the distribution of the binary phenotype is discrete (not continuous). 

We can use a generalized linear mixed model to fit the binary phenotype data using a logistic function that is more suitable for binary phenotype data. This can be implemented using GMMAT R package.

The link for the GMMAT package description is [https://github.com/hanchenphd/GMMAT/blob/master/inst/doc/GMMAT.pdf](https://github.com/hanchenphd/GMMAT/blob/master/inst/doc/GMMAT.pdf "GMMAT"). 

I tried to install this package but it was not working.

This package needs several updated R dependency packages which could not be updated in the R version (R-3.3.1) on the cluster due to admin issues. So I installed a recent version of R (R-3.5.1) at /mnt/nfs/nfs2/bickhart-users/binaries 

Now I downloaded all the required packages and GMMAT is ready to run from this folder /mnt/nfs/nfs2/bickhart-users/binaries/R-3.5.1/bin/R 

I tested the run on small number of SNPs to check if it is working. It is working but is taking a lot of time!!!

    /mnt/nfs/nfs2/bickhart-users/binaries/R-3.5.1/bin/R
    
    setwd("/mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/round2/gmmat")
    
    library(GMMAT)
    
    pheno <- read.table("merge1_pheno_cov", header = TRUE, stringsAsFactors=F)
    
    GRM <- as.matrix(read.table("/mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/round2/output/merge1.cXX.txt"))

	SeqArray::seqBED2GDS("/mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/round2/merge1.bed","/mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/round2/merge1.fam","/mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/round2/merge1.bim", "merge1.gds") 
    
    infile <- "merge1.gds"
    
    soi <- c("ARS_PIRBRIGHT_5_99094858","ARS_PIRBRIGHT_5_99190989","ARS_PIRBRIGHT_5_99194167","ARS_PIRBRIGHT_5_99203402","ARS_PIRBRIGHT_5_99412564","ARS_PIRBRIGHT_5_99437819","ARS_PIRBRIGHT_5_99749976","ARS_PIRBRIGHT_18_62572950","ARS_PIRBRIGHT_18_62665698","ARS_PIRBRIGHT_18_62716825","ARS_PIRBRIGHT_18_62766196","ARS_PIRBRIGHT_18_62994372","ARS_PIRBRIGHT_18_63020431","ARS_PIRBRIGHT_18_63036451","ARS_PIRBRIGHT_18_63067935","ARS_PIRBRIGHT_18_63082203","ARS_PIRBRIGHT_18_63084493","ARS_PIRBRIGHT_18_63114542","ARS_PIRBRIGHT_18_63122963","ARS_PIRBRIGHT_18_63131111","ARS_PIRBRIGHT_18_63141688","ARS_PIRBRIGHT_18_63156196","ARS_PIRBRIGHT_18_63185710","ARS_PIRBRIGHT_18_63186702","ARS_PIRBRIGHT_18_63284810","ARS_PIRBRIGHT_18_63294645","ARS_PIRBRIGHT_18_63340289","ARS_PIRBRIGHT_18_63385249","ARS_PIRBRIGHT_18_63399311","ARS_PIRBRIGHT_18_63417698","ARS_PIRBRIGHT_23_28648192","ARS_PIRBRIGHT_LIB14427_MHC_2500","ARS_PIRBRIGHT_LIB14427_MHC_3538","ARS_PIRBRIGHT_LIB14427_MHC_36271","ARS_PIRBRIGHT_LIB14427_MHC_45690","ARS_PIRBRIGHT_LIB14427_MHC_57505","ARS_PIRBRIGHT_LIB14427_MHC_73766","ARS_PIRBRIGHT_LIB14427_MHC_116810","ARS_PIRBRIGHT_LIB14427_MHC_122678","ARS_PIRBRIGHT_LIB14427_MHC_126886","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_120784","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_181938","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_302664","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_395846","BovineHD0700028028","BovineHD1300022450","BovineHD0300013035","BovineHD1300021146")
    
    save.session("R.3.5_GMMAT_merg1.Rda")
    
    merge1_model <- glmm.wald(fixed = bTB_status ~ breed + age + yr + season + reason + prevalance, data = pheno, kins = GRM, family = binomial(link = "logit"), infile = infile, snps = soi)


GLMM is computationally intensive than LMM and I need to queue it up on the cluster.

Now drafting an Rscript:  

    #!/mnt/nfs/nfs2/bickhart-users/binaries/R-3.5.1/bin/Rscript
    # Date: Nov 30 2018
    # This script loads a svaed R session data and runs GWAS using GMMAT package on the list of snps that the user inputs
    
    
    library(GMMAT)
    
    setwd("/mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/round2/gmmat")
    
    args = commandArgs(trailingOnly=TRUE)
    
    # test if there are at least one argument: if not, return an error
    if (length(args)==0) {
      stop("At least one arguments must be supplied (snp ids file).\n", call.=FALSE)
    } 
    
    
    load("R.3.5_GMMAT_merg1.Rda")
    
    snps <- scan(args[1], what="character")
    
    results<-glmm.wald(fixed = bTB_status ~ breed + age + yr + season + reason + prevalance, data = pheno, kins = GRM, family = binomial(link = "logit"), infile = infile, snps = snps)
    
    write.table(results, file=stdout(), row.names=F, quote=F, sep = "\t")


Now split the 187,281 snpids to 188 files each containing 1000 snps so that we can pass each file to run GWAS in parallel 

    split -d snp_ids snpset
    ls snpset* > snpset.filelist

Now a bash script (runR.sh) to start the queue up of jobs...

    #!/bin/bash
    
    #SBATCH --nodes=1
    #SBATCH --mem=30000
    #SBATCH --job-name=GWAS
    #SBATCH -o output.%j
    #SBATCH -e errLog.%j
    #SBATCH --ntasks-per-node=5
    #SBATCH --cpus-per-task=1
    
    
    cd $SLURM_SUBMIT_DIR
    echo $PWD
    
    
    echo "/mnt/nfs/nfs2/bickhart-users/binaries/R-3.5.1/bin/Rscript --vanilla GWAS.R $1"
    /mnt/nfs/nfs2/bickhart-users/binaries/R-3.5.1/bin/Rscript --vanilla GWAS.R $1
    
    
Start jobs on cluster

    cat snp.filelist | xargs -I {} sbatch runR.sh {}
    
    squeue | grep 'GWAS' | head -n 80 | tail -20 | perl -lane 'print $F[0];'  | xargs -I {} scontrol update jobid={} partition=assemble3
    
    squeue | grep 'GWAS' | head -n 20 | perl -lane 'print $F[0];'  | xargs -I {} scontrol update jobid={} partition=assemble2
    
    squeue | grep 'GWAS' | head -n 40 | tail -20 | perl -lane 'print $F[0];'  | xargs -I {} scontrol update jobid={} partition=assemble1



I made a mistake by not writing the output to a file and then printing. The data seems to be printed directly to the output. 

I stopped the pending jobs and made changes to the Rscript. Now again submitted jobs to the cluster. 

But there are around 37 files with weird output format. Lets reformat using some linux commands...

    for i in output.826*; do sort -k1 -n $i | tail -n+5 | awk 'ORS=NR%2?"\t":"\n"'| awk -v OFS='\t' '{print $2,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13}' >> gwas_output; done
    
The above command pipeline sorts the output by first column, removes the first 4 lines (program commands and header lines) which we don't need, merges every second line, prints tab-delimited columns in required order...

Now waiting to complete the association tests...and then go to plotting.


	cat round3/snpfile.list3 | xargs -I {} sbatch runR.sh {}
 
	squeue -u kiranmayee.bakshy | grep 'PD' | grep 'GWAS' | tail -20 | perl -lane 'print $F[0];'  | xargs -I {} scontrol update jobid={} partition=assemble2


Here is the plotting script which is queued up and will start once the GWAS results are ready...

    #!/usr/bin/env Rscript
    #Date:Dec 7 2018
    #Author:Kiranmayee Bakshy
    
    # A program to make QQplot and manhattan plot of GWAS results obtained using GMMAT package
    library(qqman)
    library(dplyr)
    
    gwasresults<-read.delim("gwas_output_sorted", sep="\t", header=F, stringsAsFactors=F)
    
    colnames(gwasresults)<-c("rs", "chr", "pos", "ref", "alt", "n", "af", "beta", "se", "pval", "converged")
    
    #SnpsOfInterest
    #soi1 <- c("ARS_PIRBRIGHT_5_98985645","ARS_PIRBRIGHT_5_99036250","ARS_PIRBRIGHT_5_99333595","ARS_PIRBRIGHT_5_99445508","ARS_PIRBRIGHT_5_99780983","ARS_PIRBRIGHT_18_62460291","ARS_PIRBRIGHT_18_62559417","ARS_PIRBRIGHT_18_62644240","ARS_PIRBRIGHT_18_62670367","ARS_PIRBRIGHT_18_62774779","ARS_PIRBRIGHT_18_62812297","ARS_PIRBRIGHT_18_63089154","ARS_PIRBRIGHT_23_28535354","ARS_PIRBRIGHT_23_28651088","ARS_PIRBRIGHT_CH240_391K10_KIR_41480","ARS_PIRBRIGHT_CH240_370M3_LILR_LRC_61352","ARS_PIRBRIGHT_LIB14413_LRC_65014","ARS_PIRBRIGHT_CH240_370M3_LILR_LRC_72198","ARS_PIRBRIGHT_LIB14413_LRC_81028","ARS_PIRBRIGHT_CH240_370M3_LILR_LRC_98785","ARS_PIRBRIGHT_LIB14413_LRC_106729","ARS_PIRBRIGHT_CH240_370M3_LILR_LRC_174904","ARS_PIRBRIGHT_LIB14427_MHC_9213","ARS_PIRBRIGHT_LIB14427_MHC_43656","ARS_PIRBRIGHT_LIB14427_MHC_59013","ARS_PIRBRIGHT_LIB14427_MHC_86084","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_115082","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_143922","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_154399","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_208321","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_260913","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_286137","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_317666","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_324231","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_380177","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_400205")
    soi <- c("ARS_PIRBRIGHT_5_99094858","ARS_PIRBRIGHT_5_99190989","ARS_PIRBRIGHT_5_99194167","ARS_PIRBRIGHT_5_99203402","ARS_PIRBRIGHT_5_99412564","ARS_PIRBRIGHT_5_99437819","ARS_PIRBRIGHT_5_99749976","ARS_PIRBRIGHT_18_62572950","ARS_PIRBRIGHT_18_62665698","ARS_PIRBRIGHT_18_62716825","ARS_PIRBRIGHT_18_62766196","ARS_PIRBRIGHT_18_62994372","ARS_PIRBRIGHT_18_63020431","ARS_PIRBRIGHT_18_63036451","ARS_PIRBRIGHT_18_63067935","ARS_PIRBRIGHT_18_63082203","ARS_PIRBRIGHT_18_63084493","ARS_PIRBRIGHT_18_63114542","ARS_PIRBRIGHT_18_63122963","ARS_PIRBRIGHT_18_63131111","ARS_PIRBRIGHT_18_63141688","ARS_PIRBRIGHT_18_63156196","ARS_PIRBRIGHT_18_63185710","ARS_PIRBRIGHT_18_63186702","ARS_PIRBRIGHT_18_63284810","ARS_PIRBRIGHT_18_63294645","ARS_PIRBRIGHT_18_63340289","ARS_PIRBRIGHT_18_63385249","ARS_PIRBRIGHT_18_63399311","ARS_PIRBRIGHT_18_63417698","ARS_PIRBRIGHT_23_28648192","ARS_PIRBRIGHT_LIB14427_MHC_2500","ARS_PIRBRIGHT_LIB14427_MHC_3538","ARS_PIRBRIGHT_LIB14427_MHC_36271","ARS_PIRBRIGHT_LIB14427_MHC_45690","ARS_PIRBRIGHT_LIB14427_MHC_57505","ARS_PIRBRIGHT_LIB14427_MHC_73766","ARS_PIRBRIGHT_LIB14427_MHC_116810","ARS_PIRBRIGHT_LIB14427_MHC_122678","ARS_PIRBRIGHT_LIB14427_MHC_126886","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_120784","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_181938","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_302664","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_395846")
    
    
    png("Manhattan.png", w=2300, h=1200, pointsize=30)
    manhattan(gwasresults, chr="chr", bp="ps", snp="rs", p="pval", main = "Cases vs Controls", ymax = 20, col = c("grey", "skyblue"), 
       suggestiveline = -log10(5.3e-06), genomewideline = -log10(2.7e-07), highlight=soi)
    dev.off()
    
    
    png("QQplot.png", w=2300, h=1200, pointsize=30)
    qq(gwasresults$pval, main = "Cases vs Controls", xlim = c(0, 7), ylim = c(0, 21), pch = 1, col = "blue4")
    dev.off()



Moving on to CERES as the AGIL cluster is still down. 

I transfered all the data and tried to restart the GMMAT run for rest of the snps (around 17000 SNPs to go). The R version on the CERES cluster needs to be updated for the recent version R/3.5.2

The R version update by the VRSC team is done. I have installed all the required packages and GMMAT is setup.

My script seems to fail. It seems that there are some changes to the GMMAT in R/3.5.2

All the GMMAT GWAS on AGIL cluster was run on R/3.5.1

So there are some discrepancies I guess and hence I am trying to find out where the script fails.

It seems there is an addition to the default parameters that the glmm.wald function needs in order to work ("id")

I don't have the ID column set in the pheno data of my previous saved RData.

I am rebuilding all the required data. 


    setwd("/beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/IGC/gmmat")
    
    library(GMMAT)
    
	fam<-read.table("merge1.fam", sep="\t", header = F, stringsAsFactors=F)
	colnames(fam)<-c("iid","fid","mid","pid","sex","bTB_status")
	pheno <- read.table("merg1_pheno_cov", header = TRUE, stringsAsFactors=F)
	pheno$id<-fam$iid
	    
    GRM <- as.matrix(read.table("merge1.cXX.txt"))
	
	row.names(GRM)<-pheno$id
	colnames(GRM)<-pheno$id

	SeqArray::seqBED2GDS("merge1.bed","merge1.fam","merge1.bim", "merge1.gds") 
	Mon Feb 25 11:27:51 2019
	PLINK BED to SeqArray GDS Format:
    BED file: 'merge1.bed' in the SNP-major mode (Sample X SNP)
    FAM file: 'merge1.fam' (1,797 samples)
    BIM file: 'merge1.bim' (187,281 variants)
    sample.id  [md5: 0a62883ab4e1acc8795c0d98736e4f50]
    variant.id  [md5: 0be5e97ff7bdc9f4ace3998ae7165297]
    position  [md5: 5561756550f0ba98fb12c92358e3d6dc]
    chromosome  [md5: 0bb1828d07890580f867520a5a5de9b3]
    allele  [md5: c82f567225b178e84a216ea8f466b81c]
    genotype  [md5: 96e03ebae4a884aa8051071b7c82e3a5]
    phase  [md5: cee21a1da03e91d72562aa16250cdbc6]
    annotation/id  [md5: c7a60874372350ce93791337a4efb78d]
    annotation/qual  [md5: 27e8297252c73098943854bf45d0186a]
    annotation/filter  [md5: a5abcdc8c13c31aa4b353f6d42afe8c6]
    sample.annotation
	Done.
	Mon Feb 25 11:31:11 2019
	Optimize the access efficiency ...
	Clean up the fragments of GDS file:
    open the file 'merge1.gds' (50.7M)
    # of fragments: 104
    save to 'merge1.gds.tmp'
    rename 'merge1.gds.tmp' (50.7M, reduced: 672B)
    # of fragments: 48
	Mon Feb 25 11:31:11 2019

    
    infile <- "merge1.gds"
    
    soi <- c("ARS_PIRBRIGHT_5_99094858","ARS_PIRBRIGHT_5_99190989","ARS_PIRBRIGHT_5_99194167","ARS_PIRBRIGHT_5_99203402","ARS_PIRBRIGHT_5_99412564","ARS_PIRBRIGHT_5_99437819","ARS_PIRBRIGHT_5_99749976","ARS_PIRBRIGHT_18_62572950","ARS_PIRBRIGHT_18_62665698","ARS_PIRBRIGHT_18_62716825","ARS_PIRBRIGHT_18_62766196","ARS_PIRBRIGHT_18_62994372","ARS_PIRBRIGHT_18_63020431","ARS_PIRBRIGHT_18_63036451","ARS_PIRBRIGHT_18_63067935","ARS_PIRBRIGHT_18_63082203","ARS_PIRBRIGHT_18_63084493","ARS_PIRBRIGHT_18_63114542","ARS_PIRBRIGHT_18_63122963","ARS_PIRBRIGHT_18_63131111","ARS_PIRBRIGHT_18_63141688","ARS_PIRBRIGHT_18_63156196","ARS_PIRBRIGHT_18_63185710","ARS_PIRBRIGHT_18_63186702","ARS_PIRBRIGHT_18_63284810","ARS_PIRBRIGHT_18_63294645","ARS_PIRBRIGHT_18_63340289","ARS_PIRBRIGHT_18_63385249","ARS_PIRBRIGHT_18_63399311","ARS_PIRBRIGHT_18_63417698","ARS_PIRBRIGHT_23_28648192","ARS_PIRBRIGHT_LIB14427_MHC_2500","ARS_PIRBRIGHT_LIB14427_MHC_3538","ARS_PIRBRIGHT_LIB14427_MHC_36271","ARS_PIRBRIGHT_LIB14427_MHC_45690","ARS_PIRBRIGHT_LIB14427_MHC_57505","ARS_PIRBRIGHT_LIB14427_MHC_73766","ARS_PIRBRIGHT_LIB14427_MHC_116810","ARS_PIRBRIGHT_LIB14427_MHC_122678","ARS_PIRBRIGHT_LIB14427_MHC_126886","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_120784","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_181938","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_302664","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_395846","BovineHD0700028028","BovineHD1300022450","BovineHD0300013035","BovineHD1300021146")
    
    save.image("GMMAT_merg1.Rdata")

Finally this seems to work without errors. I ran the command on a couple of SNPs to see if I am not getting any errors.

Now queuing up the GWAS jobs on the ceres cluster using the shell and R scripts (GWAS.R, runR.sh) modified according to the ceres cluster parameters...

	cat round3/snpfile.list3 | xargs -I {} sbatch runR.sh {}


When GWAS for all the SNPs complete, I need to concatenate all the results, sort and plot the Manhattan and QQ plots.

There were some output files which needs to be formatted properly, I did it using an Rscript ()

After concatenation and sorting (according to chr and position) of all the GMMAT output files, I plotted both the Manhattan and QQ plots.

Here are the plots:


![Manhattan plot](https://i.imgur.com/iBblV4W.png)

![QQplot](https://i.imgur.com/jXN6eFy.png)


**Summary of the significant SNPs:**

| ID                        | CHR | POS      | REF | ALT | N    | Alt AF      | BETA         | SE          | pval     | converged | MAF/Irish Holstein      | MAF/US Holstein (172)|
|---------------------------|-----|----------|-----|-----|------|-------------|--------------|-------------|----------|-----------|----------|---------|
| BovineHD0700028028        | 7   | 96196083 | C   | A   | 1765 | 0.837677054 | 0.571489432  | 0.113696428 | 5.00E-07 | TRUE      | 0.16 |         |
| ARS_PIRBRIGHT_18_62766196 | 18  | 63021818 | G   | A   | 1797 | 0.517529215 | -2.074083012 | 0.276753219 | 6.66E-14 | TRUE      | 0.482    | 0.66    |
| ARS_PIRBRIGHT_18_63141688 | 18  | 63397310 | C   | G   | 1797 | 0.681691708 | -0.734112587 | 0.109094639 | 1.71E-11 | TRUE      | 0.32 | 0.31    |

**The significant SNP on BTA 7 was not previously reported either by the Roslin GWAS study or in the thesis by Raphaka et al.**

LD between the 2 significant SNPs on BTA 18 has been calculated: **R2 = 0.32**

### Principle Components Analysis: ###

PCA was carried out using plink:

PLINK v1.90b4.4 64-bit (21 May 2017)
Options in effect:
  --1
  --bfile /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/round2/merge1
  --covar /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/trial2/merge1_cov.txt
  --cow
  --out pc_cov
  --pca header var-wts

Hostname: assembler3.agil.barc.ba.ars.usda.gov
Working directory: /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/round2/gmmat/pca
Start time: Mon Mar  4 16:28:08 2019


End time: Mon Mar  4 16:28:08 2019


![pc12_plot](https://i.imgur.com/QMiNNR7.png)

![pc13_plot](https://i.imgur.com/pgYGEt0.png)

![pc23_plot](https://i.imgur.com/bdk95qd.png)


I have calculated the inflation factor (lambda) using p-values to be **1.012** 
 
	# Here is my R code for calculating the lambda: 
    chisq<-qchisq(1-df$pval, 1)
    lambda<-median(chisq)/qchisq(0.5,1)
 

The genomic inflation factor is close to 1 which indicates there is no population stratification that can effect the GWAS.

It seems like the population substructure has been taken care of by the GRM in the GWAS.


Significant SNPs identified on BTA 23 as reported in the thesis by Raphaka et al

 SNP1 = ARS-BFGL-NGS-40833; SNP2= Hapmap38114-BTA-57971; SNP3 = BTA-56563-no-rs













    
    
 