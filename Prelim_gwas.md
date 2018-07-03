## Preliminary GWAS trials ##

### Preparing the data ###

The data we got from the Roslin Institute is at /mnt/nfs/nfs2/bickhart-users/natdb\_sequencing/prelim_gwas

First do the LD pruning of the case-control dataset (datConCasBed), I used  plink/1.90b4.4-2017-05-21 for this:

    [kiranmayee.bakshy@assembler2 datConCas]$ plink --bfile datConCasBed --cow --indep 10 5 4
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to plink.log.
    Options in effect:
      --bfile datConCasBed
      --cow
      --indep 10 5 4
    
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    538231 variants loaded from .bim file.
    1966 cattle (0 males, 1966 females) loaded from .fam.
    1966 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 1966 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.995131.
    538231 variants and 1966 cattle pass filters and QC.
    Among remaining phenotypes, 1407 are cases and 559 are controls.
    --indep: Ignoring 645 chromosome 0 variants.
    Pruned 22452 variants from chromosome 1, leaving 11555.
    Pruned 19216 variants from chromosome 2, leaving 9531.
    Pruned 17137 variants from chromosome 3, leaving 8780.
    Pruned 16801 variants from chromosome 4, leaving 8998.
    Pruned 16534 variants from chromosome 5, leaving 8379.
    Pruned 17823 variants from chromosome 6, leaving 9085.
    Pruned 15969 variants from chromosome 7, leaving 7811.
    Pruned 13030 variants from chromosome 8, leaving 7243.
    Pruned 15166 variants from chromosome 9, leaving 7760.
    Pruned 14974 variants from chromosome 10, leaving 8112.
    Pruned 16227 variants from chromosome 11, leaving 8069.
    Pruned 12470 variants from chromosome 12, leaving 6755.
    Pruned 9566 variants from chromosome 13, leaving 5355.
    Pruned 9887 variants from chromosome 14, leaving 5451.
    Pruned 11779 variants from chromosome 15, leaving 6686.
    Pruned 11857 variants from chromosome 16, leaving 5835.
    Pruned 10937 variants from chromosome 17, leaving 5913.
    Pruned 9622 variants from chromosome 18, leaving 5351.
    Pruned 9055 variants from chromosome 19, leaving 5338.
    Pruned 10985 variants from chromosome 20, leaving 5645.
    Pruned 10038 variants from chromosome 21, leaving 5305.
    Pruned 9246 variants from chromosome 22, leaving 4956.
    Pruned 7025 variants from chromosome 23, leaving 4604.
    Pruned 9011 variants from chromosome 24, leaving 4727.
    Pruned 6189 variants from chromosome 25, leaving 3815.
    Pruned 7571 variants from chromosome 26, leaving 4188.
    Pruned 6391 variants from chromosome 27, leaving 3984.
    Pruned 6431 variants from chromosome 28, leaving 3905.
    Pruned 6907 variants from chromosome 29, leaving 4154.
    Pruning complete.  350296 of 537586 variants removed.
    Marker lists written to plink.prune.in and plink.prune.out .
 

The above command gives 2 files, prune.in and prune.out which can be used to subset the main data

    [kiranmayee.bakshy@assembler2 datConCas]$ plink --bfile datConCasBed --extract plink.prune.in --out datConCaspr
    PLINK v1.90b4.4 64-bit (21 May 2017)
    Options in effect:
      --bfile datConCasBed
      --cow
      --extract plink.prune.in
      --make-bed
      --out datConCaspr
    
    Hostname: assembler2.agil.barc.ba.ars.usda.gov
    Working directory: /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/datConCas
    Start time: Thu Jun 28 13:42:49 2018
    
    Random number seed: 1530207769
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    538231 variants loaded from .bim file.
    1966 cattle (0 males, 1966 females) loaded from .fam.
    1966 phenotype values loaded from .fam.
    --extract: 187290 variants remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 1966 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.995255.
    187290 variants and 1966 cattle pass filters and QC.
    Among remaining phenotypes, 1407 are cases and 559 are controls.
    --make-bed to datConCaspr.bed + datConCaspr.bim + datConCaspr.fam ... done.
    
    End time: Thu Jun 28 13:42:51 2018

    
Now I used this pruned dataset to test for association using Stratified analyses explained [here](http://zzz.bwh.harvard.edu/plink/anal.shtml#strat).

> **This is an attempt to test the hint explained in this section. Not sure if it works for our analyses** 

**A trick to analyse phenotypes with more than two categories (but only with nominal, not ordinal outcomes) is to use the --mh2 option with the phenotype in the cluster file and the phenotype in the PED file set all to a single value.**

I made a cluster file manually with only the cases coded as Nvls and Vls. The rest of the files (bed, bim and fam) being same. Phenotypes coded as 1:unaffected; 2: affected. 
 
Note: This test worked in plink/1.07-2009-08-10 but failed in plink/1.90b4.4-2017-05-21
    
    [kiranmayee.bakshy@assembler2 datConCas]$ plink --bfile datConCaspr --cow --mh2 --within datConCasBed_cluster.txt --covar datConCasCovarFile.txt
    
    @----------------------------------------------------------@
    |PLINK!   | v1.07  |   10/Aug/2009 |
    |----------------------------------------------------------|
    |  (C) 2009 Shaun Purcell, GNU General Public License, v2  |
    |----------------------------------------------------------|
    |  For documentation, citation & bug-report instructions:  |
    |http://pngu.mgh.harvard.edu/purcell/plink/|
    @----------------------------------------------------------@
    
    Web-based version check ( --noweb to skip )
    Recent cached web-check found...Problem connecting to web
    
    Writing this text to log file [ plink.log ]
    Analysis started: Fri Jun 29 11:18:31 2018
    
    Options in effect:
    	--bfile datConCaspr
    	--cow
    	--mh2
    	--within datConCasBed_cluster.txt
    	--covar datConCasCovarFile.txt
    
    Reading map (extended format) from [ datConCaspr.bim ] 
    187290 markers to be included from [ datConCaspr.bim ]
    Reading pedigree information from [ datConCaspr.fam ] 
    1966 individuals read from [ datConCaspr.fam ] 
    1966 individuals with nonmissing phenotypes
    Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
    Missing phenotype value is also -9
    1407 cases, 559 controls and 0 missing
    0 males, 1966 females, and 0 of unspecified sex
    Reading genotype bitfile from [ datConCaspr.bed ] 
    Detected that binary PED file is v1.00 SNP-major mode
    Reading 5 covariates from [ datConCasCovarFile.txt ] with nonmissing values for 1966 individuals
    Reading clusters from [ datConCasBed_cluster.txt ]
    1407 of 1966 individuals assigned to 2 cluster(s)
    Before frequency and genotyping pruning, there are 187290 SNPs
    1966 founders and 0 non-founders found
    Total genotyping rate in remaining individuals is 0.995255
    0 SNPs failed missingness test ( GENO > 1 )
    0 SNPs failed frequency test ( MAF < 0 )
    After frequency and genotyping pruning, there are 187290 SNPs
    After filtering, 1407 cases, 559 controls and 0 missing
    After filtering, 0 males, 1966 females, and 0 of unspecified sex
    Cochran-Mantel-Haenszel IxJxK test, K = 2
    Testing SNP x STRATUM | DISEASE (option --mh2)
    Writing results to [ plink.cmh2 ]
    
    Analysis finished: Fri Jun 29 11:19:29 2018


OK, now lets plot the p-values and see if something is really significant. 


