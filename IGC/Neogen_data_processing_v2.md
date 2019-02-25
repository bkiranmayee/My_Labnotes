## Prepare Neogen genotyping results file for analyses##

**1a) Convert the neogen genotype file to plink format** 

This was done by a custom script: neogen\_to\_ped.pl

    #!/usr/bin/perl 
    
    ## A script to convert neogen results to ped format
    ## Date: June 29 2018
    
    use strict;
    
    my $input = "neogen_results_5_18_2018.tab";
    open (my $OUTFILE, "> neogen_updated.ped");
    
    my @newcols;
    my @new;
    
    my %varRefAlleles = (
    "18_62404864" => "G",
    "18_62438492" => "T",
    "18_62460291" => "C",
    "18_62469715" => "C",
    "18_62559417" => "A",
    "18_62575095" => "A",
    "18_62644240" => "T",
    "18_62670367" => "G",
    "18_62723255" => "C",
    "18_62774779" => "C",
    "18_62812297" => "A",
    "18_62824331" => "G",
    "18_62914384" => "C",
    "18_62985809" => "C",
    "18_63089154" => "T",
    "18_63138493" => "C",
    "18_63222615" => "C",
    "18_63269795" => "C",
    "18_63308589" => "C",
    "18_63356424" => "C",
    "23_28481840" => "T",
    "23_28518797" => "A",
    "23_28535354" => "C",
    "23_28651088" => "A",
    "5_98985645" => "A",
    "5_99036250" => "C",
    "5_99073686" => "T",
    "5_99228484" => "C",
    "5_99266853" => "T",
    "5_99333595" => "A",
    "5_99392730" => "C",
    "5_99445508" => "A",
    "5_99560651" => "C",
    "5_99706728" => "A",
    "5_99763326" => "G",
    "5_99766613" => "C",
    "5_99780983" => "C",
    "KIR_41480" => "T",
    "LRC_106729" => "G",
    "LRC_118549" => "T",
    "LRC_174904" => "G",
    "LRC_61352" => "A",
    "LRC_65014" => "G",
    "LRC_72198" => "G",
    "LRC_81028" => "T",
    "LRC_82562" => "G",
    "LRC_98357" => "G",
    "LRC_98785" => "A",
    "MHC_103858" => "C",
    "MHC_115082" => "G",
    "MHC_121069" => "G",
    "MHC_143922" => "A",
    "MHC_154399" => "A",
    "MHC_208321" => "T",
    "MHC_260913" => "A",
    "MHC_286137" => "C",
    "MHC_317666" => "A",
    "MHC_324231" => "T",
    "MHC_356873" => "T",
    "MHC_359718" => "A",
    "MHC_380177" => "T",
    "MHC_400205" => "A",
    "MHC_43656" => "G",
    "MHC_59013" => "A",
    "MHC_63411" => "A",
    "MHC_86084" => "A",
    "MHC_9213" => "A");
    
    open(my $IN, "< $input") || die "Could not open input file: $input!\n";
    
    my %varIdx;
    while(my $line = <$IN>){
    	chomp $line;
    	my @hsegs = split(/\t/, $line);
    	
    	for(my $x = 2; $x < scalar(@hsegs); $x++){
    		$varIdx{$x} = $hsegs[$x];
    		}
    	print $OUTFILE "$line\n";
    	last;
    }
    
    my @newsegs;
    while(my $line = <$IN>){
    	chomp $line;
    	my @segs = split(/\t/, $line);
    	
    	for(my $x = 2; $x < scalar(@segs); $x++){
    			my $var = $varIdx{$x};
    			my $hallele = $varRefAlleles{$var};
    			if(!defined($segs[$x]) || $segs[$x] eq ""){
    				$newsegs[$x] = join (" ", "0","0");
    			}elsif(length($segs[$x]) > 1){
    				my @newgt = split("", $segs[$x]);
    				$newsegs[$x] = join " ", @newgt;
    			}elsif($segs[$x] eq $hallele){
    				$newsegs[$x] = join (' ', $hallele, $hallele);
    			}else{
    				$newsegs[$x] = join (' ', $segs[$x], $segs[$x]);
    				}
    			
    		}
    		print {$OUTFILE} join("\t", $segs[0], @newsegs) . "\n";
    }
    			
    
    close $IN;
    close $OUTFILE;
    exit;


The output was then manually edited to include the extra columns such as Family ID, maternal ID, paternal ID, Sex and phenotype from NI\_samples\_forJHDB recoded_formatted.xlsx file. 

 A map file was then generated


**1b) Converting the plink file format to vcf (since BEAGLE imputation requires vcf as input) keeping in mind the Ref allele:**

     [kiranmayee.bakshy@assembler2 LD]$ plink --file neogen_updated --allow-extra-chr --recode vcf-iid --a2-allele Ref_allele.txt --out neogen_cleaned


**1c) Imputation of Genotypes using BEAGLE**

Input and output files for Beagle is VCF format

    [kiranmayee.bakshy@assembler2 LD]$ srun java -Xmx50g -D$JAVA_HOME -jar /mnt/nfs/nfs2/bickhart-users/binaries/beagle.03Jul18.40b.jar gt=neogen_cleaned.vcf out=neogen-imputed

The imputation seems to be successful for all the variants except for KIR_41480 because there is only one marker that belongs to that Chromosome (KIR).
Mmmm! I think I should first make arbitary UMD3 coordinates for these and then do the imputation I guess!!

OK redo the whole thing again.

I converted the map and ped files to binary plink format by giving the Ref allele text file so that it does not change the alleles based on frequency. I have also cleaned the dataset in the same step:

     [kiranmayee.bakshy@assembler2 LD]$ plink --file neogen_updated --geno 0.2 --maf 0.05 --cow --allow-extra-chr --make-bed --recode --a2-allele Ref_allele.txt --out neogen_cleaned

I cleaned the data based on call rate (reject > 80% missing per variant (--geno 0.2) and MAF (--maf 0.05). I get only 36 variants left after the cleaning step. (saved the dataset as neogen_cleaned)

Then I gave arbitrary UMD3 coordinates and updated the SNP IDs manually in an excel sheet (UMD3_arbitrary_coordinates.xlsx).

I saved the original SNP coordinates map file as neogen_cleaned_original.map. 

I then converted it to VCF format for imputation using Beagle

    [kiranmayee.bakshy@assembler2 LD]$ plink --file neogen_cleaned --geno 0.2 --maf 0.05 --cow --allow-extra-chr --recode vcf-iid --a2-allele Ref_allele.txt --out neogen_cleaned
    [kiranmayee.bakshy@assembler2 LD]$ java -Xmx50g -D$JAVA_HOME -jar /mnt/nfs/nfs2/bickhart-users/binaries/beagle.03Jul18.40b.jar gt=neogen_cleaned.vcf out=neogen-imputed

Now the imputation was successful and the process is logged in neogen-imputed.log

**Converting the imputed vcf back to plink format:**
    
     plink --vcf neogen-imputed.vcf --double-id --keep-allele-order --recode --make-bed --out neogen_cleaned_imputed

I realized that the phenotype data and sex could not be recovered from the imputed VCF file.

So manually edit the files to include the sex and phenotype in .ped and .fam files using R


## Merging neogen results markers and bovineHD markers from Roslin institute ##


**Subset the main ped and map files of Roslin data to include only the animals that were genotyped by neogen** 

    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile datConCas/datConCasBed --cow --keep /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results/LD/neogen_cleaned_imputed.fam --make-bed --out trial2/set1

**Apply filters similar to the ones mentioned in the paper with 1 extra step**

**Apply HWE filter**

**Apply LD filter**


    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile trial2/set1 --cow --allow-extra-chr --mind 0.1 --maf 0.05 --geno 0.1 --hwe 8.1e-6 --indep 10 5 4 --make-bed --out trial2/set1_clean

**Extract the SNPs that passed the LD filter**

    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile trial2/set1_clean --cow --extract trial2/set1_clean.prune.in --make-bed --out trial2/set1_ld_filtered

**Merging both the datsets**
         
    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile trial2/set1_ld_filtered --bmerge /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results/LD/neogen_cleaned_imputed --cow --out trial2/merge1
    
**Now ready to prepare the 3 pheno groups mentioned in the thesis by Raphaka from University of Edinburg**

The 3 pheno groups are manually made in excel spreadsheet NI_samples_forJHDB recoded_formatted.xlsx

- Merge1: all 1797 samples genotyped (1331 cases and 466 controls); 187237 variants which includes LD filtered Roslin data + 36 new SNPs (filtered MAF < 0.05 and missing CALL_RATE per variant > 80%) 

- Pheno1: VL+P (P in st_res, positive reactors) (525 cases, 460 controls)
 
- Pheno2: (VL+P)+(NVL+P) (1083 cases, 460 controls)

- Pheno3: (VL+P)+(NVL+P)+2.6% I+0.1% N (I : inconclusive in sv_res, N in st_res: negative reactors for skin test) (1115 cases,456 controls)


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

First load covariates and Q covariates files of all the samples (total no. of samples 1966):

    c<-read.delim("/mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/datConCas/datConCasCovarFile.txt",sep=" ", header=F, stringsAsFactors=F)  # load covariates file
    
    q<-read.delim("/mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/datConCas/datConCasQcovarFile.txt",sep=" ", header=F, stringsAsFactors=F) # load Quantitative covariate file

    colnames(c)<-c("id","V2","COV1","COV2","COV3","COV4","COV5")  # Assign column names for covariates 
    colnames(q)<-c("id","V2","QCOV")  		# and Quantitative covariates dataframes
    phe<-read.delim("pheno3.fam",sep=" ",header=F,stringsAsFactors=F)		# Load the .fam file to add the 																			# covariate columns 
    colnames(phe)<-c("id","V2","V3","V4","V5","pheno")		Assign column names for merging the dataframes
    head(phe)
    dim(phe)
    m<-merge(phe,c,by="id")
    m<-merge(m,q,by="id")
    m$COV<-rep(1,nrow(m)) 		# Add a column of 1s as the first column which acts as an intercept
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
    
    
    
    srun /mnt/nfs/nfs2/bickhart-users/binaries/gemma/bin/gemma -bfile /mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/trial2/merge1 -k output/merge1.cXX.txt -c merge1_cov.txt -lmm 1 -o merge1_gwas
    
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

