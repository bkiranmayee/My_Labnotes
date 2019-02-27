# Neogen genotyping data second round markers data processing #

**Lab notes on Nov 23 2018**

## Steps to be done ##

- Generate the initial ped and map files as done previously
- Filter based on MAF > 0.05 and genotyping call rate > 0.8
- Give arbitrary UMD3 coordinates to be included in the GWAS 
- Imputation using Beagle
- Merge this filtered data into the previous data set and run GWAS
- Plot manhattan and QQplots
- Run LD analysis and haplotype block analysis
 
**1) Convert the neogen genotype file to plink format** 

This was done by a custom script: neogen\_to\_ped\_2.pl

    #!/usr/bin/perl 
    
    ## A script to convert neogen results to ped format
    ## Date: Nov 23 2018
    
    use strict;
    
    my $input = "neogen_results_11_13_2018.tab";
    open (my $OUTFILE, "> neogen_2ndRound.ped");
    
    my @newcols;
    my @new;
    
    my %varRefAlleles = (
    	"18_62572950" => "C",
    	"18_62665698" => "A",
    	"18_62716825" => "G",
    	"18_62766196" => "G",
    	"18_62849303" => "G",
    	"18_62994372" => "C",
    	"18_63020431" => "A",
    	"18_63036451" => "G",
    	"18_63050618" => "T",
    	"18_63067935" => "C",
    	"18_63082203" => "A",
    	"18_63084493" => "C",
    	"18_63114542" => "G",
    	"18_63122963" => "C",
    	"18_63131111" => "T",
    	"18_63140201" => "G",
    	"18_63141688" => "G",
    	"18_63156196" => "C",
    	"18_63185710" => "G",
    	"18_63186702" => "T",
    	"18_63284810" => "T",
    	"18_63294645" => "G",
    	"18_63340289" => "T",
    	"18_63385249" => "C",
    	"18_63399311" => "C",
    	"18_63417698" => "T",
    	"18_63442255" => "T",
    	"23_28648192" => "T",
    	"5_99023568" => "C",
    	"5_99094858" => "A",
    	"5_99190989" => "T",
    	"5_99194167" => "T",
    	"5_99203402" => "C",
    	"5_99329862" => "T",
    	"5_99397443" => "T",
    	"5_99412564" => "G",
    	"5_99437819" => "A",
    	"5_99462993" => "A",
    	"5_99475790" => "A",
    	"5_99714227" => "A",
    	"5_99749976" => "T",
    	"LIB14427_MHC_2500" => "T",
    	"LIB14427_MHC_3538" => "G",
    	"LIB14427_MHC_36271" => "C",
    	"LIB14427_MHC_45690" => "G",
    	"LIB14427_MHC_57505" => "C",
    	"LIB14427_MHC_73766" => "G",
    	"LIB14427_MHC_116810" => "G",
    	"LIB14427_MHC_122678" => "A",
    	"LIB14427_MHC_126886" => "G",
    	"MHCclassI_MHC_120784" => "A",
    	"MHCclassI_MHC_150472" => "A",
    	"MHCclassI_MHC_181938" => "T",
    	"MHCclassI_MHC_223262" => "A",
    	"MHCclassI_MHC_302664" => "A",
    	"MHCclassI_MHC_392940" => "A",
    	"MHCclassI_MHC_395846" => "G");
    
    open(my $IN, "< $input") || die "Could not open input file: $input!\n";
    
    my %varIdx;
    while(my $line = <$IN>){
    	chomp $line;
    	my @hsegs = split(/\t/, $line);
    	
    	for(my $x = 5; $x < scalar(@hsegs); $x++){
    		$varIdx{$x} = $hsegs[$x];
    		}
    	print $OUTFILE "$line\n";
    	last;
    }
    
    my @newsegs;
    while(my $line = <$IN>){
    	chomp $line;
    	my @segs = split(/\t/, $line);
    	
    	for(my $x = 5; $x < scalar(@segs); $x++){
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

 A map file was then generated...

 grep -f neogen_marker_ids.txt marker_composition_master_file.txt > marker_comp_2ndRound.txt

I converted the map and ped files to binary plink format by giving the Ref allele text file so that it does not change the alleles based on frequency. I have also cleaned the dataset in the same step:

    working directory: /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results/round_2
    
    plink --file neogen_2ndRound --geno 0.2 --maf 0.05 --cow --allow-extra-chr --make-bed --recode --a2-allele marker_ref_alleles_2ndround.txt --out neogen_cleaned


I cleaned the data based on call rate (reject > 80% missing per variant (--geno 0.2) and MAF (--maf 0.05). I get 44 variants left after the cleaning step. (saved the dataset as neogen_cleaned)


Then I gave arbitrary UMD3 coordinates to the passed markers and updated the SNP IDs manually in an excel sheet (UMD3_arbitrary_coordinates.xlsx).

I saved the original SNP coordinates map file as neogen_cleaned_original.map. 

I also updated the IDs and Ref alleles in the ref alleles file and saved it as selected_markers_ref 

I then converted it to VCF format for imputation using Beagle.

    plink --file neogen_cleaned --geno 0.2 --maf 0.05 --cow --allow-extra-chr --recode vcf-iid --a2-allele selected_markers_ref --out neogen_cleaned
    
    java -Xmx50g -D$JAVA_HOME -jar /mnt/nfs/nfs2/bickhart-users/binaries/beagle.03Jul18.40b.jar gt=neogen_cleaned.vcf out=neogen-imputed


The imputation was successful and the process is logged in neogen-imputed.log


I realized that the phenotype data and sex could not be recovered from the imputed VCF file.

So I manually edit the files to include the sex and phenotype in .ped and .fam files


     [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile trial2/merge1 --bmerge /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results/round_2/neogen_cleaned_imputed --cow --out round2/merge1


**Now ready to prepare the 3 pheno groups mentioned in the thesis by Raphaka from University of Edinburg**

The 3 pheno groups are manually made in excel spreadsheet NI_samples_forJHDB recoded_formatted.xlsx

- Merge1: all 1797 samples genotyped (1331 cases and 466 controls); 187281 variants which includes LD filtered Roslin data + 36 first round markers + 44 second round SNPs (filtered MAF < 0.05 and missing CALL_RATE per variant > 80%) 

- Pheno1: VL+P (P in st_res, positive reactors) (525 cases, 460 controls)
 
- Pheno2: (VL+P)+(NVL+P) (1083 cases, 460 controls)

- Pheno3: (VL+P)+(NVL+P)+2.6% I+0.1% N (I : inconclusive in sv_res, N in st_res: negative reactors for skin test) (1115 cases,456 controls)


Its easy to prepare all three groups using plink then the phenotype coding can be changed in R

    
    [kiranmayee.bakshy@assembler2 trial2]$ plink --bfile merge1 --cow --keep-fam Pheno1.txt --make-bed --out pheno1
    [kiranmayee.bakshy@assembler2 trial2]$ plink --bfile merge1 --cow --keep-fam Pheno2.txt --make-bed --out pheno2
    [kiranmayee.bakshy@assembler2 trial2]$ plink --bfile merge1 --cow --keep-fam Pheno3.txt --make-bed --out pheno3

**Phenotype recoding to controls=0, cases=1**
 I manually edited the phenotype codes using vim editing tools

I used the same covariate files which I used for my previous GWAS study...(/mnt/nfs/nfs2/bickhart-users/natdb_sequencing/prelim_gwas/trial2)

Running GWAS using GEMMA within the script gemma.sh and plotting the results within the script plot.sh

Hope it works this time....

Yes! two markers seems to be associated with bTB at both genome-wide and suggestive significance threshold...

But GEMMA package is better suited for continuous phenotype while the p-values are usually inflated for binary phenotype...

I am going to try GWAS using GMMAT package in R because as opposed to GEMMA, GMMAT implements generalized linear mixed model and also we have an option to select the type of distribution of out phenotype data. 


Meanwhile, here is the **Distance vs R2 plot and haploblock view of the 43 markers selected for the GWAS study**:


![](https://i.imgur.com/bMAsZvz.png)

Some of the selected variants are in high LD. 

It is interesting to note that the 2 SNPs which were significantly associated with bTB in GWAS implemented in GEMMA have the following features: 

- R2 = 0.32
- D` = 0.77
- Distance = 375492 kb
- MAF of 18_62766196 = 0.48
- MAF of 18_63141688 = 0.33

Also, both these SNPs have HWE p-value < 10-8 

I think these variants would get filtered out if we apply LD and HWE filters!!!


![](https://i.imgur.com/4h86gVF.png)


Haplotype association analyses for binary trait can be done using Haplo.stats package in R...
	