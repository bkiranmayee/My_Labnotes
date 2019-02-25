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


The out put was then manually edited to include the extra columns such as Family ID, maternal ID, paternal ID, Sex and phenotype from NI\_samples\_forJHDB recoded_formatted.xlsx file. 

 A map file was then generated


**1b) Converting the plink file format to vcf keeping in mind the Ref allele:**

     [kiranmayee.bakshy@assembler2 LD]$ plink --file neogen_updated --allow-extra-chr --recode vcf-iid --a2-allele Ref_allele.txt --out neogen_cleaned


**1c) Imputation of Genotypes using BEAGLE**

Input and output files for Beagle is VCF format

    [kiranmayee.bakshy@assembler2 LD]$ srun java -Xmx50g -D$JAVA_HOME -jar /mnt/nfs/nfs2/bickhart-users/binaries/beagle.03Jul18.40b.jar gt=neogen_cleaned.vcf out=neogen-imputed

The imputation seems to be successful for all the variants except for KIR_41480 because there is only one marker that belongs to that Chromosome.
Mmmm! I should first make fake UMD3 coordinates and then do the imputation I guess!!

OK redo the whole thing again.

I converted the map and ped files to binary plink format by giving the Ref allele text file so that it does not change the alleles based on frequency. I have also cleaned the dataset:

     plink --file neogen_updated --geno 0.2 --maf 0.05 --cow --allow-extra-chr --make-bed --recode --a2-allele Ref_allele.txt --out neogen_cleaned

I cleaned the data based on call rate (reject > 80% missing per variant (--geno 0.2) and MAF (--maf 0.05). I get only 36 variants left. (saved the dataset as neogen_cleaned)

Then I gave arbitrary UMD3 coordinates and updated the SNP IDs manually in an excel sheet (UMD3_fake_coordinates.xlsx).

I saved the original SNP coordinates map file as neogen_cleaned_original.map. 

I then converted it to VCF format for imputation using Beagle

    plink --file neogen_cleaned --geno 0.2 --maf 0.05 --cow --allow-extra-chr --recode vcf-iid --a2-allele Ref_allele.txt --out neogen_cleaned
    java -Xmx50g -D$JAVA_HOME -jar /mnt/nfs/nfs2/bickhart-users/binaries/beagle.03Jul18.40b.jar gt=neogen_cleaned.vcf out=neogen-imputed

Now the imputation was successful and the process is logged in neogen-imputed.log

Converting the imputed vcf back to plink format:




**2) Filtration using plink 1.9**

Filters such as missing per sample, per variant, MAF, HWE and LD were applied on the data using plink 

    [kiranmayee.bakshy@assembler2 LD]$ plink --file neogen_updated --cow --allow-extra-chr --mind 0.1 --maf 0.05 --geno 0.1 --make-bed --out neogen_cleaned
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to neogen_cleaned.log.
    Options in effect:
      --allow-extra-chr
      --file neogen_updated
      --geno 0.1
      --maf 0.05
      --make-bed
      --mind 0.1
      --out neogen_cleaned
    
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    .ped scan complete (for binary autoconversion).
    Performing single-pass .bed write (67 variants, 1797 people).
    --file: neogen_cleaned-temporary.bed + neogen_cleaned-temporary.bim +
    neogen_cleaned-temporary.fam written.
    67 variants loaded from .bim file.
    1797 people (0 males, 1797 females) loaded from .fam.
    1797 phenotype values loaded from .fam.
    945 people removed due to missing genotype data (--mind).
    IDs written to neogen_cleaned.irem .
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 852 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate in remaining samples is 0.916176.
    6 variants removed due to missing genotype data (--geno).
    25 variants removed due to minor allele threshold(s)
    (--maf/--max-maf/--mac/--max-mac).
    36 variants and 852 people pass filters and QC.
    Among remaining phenotypes, 643 are cases and 209 are controls.
    --make-bed to neogen_cleaned.bed + neogen_cleaned.bim + neogen_cleaned.fam ...
    done.

**Calculate Hardy-Weinberg exact p-value using plink**
    
    [kiranmayee.bakshy@assembler2 LD]$ plink --bfile neogen_cleaned  --cow --allow-extra-chr --hardy midp --out hardy
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to hardy.log.
    Options in effect:
      --allow-extra-chr
      --bfile neogen_cleaned
      --hardy midp
      --out hardy
    
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    36 variants loaded from .bim file.
    852 people (0 males, 852 females) loaded from .fam.
    852 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 852 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.991817.
    --hardy: Writing Hardy-Weinberg report (founders only) to hardy.hwe ... done.

**The following are the SNPs removed based on HWE p-values:**

**CHR** |**SNP**|**TEST**|**A1**|**A2**|**GENO**|**O(HET)**|**E(HET)**|**P**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
5|5\_98985645|UNAFF|A|G|0/209/0|1|0.5|1.56E-62
MHC|MHC\_380177|UNAFF|T|C|77/0/128|0|0.4691|8.36E-60
MHC|MHC\_317666|UNAFF|A|G|7/170/30|0.8213|0.4938|8.21E-24
MHC|MHC\_260913|UNAFF|G|A|57/30/119|0.1456|0.4547|2.09E-23
18|18\_62812297|UNAFF|G|A|12/166/31|0.7943|0.4959|2.14E-19
18|18\_62774779|UNAFF|A|C|43/31/135|0.1483|0.4031|2.39E-19
MHC|MHC\_154399|UNAFF|C|A|52/35/117|0.1716|0.4492|3.91E-19
MHC|MHC\_86084|UNAFF|A|G|19/14/176|0.06699|0.2179|2.23E-16
LRC|LRC\_81028|UNAFF|T|C|58/49/102|0.2344|0.4778|6.18E-14
MHC|MHC\_121069|UNAFF|A|G|11/17/174|0.08416|0.1744|7.47E-09
MHC|MHC\_43656|UNAFF|G|T|13/22/174|0.1053|0.2033|8.09E-09


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    [kiranmayee.bakshy@assembler2 LD]$ plink --bfile neogen_cleaned --cow --allow-extra-chr --hwe 8.1e-6 --recode --out neogen_hwe_filtered
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to neogen_hwe_filtered.log.
    Options in effect:
      --allow-extra-chr
      --bfile neogen_cleaned
      --hwe 8.1e-6
      --out neogen_hwe_filtered
      --recode
      
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    36 variants loaded from .bim file.
    852 people (0 males, 852 females) loaded from .fam.
    852 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 852 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.991817.
    --hwe: 11 variants removed due to Hardy-Weinberg exact test.
    25 variants and 852 people pass filters and QC.
    Among remaining phenotypes, 643 are cases and 209 are controls.
    --recode ped to neogen_hwe_filtered.ped + neogen_hwe_filtered.map ... done.


**Calculate LD for the remaining 25 variants**


        [kiranmayee.bakshy@assembler2 LD]$ plink --file neogen_hwe_filtered --allow-extra-chr --r2 inter-chr dprime with-freqs
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to plink.log.
    Options in effect:
      --allow-extra-chr
      --file neogen_hwe_filtered
      --r2 inter-chr dprime with-freqs
    
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    .ped scan complete (for binary autoconversion).
    Performing single-pass .bed write (25 variants, 852 people).
    --file: plink-temporary.bed + plink-temporary.bim + plink-temporary.fam
    written.
    25 variants loaded from .bim file.
    852 people (0 males, 852 females) loaded from .fam.
    852 phenotype values loaded from .fam.
    Using up to 63 threads (change this with --threads).
    Before main variant filters, 852 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.991831.
    25 variants and 852 people pass filters and QC.
    Among remaining phenotypes, 643 are cases and 209 are controls.
    --r2 inter-chr dprime with-freqs... done.
    Results written to plink.ld .



**CHR\_A**|**BP\_A**|**SNP\_A**|**MAF\_A**|**CHR\_B**|**BP\_B**|**SNP\_B**|**MAF\_B**|**R2**|**DP**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
5|99036250|5\_99036250|0.265492|5|99333595|5\_99333595|0.241197|0.228196|0.519767
23|28535354|23\_28535354|0.272727|MHC|115082|MHC\_115082|0.226148|0.25692|0.574356
23|28535354|23\_28535354|0.272727|MHC|286137|MHC\_286137|0.0766944|0.223761|1
LRC|65014|LRC\_65014|0.493381|LRC|98785|LRC\_98785|0.288732|0.251829|0.796809
LRC|65014|LRC\_65014|0.493381|LRC|174904|LRC\_174904|0.176887|0.202538|0.989097
LRC|98785|LRC\_98785|0.288732|LRC|174904|LRC\_174904|0.176887|0.420333|0.888905
MHC|9213|MHC\_9213|0.10446|MHC|115082|MHC\_115082|0.226148|0.23176|0.7605


**The distance between the SNPs was calculated and R2 was plotted against the distance (bp) in R**

	R    
	df<-read.delim("plink_lifted.ld", header=T, sep="", stringsAsFactors = F)
	df[3,]<-NULL  ## remove the SNP pairs which are not on the same Chromosome
    df$dist<-df$BP_B-df$BP_A ## add new column, dist which is the difference between the SNPs positions
    library(ggplot2)
    p<-ggplot(df, aes(dist,R2)) + geom_point(aes(size=R2))+labs(x="Distance (bp)",y="LD (R squared)") + geom_text(label=paste(df$SNP_A,"-", df$SNP_B), hjust="inward", vjust="inward")
	dev.copy(p, 'R2_distance_plot.jpg')
	dev.off()
    write.table(df, file="r2_dist.txt", row.names=F, quote=F, sep="\t")



**CHR\_A**|**BP\_A**|**SNP\_A**|**MAF\_A**|**CHR\_B**|**BP\_B**|**SNP\_B**|**MAF\_B**|**R2**|**DP**|**dist**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
5|99036250|5\_99036250|0.265492|5|99333595|5\_99333595|0.241197|0.228196|0.519767|297345
23|28535354|23\_28535354|0.272727|23|28631651|*MHC\_115082(lifted_coordinate)|0.226148|0.25692|0.574356|96297
LRC|65014|LRC\_65014|0.493381|LRC|98785|LRC\_98785|0.288732|0.251829|0.796809|33771
LRC|65014|LRC\_65014|0.493381|LRC|174904|LRC\_174904|0.176887|0.202538|0.989097|109890
LRC|98785|LRC\_98785|0.288732|LRC|174904|LRC\_174904|0.176887|0.420333|0.888905|76119
MHC|9213|MHC\_9213|0.10446|MHC|115082|MHC\_115082|0.226148|0.23176|0.7605|105869




![R2-dist_plot](https://i.imgur.com/Qr1yp8x.jpg)


**The following LD-plot for all 25 SNPs was created using Haploview:**

![LD-plot](https://i.imgur.com/Hp6UkdE.png)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    [kiranmayee.bakshy@assembler2 LD]$ plink --file neogen_hwe_filtered --allow-extra-chr --indep-pairwise 10 5 0.2 --recode --out neogen_ld_filtered
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to neogen_ld_filtered.log.
    Options in effect:
      --allow-extra-chr
      --file neogen_hwe_filtered
      --indep-pairwise 10 5 0.2
      --out neogen_ld_filtered
      --recode
    
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    .ped scan complete (for binary autoconversion).
    Performing single-pass .bed write (25 variants, 852 people).
    --file: neogen_ld_filtered-temporary.bed + neogen_ld_filtered-temporary.bim +
    neogen_ld_filtered-temporary.fam written.
    25 variants loaded from .bim file.
    852 people (0 males, 852 females) loaded from .fam.
    852 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 852 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.991831.
    25 variants and 852 people pass filters and QC.
    Among remaining phenotypes, 643 are cases and 209 are controls.
    --recode ped to neogen_ld_filtered.ped + neogen_ld_filtered.map ... done.
    Pruned 1 variant from chromosome 5, leaving 3.
    Pruned 0 variants from chromosome 18, leaving 5.
    Pruned 0 variants from chromosome 23, leaving 2.
    Pruned 0 variants from chromosome 27, leaving 1.
    Pruned 2 variants from chromosome 28, leaving 4.
    Pruned 2 variants from chromosome 29, leaving 5.
    Pruning complete.  5 of 25 variants removed.
    Marker lists written to neogen_ld_filtered.prune.in and
    neogen_ld_filtered.prune.out .



**List of SNPs removed based on LD calculations:**

 * 5\_99333595
 * LRC\_61352
 * LRC\_174904
 * MHC\_9213
 * MHC\_324231
    

**Now creating the final filtered files using plink** 

    [kiranmayee.bakshy@assembler2 LD]$ plink --file  neogen_hwe_filtered --allow-extra-chr --extract neogen_ld_filtered.prune.in --recode --make-bed --out neogen_ld_filtered
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to neogen_ld_filtered.log.
    Options in effect:
      --allow-extra-chr
      --extract neogen_ld_filtered.prune.in
      --file neogen_hwe_filtered
      --make-bed
      --out neogen_ld_filtered
      --recode
    
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    .ped scan complete (for binary autoconversion).
    Performing single-pass .bed write (25 variants, 852 people).
    --file: neogen_ld_filtered-temporary.bed + neogen_ld_filtered-temporary.bim +
    neogen_ld_filtered-temporary.fam written.
    25 variants loaded from .bim file.
    852 people (0 males, 852 females) loaded from .fam.
    852 phenotype values loaded from .fam.
    --extract: 20 variants remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 852 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.991784.
    20 variants and 852 people pass filters and QC.
    Among remaining phenotypes, 643 are cases and 209 are controls.
    --make-bed to neogen_ld_filtered.bed + neogen_ld_filtered.bim +
    neogen_ld_filtered.fam ... done.
    --recode ped to neogen_ld_filtered.ped + neogen_ld_filtered.map ... done.
