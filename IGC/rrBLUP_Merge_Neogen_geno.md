## Merging neogen results markers and HD bovine markers from Roslin institute ##


**Subset the main ped and map files to include only the animals that were genotyped by neogen** 

    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile datConCas/datConCasBed --cow --keep /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results/LD/neogen_ld_filtered.fam --make-bed --out trial_1/set1
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to trial_1/set1.log.
    Options in effect:
      --bfile datConCas/datConCasBed
      --cow
      --keep /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results/LD/neogen_ld_filtered.fam
      --make-bed
      --out trial_1/set1
    
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    538231 variants loaded from .bim file.
    1966 cattle (0 males, 1966 females) loaded from .fam.
    1966 phenotype values loaded from .fam.
    --keep: 852 cattle remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 852 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate in remaining samples is 0.99513.
    538231 variants and 852 cattle pass filters and QC.
    Among remaining phenotypes, 643 are cases and 209 are controls.
    --make-bed to trial_1/set1.bed + trial_1/set1.bim + trial_1/set1.fam ... done.




**Apply plink filters similar to the neogen markers cleaning**

    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile trial_1/set1 --cow --allow-extra-chr --mind 0.1 --maf 0.05 --geno 0.1 --make-bed --out trial_1/set1_cleaned
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to trial_1/set1_cleaned.log.
    Options in effect:
      --allow-extra-chr
      --bfile trial_1/set1
      --cow
      --geno 0.1
      --maf 0.05
      --make-bed
      --mind 0.1
      --out trial_1/set1_cleaned
    
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    538231 variants loaded from .bim file.
    852 cattle (0 males, 852 females) loaded from .fam.
    852 phenotype values loaded from .fam.
    0 cattle removed due to missing genotype data (--mind).
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 852 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.99513.
    13 variants removed due to missing genotype data (--geno).
    833 variants removed due to minor allele threshold(s)
    (--maf/--max-maf/--mac/--max-mac).
    537385 variants and 852 cattle pass filters and QC.
    Among remaining phenotypes, 643 are cases and 209 are controls.
    --make-bed to trial_1/set1_cleaned.bed + trial_1/set1_cleaned.bim +
    trial_1/set1_cleaned.fam ... done.



**Apply HWE filter**

    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile trial_1/set1_cleaned --cow --allow-extra-chr --hwe 8.1e-6 --make-bed --out trial_1/set1_hwe_filtered
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to trial_1/set1_hwe_filtered.log.
    Options in effect:
      --allow-extra-chr
      --bfile trial_1/set1_cleaned
      --cow
      --hwe 8.1e-6
      --make-bed
      --out trial_1/set1_hwe_filtered
    
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    537385 variants loaded from .bim file.
    852 cattle (0 males, 852 females) loaded from .fam.
    852 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 852 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.99513.
    Warning: --hwe observation counts vary by more than 10%.  Consider using
    --geno, and/or applying different p-value thresholds to distinct subsets of
    your data.
    --hwe: 10 variants removed due to Hardy-Weinberg exact test.
    537375 variants and 852 cattle pass filters and QC.
    Among remaining phenotypes, 643 are cases and 209 are controls.
    --make-bed to trial_1/set1_hwe_filtered.bed + trial_1/set1_hwe_filtered.bim +
    trial_1/set1_hwe_filtered.fam ... done.



**Apply LD filter**

    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile trial_1/set1_hwe_filtered --cow --indep 10 5 4 --out trial_1/ld_filter
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to trial_1/ld_filter.log.
    Options in effect:
      --bfile trial_1/set1_hwe_filtered
      --cow
      --indep 10 5 4
      --out trial_1/ld_filter
    
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    537375 variants loaded from .bim file.
    852 cattle (0 males, 852 females) loaded from .fam.
    852 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 852 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.99513.
    537375 variants and 852 cattle pass filters and QC.
    Among remaining phenotypes, 643 are cases and 209 are controls.
    --indep: Ignoring 645 chromosome 0 variants.
    Pruned 22367 variants from chromosome 1, leaving 11576.
    Pruned 19174 variants from chromosome 2, leaving 9516.
    Pruned 17187 variants from chromosome 3, leaving 8689.
    Pruned 16842 variants from chromosome 4, leaving 8917.
    Pruned 16466 variants from chromosome 5, leaving 8414.
    Pruned 17783 variants from chromosome 6, leaving 9085.
    Pruned 15909 variants from chromosome 7, leaving 7837.
    Pruned 12994 variants from chromosome 8, leaving 7250.
    Pruned 15148 variants from chromosome 9, leaving 7749.
    Pruned 14976 variants from chromosome 10, leaving 8081.
    Pruned 16208 variants from chromosome 11, leaving 8059.
    Pruned 12427 variants from chromosome 12, leaving 6767.
    Pruned 9576 variants from chromosome 13, leaving 5314.
    Pruned 9864 variants from chromosome 14, leaving 5438.
    Pruned 11779 variants from chromosome 15, leaving 6649.
    Pruned 11864 variants from chromosome 16, leaving 5804.
    Pruned 10902 variants from chromosome 17, leaving 5932.
    Pruned 9607 variants from chromosome 18, leaving 5344.
    Pruned 9087 variants from chromosome 19, leaving 5297.
    Pruned 10947 variants from chromosome 20, leaving 5647.
    Pruned 10050 variants from chromosome 21, leaving 5267.
    Pruned 9187 variants from chromosome 22, leaving 4992.
    Pruned 7035 variants from chromosome 23, leaving 4563.
    Pruned 9011 variants from chromosome 24, leaving 4701.
    Pruned 6163 variants from chromosome 25, leaving 3817.
    Pruned 7545 variants from chromosome 26, leaving 4201.
    Pruned 6400 variants from chromosome 27, leaving 3952.
    Pruned 6429 variants from chromosome 28, leaving 3903.
    Pruned 6884 variants from chromosome 29, leaving 4158.
    Pruning complete.  349811 of 536730 variants removed.
    Marker lists written to trial_1/ld_filter.prune.in and
    trial_1/ld_filter.prune.out .


**Extract the SNPs that passed the LD filter**

    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile trial_1/set1_hwe_filtered --cow --extract trial_1/ld_filter.prune.in --make-bed --out trial_1/set1_ld_filtered
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to trial_1/set1_ld_filtered.log.
    Options in effect:
      --bfile trial_1/set1_hwe_filtered
      --cow
      --extract trial_1/ld_filter.prune.in
      --make-bed
      --out trial_1/set1_ld_filtered
    
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    537375 variants loaded from .bim file.
    852 cattle (0 males, 852 females) loaded from .fam.
    852 phenotype values loaded from .fam.
    --extract: 186919 variants remaining.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 852 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.995283.
    186919 variants and 852 cattle pass filters and QC.
    Among remaining phenotypes, 643 are cases and 209 are controls.
    --make-bed to trial_1/set1_ld_filtered.bed + trial_1/set1_ld_filtered.bim +
    trial_1/set1_ld_filtered.fam ... done.


**Merging both the datsets**

    [kiranmayee.bakshy@assembler2 prelim_gwas]$ plink --bfile trial_1/set1_ld_filtered --bmerge /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results/LD/neogen_ld_filtered --cow --out trial_1/merge1
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Logging to trial_1/merge1.log.
    Options in effect:
      --bfile trial_1/set1_ld_filtered
      --bmerge /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results/LD/neogen_ld_filtered
      --cow
      --out trial_1/merge1
    
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    852 cattle loaded from trial_1/set1_ld_filtered.fam.
    852 cattle to be merged from
    /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results/LD/neogen_ld_filtered.fam.
    Of these, 0 are new, while 852 are present in the base dataset.
    186919 markers loaded from trial_1/set1_ld_filtered.bim.
    20 markers to be merged from
    /mnt/nfs/nfs2/bickhart-users/cattle_asms/ars_ucd_114_igc/neogen_results/LD/neogen_ld_filtered.bim.
    Of these, 20 are new, while 0 are present in the base dataset.
    Warning: Variants 'BovineHD0400001282' and 'BTA-98753-no-rs' have the same
    position.
    Warning: Variants 'BovineHD4100003297' and 'ARS-BFGL-NGS-119660' have the same
    position.
    Warning: Variants 'Hapmap53362-rs29013727' and
    'ARS-USMARC-Parent-DQ990834-rs29013727' have the same position.
    Performing single-pass merge (852 cattle, 186939 variants).
    Merged fileset written to trial_1/merge1.bed + trial_1/merge1.bim +
    trial_1/merge1.fam .


**Step1:Recoding the genotypes to get a partially compatible rrBLUP geno format**

    [kiranmayee.bakshy@assembler2 trial_1]$ plink --bfile merge1 --cow --recode12 tab --transpose --out merge1_recoded
    PLINK v1.90b4.4 64-bit (21 May 2017)   www.cog-genomics.org/plink/1.9/
    (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    Note: --recode12 flag deprecated.  Use 'recode 12 ...'.
    Logging to merge1_recoded.log.
    Options in effect:
      --bfile merge1
      --cow
      --out merge1_recoded
      --recode 12 tab
      --transpose
    
    Note: --transpose flag deprecated.  Use '--recode transpose ...'.
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    186939 variants loaded from .bim file.
    852 cattle (0 males, 852 females) loaded from .fam.
    852 phenotype values loaded from .fam.
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 852 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.995283.
    186939 variants and 852 cattle pass filters and QC.
    Among remaining phenotypes, 643 are cases and 209 are controls.
    --recode transpose to merge1_recoded.tped + merge1_recoded.tfam ... done.
 
**Plink --recode12 recodes the alleles; minor allele A1(a)=1, Major allele A2(A)=2, missing=0 0**

**Genotypes A A = 2 2; A a = 1 2 or 2 1; a a = 1 1**

**Step2: Custom perl script (tped_2_geno.pl) to achieve the genotypes coded as -1, 0 and 1**

    #!/usr/bin/perl 
    
    ## A script to recode tped format file to geno compatible with rrBLUP for GWAS
    ## Date: Jul 5 2018
    
    use strict;
    
    # make sure the user passed the required arguments
    if (scalar @ARGV != 2 ) {
       print STDERR "Usage: tped_2geno.pl <tped file> <output file>\n";
     	  exit(1);
    }
    chomp(@ARGV);
    my ($tpedfile, $outfile) = @ARGV;
    
    open(my $ifh, "<", $tpedfile) || die "ERROR: failed to read tped file: $!";
    
    open (my $OUT, "> $outfile");
    
    my @newsegs;
    
    while(my $line = <$ifh>){
    	chomp $line;
    	my @segs = split(/\t/, $line);
    	
    			for(my $x = 0; $x < scalar(@segs); $x++){
    					if($segs[$x] eq "0 0"){
    						$newsegs[$x] = "NA";
    					}
    					elsif($segs[$x] eq "1 1"){
    						$newsegs[$x] = "-1";
    					}
    					elsif($segs[$x] eq "1 2" || $segs[$x] eq "2 1"){
    						$newsegs[$x] = "0";
    					}
    					elsif($segs[$x] eq "2 2"){
    						$newsegs[$x] = "1";
    					}else{$newsegs[$x] = $segs[$x];
    						}
    		}
    		print {$OUT} join("\t", @newsegs) . "\n";
    }
    			
    
    close $ifh;
    close $OUT;
    exit;

**Using the above script I have coded the genotypes AA=1, Aa=0, aa=-1 and missing(0 in tped)=NA**

