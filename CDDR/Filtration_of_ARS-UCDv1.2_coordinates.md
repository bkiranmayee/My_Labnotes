# Variants filtration using the coordinates directly lifted from ARS-UCDv14 to ARS-UCDv1.2

## This is a continuation from the assembly liftover notes [here](https://github.com/bkiranmayee/CDDR-Assembly-liftover-Project/blob/master/ARS-UCDV14%20to%20ARS-UCD.v1.2.md)

Prepared lifted over regions-sample file for progressive filtration

Out of 1784 haplotype groups only 1690 groups have been lifted over to the current assembly...

Working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration 
    
    [kiranmayee.bakshy@assembler2 filtration]$ bcftools stats f1.vcf.gz | grep -P "SN\t"
    # SN[2]id   [3]key  [4]value
    SN  0   number of samples:  172
    SN  0   number of records:  17883218
    SN  0   number of no-ALTs:  0
    SN  0   number of SNPs: 17883218
    SN  0   number of MNPs: 0
    SN  0   number of indels:   0
    SN  0   number of others:   0
    SN  0   number of multiallelic sites:   0
    SN  0   number of multiallelic SNP sites:   0
    [kiranmayee.bakshy@assembler2 filtration]$ bcftools stats f2.vcf.gz | grep -P "SN\t"
    # SN[2]id   [3]key  [4]value
    SN  0   number of samples:  172
    SN  0   number of records:  17883171
    SN  0   number of no-ALTs:  0
    SN  0   number of SNPs: 17883171
    SN  0   number of MNPs: 0
    SN  0   number of indels:   0
    SN  0   number of others:   0
    SN  0   number of multiallelic sites:   0
    SN  0   number of multiallelic SNP sites:   0
    [kiranmayee.bakshy@assembler2 filtration]$ bcftools annotate --rename-chrs ARS-UCDv1.2_NKLS_chr_num.tab -Oz -o f2.1.vcf.gz f2.vcf.gz
    

First getting rid of singletons and selecting only the variants which are homozygous in the target haplotype groups.

*Working directory:/mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/*

NoSingletons_OnlyHomo.vcf.gz

Summary:
 
	Kept 3440943 out of 17883171 variant sites
	Removed 12644107 not valid SNPs
	Removed 1652796 singletons
	Removed 145325 non-polymorphic SNPs

I tried to generate the initial truth and false datasets using the scripts progressiveSelection_1.pl and progressiveSelection_2.pl through pipeline.sh and pipeline2.sh

initial_truthset: all the sites homozygous haplotypes variants and no heterozygous sites at all
	
Summary: 

	Kept 2279176 out of 10418654 variant sites
	Removed 5 not valid sites
	Removed 1652796 singletons
	Removed 145325 non-polymorphic SNPs
	Removed 327885 that were heterozygous in all other regions
	Removed 471868 sites with AN < 337
	Removed 273314 sites with MQ0F >= 0.1
	Removed 5110651 sites with MQSB <= 0.95
	Removed 157634 328> DP >9002

initial_falseset: just the opposite of the above, select all the heterozygous sites in the homozygous haplotype regions...

Summary:

	Kept 36627 out of 23605381 variant sites
	Removed 5 not valid sites
	Removed 1652796 singletons
	Removed 145325 non-polymorphic SNPs
	Removed 15757160 that were not heterozygous in all other regions
	Removed 471868 sites with AN < 337
	Removed 273314 sites with MQ0F >= 0.1
	Removed 5110652 sites with MQSB <= 0.95
	Removed 157634328> DP >9002
	

I guess the falseset is not very useful...

Going on with further filtration based on LD and HWE:

*working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/initial_truthset*

	PLINK v1.90b4.4 64-bit (21 May 2017)
	Options in effect:
	  --cow
	  --double-id
	  --hwe 1e-05 midp
	  --indep-pairwise 1 kb 5 0.9
	  --keep-allele-order
	  --keep-autoconv
	  --out initial_truthset_hwe_ld
	  --vcf initial_truthset.vcf.gz
	  --write-snplist
	
	Hostname: assembler2.agil.barc.ba.ars.usda.gov
	Working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/initial_truthset
	Start time: Mon Oct 29 12:53:26 2018
	
	Random number seed: 1540832006
	515987 MB RAM detected; reserving 257993 MB for main workspace.
	--vcf: initial_truthset_hwe_ld.bed + initial_truthset_hwe_ld.bim +
	initial_truthset_hwe_ld.fam written.
	2279176 variants loaded from .bim file.
	172 cattle (0 males, 0 females, 172 ambiguous) loaded from .fam.
	Ambiguous sex IDs written to initial_truthset_hwe_ld.nosex .
	Using 1 thread (no multithreaded calculations invoked).
	Before main variant filters, 172 founders and 0 nonfounders present.
	Calculating allele frequencies... done.
	Total genotyping rate is 0.999748.
	--hwe: 6913 variants removed due to Hardy-Weinberg exact test.
	2272263 variants and 172 cattle pass filters and QC.
	Note: No phenotypes present.
	List of variant IDs written to initial_truthset_hwe_ld.snplist .
	Pruned 71638 variants from chromosome 1, leaving 99050.
	Pruned 56848 variants from chromosome 2, leaving 81413.
	Pruned 40010 variants from chromosome 3, leaving 55600.
	Pruned 28189 variants from chromosome 4, leaving 38360.
	Pruned 39123 variants from chromosome 5, leaving 54709.
	Pruned 35286 variants from chromosome 6, leaving 46693.
	Pruned 69269 variants from chromosome 7, leaving 94129.
	Pruned 31986 variants from chromosome 8, leaving 43737.
	Pruned 15531 variants from chromosome 9, leaving 21122.
	Pruned 41166 variants from chromosome 10, leaving 61391.
	Pruned 72539 variants from chromosome 11, leaving 105381.
	Pruned 9498 variants from chromosome 12, leaving 16400.
	Pruned 38633 variants from chromosome 13, leaving 60851.
	Pruned 61508 variants from chromosome 14, leaving 73642.
	Pruned 5140 variants from chromosome 15, leaving 7225.
	Pruned 35164 variants from chromosome 16, leaving 38603.
	Pruned 40904 variants from chromosome 17, leaving 49834.
	Pruned 14012 variants from chromosome 18, leaving 16430.
	Pruned 35051 variants from chromosome 19, leaving 47473.
	Pruned 49486 variants from chromosome 20, leaving 60623.
	Pruned 60613 variants from chromosome 21, leaving 80981.
	Pruned 13638 variants from chromosome 22, leaving 17145.
	Pruned 3257 variants from chromosome 23, leaving 4879.
	Pruned 27185 variants from chromosome 24, leaving 33462.
	Pruned 2696 variants from chromosome 25, leaving 5311.
	Pruned 31851 variants from chromosome 26, leaving 39014.
	Pruned 7796 variants from chromosome 27, leaving 10591.
	Pruned 3089 variants from chromosome 28, leaving 5247.
	Pruned 27235 variants from chromosome 29, leaving 34626.
	Pruning complete.  968341 of 2272263 variants removed.
	Marker lists written to initial_truthset_hwe_ld.prune.in and
	initial_truthset_hwe_ld.prune.out .
	
	End time: Mon Oct 29 12:54:11 2018

	
*The prunned in variants were then used to do the overlap analyses in the following folder: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/venn/ARS-UCDv1.2*
	
	
Link to the overlap [analyses](https://github.com/bkiranmayee/My_Labnotes/blob/master/CDDR/prep_dbsnp_for_vqsr.md)...
	