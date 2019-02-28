# Data generation and pre-processing for NatDb imputation #
##  ##

### Non-IGC variants processing ##

working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3

- The high quality filtered data set that was generated in the UMD.3.1.1 coordinates system consists of 1,288,386 variants. (filtered_ld.vcf.gz)

   `bcftools view -R selected_snps.sorted -O z -o selected.snps.vcf.gz filtered_ld.vcf.gz`

- After removing the variants which are redundant in the bovineHD and 1000bulls dataset, we are left with a total of 218714 variants. (selected.snps.vcf.gz)

Step 1: Further reduce the no. of variants by removing the markers which are in LD with R2 > 0.5

    module load plink/1.90b4.4-2017-05-21
    
    plink --vcf selected.snps.vcf.gz --cow --double-id --indep-pairwise 1 kb 5 0.5
    
    PLINK v1.90b4.4 64-bit (21 May 2017)
    Options in effect:
      --cow
      --double-id
      --indep-pairwise 1 kb 5 0.5
      --out ld.filtration.2
      --vcf selected.snps.vcf.gz
    
    Hostname: assembler2.agil.barc.ba.ars.usda.gov
    Working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3
    Start time: Fri Sep 14 10:50:08 2018
    
    Random number seed: 1536936608
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    --vcf: ld.filtration.2-temporary.bed + ld.filtration.2-temporary.bim +
    ld.filtration.2-temporary.fam written.
    218714 variants loaded from .bim file.
    172 cattle (0 males, 0 females, 172 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to ld.filtration.2.nosex .
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 172 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.999436.
    218714 variants and 172 cattle pass filters and QC.
    Note: No phenotypes present.
    Pruned 319 variants from chromosome 1, leaving 13495.
    Pruned 298 variants from chromosome 2, leaving 11879.
    Pruned 217 variants from chromosome 3, leaving 8194.
    Pruned 75 variants from chromosome 4, leaving 3916.
    Pruned 172 variants from chromosome 5, leaving 7219.
    Pruned 136 variants from chromosome 6, leaving 4896.
    Pruned 238 variants from chromosome 7, leaving 11078.
    Pruned 89 variants from chromosome 8, leaving 4280.
    Pruned 37 variants from chromosome 9, leaving 2305.
    Pruned 155 variants from chromosome 10, leaving 7434.
    Pruned 284 variants from chromosome 11, leaving 12907.
    Pruned 35 variants from chromosome 12, leaving 1996.
    Pruned 143 variants from chromosome 13, leaving 8123.
    Pruned 175 variants from chromosome 14, leaving 7965.
    Pruned 14 variants from chromosome 15, leaving 1083.
    Pruned 158 variants from chromosome 16, leaving 5336.
    Pruned 123 variants from chromosome 17, leaving 5617.
    Pruned 47 variants from chromosome 18, leaving 1694.
    Pruned 83 variants from chromosome 19, leaving 4517.
    Pruned 263 variants from chromosome 20, leaving 9505.
    Pruned 257 variants from chromosome 21, leaving 11013.
    Pruned 25 variants from chromosome 22, leaving 1412.
    Pruned 15 variants from chromosome 23, leaving 432.
    Pruned 93 variants from chromosome 24, leaving 4471.
    Pruned 2 variants from chromosome 25, leaving 331.
    Pruned 171 variants from chromosome 26, leaving 6072.
    Pruned 7 variants from chromosome 27, leaving 797.
    Pruned 5 variants from chromosome 28, leaving 587.
    Pruned 52 variants from chromosome 29, leaving 2798.
    Pruned 3933 variants from chromosome 30, leaving 49741.
    Pruning complete.  7621 of 218714 variants removed.
    Marker lists written to ld.filtration.2.prune.in and ld.filtration.2.prune.out
    .
    
    End time: Fri Sep 14 10:50:13 2018

    gunzip -c selected.snps.vcf.gz | grep -v -f <(cat ld.filtration.2.prune.out) > selected.ld.prunned.vcf

Oh! Ok this removed only 7621 variants...we are left with 211093 variants (selected.ld.prunned.vcf) which is still out of our expected number.


Lifting over these coordinates to ARS-UCD.v1.2 assembly version using Liftovervcf 

    java -jar /mnt/nfs/nfs1/kiranmayee.bakshy/picard/build/libs/picard.jar LiftoverVcf I=selected.ld.prunned.vcf O=selected.ARS-UCD.v1.2.vcf CHAIN=/mnt/nfs/nfs2/bickhart-users/cattle_asms/liftovers/UMD3.11_to_ARS-UCD1.2/umd3_kary_unmask_ngap_to_ARS-UCD1.2.mmap.liftover.chain REJECT=selected.r.vcf R=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa WRITE_ORIGINAL_POSITION=true
    12:28:08.305 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/mnt/nfs/nfs1/kiranmayee.bakshy/picard/build/libs/picard.jar!/com/intel/gkl/native/libgkl_compression.so
    [Fri Sep 14 12:28:08 EDT 2018] LiftoverVcf INPUT=selected.ld.prunned.vcf OUTPUT=selected.ARS-UCD.v1.2.vcf CHAIN=/mnt/nfs/nfs2/bickhart-users/cattle_asms/liftovers/UMD3.11_to_ARS-UCD1.2/umd3_kary_unmask_ngap_to_ARS-UCD1.2.mmap.liftover.chain REJECT=selected.r.vcf WRITE_ORIGINAL_POSITION=true REFERENCE_SEQUENCE=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.faWARN_ON_MISSING_CONTIG=false LOG_FAILED_INTERVALS=true WRITE_ORIGINAL_ALLELES=false LIFTOVER_MIN_MATCH=1.0 ALLOW_MISSING_FIELDS_IN_HEADER=false RECOVER_SWAPPED_REF_ALT=false TAGS_TO_REVERSE=[AF] TAGS_TO_DROP=[MAX_AF] VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false
    [Fri Sep 14 12:28:08 EDT 2018] Executing as kiranmayee.bakshy@assembler2.agil.barc.ba.ars.usda.gov on Linux 4.8.13-100.fc23.x86_64 amd64; OpenJDK 64-Bit Server VM 1.8.0_111-b16; Deflater: Intel; Inflater: Intel; Provider GCS is not available; Picard version: 2.18.4-SNAPSHOT
    INFO2018-09-14 12:28:10 LiftoverVcf Loading up the target reference genome.
    INFO2018-09-14 12:28:37 LiftoverVcf Lifting variants over and sorting (not yet writing the output file.)
    INFO2018-09-14 12:29:34 LiftoverVcf Processed 211093 variants.
    INFO2018-09-14 12:29:34 LiftoverVcf 605 variants failed to liftover.
    INFO2018-09-14 12:29:34 LiftoverVcf 602 variants lifted over but had mismatching reference alleles after lift over.
    INFO2018-09-14 12:29:34 LiftoverVcf 0.5718% of variants were not successfully lifted over and written to the output.
    INFO2018-09-14 12:29:34 LiftoverVcf liftover success by source contig:
    INFO2018-09-14 12:29:34 LiftoverVcf chr1: 13427 / 13495 (99.4961%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr10: 7401 / 7434 (99.5561%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr11: 12880 / 12907 (99.7908%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr12: 1976 / 1996 (98.9980%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr13: 8093 / 8123 (99.6307%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr14: 7890 / 7965 (99.0584%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr15: 1070 / 1083 (98.7996%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr16: 5293 / 5336 (99.1942%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr17: 5593 / 5617 (99.5727%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr18: 1648 / 1694 (97.2845%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr19: 4487 / 4517 (99.3358%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr2: 11853 / 11879 (99.7811%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr20: 9440 / 9505 (99.3161%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr21: 10934 / 11013 (99.2827%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr22: 1395 / 1412 (98.7960%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr23: 430 / 432 (99.5370%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr24: 4463 / 4471 (99.8211%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr25: 327 / 331 (98.7915%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr26: 5984 / 6072 (98.5507%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr27: 786 / 797 (98.6198%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr28: 583 / 587 (99.3186%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr29: 2761 / 2798 (98.6776%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr3: 8159 / 8194 (99.5729%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr4: 3883 / 3916 (99.1573%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr5: 7189 / 7219 (99.5844%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr6: 4807 / 4896 (98.1822%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr7: 10963 / 11078 (98.9619%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr8: 4245 / 4280 (99.1822%)
    INFO2018-09-14 12:29:34 LiftoverVcf chr9: 2288 / 2305 (99.2625%)
    INFO2018-09-14 12:29:34 LiftoverVcf chrX: 49638 / 49741 (99.7929%)
    INFO2018-09-14 12:29:34 LiftoverVcf lifted variants by target contig:
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000001.1: 13496
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000002.1: 11855
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000003.1: 8153
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000004.1: 3874
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000005.1: 7131
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000006.1: 4801
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000007.1: 10917
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000008.1: 4202
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000009.1: 2294
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000010.1: 7404
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000011.1: 12872
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000012.1: 2040
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000013.1: 8029
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000014.1: 7871
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000015.1: 1068
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000016.1: 5292
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000017.1: 5450
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000018.1: 1652
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000019.1: 4484
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000020.1: 9457
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000021.1: 10931
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000022.1: 1395
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000023.1: 486
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000024.1: 4463
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000025.1: 327
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000026.1: 5981
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000027.1: 787
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000028.1: 597
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000029.1: 2764
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000030.1: 49798
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000412.1: 1
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000591.1: 2
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02000661.1: 1
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02001064.1: 1
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02001545.1: 1
    INFO2018-09-14 12:29:34 LiftoverVcf NKLS02002204.1: 9
    WARNING 2018-09-14 12:29:34 LiftoverVcf 597 variants with a swapped REF/ALT were identified, but were not recovered.  See RECOVER_SWAPPED_REF_ALT and associated caveats.
    INFO2018-09-14 12:29:34 LiftoverVcf Writing out sorted records to final VCF.
    [Fri Sep 14 12:29:50 EDT 2018] picard.vcf.LiftoverVcf done. Elapsed time: 1.71 minutes.
    Runtime.totalMemory()=8853651456

 
Liftover was around 99% successful. 209886 out of 211093 variants were lifted over to the new assembly.


But the chromosome codes have to be changed from the NKLS02 to NCBI codes...

    module load bcftools
    bcftools annotate --rename-chrs /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/ARS-UCD_NKLS_to_ncbi.tab -O z -o selected.ARS-UCDv1.2.NCBI.vcf.gz selected.ARS-UCD.v1.2.vcf
    bcftools annotate -I 'bakshy\_%CHROM\_%POS\_%REF\_%ALT' -O z -o selected.ARS-UCDv1.2.NCBI_ID.vcf.gz selected.ARS-UCDv1.2.NCBI.vcf.gz
    
There are some contig ID lines in the vcf file which do not belong to the NCBI chr codes. So removing those lines
    `gunzip -c selected.ARS-UCDv1.2.NCBI_ID.vcf.gz | grep -v 'NKLS02' | bgzip > selected.ARS-UCDv1.2.NCBI_ID_trimmed.vcf.gz`

But this step does not affect the no. of variants in the vcf file.

So better remove redundant files...

`rm selected.ARS-UCDv1.2.NCBI_ID.vcf.gz`

Our goal is around 80-120k final variants, so we should prune the dataset further based on LD.

This time I prunned the dataset based on Variance Inflation Factor i.e. --indep which prunes based on the variance inflation factor (VIF), which recursively removes SNPs within a sliding window.

The following are the lines from the Plink/1.7 [website](http://zzz.bwh.harvard.edu/plink/summary.shtml#prune) 

- The parameters for --indep are: window size in SNPs (e.g. 50), the number of SNPs to shift the window at each step (e.g. 5), the VIF threshold. 
- The VIF is 1/(1-R^2) where R^2 is the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously. That is, this considers the correlations between SNPs but also between linear combinations of SNPs.
- A VIF of 10 is often taken to represent near collinearity problems in standard multiple regression analyses (i.e. implies R^2 of 0.9). 
- A VIF of 1 would imply that the SNP is completely independent of all other SNPs. Practically, values between 1.5 and 2 should probably be used; particularly in small samples, if this threshold is too low and/or the window size is too large, too many SNPs may be removed.

    
     plink --vcf selected.ARS-UCDv1.2.NCBI_ID_trimmed.vcf.gz --cow --allow-extra-chr --double-id --indep 5 kb 5 1 --out ld.filter3 
    
    PLINK v1.90b4.4 64-bit (21 May 2017)
    Options in effect:
      --allow-extra-chr
      --cow
      --double-id
      --indep 5 kb 5 1
      --out ld.filter3
      --vcf selected.ARS-UCDv1.2.NCBI_ID_trimmed.vcf.gz
    
    Hostname: assembler2.agil.barc.ba.ars.usda.gov
    Working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3
    Start time: Fri Sep 14 19:10:49 2018
    
    Random number seed: 1536966649
    515987 MB RAM detected; reserving 257993 MB for main workspace.
    --vcf: ld.filter3-temporary.bed + ld.filter3-temporary.bim +
    ld.filter3-temporary.fam written.
    209886 variants loaded from .bim file.
    172 cattle (0 males, 0 females, 172 ambiguous) loaded from .fam.
    Ambiguous sex IDs written to ld.filter3.nosex .
    Using 1 thread (no multithreaded calculations invoked).
    Before main variant filters, 172 founders and 0 nonfounders present.
    Calculating allele frequencies... done.
    Total genotyping rate is 0.999446.
    209886 variants and 172 cattle pass filters and QC.
    Note: No phenotypes present.
    Pruned 5958 variants from chromosome 34, leaving 7538.
    Pruned 5332 variants from chromosome 35, leaving 6523.
    Pruned 3738 variants from chromosome 36, leaving 4415.
    Pruned 1162 variants from chromosome 37, leaving 2712.
    Pruned 2237 variants from chromosome 38, leaving 4894.
    Pruned 1787 variants from chromosome 39, leaving 3014.
    Pruned 4131 variants from chromosome 40, leaving 6786.
    Pruned 1206 variants from chromosome 41, leaving 2996.
    Pruned 326 variants from chromosome 42, leaving 1968.
    Pruned 2632 variants from chromosome 43, leaving 4772.
    Pruned 5394 variants from chromosome 44, leaving 7478.
    Pruned 328 variants from chromosome 45, leaving 1712.
    Pruned 2555 variants from chromosome 46, leaving 5474.
    Pruned 2736 variants from chromosome 47, leaving 5135.
    Pruned 198 variants from chromosome 48, leaving 870.
    Pruned 2322 variants from chromosome 49, leaving 2970.
    Pruned 1912 variants from chromosome 50, leaving 3538.
    Pruned 315 variants from chromosome 51, leaving 1337.
    Pruned 1394 variants from chromosome 52, leaving 3090.
    Pruned 4811 variants from chromosome 53, leaving 4646.
    Pruned 4646 variants from chromosome 54, leaving 6285.
    Pruned 252 variants from chromosome 55, leaving 1143.
    Pruned 51 variants from chromosome 56, leaving 435.
    Pruned 1846 variants from chromosome 57, leaving 2617.
    Pruned 18 variants from chromosome 58, leaving 309.
    Pruned 3133 variants from chromosome 59, leaving 2848.
    Pruned 76 variants from chromosome 60, leaving 711.
    Pruned 86 variants from chromosome 61, leaving 511.
    Pruned 499 variants from chromosome 62, leaving 2265.
    Pruned 43062 variants from chromosome 63, leaving 6736.
    Pruned 0 variants from chromosome 64, leaving 1.
    Pruned 0 variants from chromosome 65, leaving 2.
    Pruned 0 variants from chromosome 66, leaving 1.
    Pruned 0 variants from chromosome 67, leaving 1.
    Pruned 0 variants from chromosome 68, leaving 1.
    Pruned 1 variant from chromosome 69, leaving 8.
    Pruning complete.  104144 of 209886 variants removed.
    Marker lists written to ld.filter3.prune.in and ld.filter3.prune.out .
    
    End time: Fri Sep 14 19:10:55 2018

Now create the final vcf file by removing the variants which are in the prune.out file

    gunzip -c selected.ARS-UCDv1.2.NCBI_ID_trimmed.vcf.gz | grep -v -f <(cat ld.filter3.prune.out) > final.list1.vcf
    module load samtools htslib
    bgzip final.list1.vcf
    tabix -C final.list1.vcf.gz

Now process the IGC variants to merge with this final dataset...

### IGC variants filtering and pre-processing ###

working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs

First concatenate all the IGC haplotype coordinates to get a single file

    bcftools concat -a CH240.vcf.gz Domino.vcf.gz TPI4222.vcf.gz HF.vcf.gz LIB14427.vcf.gz LIB14413.vcf.gz -O z -o igc.vcf.gz

|Haplotype                           | Count|
|:--------------------------------|-----:|
|CH240_370M3_LILR_LRC             |  1500|
|CH240_391K10_KIR                 |   695|
|Domino_MHCclassI_gene2-5hap1_MHC |  2206|
|HF_LRC_hap1_KIR_LRC              |  4325|
|LIB14413_LRC                     |  2985|
|LIB14427_MHC                     |  5924|
|TPI4222_A14_MHCclassI_MHC        | 16299|


Now start filtering for high quality SNPs... 


    bcftools view -m2 -M2 -v snps igc.vcf.gz | bcftools view -q 0.25 | bcftools view -i 'QUAL == 999 && AN > 258' | bgzip > igc.filtered2.vcf.gz

- -m2 -M2 -v snps to only view biallelic SNPs
- -q 0.25 minimum allele frequency
- -i 'QUAL == 999 && AN > 258' include only sites with QUAL == 999 and total allele number > 258

This is how the data looks like now:

|Haplotype                          | Count|
|:--------------------------------|-----:|
|CH240_370M3_LILR_LRC             |   180|
|CH240_391K10_KIR                 |    66|
|Domino_MHCclassI_gene2-5hap1_MHC |   330|
|HF_LRC_hap1_KIR_LRC              |   664|
|LIB14413_LRC                     |   209|
|LIB14427_MHC                     |   939|
|TPI4222_A14_MHCclassI_MHC        |  2617|


I have also generated some stats for comparison in the following folders: igc.filtered and igc


####  Now liftover these haplotype coordinates to the ARS-UCDv1.2 assembly ####

    module load java
    module load picard
    
    java -jar $PICARD CreateSequenceDictionary R=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD.v1.2_NCBI_IDs.fa O=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD.v1.2_NCBI_IDs.fa.dict
    
    # check for the liftover filtered IGC VCF file
    [kiranmayee.bakshy@assembler2 condensed_vcfs]$ java -jar $PICARD LiftoverVcf I=igc.filtered2.vcf.gz O=liftover_to_v1.2/igc.vcf CHAIN=/mnt/nfs/nfs2/bickhart-users/cattle_asms/liftovers/pirbright_to_ARS-UCDv1.2/pirbright_to_ARS-UCDv1.2.liftover.chain REJECT=liftover_to_v1.2/igc.r.vcf R=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD.v1.2_NCBI_IDs.fa WRITE_ORIGINAL_POSITION=true


- VCF Liftover was 16.3% successful. Out of 5005 variants only 816 have lifted over to the new assembly.

- So the remaining rejected variants have to be arbitrarily placed into the ARS-UCDv1.2 system based on their relative positions and then they have to checked for potential overlap with the actual variants in the main variants file (final.set1.vcf.gz)

- Placing these markers was done in an excel spreadsheet using simple formulae and double checking for overlap with the actual variants coordinates. The spreadsheet is on my user server at U:\IGC_Project\ars-ucdv1.2\Arbitrary_ARS-UCDv1.2_coords.xlsx


**Now merge all the three vcf files**

1. High quality filtered variants 
2. IGC lifted variants and
3. IGC liftover rejected arbitrary variants.

This step can be easily done using bcftools but for this step to work the samples in the 3 VCF files should be in the same order.

Before that variants in all the 3 files should be given their SNP IDs separately and then the NCBI chr coordinates have to be converted to Chr numbers.

- SNP IDs for non-IGC variants: bakshy_chr_pos_ref_alt
- SNP IDs for lifted over IGC variants: ARS_PIRBRIGHT_chr_pos_ref_alt
- SNP IDs for liftover rejected variants which were given arbitrary ARS-UCD.v1.2 coordinates: ARS_PIRBRIGHT_Haplotype_pos_ref_alt

working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2

    bcftools view -S /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/sample_list_final -O z -o igc_filtered_rejected_arb.sampled.vcf.gz igc_filtered_rejected_arb.vcf
    
    bcftools annotate --rename-chr ARS-UCD_ncbi_to_chrnum.tab -O z -o igc.filtered.rejected.arb.sampled.chrnum.vcf.gz igc.filtered.rejected.arb.sampled.vcf.gz 
    
    bcftools view -S /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/sample_list_final -O z -o igc.id.sampled.vcf.gz igc.id.vcf.gz
    
    bcftools annotate --rename-chr ARS-UCD_ncbi_to_chrnum.tab -O z -o igc.id.sampled.chrnum.vcf.gz igc.id.sampled.vcf.gz
    
    bcftools annotate -I 'ARS\_PIRBRIGHT\_%CHROM\_%POS\_%REF\_%ALT' -O z -o igc.id.sampled.chrnum.ID.vcf.gz igc.id.sampled.chrnum.vcf.gz
    
    bcftools annotate --rename-chr ARS-UCD_ncbi_to_chrnum.tab -O z -o  /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/final.list1.chrnum.vcf.gz /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/final.list1.vcf.gz 
    
    bcftools annotate -I 'bakshy\_%CHROM\_%POS\_%REF\_%ALT' -O z -o /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/final.list1.chrnum.id.vcf.gz /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/final.list1.chrnum.vcf.gz

    bcftools concat -a igc.id.sampled.chrnum.ID.vcf.gz igc.filtered.rejected.arb.sampled.chrnum.vcf.gz /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/final.list1.chrnum.id.vcf.gz | bgzip > final_set1.vcf.gz

    
  The final dataset consists of a total of 110773 variants. This includes all the 5005 filtered IGC variants. The following is the per chromosome count: 


| Chr | Count |
|-----|-------|
| 1   | 7538  |
| 2   | 6548  |
| 3   | 4438  |
| 4   | 2712  |
| 5   | 4896  |
| 6   | 3014  |
| 7   | 6786  |
| 8   | 2996  |
| 9   | 1968  |
| 10  | 4774  |
| 11  | 7480  |
| 12  | 1713  |
| 13  | 5474  |
| 14  | 5136  |
| 15  | 870   |
| 16  | 2970  |
| 17  | 3538  |
| 18  | 2450  |
| 19  | 3090  |
| 20  | 4646  |
| 21  | 6285  |
| 22  | 1151  |
| 23  | 4239  |
| 24  | 2617  |
| 25  | 321   |
| 26  | 2860  |
| 27  | 711   |
| 28  | 511   |
| 29  | 2265  |
| X   | 6736  |


Now converting the final variants vcf file to AIPL format is simple and can be done using Plink and R.

Here is the trick: 
Plink2 --export A-transpose counts for the number of Ref alleles and codes as 0 = HomoAlt, 1 = Hetero, 2 = HomoRef and NA = missing.
AIPL format counts for number of ALT alleles. So we just need to exchange 0 and 2; recode NA as 5.

Here are my bash and R scripts for the format conversion (vcfToAipl.sh):

    #!/bin/bash
     
    #SBATCH -p assemble2
    #SBATCH --mem 50000 
    #SBATCH --nodes=1 
    #SBATCH --ntasks-per-node=1
    #SBATCH --cpus-per-task=1
    
	# A simple bash script which takes a vcf file as an input and outputs ped and map files in AIPL format
    
    module load plink/2.00alM-2017-05-22
    
    vcffile=$1
    outfile=$2
    
    plink2 --vcf $vcffile --const-fid --cow --export A-transpose --out $outfile
    
    Rscript --vanilla vcfToAipl.R $outfile.traw $outfile
    
    wait


The Rscript is as follows: 

    #!/usr/bin/env Rscript
    #Date:Sep 21 2018
    #Author:Kiranmayee Bakshy
    
    # A program to convert the plink genotype recodings (A-transpose) output to AIPL format ped and map files
    
    args = commandArgs(trailingOnly=TRUE)
    
    # test if there is at least one argument: if not, return an error
    if (length(args)==0) {
      stop("Usage: vcfToAipl.R plink.traw output filename.\n.At least one argument must be supplied (plink.traw).\n", call.=FALSE)
    } else if (length(args)==1) {
      # default output file
      args[2] = "out.aipl"
    }
    
    
    dat <- read.delim(args[1], header=T, stringsAsFactors=FALSE)
    dat.map<-dat[,c(1:6)]
    names(dat.map)<-c("CHR","SNP","Dist","POS","REF","ALT")
    write.table(x = dat.map[,c(1:4)], file = paste(args[2],"aipl.map", sep="."), sep=" ", row.names = F, col.names = F, quote=F)
    samples<-names(dat[,c(7:ncol(dat))])
    samples.list<-gsub("X0_","",samples)
    names(dat[,c(7:ncol(dat))])<-samples.list
    dat.ped<-t(dat[,c(7:ncol(dat))])
    dat1.ped<-as.data.frame(cbind(samples.list,dat.ped),stringsAsFactors = F)
    row.names(dat1.ped)<-samples.list
    dat.ped<-dat1.ped
    dat.ped[2:ncol(dat.ped)]<-lapply(dat1.ped[2:ncol(dat1.ped)], function(x) ifelse(x == 0, 2, ifelse(x == 2, 0, x)))
    dat.ped[2:ncol(dat.ped)]<-lapply(dat.ped[2:ncol(dat.ped)], function(x) ifelse(x == "NA", 5, x))
    write.table(dat.ped, file = paste(args[2],"aipl.ped", sep="."), sep=" ",col.names=F, quote=F)



**Next step is resolving the ChrX coordinates into Chr30 (PAR=Pseudo Autosomal Regions) and Chr31 (ChrX).**

The PAR region for ARS-UCDv1.2 is approximately at X:1-5700000.

There is only one SNP in this region, so I manually changed the chromosome coordinates to 30.


AIPL uses fixed width format files. 

Example format for reference: 

	AnKey Gender      NumSnps     Genotype  (this header is not part of the file)
	1182835 1       641459  002224130020…
	1496541 1       641459  002000100020…
	1857310 1       641459  002222200020…
	1979779 17      641459  002111111021…
	2117682 1       641459  002111111021…
	2160886 1       641459  002221111021…
	2214378 1       641459  111110111111…
	2409276 17      641459  002222200020…
	2446336 1       641459  002111111021…
	2513929 2       641459  002111100020…


The “NumSNPs” is the total number of present or missing SNPs in the genotype array. It is a the sum total of all characters in the genotype string.

I can do this easily in R:
	
	
	library(tidyr)
	library(gdata)
	setwd("U:\\IGC_Project\\ars-ucdv1.2\\ARS-UCDv1.2_aipl_formatting")
	# read in the current format
	dat<-read.delim("final_set.aipl.ped", sep=" ", header=F, stringsAsFactors=F)
	#drop second column to remove duplicate column of animal id
	dat$V2<-NULL
	# concatenate all the columns except the first without space
	aipl<-unite(aipl, genotype, sep="", -1)
	# give names to the existing columns for ref
	names(aipl2)<-c("AnKey", "genotype")
	# add number of snps column
	aipl2$NumSnps<-110733
	# Add gender column, in this case all are bulls might be 1 or 2
	aipl2$Gender<-1
	# rearrange columns in required order
	aipl3<-aipl2[,c("AnKey","Gender","NumSnps","genotype")]
	# check the structure of the dataframe
	str(aipl3)
	# write as a fixed width file
	write.fwf(aipl3, "reformatted.aipl.ped", quote=F, sep=" ", rownames=F, colnames=F)
	
I am saving the R work space to quickly convert the gender coding in case my guess is wrong.

	save.image("ped_to_aipl.RData")
	
The final data set in AIPL format are at the following location in my [Google drive](https://drive.google.com/drive/folders/1S2zcb8ZTxw638Wod_I6okaO9tlHln-aq?usp=sharing "Google drive")

	