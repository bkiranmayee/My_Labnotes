# Data preparation for GATK VQSR Run3 #

##  ##
### Processing the dbSNP markers for use in the GATK VQSR analysis

The dbsnp dataset given by Derek was a snapshot from the previous analysis done in 2013

It is a UMD3.1.1 filtered version of dbSNP  


working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/resources

	# Liftover of filtered dbsnp markers from umd3.1.1 to ARS-UCDv1.2
    java -jar $PICARD LiftoverVcf I=dbsnp_filtered_umd3_snps.vcf O=dbsnp_filtered_ARS-UCDv1.2_snps.vcf CHAIN=/mnt/nfs/nfs2/bickhart-users/cattle_asms/liftovers/UMD3.11_to_ARS-UCD1.2/umd3_kary_unmask_ngap_to_ARS-UCD1.2.mmap.liftover.chain REJECT=rejected.vcf R=/mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARS-UCD1.2.PlusY.fa
	# Rename chromosome names from the NKLS02 to numeric
    bcftools annotate --rename-chrs ARS-UCDv1.2_NKLS_to_chrnum dbsnp_filtered_ARS-UCDv1.2_snps.vcf -O z -o dbsnp_filtered_ARS-UCDv1.2_snps_chrnum.vcf.gz
    # remove extra NKLS02 IDs of all the leftover contigs
    gunzip -c dbsnp_filtered_ARS-UCDv1.2_snps_chrnum.vcf.gz | grep -v 'NKLS02' > dbsnp_filtered_ARS-UCDv1.2_snps.vcf
    module load samtools htslib
    bgzip dbsnp_filtered_ARS-UCDv1.2_snps.vcf
    tabix dbsnp_filtered_ARS-UCDv1.2_snps.vcf.gz
	# Now select only SNPs whose Ref allele matches with the new reference genome assembly
    perl select_snps.pl dbsnp_filtered_ARS-UCDv1.2_snps.vcf.gz dbsnp_filtered_ARS-UCDv1.2_selected_snps.vcf
    # Remove extra information from the ID field
	sed -i 's/\;Reference_seq\=[ACGT]//' dbsnp_filtered_ARS-UCDv1.2_selected_1_snps.vcf
	# Remove duplicate entries
	/mnt/nfs/nfs2/bickhart-users/binaries/vcflib/bin/vcfuniq dbsnp_filtered_ARS-UCDv1.2_selected_snps.vcf.gz >  dbsnp_filtered_ARS-UCDv1.2_selected_uniq_snps.vcf
    bgzip dbsnp_filtered_ARS-UCDv1.2_selected_uniq_snps.vcf
    tabix dbsnp_filtered_ARS-UCDv1.2_selected_uniq_snps.vcf.gz
	# Print selected dbsnp IDs as chrom.pos.ref.alt to find overlap with other datasets 
	gunzip -c dbsnp_filtered_ARS-UCDv1.2_selected_uniq_snps.vcf.gz | grep -v '#' | awk '{print $1"."$2"."$4"."$5}' > dbsnp_selected_ids
	gunzip -c dbsnp_filtered_ARS-UCDv1.2_selected_uniq_snps.vcf.gz | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -o v_per_chr_m -e '#'
	sort -k1,1 -V -s v_per_chr_m > sorted
	


**CHROM**|**COUNT**
:-----:|:-----:
1|603428
2|532677
3|490337
4|508781
5|524150
6|461740
7|433389
8|404725
9|408182
10|421991
11|424956
12|353358
13|341078
14|347560
15|359097
16|309624
17|320578
18|265525
19|280289
20|297531
21|262898
22|258018
23|278132
24|278039
25|192686
26|205396
27|196744
28|200462
29|239433
X|247955


Now do the overlap analysis with other datasets such as 1KBulls, bovineHD and the Holstein specific markers set
 
    [kiranmayee.bakshy@assembler2 ARS-UCDv1.2]$ wc -l /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/resources/dbsnp_selected_ids filtered_ids combined_ids bovineHD_ids 1kbulls_ids
     10451347 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/resources/dbsnp_selected_ids
      1303922 filtered_ids
     17883171 combined_ids
       771291 bovineHD_ids
     22761723 1kbulls_ids
     53171753 total
    

    [kiranmayee.bakshy@assembler2 ARS-UCDv1.2]$ perl /mnt/nfs/nfs2/bickhart-users/binaries/perl_toolchain/bed_cnv_fig_table_pipeline/nameListVennCount.pl -o filtered_ids combined_ids bovineHD_ids 1kbulls_ids /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/resources/dbsnp_selected_ids
    File Number 1: filtered_ids
    File Number 2: combined_ids
    File Number 3: bovineHD_ids
    File Number 4: 1kbulls_ids
    File Number 5: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/resources/dbsnp_selected_ids
    Set Count
    1;2 248793
    1;2;3   3500
    1;2;3;4 2006
    1;2;3;4;5   15794
    1;2;3;5 892
    1;2;4   402375
    1;2;4;5 594969
    1;2;5   35593
    2   5902346
    2;3 19733
    2;3;4   11158
    2;3;4;5 93017
    2;3;5   6064
    2;4 5089850
    2;4;5   4977725
    2;5 479356
    3   545334
    3;4 41525
    3;4;5   12872
    3;5 19350
    4   10485070
    4;5 1032271
    5   3183444
    
    Group: 5 in output file: group_5.txt
    Group: 1;2;3 in output file: group_1_2_3.txt
    Group: 2;3;4 in output file: group_2_3_4.txt
    Group: 3;5 in output file: group_3_5.txt
    Group: 1;2 in output file: group_1_2.txt
    Group: 1;2;3;5 in output file: group_1_2_3_5.txt
    Group: 1;2;4;5 in output file: group_1_2_4_5.txt
    Group: 2;3;5 in output file: group_2_3_5.txt
    Group: 2;3 in output file: group_2_3.txt
    Group: 2 in output file: group_2.txt
    Group: 3;4;5 in output file: group_3_4_5.txt
    Group: 4 in output file: group_4.txt
    Group: 1;2;3;4;5 in output file: group_1_2_3_4_5.txt
    Group: 1;2;5 in output file: group_1_2_5.txt
    Group: 2;5 in output file: group_2_5.txt
    Group: 1;2;3;4 in output file: group_1_2_3_4.txt
    Group: 2;3;4;5 in output file: group_2_3_4_5.txt
    Group: 2;4;5 in output file: group_2_4_5.txt
    Group: 4;5 in output file: group_4_5.txt
    Group: 3 in output file: group_3.txt
    Group: 3;4 in output file: group_3_4.txt
    Group: 1;2;4 in output file: group_1_2_4.txt
    Group: 2;4 in output file: group_2_4.txt



- Uniq to dbSNP : group_5.txt **~3M**
- NonHolstein-1kbull: group_2_4.txt **~5M** => needs further filtering on QUAL
- Highly-filtered-uniq-ARS: group_1_2.txt **~200K**
- 1kbull: group_1_2_4.txt **~400K** => needs further filtering on QUAL


    	sed -i 's/\./\t/g' group_1_2.txt
    	sed -i 's/\./\t/g' group_5.txt
   	 	sed -i 's/\./\t/g' group_2_4.txt
    	sed -i 's/\./\t/g' group_1_2_4.txt
    	sort -k 1,1 -k2,2n -V -s group_5.txt > bed/group_5.txt
    	sort -k 1,1 -k2,2n -V -s group_1_2.txt > bed/group_1_2.txt
    	sort -k 1,1 -k2,2n -V -s group_2_4.txt > bed/group_2_4.txt
    	sort -k 1,1 -k2,2n -V -s group_1_2_4.txt > bed/group_1_2_4.txt
    	wc -l bed/group_1_2.txt bed/group_2_4.txt bed/group_5.txt bed/group_1_2_4.txt
	       248793 bed/group_1_2.txt
	      5089850 bed/group_2_4.txt
	      3183444 bed/group_5.txt
	       402375 bed/group_1_2_4.txt
	      8924462 total
    

Now subsetting......

	module load bcftools samtools htslib
	# Fill in MAF and HWE tags and tabix using bcftools and tabix
	bcftools plugin fill-tags /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/f2.vcf.gz -O z -o f2.tagged.vcf.gz -- -t MAF,HWE
	tabix -p vcf f2.tagged.vcf.gz
	#Subset high quality 1000 bulls run 5 data set	
	bcftools view -i 'QUAL>=999' -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/venn/ARS-UCDv1.2/bed/group_2_4.txt -O z -o 1kbull_sites_q999.vcf.gz f2.tagged.vcf.gz
	
	#Subset high quality Holstein specific 1000 bulls dataset
	bcftools view -i 'QUAL>=999' -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/venn/ARS-UCDv1.2/bed/group_1_2_4.txt -O z -o 1kbull_Hol_sites_q999.vcf.gz f2.tagged.vcf.gz
	#Subset high quality Holstein specific ARS sites
	bcftools view -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/venn/ARS-UCDv1.2/bed/group_1_2.txt -O z -o Hol_filtered.vcf.gz f2.tagged.vcf.gz
	# Subset dbSNP sites
	bcftools view -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/venn/ARS-UCDv1.2/bed/group_5.txt -O z -o dbSNP.vcf.gz /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/resources/dbsnp_filtered_ARS-UCDv1.2_selected_uniq_snps.vcf.gz
	# Subset bovineHD sites
	bcftools view -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/venn/ARS-UCDv1.2/bovineHDmapped -O z -o bovineHD.vcf.gz f2.tagged.vcf.gz
	
	tabix -p vcf bovineHD.vcf.gz
	tabix -p vcf 1kbull_sites_q999.vcf.gz
	tabix -p vcf 1kbull_Hol_sites_q999.vcf.gz
	tabix -p vcf Hol_filtered.vcf.gz
	tabix -p vcf dbSNP.vcf.gz
	
	gunzip -c 1kbull_sites_q999.vcf.gz | grep -v '#' | wc -l
	3583444
	gunzip -c 1kbull_Hol_sites_q999.vcf.gz | grep -v '#' | wc -l
	332403
	gunzip -c Hol_filtered.vcf.gz | grep -v '#' | wc -l
	248793
	gunzip -c dbSNP.vcf.gz | grep -v '#' | wc -l
	3180858
	gunzip -c bovineHD.vcf.gz | grep -v '#' | wc -l
	584638
	
	gunzip -c 1kbull_sites_q999.vcf.gz | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -o 1kbull_v_per_chr -e '#'
	sort -k1,1 -V -s 1kbull_v_per_chr > 1kbull_sorted.v_per_chr
	
	gunzip -c 1kbull_Hol_sites_q999.vcf.gz | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -o 1kbull_hol_v_per_chr -e '#'
	sort -k1,1 -V -s 1kbull_hol_v_per_chr > 1kbull_hol_sorted.v_per_chr
	
	gunzip -c Hol_filtered.vcf.gz | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -o filtered_v_per_chr -e '#'
	sort -k1,1 -V -s filtered_v_per_chr > filtered_sorted.v_per_chr
	
	gunzip -c dbSNP.vcf.gz | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -o dbsnp_v_per_chr -e '#'
	sort -k1,1 -V -s dbsnp_v_per_chr > dbsnp_sorted.v_per_chr

	gunzip -c bovineHD.vcf.gz | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -o bovineHD_v_per_chr -e '#'
	sort -k1,1 -V -s bovineHD_v_per_chr > bovineHD_sorted.v_per_chr


The bovineHD sites that I selected is based on just the lifted over coordinates bed file...

So, I subsetted the bovineHD sites that match both the ref and alt allele as follows:
 
	bcftools view -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/venn/ARS-UCDv1.2/overlapped_bovineHD_bed.txt -O z -o bovineHD2.vcf.gz f2.tagged.vcf.gz 

Link to the VQSR trial [runs](https://github.com/bkiranmayee/My_Labnotes/blob/master/CDDR/initial_vqsr_trials.md)