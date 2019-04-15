# Preparation of test and training datasets for GATK VQSR using the new assembly coordinates

Date: Oct 31 2018

### Test and training of the highly filtered dataset ###

working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration

    mkdir gatk
    
    cd gatk
    
    # using vcflib/vcfrandomsample for randomly selecting variants for test set 
    /mnt/nfs/nfs2/bickhart-users/binaries/vcflib/bin/vcfrandomsample -r 0.5 -p 1 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/initial_truthset/initial_truthset_hwe_ld.vcf.gz > filtered.testset.vcf 
	bgzip filtered.testset.vcf
	tabix -p vcf filtered.testset.vcf.gz
    
    # now select the other half of the vcf as training dataset
    module load bcftools
    bcftools isec -p diff filtered.testset.vcf.gz /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/initial_truthset/initial_truthset_hwe_ld.vcf.gz
	
	mv diff/0001.vcf ./filtered.trainingset.vcf
	bgzip filtered.trainingset.vcf
	tabix -p vcf filtered.trainingset.vcf.gz
    
    # Comparing the number of variants per chromosome in test and training sets
    gunzip -c filtered.testset.vcf.gz | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -o filtered.testset.v_per_chr_m -e '#'
	sort -k1,1 -V -s filtered.testset.v_per_chr > filtered.testset._sorted.v_per_chr
	gunzip -c filtered.trainingset.vcf.gz | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -o filtered.trainingset.v_per_chr_m -e '#' 
	sort -k1,1 -V -s filtered.trainingset.v_per_chr > filtered.trainingset.sorted.v_per_chr



**CHROM**|**TEST**|**TRAIN**
:-----:|:-----:|:-----:
1|49515|49535
2|40920|40493
3|27642|27958
4|19079|19281
5|27109|27600
6|23317|23376
7|46941|47188
8|22083|21654
9|10445|10677
10|30748|30643
11|52801|52580
12|8174|8226
13|30388|30463
14|36678|36964
15|3552|3673
16|19200|19403
17|25085|24749
18|8244|8186
19|23661|23812
20|30363|30260
21|40306|40675
22|8570|8575
23|2450|2429
24|16731|16731
25|2609|2702
26|19526|19488
27|5315|5276
28|2604|2643
29|17145|17481

The number of variants per chromosome seems to be comparable.

### Now do similar processing for the test and training set from the 1kbulls variants present in Holstein with QUAL >= 999 ###

First subset the 1kbulls variants and filter by QUAL
	
	# get the SNP ids that overlap with 1kbulls run5 data
    sed -i 's/\./\t/g' /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/venn/ARS-UCDv1.2/group_2_4.txt
	# sort the bed file
    sort -k1,1 -V -s /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/venn/ARS-UCDv1.2/group_2_4.txt > /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/venn/ARS-UCDv1.2/group_2_4.bed
	# Now subset using bcftools
	# select only variants which overlap with 1kbulls run5 data with QUAL > = 999
    bcftools view -i 'QUAL>=999' -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/venn/ARS-UCDv1.2/group_2_4.bed.txt -O z -o 1kbull_sites_q999.vcf.gz /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/f2.vcf.gz

	tabix -p vcf 1kbull_sites_q999.vcf.gz
	# Exclude all non-informative sites (Alt allele count = total allele number)
	bcftools view -e 'AC==AN' -O z -o 1kbull_sites_Hol_q999.vcf.gz 1kbull_sites_q999.vcf.gz

	tabix -p vcf 1kbull_sites_Hol_q999.vcf.gz
	# Now fill in MAF AND HWE tags and filter based on these values
    bcftools plugin fill-tags 1kbull_sites_Hol_q999.vcf.gz -O z -o 1kbull_sites_Hol_q999_tagged.vcf.gz -- -t MAF,HWE
    tabix -p vcf 1kbull_sites_Hol_q999_tagged.vcf.gz
	# Select only variants with MAF >= 0.05 and HWE p-value > 1e-05
    bcftools view -i 'INFO/MAF>=0.05 & INFO/HWE>1e-05' -O z -o 1kbull_sites_Hol_q999_maf_hwe.vcf.gz 1kbull_sites_Hol_q999_tagged.vcf.gz
    tabix -p vcf 1kbull_sites_Hol_q999_maf_hwe.vcf.gz
    # Check the number of variants in each dataset 
    wc -l /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/venn/ARS-UCDv1.2/group_2_4.bed.txt
	10067575 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/venn/ARS-UCDv1.2/group_2_4.bed.txt
    gunzip -c 1kbull_sites_Hol_q999.vcf.gz | grep -v '#' | wc -l
	8162017
    gunzip -c 1kbull_sites_Hol_q999_maf_hwe.vcf.gz | grep -v '#' | wc -l
	6718071
    
    # using vcflib/vcfrandomsample for randomly selecting variants for test dataset
    /mnt/nfs/nfs2/bickhart-users/binaries/vcflib/bin/vcfrandomsample -r 0.3 -p 1  1kbull_sites_Hol_q999_maf_hwe.vcf.gz > 1kbull.testset.vcf 
    bgzip 1kbull.testset.vcf
    tabix -p vcf 1kbull.testset.vcf.gz 
    # using vcflib/vcfrandomsample for randomly selecting variants for training dataset
    /mnt/nfs/nfs2/bickhart-users/binaries/vcflib/bin/vcfrandomsample -r 0.4 -p 2  1kbull_sites_Hol_q999_maf_hwe.vcf.gz > 1kbull.trainingset.vcf 
    bgzip 1kbull.trainingset.vcf
    tabix -p vcf 1kbull.trainingset.vcf.gz
    
    bcftools isec -p diff/kbull 1kbull.testset.vcf.gz 1kbull.trainingset.vcf.gz
    
Looks like there are 807086 sites that overlapped in training and test sets

So I should be selecting the unique sites from both sets

    mv diff/kbull/0000.vcf ./kbull_testset_unique.vcf 
    mv diff/kbull/0001.vcf ./kbull_trainingset_unique.vcf
    bgzip kbull_testset_unique.vcf
    tabix -p vcf kbull_testset_unique.vcf.gz
    bgzip kbull_trainingset_unique.vcf
    tabix -p vcf kbull_trainingset_unique.vcf.gz


    # Comparing the number of variants per chromosome in test and training sets
    gunzip -c kbull_testset_unique.vcf.gz | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -o 1kbull.testset.unique_v_per_chr -e '#'
    sort -k1,1 -V -s 1kbull.testset.unique_v_per_chr > 1kbull.testset_unique_sorted.v_per_chr
    gunzip -c kbull_trainingset_unique.vcf.gz | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -o 1kbull.trainingset_unique.v_per_chr -e '#' 
    sort -k1,1 -V -s 1kbull.trainingset_unique.v_per_chr > 1kbull.trainingset_unique.sorted.v_per_chr

	

**CHROM**|**training**|**test**
:-----:|:-----:|:-----:
1|120556|77060
2|86393|55125
3|96163|61512
4|101546|65031
5|98933|63244
6|93044|59490
7|64841|42313
8|79144|51258
9|93263|60109
10|63295|40555
11|71343|46054
12|82974|53574
13|51171|32603
14|55636|35920
15|77892|49708
16|58607|37813
17|62976|40766
18|46532|30136
19|34142|21380
20|49114|31488
21|34258|22048
22|43876|28413
23|61840|39905
24|53861|35116
25|37361|24080
26|39878|25301
27|42131|26905
28|38544|24736
29|41602|26713
Total|1880916|1208356

All the chromosomes are sufficiently represented in these datasets.

Now run VQSR with these datasets...

I realized that this type of modelling and cross validation is not suitable at this moment...

Changing the approach...

 