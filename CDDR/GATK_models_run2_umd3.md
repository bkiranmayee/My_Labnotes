# GATK VQSR training models using CDDR gold std variants set using UMD3.1.1 filtered coordinates#

## Pre-processing data sets for model training and parameter optimization ##

### Training truth set ###

The training truth set consists of a subset of (chr3-chr29) gold standard highly filtered homozygous variants in UMD3.1.1 version coordinate system

working Directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/gatk_datasets/truthset
	
	bcftools view -r chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29 -O z -o truthset.vcf.gz /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/selected.snps.vcf.gz
	# collect these IDs and label them as true sites in the training_falseset
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID:true\n' truthset/truthset.vcf.gz | bgzip -c > truthset/truthset.tab.gz
 	tabix -s1 -b2 -e2 truthset/truthset.tab.gz


### Training false set ###

This data set consists of high quality variants which includes false positives (heterozygous calls in non-target animals)

Its a merge of two datasets: initial_truthset.vcf.gz and initial_falseset.vcf.gz which were generated using the progressiveSelection.pl scripts using the specific filters

initial_truthset.vcf.gz: 

initial_falseset.vcf.gz:  


working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/gatk_datasets/falseset
 
     plink --vcf falseset.vcf.gz --cow --double-id --keep-autoconv --maf 0.05 --hwe 1e-05 midp --out falseset.maf.hwe --write-snplist
     gunzip -c falseset.vcf.gz | grep '#' > falseset1.vcf
     gunzip -c falseset.vcf.gz | grep -Fwf falseset.maf.hwe.snplist >> falseset1.vcf
     bgzip falseset1.vcf
     tabix -p vcf falseset1.vcf.gz
	# validation set
     bcftools view -r chr1 -O z -o val.set.vcf.gz falseset1.vcf.gz 
	# test set
     bcftools view -r chr2 -O z -o test.set.vcf.gz falseset1.vcf.gz 
	# training set     
	bcftools view -r chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29 -O z -o trainingset.false.vcf.gz falseset1.vcf.gz


Labelling the true sites in the falseset

	bcftools annotate -a truthset/truthset.tab.gz -m +ID -c CHROM,POS,-,-,ID -O z -o falseset/trainingset.false.tagged.vcf.gz falseset/trainingset.false.vcf.gz

	tabix -p vcf falseset/trainingset.false.tagged.vcf.gz
    
    gunzip -c falseset/trainingset.false.tagged.vcf.gz | grep -v '#' | grep -c ':true'
	56257

Only 56257 sites are true in this set out of 139049. These must have filtered out in the MAF threshold. 

Let's see the overlap of the original falseset.vcf.gz from which I have filtered out using plink hwe and maf filters.
	
	bcftools isec -p dir falseset/falseset.vcf.gz /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/selected.snps.vcf.gz

This has a full overlap of all the gold std variant sites. How about the falseset1.vcf.gz which was filtered form the falseset.vcf.gz using the maf and hwe filters?

	bcftools isec -p dir falseset/falseset1.vcf.gz /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/selected.snps.vcf.gz

OK, only 64450 overlap among 218714 sites. The rest of the 154264 sites must have filtered out by MAF threshold because we used the HWE filter for the gold std but not the MAF!!! We used the total allele number (AN in the vcf info) and Alt allele count (AC in the vcf info) instead. 

Lets see how these data sets work as training sets for GATK VQSR recalibration.

## GATK VQSR parameter optimization ##

### Trial 1 ###

	GenomeAnalysisTK -T VariantRecalibrator \
					 -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/test_files/umd3_kary_unmask_ngap.fa \ 
					 -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/f2.vcf.gz \
					 -resource:HQF,known=false,training=true,truth=true,prior=10.0 truthset/truthset.vcf.gz  \
					 -resource:HQFHET,known=false,training=true,truth=false,prior=5.0 falseset/trainingset.false.tagged.vcf.gz \
					 -an DP -an MQSB -an BQB -an RPB -an MQB -an ICB
					 -mode SNP --maxGaussians 4 \
					 -recalFile Run2/trial_1_output.recal 
					 -tranchesFile Run2/trial_1_output.tranches 
					 -rscriptFile Run2/trial_1_output.plots.R



### Trial 2 ###




GenomeAnalysisTK -T VariantRecalibrator -R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/test_files/umd3_kary_unmask_ngap.fa -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/f2.vcf.gz -resource:HQF,known=false,training=true,truth=true,prior=10.0 truthset/truthset.vcf.gz -resource:HQFHET,known=false,training=true,truth=false,prior=5.0 falseset/trainingset.false.tagged.vcf.gz -an DP -an MQSB -an BQB -an RPB -an MQB -an ICB -mode SNP --maxGaussians 4 -recalFile Run2/trial_1_output.recal -tranchesFile Run2/trial_1_output.tranches -rscriptFile Run2/trial_1_output.plots.R
 



  
