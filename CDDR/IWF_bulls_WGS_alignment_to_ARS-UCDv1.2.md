# WGS Realignment to the new assembly ARS.UCD-v1.2 #

This documents all the steps and files generated during realignment of WGS bams (aligned to unpublished version ARS-UCDv14) of the 172 CDDR bulls selected using IWF.

 
The aligned bams (to ARS-UCDv14) are on AGIL cluster (/mnt/nfs/nfs1/derek.bickhart/CDDR-Project/??????/*/*.sorted.merged.bam) and were transferred to the CERES cluster by Steve to the following folder: 

*/beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR/iwf_prioritized_bulls_arsucdv14_bams*


	# Check the md5sum to ensure that the transferred files are intact
	screen -S checksum
	md5sum ??????/*sorted.merged.bam > checksum_bam_ceres
	
The bams are intact, now start alignment to the new reference genome assembly using BWA

	ls $PWD/iwf_prioritized_bulls_arsucdv14_bams/aligns/*bam > iwf_prioritized_bulls_arsucdv14_bams/aligns.bulls.tab.
	
I manually edited the tab file created above using vim to include the sample and library names

	perl /beegfs/project/rumen_longread_metagenome_assembly/binaries/perl_toolchain/sequence_data_pipeline/alignBamReadsToNewAssemSlurm.pl -b iwf_prioritized_bulls_arsucdv14_bams/aligns_v1.2 -t iwf_prioritized_bulls_arsucdv14_bams/aligns.bulls.tab -f ARSUCD1.2.current_ref.fa -m true
	
	ls $PWD/round2/*bam > round2.bulls.tab
	
	perl /beegfs/project/rumen_longread_metagenome_assembly/binaries/perl_toolchain/sequence_data_pipeline/alignBamReadsToNewAssemSlurm.pl -b iwf_prioritized_bulls_arsucdv14_bams/round2_v1.2 -t iwf_prioritized_bulls_arsucdv14_bams/round2.bulls.tab -f ARSUCD1.2.current_ref.fa -m true
	
	
I had to cancel all the jobs and did the fresh start because it takes more than 2 days. After submitting I updated all the alignment jobs to medium partition and let the samMerger jobs to stay on short partition.

But as soon as I submitted the jobs some of them already started running on the short partition. I was not able to update the partition using the scontrol command for running jobs.

I let them run and timeout after 2 days and resubmitted these scripts manually on medium partition.

I had to resubmit around 21 (out of 26) scripts for the bams in aligns folder and 11 (out of 35) scripts for the bams in round2 folder.

I am waiting for these jobs to complete so that I can start 3 variant calling pipelines: 

- SAMtools
- GATK-HC
- Freebayes

Meanwhile I can run BQSR on one of the bams after indel realignment and marking duplicates.

So submit a job for one of the sample and check if the BQSR really makes a huge difference.


Markduplicates and INDEL realign for one of the bams:

	sbatch indel_realign_markdups_kb.sh iwf_prioritized_bulls_arsucdv14_bams/aligns_v1.2/014HO06441-6441_qia_test/014HO06441-6441_qia_test.sorted.merged.bam ARSUCD1.2.current_ref.fa

Running GATK BQSR for the above deduplicated and indel realigned bam...
	
	module load gatk/3.6

	java -jar /software/7/apps/gatk/3.6/GenomeAnalysisTK.jar -T BaseRecalibrator -R ARSUCD1.2.current_ref.fa -I iwf_prioritized_bulls_arsucdv14_bams/aligns_v1.2/014HO06441-6441_qia_test/014HO06441-6441_qia_test.dedup.realn.bam -o bqsr/recal_data.table

### Variant calling using SAMtools ###


### Variant calling using GATK HaplotypeCaller ###