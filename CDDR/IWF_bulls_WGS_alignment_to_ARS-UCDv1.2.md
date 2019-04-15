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


Markduplicates and INDEL realign:
---------------------------------

	sbatch indel_realign_markdups_kb.sh iwf_prioritized_bulls_arsucdv14_bams/aligns_v1.2/014HO06441-6441_qia_test/014HO06441-6441_qia_test.sorted.merged.bam ARSUCD1.2.current_ref.fa

	ls iwf_prioritized_bulls_arsucdv14_bams/aligns_v1.2/*/*sorted.merged.bam > new.bamlist
	ls iwf_prioritized_bulls_arsucdv14_bams/round2_v1.2/*/*sorted.merged.bam >> new.bamlist
	
	# some bam files were not created or I suspect got accidentally deleted (on March 22 I guess?) 
	007HO092426, 007HO10687, 007HO09221, 007HO7004, 001HO11389, 200HO03453

	# Queueing up the alignemnt of these 6 bams again after transferring the arsucdv14 bams again from AGIL to ceres
	perl /beegfs/project/rumen_longread_metagenome_assembly/binaries/perl_toolchain/sequence_data_pipeline/alignBamReadsToNewAssemSlurm.pl -b iwf_prioritized_bulls_arsucdv14_bams/aligns_v1.2 -t iwf_prioritized_bulls_arsucdv14_bams/aligns.repeat.tab -f ARSUCD1.2.current_ref.fa -p msn -m true
	 
	# Calculating some stats 
	perl /beegfs/project/rumen_longread_metagenome_assembly/binaries/perl_toolchain/sequence_data_scripts/getBamStats.pl -n new.bamlist -o new.bam.stats
	
	
	cat new.bamlist | xargs -I {} sbatch indel_realign_markdups_kb.sh {} ARSUCD1.2.current_ref.fa

	ls iwf_prioritized_bulls_arsucdv14_bams/aligns_v1.2/*/*dedup.realn.bam > new.dedup.realn.bamlist
	ls iwf_prioritized_bulls_arsucdv14_bams/round2_v1.2/*/*dedup.realn.bam >> new.dedup.realn.bamlist

	perl /beegfs/project/rumen_longread_metagenome_assembly/binaries/perl_toolchain/sequence_data_scripts/getBamStats.pl -n new.dedup.realn.bamlist -o new.dedup.realn.bam.stats
	
	# some dedup.realn.bam files in aligns_v1.2 folder are truncated for some reason
	011HO10331, 011HO11231, 014HO05680
	
	# Doing indel realignment for these 3 bams again
	 sbatch indel_realign.sh /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR/iwf_prioritized_bulls_arsucdv14_bams/aligns_v1.2/011HO10331/011HO10331.sorted.merged.bam.dedup.bam ARSUCD1.2.current_ref.fa /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR/iwf_prioritized_bulls_arsucdv14_bams/aligns_v1.2/011HO10331/011HO10331.sorted.merged.bam.intervals -partition msn

	 sbatch indel_realign.sh /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR/iwf_prioritized_bulls_arsucdv14_bams/aligns_v1.2/011HO11231/011HO11231.sorted.merged.bam.dedup.bam ARSUCD1.2.current_ref.fa /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR/iwf_prioritized_bulls_arsucdv14_bams/aligns_v1.2/011HO11231/011HO11231.sorted.merged.bam.intervals -partition msn

	 sbatch indel_realign.sh /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR/iwf_prioritized_bulls_arsucdv14_bams/aligns_v1.2/014HO05680/014HO05680.sorted.merged.bam.dedup.bam ARSUCD1.2.current_ref.fa /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR/iwf_prioritized_bulls_arsucdv14_bams/aligns_v1.2/014HO05680/014HO05680.sorted.merged.bam.intervals -partition msn

	INDEL realignment complete! Hopefully this time we get no truncated verisons...Will check these after the other 6 bams are ready...
	
	# Let's see how many bams are under 5X coverage
	perl -lane 'if($F[4] < 6){print $_;}' < new.dedup.realn.bam.stats | wc -l
	9
	# Not too many. Now let's check high deviation from the mapping percentages
	perl -lane 'if($F[5] == 0){next;} if($F[2]/$F[1] < 0.97){print $_;}' < new.dedup.realn.bam.stats | wc -l
	4
	
	
	# So I will drop all bams under 5X coverage and the < 97% mapping rate bams
	perl -lane 'if($F[5] == 0){next;} if($F[4] > 6 & $F[2]/$F[1] > 0.97){print $F[0];}' < new.dedup.realn.bam.stats > dedup.realn.filtered.final.bam.list

	perl -ne '$_ =~ s/^\./\/beegfs\/project\/rumen_longread_metagenome_assembly\/kiranmayee\/CDDR/; print $_;' < dedup.realn.filtered.final.bam.list > dedup.realn.filtered.final.bam.fullpath.list


Running GATK BQSR for the above deduplicated and indel realigned bam...
	
	module load gatk/3.6

	java -jar /software/7/apps/gatk/3.6/GenomeAnalysisTK.jar -T BaseRecalibrator -R ARSUCD1.2.current_ref.fa -I iwf_prioritized_bulls_arsucdv14_bams/
	
	aligns_v1.2/014HO06441-6441_qia_test/014HO06441-6441_qia_test.dedup.realn.bam -o bqsr/recal_data.table

Ok looks like it needs a known sites file, lets use bovineHD SNPs called by SAMtools using the arsucdv14 coordinates and lifted over to arsucdv1.2 assembly

Ok now it throws another error: the input bam file needs PL tags filled in by GATK... 


### Variant calling using SAMtools ###


### Variant calling using GATK HaplotypeCaller ###