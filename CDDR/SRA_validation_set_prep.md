# Preparation of SRA Holstein WGS dataset for GATK cross validation #

*3/7/2019*

These are my notes on selecting, downloading and downstream processing of the SRA datasets specially chosen from Holstein WGS at high coverage.

This dataset after alignment to the ARSUCDv1.2, variant calling using SAMtools and filtration will be used for validating the GATK VQSR model.

*working directory is on CERES cluster: /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR/sra_data*
  
SRA search term: •	(cattle holstein) AND "Bos taurus"[orgn:__txid9913]

All the 100 accessions selected were saved in the file at /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR/sra_data/SraAccList_all100.txt

The experimental details are saved in the same directory as sra_results.csv

### Download ###

The following script (download-sra.sh) was used to download all the accessions from SRA:

	#!/usr/bin/sh
	#SBATCH --nodes=1
	#SBATCH --mem=80G
	#SBATCH --partition=medium
	#SBATCH --ntasks-per-node=10
	
	module load sratoolkit/2.9.0
	
	START=$(date +%s.%N)
	
	echo $1
	mkdir $1
	fastq-dump -A $1 -F --split-files --gzip -O $1
	
	
	END=$(date +%s.%N)
	DIFF=$(echo "$END - $START" | bc)
	echo "Executed in $DIFF time from $START and $END"


	# The commands used to queue up the jobs:
	cat SraAccList_all100.txt | xargs -I {} sbatch download-sra.sh {}

	
I had issues with temporary directory getting filled because the default tmp directory is in the home folder.

I created a file which will redirect ncbi default tmp to the specified folder which has lot of memory for the sra.cache files.

Basically I created a file in my home folder: ".ncbi/user-settings.mkfg"
with the following line: /repository/user/main/public/root = "/beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR/sra_data/tmp"
 

If the download was successful, the fastq-dump would give us some stats. 

After checking the slurm outputs, I found that 2 of the files have not downloaded properly: ERR2694977, ERR2694981

I queued these up again manually. I removed all the sra.cache files from the temporary directory to avoid filling up the memory.

  
### Quality check ###
 

I have to run fastqc to check for the quality and get a consolidated report using MultiQC

Here is the script (runFastQC.sh) to run fastqc:

	#!/bin/sh
	#SBATCH -p short
	#SBATCH --nodes=1
	#SBATCH --mem=3000
	#SBATCH --ntasks-per-node=5
	
	module load fastqc/0.11.5
	
	fastqc -o /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR/sra_data/qc --noextract -t 4 -d /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR/sra_data/tmp/fqc $1 


	# Now queuing up the fastqc jobs for all except the 2 unsuccessful samples:

	ls ./*/*.fastq.gz | grep -v 'ERR2694977' | grep -v 'ERR2694981' | xargs -I {} sbatch runFastQC.sh {}