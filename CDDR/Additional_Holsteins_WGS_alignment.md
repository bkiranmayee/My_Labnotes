# Alignment of additional Holstein WGS to the new assembly version ARS-UCD.v1.2
---
*2/6/2019*

These are my notes and commands for alignment and variant calling of additional Holsteins that serve as a validation dataset for GATK variant quality score recalibration modelling to identify tru variant locations.

The additional Holstein WGS data was already aligned to UMD3.1.1 and we need to re-align these aligned bams to the current assembly verion and perform variant calling.

Location of the aligned bams: AGIL cluster : /mnt/cifs/bickhart-qnap/1000_bulls_bams

## Table of Contents
* [Initial alignment runs](#testone)
* [Actual alignment runs]
* [Indel realignment and mark duplicates]
* [Variant calling]

<a name="testone"></a>
## Testing the alignment run

This is the first time I am going to align the WGS data to a reference genome.
The AGIL cluster is down, so I had to transfer the bams to the SCINET and start the work.

There are 37 bam files and 2 of them seem to have no data based on the file size.
The file sizes differ and this may because they were sequenced at different coverages. It sould be confirmed by running bamstats...

working directory: /beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR

```bash
# data transfer
scp -pr /mnt/cifs/bickhart-qnap/1000_bulls_bams/HO*reformatted.sorted.bam kiranmayee.bakshy@ceres:/beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR/Additional_Holstein_bams

scp -pr /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARSUCD1.2.current_ref.fa kiranmayee.bakshy@ceres:/beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR

scp -pr /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/indel_realign_markdups_kb.sh kiranmayee.bakshy@ceres:/beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/CDDR

# list all the bams, generate alignment scripts and queue them all up on the cluster
ls aligns/HO*/HO*sorted.merged.bam > sorted.merged.bamlist
 
perl /beegfs/project/rumen_longread_metagenome_assembly/binaries/perl_toolchain/sequence_data_pipeline/alignBamReadsToNewAssemSlurm.pl - b aligns -t bulls.tab -f ARSUCD1.2.current_ref.fa -m true

# get bam stats
perl ~/perl_toolchain/sequence_data_scripts/getBamStats.pl -b aligns/HOCAN000006193092/HOCAN000006193092.sorted.merged.bam,aligns/HOCAN000006229227/HOCAN000006229227.sorted.merged.bam,aligns/HOCAN000008432142/HOCAN000008432142.sorted.merged.bam,aligns/HODEU000000253642/HODEU000000253642.sorted.merged.bam,aligns/HODEU000341037501/HODEU000341037501.sorted.merged.bam,aligns/HOGBR000000598172/HOGBR000000598172.sorted.merged.bam,aligns/HOUSA000001244845/HOUSA000001244845.sorted.merged.bam,aligns/HOUSA000001417390/HOUSA000001417390.sorted.merged.bam,aligns/HOUSA000001427381/HOUSA000001427381.sorted.merged.bam,aligns/HOUSA000001447141/HOUSA000001447141.sorted.merged.bam,aligns/HOUSA000001537060/HOUSA000001537060.sorted.merged.bam,aligns/HOUSA000001556373/HOUSA000001556373.sorted.merged.bam,aligns/HOUSA000001563453/HOUSA000001563453.sorted.merged.bam,aligns/HOUSA000001667366/HOUSA000001667366.sorted.merged.bam,aligns/HOUSA000001672325/HOUSA000001672325.sorted.merged.bam,aligns/HOUSA000001682485/HOUSA000001682485.sorted.merged.bam,aligns/HOUSA000001697572/HOUSA000001697572.sorted.merged.bam,aligns/HOUSA000001721881/HOUSA000001721881.sorted.merged.bam,aligns/HOUSA000001810969/HOUSA000001810969.sorted.merged.bam,aligns/HOUSA000001879149/HOUSA000001879149.sorted.merged.bam,aligns/HOUSA000001903604/HOUSA000001903604.sorted.merged.bam,aligns/HOUSA000002026215/HOUSA000002026215.sorted.merged.bam,aligns/HOUSA000002030882/HOUSA000002030882.sorted.merged.bam,aligns/HOUSA000002040728/HOUSA000002040728.sorted.merged.bam,aligns/HOUSA000002041271/HOUSA000002041271.sorted.merged.bam,aligns/HOUSA000002064459/HOUSA000002064459.sorted.merged.bam,aligns/HOUSA000002103297/HOUSA000002103297.sorted.merged.bam,aligns/HOUSA000002125714/HOUSA000002125714.sorted.merged.bam,aligns/HOUSA000002147486/HOUSA000002147486.sorted.merged.bam,aligns/HOUSA000002266677/HOUSA000002266677.sorted.merged.bam,aligns/HOUSA000002284985/HOUSA000002284985.sorted.merged.bam,aligns/HOUSA000002290977/HOUSA000002290977.sorted.merged.bam,aligns/HOUSA000017056520/HOUSA000017056520.sorted.merged.bam,aligns/HOUSA000017062963/HOUSA000017062963.sorted.merged.bam,aligns/HOUSA000017349617/HOUSA000017349617.sorted.merged.bam,aligns/HOUSA000120754720/HOUSA000120754720.sorted.merged.bam,aligns/HOUSA000122358313/HOUSA000122358313.sorted.merged.bam -o bamStats.tab
```

## Realign indel and mark duplicates using picard

```bash
cat sorted.merged.bamlist | xargs -I {} sbatch indel_realign_markdups_kb.sh {} ARSUCD1.2.current_ref.fa
```

This script needs to be corrected because it was designed for AGIL cluster.

Ok now after correcting, it again failed because the sequence dictionary file for the ref genome fasta is absent.

Now generating one...

```bash
java -jar /software/7/apps/picard/64/2.9.2/picard.jar CreateSequenceDictionary R=ARSUCD1.2.current_ref.fa O=ARSUCD1.2.current_ref.dict
```
Now everything seems to work... 


## Variant calling

I need to select BAMS with atleast 20-25X coverage for variant calling, as this is one of the requirement to improve signal to noise ratio. 
So run the bam stats again on all the dedup bams now and tabulate.

```bash
# list all the dedup bams
ls aligns/*/*.dedup.bam > dedup.bam.list
perl /beegfs/project/rumen_longread_metagenome_assembly/binaries/perl_toolchain/sequence_data_scripts/getBamStats.pl -n dedup.bam.list
```

I realized that 10 samples are missing sorted.merged.dedup.bam.sorted.bam files.

I wonder if all the jobs ran to completion.

sacct -S 2019-01-01 -u kiranmayee.bakshy --format=User,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss,MaxVMSize,nnodes,ncpus,nodelist

Okay there are only 14 indel realignemnt jobs that exited with zero status. What happened to the rest?

There are 10 jobs which show the following error message: 

```bash
[kiranmayee.bakshy@sn-cn-8-1 CDDR]$ ls slurm* | xargs grep 'failed to open'
slurm-505033.out:samtools index: failed to open "aligns/HOCAN000006193092/HOCAN000006193092.sorted.merged.bam.dedup.bam": No such file or directory
slurm-505034.out:samtools index: failed to open "aligns/HOCAN000006229227/HOCAN000006229227.sorted.merged.bam.dedup.bam": No such file or directory
slurm-505035.out:samtools index: failed to open "aligns/HOCAN000008432142/HOCAN000008432142.sorted.merged.bam.dedup.bam": No such file or directory
slurm-505041.out:samtools index: failed to open "aligns/HOUSA000001427381/HOUSA000001427381.sorted.merged.bam.dedup.bam": No such file or directory
slurm-505044.out:samtools index: failed to open "aligns/HOUSA000001556373/HOUSA000001556373.sorted.merged.bam.dedup.bam": No such file or directory
slurm-505056.out:samtools index: failed to open "aligns/HOUSA000002040728/HOUSA000002040728.sorted.merged.bam.dedup.bam": No such file or directory
slurm-505059.out:samtools index: failed to open "aligns/HOUSA000002103297/HOUSA000002103297.sorted.merged.bam.dedup.bam": No such file or directory
slurm-505060.out:samtools index: failed to open "aligns/HOUSA000002125714/HOUSA000002125714.sorted.merged.bam.dedup.bam": No such file or directory
slurm-505064.out:samtools index: failed to open "aligns/HOUSA000002290977/HOUSA000002290977.sorted.merged.bam.dedup.bam": No such file or directory
slurm-505067.out:samtools index: failed to open "aligns/HOUSA000017349617/HOUSA000017349617.sorted.merged.bam.dedup.bam": No such file or directory
```

Now to restart indel realignment jobs that failed...

```bash
[kiranmayee.bakshy@sn-cn-8-1 CDDR]$ grep -Fvwf dedup.final.bam.list sorted.merged.bamlist > restart_indel_realignment
[kiranmayee.bakshy@sn-cn-8-1 CDDR]$ cat restart_indel_realignment.bamlist | xargs -I {} sbatch indel_realign_markdups_kb.sh {} ARSUCD1.2.current_ref.fa
```

The problem still exists, upon careful inspection I realized that the mark duplicates step is not running to completion.

The output.slurm records that the disk quota is full. I tried to find the solution via stackoverflow:

I need to set the tmp dir for the picard mark duplicates so that it doesn't use the default home directory.

Restarting the Markduplicates step for the 10 bams which failed at this step...

	# The modification for my script (markdups_kb.sh) was the following:
	java -jar /software/7/apps/picard/64/2.9.2/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT I=$1 O=$dedup_bam M=$duplicates_metrics TMP_DIR=`pwd`/tmp
 

### Preparing for variant calling ###

	# Select samples for variant calling
	ls ./*/*/*.dedup.bam > dedup.final.bam.list
	module load samtools; 
	perl /beegfs/project/rumen_longread_metagenome_assembly/binaries/perl_toolchain/sequence_data_scripts/getBamStats.pl -n dedup.final.bam.list -o summary_stats_dedup_bams.2019.02.20.tab
	
	
	# Let's see how many bams are under 5X coverage
	perl -lane 'if($F[4] < 6){print $_;}' < summary_stats_dedup_bams.2019.02.20.tab | wc -l
	24
	
	
	# Too many. Now let's check high deviation from the mapping percentages
	perl -lane 'if($F[5] == 0){next;} if($F[2]/$F[1] < 0.97){print $_;}' < summary_stats_dedup_bams.2019.02.20.tab | wc -l
	6
	
	# So I will drop all bams under 5X coverage and the < 97% mapping rate bams
	perl -lane 'if($F[5] == 0){next;} if($F[4] > 6 & $F[2]/$F[1] > 0.97){print $F[0];}' < summary_stats_dedup_bams.2019.02.20.tab > dedup.filtered.final.bam.list
	
	# The filtered list of bams is dedup.filtered.final.bam.list
	# Now to queue up the samtools mpileup scripts
	perl -ne '$_ =~ s/^\./\/beegfs\/project\/rumen_longread_metagenome_assembly\/kiranmayee\/CDDR\/Additional_Holstein/; print $_;' < dedup.filtered.final.bam.list > dedup.filtered.final.bam.fullpath.list

	perl /beegfs/project/rumen_longread_metagenome_assembly/binaries/perl_toolchain/sequence_data_pipeline/samtoolsMpileupBamsSlurm_kb.pl -b Additional_Holstein/vcfs -s ARSUCD1.2.current_ref.fa.samtools.mpileup.wins -t Additional_Holstein/dedup.filtered.final.bam.fullpath.list -f ARSUCD1.2.current_ref.fa -m true
	
For some reason that I can't figure out the above command doesn't work. The script fails to read the segments file...

I hard coded in the script samtoolsMpileupBamsSlurm_kb.pl to generate 10 Mb segments to generate the scripts and queue up. This seems to work...

	
	perl /beegfs/project/rumen_longread_metagenome_assembly/binaries/perl_toolchain/sequence_data_pipeline/samtoolsMpileupBamsSlurm_kb.pl -b Additional_Holstein/vcfs -t Additional_Holstein/dedup.filtered.final.bam.fullpath.list -f ARSUCD1.2.current_ref.fa -m true

The above command has generated 282 segments and queued up around 282 samtools mpileup scripts and 282 bcf call scripts and 1 concat script

It seems all the call scripts have mpileup command and not the call command.

I need to cancel all the call scripts and process them separately...

	squeue -u kiranmayee.bakshy | grep "call_" | awk '{print $1}' | xargs -I {} scancel {}

Now wait until the mpileup bcf files are ready and then process them to call the variants using process_bcfs.sh script. 








