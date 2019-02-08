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




