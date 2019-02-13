# Calculation of Genomic coverage by a feature

Input files needed are gff3 annotation file and genomic sizes file

We can use awk and BEDOPS for calculating the percentage of genome covered by each feature (exon, mRNA, intron, intergenic and so on).

Working directory: /mnt/nfs/nfs1/kiranmayee.bakshy/snpEff/data/genomes

```bash
# First prepare the bed file for the required feature:
awk '$3 == "exon"' original.genes.gff | /mnt/nfs/nfs2/bickhart-users/binaries/bin/convert2bed -i gff - -d > exons.bed
# sort the bed file
/mnt/nfs/nfs2/bickhart-users/binaries/bin/sort-bed exons.bed > exons.sorted.bed
# Prepare the genome size file 
cut -f1,2 /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARSUCD1.2.current_ref.fa.fai > genomeSize.bed
# Now calculate the total genome size
cowSize=`awk '{s+=$2;}END{print s;}' genomeSize.bed` = 2628411261
# Use BEDOPS bedmap to map the sizes of merged (overlapping) exons over the current genomic space. 
# Then calculate the fraction (percentage/100) of space taken up by exons by dividing by the size of the genome build
/mnt/nfs/nfs2/bickhart-users/binaries/bin/bedops --merge exons.sorted.bed | /mnt/nfs/nfs2/bickhart-users/binaries/bin/bedmap --echo-ref-size - genomeSize.bed | awk '{s+=$1;}END{print s/2628411261;}'
0.0208454
# What would the coverage be if we don't merge the bed coordinates of the exons
/mnt/nfs/nfs2/bickhart-users/binaries/bin/bedmap --echo-ref-size exons.sorted.bed genomeSize.bed | awk '{s+=$1;}END{print s/2628411261;}'
0.0364536
```

The exons from the Ensembl gff3 file cover around 2% of the current bovine reference genome ARS-UCDv1.2

```bash
awk '$3 == "mRNA"' genes.gff | /mnt/nfs/nfs2/bickhart-users/binaries/bin/convert2bed -i gff - -d > mRNA.bed
/mnt/nfs/nfs2/bickhart-users/binaries/bin/sort-bed mRNA.bed > mRNA.sorted.bed
/mnt/nfs/nfs2/bickhart-users/binaries/bin/bedops --merge mRNA.sorted.bed | /mnt/nfs/nfs2/bickhart-users/binaries/bin/bedmap --echo-ref-size - genomeSize.bed | awk '{s+=$1;}END{print s/2628411261;}'
0.360348
/mnt/nfs/nfs2/bickhart-users/binaries/bin/bedmap --echo-ref-size mRNA.sorted.bed genomeSize.bed | awk '{s+=$1;}END{print s/2628411261;}' 0.832335
```

```bash
wc -l exons.sorted.bed
433783
/mnt/nfs/nfs2/bickhart-users/binaries/bin/bedops --merge exons.sorted.bed | wc -l
215169
awk '{print $4}' mRNA.sorted.bed | uniq | wc -l
37506
/mnt/nfs/nfs2/bickhart-users/binaries/bin/bedops --merge mRNA.sorted.bed | wc -l
19813
```
Now to make a list of features and calculate the genome coverage by each feature, I will automate by wrinting a small bash script.

Here is the bash script to do this

```bash
#!/usr/bin/sh

convert2bed=/mnt/nfs/nfs2/bickhart-users/binaries/bin/convert2bed
sortbed=/mnt/nfs/nfs2/bickhart-users/binaries/bin/sort-bed
bedops=/mnt/nfs/nfs2/bickhart-users/binaries/bin/bedops
bedmap=/mnt/nfs/nfs2/bickhart-users/binaries/bin/bedmap
#cowGSize=awk '{s+=$2;}END{print s;}' genomeSize.bed

#$1 = gff file
#$2 = genome_size.bed

grep -v '#' $1 | awk '{print $3}' | sort | uniq  > features_all

while read f; do
	gen_cov=`grep "$f" $1 | $convert2bed -i gff - -d | $sortbed - | $bedops --merge - | $bedmap --echo-ref-size - $2 | awk '{s+=$1;}END{print s/2628411261;}'`
	echo "$f"	"$gen_cov"
done <features_all
```


```bash
sh genome_cov.sh genes.gff genomeSize.bed
```



**Feature**|**Fraction of Genome coverage**
:-----:|:-----:
J\_gene\_segment|1.27E-07
scRNA|1.96E-06
C\_gene\_segment|1.34E-05
V\_gene\_segment|1.85E-05
miRNA|2.84E-05
snoRNA|3.23E-05
snRNA|6.07E-05
rRNA|9.15E-05
pseudogenic\_transcript|0.000701685
pseudogene|0.000745154
five\_prime\_UTR|0.00184731
three\_prime\_UTR|0.00482611
lnc\_RNA|0.0127188
ncRNA|0.0128723
ncRNA\_gene|0.0128723
CDS|0.0133584
biological\_region|0.0166946
exon|0.0210438
mRNA|0.360348
gene|0.373804
region|1.03326


The above table is sorted by increasing order of genome coverage by the feature. Region indicates the chromosome, so it is represented as 100%.

The highest coverage is by mRNA and gene. So it is better we cut off the mRNA from the gff file to re-build annotation database. 



