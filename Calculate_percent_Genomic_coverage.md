# Calculation of Genomic coverage by a feature

Input files needed are gff3 annotation file and genomic sizes file

We can use awk and BEDOPS for calculating the percentage of genome covered by each feature (exon, mRNA, intron, intergenic and so on).

Working directory: /mnt/nfs/nfs1/kiranmayee.bakshy/snpEff/data/GCA_002263795.2_ARS-UCD1.2

```bash
# First prepare the bed file for the required feature:
awk '$3 == "exon"' genes.gff | /mnt/nfs/nfs2/bickhart-users/binaries/bin/convert2bed -i gff - -d > exons.bed
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





