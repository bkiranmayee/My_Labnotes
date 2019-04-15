#Construction of New bovine assembly (ARS-UCD.v.1.2) annotation database


## Building ANNOVAR database for the new bovine reference assembly

The GFF3 file downloaded from NCBI consists of contig lines (feature field=region) which are not recognized by gff3_to_gtf converter programs such as gffread, gt gff3_to_gtf and so on...They just ignore these region lines but still give an output gtf which has no mRNA features in it. It only contains CDS and exon records...

I am trying to build ANNOVAR database which is more flexible for other file formats specially GFF3

Building a database in ANNOVAR is explained [here](http://annovar.openbioinformatics.org/en/latest/user-guide/gene/) in a section called "What about GFF3 file for new species?"

Trying to convert GFF3 to GenePred format...

I first removed the lines which are not recognized by the converter

    [kiranmayee.bakshy@assembler2 ARS-UCD.v1.2]$ cat genes_chr1_29_x.gff | grep -v 'tissue-type' | grep -v 'Is_circular' > genes_chr_v2.gff 
	
	[kiranmayee.bakshy@assembler2 ARS-UCD.v1.2]$ /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/gff3ToGenePred -maxParseErrors=-1 -maxConvertErrors=-1 genes_chr_v2.gff Bt_refGene.txt
    Error: no exon in rna30616 contains CDS 45375559-45375640
    Error: no exon in rna30615 contains CDS 45375559-45375640
    Error: no exon in rna22183 contains CDS 21284669-21284727
    Error: no exon in rna11247 contains CDS 12071563-12072733
    Error: no exon in rna11245 contains CDS 12071563-12072733
    Error: no exon in rna7688 contains CDS 19018057-19018121
    6 errors converting GFF3 file: genes_chr_v2.gff


Ok now this seems to work...

Next step is making the master mRNA fasta file

For this step we need to have a working version of ANNOVAR.

ANNOVAR is freely available for non-profit research institutions. I got a copy by registering at the [ANNOVAR](http://www.openbioinformatics.org/annovar/annovar_download_form.php) website.

Just download the tarball and do the following:

    tar -xvf annovar.latest.tar.gz
	mkdir annovar/ARS-UCD.v1.2_db
    	
    # download a copy of reference fasta and the genes from the gff3 file to the database folder
	cp snpEff/data/ARS-UCD.v1.2/Bt_refGene.txt annovar/ARS-UCD.v1.2_db/
	cp snpEff/data/genomes/ARS-UCD.v1.2.fa annovar/ARS-UCD.v1.2_db/
	# Now create the mRNA fasta file
    perl annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile annovar/ARS-UCD.v1.2_db/ARS-UCD.v1.2.fa annovar/ARS-UCD.v1.2_db/Bt_refGene.txt --out Bt_refGeneMrna.fa
    NOTICE: Finished writting FASTA for 65980 genomic regions to Bt_refGeneMrna.fa
    WARNING: 342 gene regions do not have complete ORF (for example, rna48852NC_037344.1:69961538, rna64507NC_037350.1:27782137, rna43899NC_037342.1:75880413, rna2094NC_037328.1:117507279, rna53460NC_037346.1:9130187)

The ANNOVAR annotation database ARS-UCD.v1.2_db is ready for gene-based annotations.

However, there were several warnings during the construction of this database all of which were due to the unplaced contig information in the GFF3 file. So there shouldn't be any further problems I guess.

For example, 

    WARNING: Cannot identify sequence for rna78146 (starting from NW_020190140.1:5816)
    WARNING: Cannot identify sequence for rna78145 (starting from NW_020190140.1:3960)
    WARNING: Cannot identify sequence for rna78144 (starting from NW_020190140.1:3960)
    WARNING: Cannot identify sequence for rna78143 (starting from NW_020190140.1:3960)
    WARNING: Cannot identify sequence for rna78142 (starting from NW_020190140.1:3960)
    WARNING: Cannot identify sequence for rna78141 (starting from NW_020190140.1:1638)
    WARNING: Cannot identify sequence for rna78140 (starting from NW_020190140.1:1638)

    
The NW_ tags all belong to the unplaced sequences...

This is how we can use this database to annotate the vcf files:

    perl /mnt/nfs/nfs1/kiranmayee.bakshy/annovar/table_annovar.pl <(gunzip -c 1.ncbi.vcf.gz) /mnt/nfs/nfs1/kiranmayee.bakshy/annovar/ARS-UCD.v1.2_db/ -buildver Bt -out 1.ncbi.vcf.ann -protocol refGene -operation g -nastring . -vcfinput
    


Meanwhile, I have to convert the NKLS02 chromosome codes in the lifted over vcf files to the NCBI codes 

    bcftools annotate --rename-chrs ARS-UCD_ncbi_to_chrnum.tab combined.vcf.gz | bgzip > combined.ncbi.vcf.gz
	
	# to test the annotation database on a smaller file...
	bcftools annotate --rename-chrs ARS-UCD_ncbi_to_chrnum.tab 1.vcf.gz | bgzip > 1.ncbi.vcf.gz
 	
	# also change the snp ids
	 bcftools annotate -I 'bakshy_%CHROM\_%POS\_%REF\_%FIRST_ALT' combined.ncbi.vcf.gz |  bgzip  > combined.ncbi.id.vcf.gz
    