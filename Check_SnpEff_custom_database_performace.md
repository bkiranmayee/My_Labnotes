# Check performance of custom build SnpEff ARS-UCDv1.2 annotation database

These are my notes on how to check the performance of the custom annotation database that I have built for ARS-UCDv1.2 bovine reference genome.

Main idea is to compare the number and type of SNP predictions of 
  a) the ARS-UCDv1.2 SnpEff summary of ARS-UCDv1.2 variant coordinates and 
  b) the UMD3.1.86 SnpEff summary of the same coordinates lifted to UMD3.1.1 
  
  Now first the liftover...
  
  working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2
  
  ```bash
[kiranmayee.bakshy@assembler3 GCA_002263795.2_ARS-UCD1.2]$ cd /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2
[kiranmayee.bakshy@assembler3 liftover_to_v1.2]$ module load bcftools
[kiranmayee.bakshy@assembler3 liftover_to_v1.2]$ bcftools query -f'%CHROM\t%POS0\t%END\t%ID\n' combined.vcf.gz > arsucdv1.2_combined.bed
[kiranmayee.bakshy@assembler3 liftover_to_v1.2]$ /mnt/nfs/nfs2/bickhart-users/binaries/kentUtils/bin/linux.x86_64/liftOver arsucdv1.2_combined.bed /mnt/nfs/nfs2/bickhart-users/cattle_asms/liftovers/ARS-UCD1.2_to_UMD3.11/ARS-UCD1.2_to_umd3_kary_unmask_ngap.mmap.liftover.chain umd3/mapped umd3/unmapped
Reading liftover chains
Mapping coordinates
[kiranmayee.bakshy@assembler3 liftover_to_v1.2]$ wc -l arsucdv1.2_combined.bed umd3/mapped umd3/unmapped
  23484333 arsucdv1.2_combined.bed
  22867323 umd3/mapped
   1234020 umd3/unmapped
  47585676 total
[kiranmayee.bakshy@assembler3 liftover_to_v1.2]$ /mnt/nfs/nfs2/bickhart-users/binaries/bin/sort-bed umd3/mapped > umd3/mapped.sorted.bed
```
  
