# Variants filtration using the coordinates directly lifted from ARS-UCDv14 to ARS-UCDv1.2

## This is a continuation from the assembly liftover notes file...

Prepared lifted over regions-sample file for progressive filtration

Out of 1784 haplotype groups only 1690 groups have been lifted over to the current assembly...

Working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration 
    
    [kiranmayee.bakshy@assembler2 filtration]$ bcftools stats f1.vcf.gz | grep -P "SN\t"
    # SN[2]id   [3]key  [4]value
    SN  0   number of samples:  172
    SN  0   number of records:  17883218
    SN  0   number of no-ALTs:  0
    SN  0   number of SNPs: 17883218
    SN  0   number of MNPs: 0
    SN  0   number of indels:   0
    SN  0   number of others:   0
    SN  0   number of multiallelic sites:   0
    SN  0   number of multiallelic SNP sites:   0
    [kiranmayee.bakshy@assembler2 filtration]$ bcftools stats f2.vcf.gz | grep -P "SN\t"
    # SN[2]id   [3]key  [4]value
    SN  0   number of samples:  172
    SN  0   number of records:  17883171
    SN  0   number of no-ALTs:  0
    SN  0   number of SNPs: 17883171
    SN  0   number of MNPs: 0
    SN  0   number of indels:   0
    SN  0   number of others:   0
    SN  0   number of multiallelic sites:   0
    SN  0   number of multiallelic SNP sites:   0
    [kiranmayee.bakshy@assembler2 filtration]$ bcftools annotate --rename-chrs ARS-UCDv1.2_NKLS_chr_num.tab -Oz -o f2.1.vcf.gz f2.vcf.gz
    

First getting rid of singletons and selecting only the variants which are homozygous in the target haplotype groups.

