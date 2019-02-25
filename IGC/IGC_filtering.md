# IGC markers filtering and formatting for Natdb imputation #

#### First identifying and filtering ####

Working directory: /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs

    module load bcftools
    
    [kiranmayee.bakshy@assembler2 condensed_vcfs]$ bcftools concat -a -O z -o igc.vcf.gz Domino.vcf.gz HF.vcf.gz LIB14413.vcf.gz LIB14427.vcf.gz TPI4222.vcf.gz CH240.vcf.gz
	[kiranmayee.bakshy@assembler2 condensed_vcfs]$ gunzip -c igc.vcf.gz | grep -v '#' | wc -l
	33934
    [kiranmayee.bakshy@assembler2 condensed_vcfs]$ gunzip -c igc.vcf.gz | grep -v '#' | grep -v 'INDEL' | perl -lane '($ac, $an) = $F[7] =~ /AC=(\d{1,3}).+AN=(\d{1,3})/; if($F[5] == 999 && $ac / $an > 0.25 && $an > 187){print $F[5];}' | wc -l
    6202
    
    [kiranmayee.bakshy@assembler2 condensed_vcfs]$ gunzip -c igc.vcf.gz | grep -v 'INDEL' | perl -lane '($ac, $an) = $F[7] =~ /AC=(\d{1,3}).+AN=(\d{1,3})/; if($F[5] == 999 && $ac / $an > 0.25 && $an > 187){print $F[5];}' > igc.filtered.vcf
    [kiranmayee.bakshy@assembler2 condensed_vcfs]$ cat igc.filtered.vcf | grep -v '#' | wc -l
    6202



     bcftools annotate -I 'ARS\_PIRBRIGHT\_%CHROM\_%POS\_%REF\_%ALT' igc.vcf -O z -o igc.id.vcf.gz
    
    
    
    [kiranmayee.bakshy@assembler2 liftover_to_v1.2]$ bgzip /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/selected.ARS-UCDv1.2.NCBI_ID_trimmed.vcf
    
    [kiranmayee.bakshy@assembler2 liftover_to_v1.2]$ bcftools index -c /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/selected.ARS-UCDv1.2.NCBI_ID_trimmed.vcf.gz
    
    [kiranmayee.bakshy@assembler2 liftover_to_v1.2]$ bcftools concat -a igc.id.sampled.vcf.gz /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_umd3/annotated/filtration/filtered_run3/selected.ARS-UCDv1.2.NCBI_ID_trimmed.vcf.gz | bgzip > final_filtered_set.vcf.gz
    
    [kiranmayee.bakshy@assembler2 liftover_to_v1.2]$ gunzip -c final_filtered_set.vcf.gz | grep -v '#' | wc -l
    

