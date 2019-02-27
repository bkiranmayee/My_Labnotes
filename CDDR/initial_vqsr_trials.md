Labnotes on VQSR trails
=======================

All datasets used are the new assembly coordinates lifted 




###run1_trial1
    
    java -Xmx8g -jar $GATK_JAR \
    -R /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARSUCD1.2.current_ref.fa \
    -T VariantRecalibrator \
    -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/f2.vcf.gz \
    -resource:1kbull_HQF,known=false,training=true,truth=false,prior=10.0 1kbull_sites_q999.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/resources/dbsnp_filtered_ARS-UCDv1.2_selected_uniq_snps.vcf.gz \
    -an DP \
    -an MQSB \
    -an BQB \
    -an RPB \
    -mode SNP \
    -recalFile modelling/run1/trial1.output.recal \
    -tranchesFile modelling/run1/trial1.output.tranches \
    -rscriptFile modelling/run1/trial1.output.plots.R \
    

###run2_trial1

    java -Xmx8g -jar $GATK_JAR \
    -R /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARSUCD1.2.current_ref.fa \
    -T VariantRecalibrator \
    -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/f2.vcf.gz \
    -resource:1kbull_Hol_HQF,known=false,training=true,truth=false,prior=10.0 1kbull_Hol_sites_q999.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbSNP.vcf.gz \
    -an DP \
    -an MQSB \
    -an BQB \
    -an RPB \
    -mode SNP \
    -recalFile modelling/run2/trial1.output.recal \
    -tranchesFile modelling/run2/trial1.output.tranches \
    -rscriptFile modelling/run2/trial1.output.plots.R \
    
### run3_trial1

    java -Xmx8g -jar $GATK_JAR \
    -R /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARSUCD1.2.current_ref.fa \
    -T VariantRecalibrator \
    -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/f2.vcf.gz \
    -resource:1kbull_HQF,known=false,training=true,truth=false,prior=10.0 1kbull_sites_q999.vcf.gz \
    -resource:Hol_HQF,known=false,training=true,truth=false,prior=8.0 Hol_filtered.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/resources/dbsnp_filtered_ARS-UCDv1.2_selected_uniq_snps.vcf.gz \
    -an DP \
    -an MQSB \
    -an BQB \
    -an RPB \
    -mode SNP \
    -recalFile modelling/run3/trial1.output.recal \
    -tranchesFile modelling/run3/trial1.output.tranches \
    -rscriptFile modelling/run3/trial1.output.plots.R \


### run1_trial2

    java -Xmx8g -jar $GATK_JAR \
    -R /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARSUCD1.2.current_ref.fa \
    -T VariantRecalibrator \
    -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/f2.vcf.gz \
    -resource:bovineHD,known=false,training=true,truth=true,prior=15.0 bovineHD.vcf.gz \
    -resource:1kbull_HQF,known=false,training=true,truth=false,prior=10.0 1kbull_sites_q999.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/resources/dbsnp_filtered_ARS-UCDv1.2_selected_uniq_snps.vcf.gz \
    -an DP \
    -an MQSB \
    -an BQB \
    -an RPB \
    -mode SNP \
    -recalFile modelling/run1/trial2.output.recal \
    -tranchesFile modelling/run1/trial2.output.tranches \
    -rscriptFile modelling/run1/trial2.output.plots.R \
 

### run3_trial2 ###
   
    java -Xmx8g -jar $GATK_JAR \
    -R /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARSUCD1.2.current_ref.fa \
    -T VariantRecalibrator \
    -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/f2.vcf.gz \
    -resource:bovineHD,known=false,training=true,truth=true,prior=15.0 bovineHD.vcf.gz \
    -resource:1kbull_HQF,known=false,training=true,truth=false,prior=10.0 1kbull_sites_q999.vcf.gz \
    -resource:Hol_HQF,known=false,training=true,truth=false,prior=8.0 Hol_filtered.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/resources/dbsnp_filtered_ARS-UCDv1.2_selected_uniq_snps.vcf.gz \
    -an DP \
    -an MQSB \
    -an BQB \
    -an RPB \
    -mode SNP \
    -recalFile modelling/run3/trial2.output.recal \
    -tranchesFile modelling/run3/trial2.output.tranches \
    -rscriptFile modelling/run3/trial2.output.plots.R \


### run2_trial2 ###

    java -Xmx8g -jar $GATK_JAR \
    -R /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARSUCD1.2.current_ref.fa \
    -T VariantRecalibrator \
    -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/f2.vcf.gz \
    -resource:bovineHD,known=false,training=true,truth=true,prior=15.0 bovineHD.vcf.gz \
    -resource:1kbull_Hol_HQF,known=false,training=true,truth=false,prior=10.0 1kbull_Hol_sites_q999.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/resources/dbsnp_filtered_ARS-UCDv1.2_selected_uniq_snps.vcf.gz \
    -an DP \
    -an MQSB \
    -an BQB \
    -an RPB \
    -mode SNP \
    -recalFile modelling/run2/trial2.output.recal \
    -tranchesFile modelling/run2/trial2.output.tranches \
    -rscriptFile modelling/run2/trial2.output.plots.R \
    

### run3_trial3 ###

    java -Xmx8g -jar $GATK_JAR \
    -R /mnt/nfs/nfs2/bickhart-users/cattle_asms/ncbi/ARSUCD1.2.current_ref.fa \
    -T VariantRecalibrator \
    -input /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover_to_v1.2/filtration/f2.vcf.gz \
    -resource:bovineHD,known=false,training=true,truth=true,prior=15.0 bovineHD.vcf.gz \
    -resource:1kbull_Hol_HQF,known=false,training=true,truth=false,prior=10.0 1kbull_Hol_sites_q999.vcf.gz \
    -resource:Hol_HQF,known=false,training=true,truth=false,prior=8.0 Hol_filtered.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/resources/dbsnp_filtered_ARS-UCDv1.2_selected_uniq_snps.vcf.gz \
    -an DP \
    -an MQSB \
    -an BQB \
    -an RPB \
    -mode SNP \
    -recalFile modelling/run3/trial3.output.recal \
    -tranchesFile modelling/run3/trial3.output.tranches \
    -rscriptFile modelling/run3/trial3.output.plots.R \
    