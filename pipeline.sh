#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=20000
#SBATCH --partition=assemble2
#SBATCH --dependency=afterok:818517
#SBATCH --job-name=Run1
#SBATCH -o run1
#SBATCH -e run1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

  
cd $SLURM_SUBMIT_DIR
echo $PWD

Rscript --vanilla extractvcfinfo.R filtered_run1/filtered_run1.vcf.gz

echo extracted VCF Info 

module load  plink/1.90b4.4-2017-05-21

plink --vcf filtered_run1/filtered_run1.vcf.gz --cow --double-ids --nonfounders --freq --missing --hardy --keep-autoconv --out filtered_run1/filtered_run1

module unload plink/1.90b4.4-2017-05-21

module load plink/2.00alM-2017-05-22

plink2 --bfile filtered_run1/filtered_run1 --cow --nonfounders --freq --out filtered_run1/filtered_run1

echo finished calculating plink stats

Rscript --vanilla convertRDS.R filtered_run1/filtered_run1.frq filtered_run1/filtered_run1.afreq

Rscript --vanilla convertRDS_fmiss.R filtered_run1/filtered_run1.lmiss filtered_run1/filtered_run1.imiss

echo Completed exporting data to rds, now ready to plot...

Rscript --vanilla plot_MAF.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.frq.rds filtered_run1/filtered_run1.frq.rds filtered_run1/MAF.pdf

Rscript --vanilla plot_MAF.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.lmiss.rds filtered_run1/filtered_run1.lmiss.rds filtered_run1/Fraction_missing_per_variant.pdf

Rscript --vanilla print_densities.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.ann.bgzip.vcf.gz.rds filtered_run1/filtered_run1.vcf.gz.rds filtered_run1/InfoDensities.pdf

Rscript --vanilla plot_Alt_Allele_Freq.R filtered_run1/filtered_run1.afreq /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.afreq filtered_run1/AAF_run1.pdf

Rscript --vanilla plot_fmiss_per_sample.R filtered_run1/filtered_run1.imiss.rds /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.imiss.rds filtered_run1/Fraction_missing_per_sample.pdf

Rscript --vanilla plot_cor_2.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.ann.bgzip.vcf.gz.rds filtered_run1/filtered_run1.vcf.gz.rds filtered_run1/Corrplot.pdf

R -e "rmarkdown::render('summary.Rmd',output_file='filtered_run1/filtered_run1.html')"