#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=20000
#SBATCH --partition=assemble2
<<<<<<< HEAD
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
=======
#SBATCH --job-name=Run2
#SBATCH -o filtered_run2/run2
#SBATCH -e filtered_run2/run2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

#srun --nodes 2 --tasks 2 mkdir -p /tmp/$USER/$SLURM_JOBID
  
cd $SLURM_SUBMIT_DIR

echo $PWD

srun perl progressiveSelectionKBv2.pl f2.vcf.gz region_sample.txt

module load samtools htslib plink/1.90b4.4-2017-05-21

srun bgzip filtered_run2/filtered_run2.vcf

srun tabix filtered_run2/filtered_run2.vcf.gz

#zgrep -v -E '^#' filtered_run2/filtered_run2.vcf.gz | cut -f 1 | sort | uniq -c > filtered_run2/v_per_chr

srun gunzip -c filtered_run2/filtered_run2.vcf.gz | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -o filtered_run2/v_per_chr -e '#' -m

srun sort -k1,1 -V -s filtered_run2/v_per_chr > filtered_run2/v_sorted_per_chr

srun Rscript --vanilla extract_vcfinfo.R filtered_run2/filtered_run2.vcf.gz &

echo extracted VCF Info 

srun plink --vcf filtered_run2/filtered_run2.vcf.gz --cow --double-id --nonfounders --freq --missing --hardy --keep-autoconv --out filtered_run2/filtered_run2 &
>>>>>>> f6a7bfd0a580a220a0a4eaad360322eb9c9c1143

module unload plink/1.90b4.4-2017-05-21

module load plink/2.00alM-2017-05-22

<<<<<<< HEAD
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
=======
srun plink2 --bfile filtered_run2/filtered_run2 --cow --nonfounders --freq --out filtered_run2/filtered_run2

echo finished calculating plink stats

srun Rscript --vanilla convertToRDS.R filtered_run2/filtered_run2.frq 

srun Rscript --vanilla convertToRDS.R filtered_run2/filtered_run2.afreq 

srun Rscript --vanilla convertToRDS.R filtered_run2/filtered_run2.lmiss

echo Completed exporting data to rds, now ready to plot...

srun Rscript --vanilla plot_MAF.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.frq.rds filtered_run2/filtered_run2.frq.rds filtered_run2/MAF_run2.pdf

srun Rscript --vanilla plot_MAF.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.lmiss.rds filtered_run2/filtered_run2.lmiss.rds filtered_run2/Fraction_missing_per_variant_run2.pdf

srun Rscript --vanilla print_densities.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.ann.bgzip.vcf.gz.RData filtered_run2/filtered_run2.vcf.gz.rds filtered_run2/InfoDensities_run2.pdf

srun Rscript --vanilla plot_Alt_Allele_Freq.R filtered_run2/filtered_run2.afreq /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.afreq filtered_run2/AAF_run2.pdf

srun Rscript --vanilla plot_fmiss_per_sample.R filtered_run2/filtered_run2.imiss /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.imiss filtered_run2/Fraction_missing_per_sample_run2.pdf

srun Rscript --vanilla plot_cor_2.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.ann.bgzip.vcf.gz.RData filtered_run2/filtered_run2.vcf.gz.rds filtered_run2/Corrplot_run2.pdf

srun R -e "rmarkdown::render('summary.Rmd',output_file='filtered_run2/filtered_run2.html')"

wait
>>>>>>> f6a7bfd0a580a220a0a4eaad360322eb9c9c1143
