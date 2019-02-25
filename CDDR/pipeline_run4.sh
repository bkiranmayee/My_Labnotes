#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=20000
#SBATCH --partition=assemble2
#SBATCH --job-name=Run3
#SBATCH -o filtered_run3/run3
#SBATCH -e filtered_run3/run3
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

#srun --nodes 2 --tasks 2 mkdir -p /tmp/$USER/$SLURM_JOBID
  
cd $SLURM_SUBMIT_DIR

echo $PWD

srun perl progressiveSelectionKBv2.pl f2.vcf.gz region_sample.txt

module load samtools htslib plink/1.90b4.4-2017-05-21

srun bgzip filtered_run3/filtered_run3.vcf

srun tabix filtered_run3/filtered_run3.vcf.gz

#zgrep -v -E '^#' filtered_run3/filtered_run3.vcf.gz | cut -f 1 | sort | uniq -c > filtered_run3/v_per_chr

srun gunzip -c filtered_run3/filtered_run3.vcf.gz | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -o filtered_run3/v_per_chr -e '#'

srun gunzip -c filtered_run3/filtered_run3.vcf.gz | perl ~/perl_toolchain/bed_cnv_fig_table_pipeline/tabFileColumnCounter.pl -f stdin -c 0 -o filtered_run3/v_per_chr_m -e '#' -m

srun sort -k1,1 -V -s filtered_run3/v_per_chr > filtered_run3/v_sorted_per_chr

srun sort -k1,1 -V -s filtered_run3/v_per_chr_m > filtered_run3/v_sorted_per_chr_m

srun Rscript --vanilla extract_vcfinfo.R filtered_run3/filtered_run3.vcf.gz &

echo extracted VCF Info 

srun plink --vcf filtered_run3/filtered_run3.vcf.gz --cow --double-id --nonfounders --freq --missing --hardy --keep-autoconv --out filtered_run3/filtered_run3 &

module unload plink/1.90b4.4-2017-05-21

module load plink/2.00alM-2017-05-22

srun plink2 --bfile filtered_run3/filtered_run3 --cow --nonfounders --freq --out filtered_run3/filtered_run3

echo finished calculating plink stats

srun Rscript --vanilla convertToRDS.R filtered_run3/filtered_run3.frq 

srun Rscript --vanilla convertToRDS.R filtered_run3/filtered_run3.afreq 

srun Rscript --vanilla convertToRDS.R filtered_run3/filtered_run3.lmiss

echo Completed exporting data to rds, now ready to plot...

srun Rscript --vanilla plot_MAF.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.frq.rds filtered_run3/filtered_run3.frq.rds filtered_run3/MAF_run3.pdf

srun Rscript --vanilla plot_MAF.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.lmiss.rds filtered_run3/filtered_run3.lmiss.rds filtered_run3/Fraction_missing_per_variant_run3.pdf

srun Rscript --vanilla print_densities.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.ann.bgzip.vcf.gz.RData filtered_run3/filtered_run3.vcf.gz.rds filtered_run3/InfoDensities_run3.pdf

srun Rscript --vanilla plot_Alt_Allele_Freq.R filtered_run3/filtered_run3.afreq /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.afreq filtered_run3/AAF_run3.pdf

srun Rscript --vanilla plot_fmiss_per_sample.R filtered_run3/filtered_run3.imiss /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.imiss filtered_run3/Fraction_missing_per_sample_run3.pdf

srun Rscript --vanilla plot_cor_2.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined.ann.bgzip.vcf.gz.RData filtered_run3/filtered_run3.vcf.gz.rds filtered_run3/Corrplot_run3.pdf

srun Rscript --vanilla per_chr_plot.R /mnt/nfs/nfs1/derek.bickhart/CDDR-Project/vcfs/condensed_vcfs/liftover/annotated/combined_sorted_per_chr filtered_run3/v_sorted_per_chr filtered_run3/Variants_per_chrom_plot_run3.pdf

srun R -e "rmarkdown::render('summary_3.Rmd',output_file='filtered_run3/filtered_run3.html')"

wait
