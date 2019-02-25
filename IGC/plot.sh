#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=20000
#SBATCH --partition=assemble2
#SBATCH --dependency=afterok:820859 
#SBATCH --job-name=plot_gemma
#SBATCH -o plot_gemma
#SBATCH -e plot_gemma
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

cd $SLURM_SUBMIT_DIR
echo $PWD
 
Rscript --vanilla gemma_results_plotting.R