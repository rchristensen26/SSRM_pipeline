#!/bin/bash
#
#
#SBATCH --job-name=2_smk_hmmer
#
#SBATCH --output=2_smk_hmmer.out
#SBATCH --error=2_smk_hmmer.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=150G
#SBATCH --partition=relman

# activate conda in general
#source /home/users/reb26/.bashrc # works because the conda init setting is here
source /home/users/jgrembi/.bashrc # works because the conda init setting is here

# activate conda environment for snakemake
conda activate snakemake

# load python and biopython
module load python/3.9
module load biology
module load py-biopython/1.79_py39
module load numpy/1.26

# go to the main project directory
#cd /home/users/reb26/SRM22/main
cd /oak/stanford/groups/relman/users/jgrembi/srb/SSRM_pipeline


snakemake --cores 24 --use-conda -j 24 --max-jobs-per-second 3 --max-status-checks-per-second 3
