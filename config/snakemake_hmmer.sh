#!/bin/bash
#
#
#SBATCH --job-name=2_smk_hmmer
#
#SBATCH --output=2_smk_hmmer.out
#SBATCH --error=2_smk_hmmer.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=150G

# activate conda in general
source /home/users/reb26/.bashrc # works because the conda init setting is here

# activate conda environment for snakemake
conda activate snakemake

# load python and biopython
module load python/3.9
module load biology
module load py-biopython/1.79_py39

# go to the main project directory
cd /home/users/reb26/SRM22/main

snakemake --cores 24 --use-conda -j 24 --max-jobs-per-second 3 --max-status-checks-per-second 3
