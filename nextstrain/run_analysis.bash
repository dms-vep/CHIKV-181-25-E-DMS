#!/bin/bash
#
#SBATCH -c 1
#SBATCH --time 1-0

# stop on errors
set -e

# Capture the Snakemake version
SNAKEMAKE_VERSION=$(snakemake --version)
echo "Running Snakemake version ${SNAKEMAKE_VERSION}..."

# Run the main analysis on `slurm` cluster
snakemake \
    -R filter_records \
    --software-deployment-method conda \
    --conda-frontend mamba \
    --rerun-incomplete