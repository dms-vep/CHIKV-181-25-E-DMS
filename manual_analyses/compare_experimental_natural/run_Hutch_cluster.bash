#!/bin/bash
#
#SBATCH -c 32

snakemake \
    -j 32 \
    --software-deployment-method conda \
    --rerun-incomplete \
    --keep-going
