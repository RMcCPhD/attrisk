#!/bin/bash

############# SLURM SETTINGS #############
#SBATCH --account=none
#SBATCH --job-name=run_dattr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=nodes
#SBATCH --time=0-00:59:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=1

############# LOADING MODULES (optional) #############
module load apps/R

############# MY CODE #############
Rscript 03_obj2_test1.R
