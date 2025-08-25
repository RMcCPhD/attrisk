#!/bin/bash

############# SLURM SETTINGS #############
#SBATCH --account=none
#SBATCH --job-name=fit_check_stan
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=nodes
#SBATCH --time=0-03:59:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=1

############# LOADING MODULES (optional) #############
module load apps/R

############# MY CODE #############
Rscript conduct.R
