#!/bin/bash

#SBATCH --job-name=delta_vs_sigma_vs_power_%j
#SBATCH --output=delta_vs_sigma_vs_power_%j.out
#SBATCH --error=delta_vs_sigma_vs_power_%j.err
#SBATCH --account=def-wjmarsha
#SBATCH --time=1-00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G

module load matlab/2019a

matlab -batch -nodisplay "delta_vs_sigma_vs_power"
