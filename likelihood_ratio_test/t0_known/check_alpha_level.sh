#!/bin/bash

#SBATCH --job-name=check_alpha_level_%j
#SBATCH --output=check_alpha_level_%j.out
#SBATCH --error=check_alpha_level_%j.err
#SBATCH --account=def-wjmarsha
#SBATCH --time=3-00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=30G

module load matlab/2019a

matlab -batch -nodisplay "check_alpha_level"
