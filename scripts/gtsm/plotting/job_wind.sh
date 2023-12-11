#!/bin/bash
#Set job requirements
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --partition=thin
#SBATCH --time=00:10:00
#SBATCH --mem=40G

module load 2021
module load Anaconda3/2021.05
source activate gtsmvalidation

python gtsm_output_windmap.py 
echo "done"
