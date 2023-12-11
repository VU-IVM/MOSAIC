#!/bin/bash
#SBATCH -p cbuild
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH --mem=40G


#module purge
#module load 2022 Cartopy/0.20.3-foss-2022a
source activate hydromt-sfincs_latest

python sfincs_hmax_plot_v2.py
echo "done"