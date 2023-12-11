#!/bin/bash
#Set job requirements
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --partition=thin
#SBATCH --time=00:10:00
#SBATCH --mem=40G
#module purge
#module load 2022
#module load Cartopy/0.20.3-foss-2022a Anaconda3/2022.05
#module load Anaconda3/2022.05
#source activate gtsmvalidation
#module load Python/3.9.5-GCCcore-10.3.0
#module load SciPy-bundle/2021.05-foss-2021a
#module load Cartopy/0.20.0-foss-2021a
#module load matplotlib/3.4.2-foss-2021a
#source ~/JHS_installations/venvs/plot_p1/bin/activate

module purge
module load 2022 Cartopy/0.20.3-foss-2022a


#python gtsm_test_partitions.py
python gtsm_output_v6.py 
echo "done"
