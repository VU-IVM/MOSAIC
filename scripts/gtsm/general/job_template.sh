#!/bin/bash
#Set job requirements
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --partition=thin
#SBATCH --time=01:00:00
#SBATCH -o log.%j.o
#SBATCH -e log.%j.e

# note the ntasks=64 above is too much and should be maybe set to 32

echo $OMP_NUM_THREADS

# modules!
source /gpfs/work2/0/einf2224/code/setdev

# input = era5 or holland
start_date="20170910"            # Default XX period that is run, and saved from the XXX day. Changes of these can be made in the gtsm_template.py script
meteo_forcing="ERA5"            # ERA5 will be used as input
#meteo_forcing="track data"      # track data will be used as input and the holland model will be executed
track_file="/gpfs/home4/benitoli/papers/paper1/scripts/gtsm/general/fine_Holland/bal112017.atcf"
bbox_gtsm=-84.,-78.,23.,34.      # bounding box where the map data will be saved at: [lon min, lon max, lat min, lat max]
# maybe by default now I would do that I have a certain period for which I run tides, and also for which I save data. To make it simple.

 
if [ "$meteo_forcing" == "ERA5" ]; then
    echo "Running GTSM with ERA5"
    mpirun --report-bindings -np 1 python /gpfs/home4/benitoli/papers/paper1/scripts/gtsm/general/gtsm_era5.py "$start_date" "$bbox_gtsm"
elif [ "$meteo_forcing" == "track data" ]; then
    echo "Running GTSM with the Holland model"
    mpirun --report-bindings -np 1 python /gpfs/home4/benitoli/papers/paper1/scripts/gtsm/general/gtsm_holland.py "$start_date" "$bbox_gtsm" "$track_file" 
else
    echo "Meteorological forcing not defined"
fi



# add in this same bash file hydromt and sfincs!
