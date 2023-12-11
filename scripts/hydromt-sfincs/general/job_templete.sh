#!/bin/bash
#SBATCH -p cbuild
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH --mem=40G

module load 2021
module load Anaconda3/2021.05
source activate hydromt-sfincs

event="xynthia"
start_date="20100215"            # Default 7 days before landfall. Results are saved from 2 days before landfall until 5 days after. Changes of these settings can be made in the gtsm_template.py script
end_date="20100218"
meteo_forcing="ERA5"            # ERA5 will be used as input
bbox_sfincs=-1.869, 45.442, -0.755, 46.555 # min lon, min lat, max lon, max lat  


echo "Setting up SFINCS for $event $meteo_forcing"
if [ "$event" == "xynthia" ]; then
    #echo "Running GTSM with ERA5"
    python ../../../hydromt-sfincs/general/build_coastal_model.py "./$event""_$meteo_forcing" "$start_date 000000" "$end_date 000000" "{'bbox': [$bbox_sfincs]}" "fabdem_xynthia"
elif [ "$event" == "irma" ]; then
    #echo "Running GTSM with the Holland model"
    python ../../../hydromt-sfincs/general/build_coastal_model.py "./$event""_$meteo_forcing" "$start_date 000000" "$end_date 000000" "{'bbox': [$bbox_sfincs]}" "fabdem_irma" 
elif [ "$event" == "haiyan" ]; then
    #echo "Running GTSM with the Holland model"
    python ../../../hydromt-sfincs/general/build_coastal_model.py "./$event""_$meteo_forcing" "$start_date 000000" "$end_date 000000" "{'bbox': [$bbox_sfincs]}" "fabdem_haiyan" 
else
    echo "No event available"
fi
echo "HydroMT done"
#cd "$event_$meteo_forcing" 
#singularity run ../sfincs-cpu_latest.sif
#echo "SFINCS done"