#!/bin/bash
#Set job requirements
#SBATCH --nodes=1
#SBATCH --ntasks=33
#SBATCH --cpus-per-task=1
#SBATCH --partition=thin
#SBATCH --time=01:00:00
#SBATCH -o log.%j.o
#SBATCH -e log.%j.e

# note the ntasks=64 above is too much and should be maybe set to 32

echo $OMP_NUM_THREADS

source /gpfs/work2/0/einf2224/code/setdev


#------------------------------------------------------------------------------------------------------------------------------------------------------
## models setup
#------------------------------------------------------------------------------------------------------------------------------------------------------
# inputs
start_date="20170910"                              # Default XX period that is run, and saved from the XXX day

# case studies with their bboxes for which gridded data will be saved at
cases=('irma', 'haiyan', 'xynthia')
declare -A bbox                                    # lon min, lon max, lat min, lat max  
bbox['irma']=(-82.2 -80.2 30.0 32.0)
bbox['haiyan']=(124.8 125.3 11.0 11.5)
bbox['xynthia']=(-1.8 -0.8 45.5 46.5)

# model configurations 
model_configs=('BS' 'TR' 'OR' 'IB')                #BS (Baseline scenario), TR (Temporal resolution Refined), OR (Output resolution Refined), IB (Improved Bathymetry)

#------------------------------------------------------------------------------------------------------------------------------------------------------
## model execution
#------------------------------------------------------------------------------------------------------------------------------------------------------

# Run GTSM for each case study and each model configuration
for case_name in "${!bbox[@]}"; do
    case_folder="$case_name"
    mkdir -p "$case_folder"                        # make a directory for the case study
    
    # extract bounding box for the current case
    bbox_gtsm=("${bbox[$case_name]}")
    
    # iterate through model_configs for each case
    for model_config in "${model_configs[@]}"; do
        model_folder="$case_folder/$model_config"
        mkdir -p "$model_folder"                   # make a directory for the model_config
        cd "$model_folder" || exit                 # move into the model_config directory
        
        # run GTSM
        if [ "$case_name" == "xynthia" ]; then
            mpirun --report-bindings -np 1 python /gpfs/home4/benitoli/papers/paper1/scripts/gtsm/general/gtsm_era5_template.py "$start_date" "${bbox_gtsm[@]}" "${model_configs[@]}"
        else
            mpirun --report-bindings -np 1 python /gpfs/home4/benitoli/papers/paper1/scripts/gtsm/general/gtsm_era5holland_template.py "$start_date" "${bbox_gtsm[@]}" "${model_configs[@]}" # - TO DO, add path to the track somewhere?
        fi
    done
done