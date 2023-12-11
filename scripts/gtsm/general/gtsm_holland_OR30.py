from datetime import datetime, timedelta

from omuse.community.delft3d.interface import DFlowFM
from omuse.ext.hurricane_models import HollandHurricane

from omuse.units import units
from matplotlib import pyplot

import sys
import os
import numpy as np
import pandas as pd
from netCDF4 import num2date, Dataset, date2num

# inputs
start_date=sys.argv[1]
start_datetime=datetime.strptime(start_date,"%Y%m%d")
tend=9. | units.day                                        # end period, I should do 9 afterwards!
dt=(30.| units.minute)                                        # timesteps
# input file for the TC track:
tc_track = sys.argv[3]
# create a file to save outputs in map format
try:
  os.mkdir('output')
except FileExistsError:
  print('Output directory already exists')
dest=Dataset(r'output/gtsm_map.nc', 'w', format='NETCDF4')  
# bbox to save files to [lon min, lon max, lat min, lat max]
#coordinates=[-84., -78., 23., 34.]                                                                                                     
coordinates=sys.argv[2].split(',')                                         

def my_plot(x,y,var,var_min,var_max,name):
    pyplot.scatter(x,y,c=var, cmap="jet", s=1., vmin=var_min, vmax=var_max)
    pyplot.xlim(-110, -50)
    pyplot.ylim(15,50)
    #pyplot.xlim(-110, -35)
    #pyplot.ylim(0,75)
    pyplot.xlabel("lon (deg)")
    pyplot.ylabel("lat (deg)")
    pyplot.colorbar()
    #pyplot.title(d.model_time)
    pyplot.savefig(name,format='png')
    pyplot.clf()

# create an array of dates
dates = np.arange(start_datetime + timedelta(days=2) + timedelta(hours=dt.value_in(units.hour)), start_datetime + timedelta(days=6), timedelta(hours=dt.value_in(units.hour)))

#d=DFlowFM( ini_file="gtsm_coarse.mdu", coordinates="spherical", redirection="none",channel_type="sockets") I commented this line because it was giving an error. I added the line below instead
d=DFlowFM(number_of_workers=32, ini_file="/gpfs/home4/benitoli/papers/paper1/scripts/gtsm/general/model_template/gtsm_fine.mdu", coordinates="spherical")
print(d.parameter_set_names())

#Use internal omuse forcings
d.parameters.use_interface_wind=True
d.parameters.use_interface_patm=True
d.ini_time.RefDate=start_date
d.ini_physics.TidalForcing=1
d.ini_output.MapInterval=0.
d.ini_output.HisInterval="1800.  172800.  777600."
d.ini_output.ObsFile="/gpfs/home4/benitoli/papers/paper1/scripts/gtsm/output_resolution/observation_locations_snapped_1p25eu/selected_output_OR_haiyan_snapped_1p25eu_unique_obs.xyn"
#d.ini_output.ObsFile="/gpfs/home4/benitoli/papers/paper1/scripts/gtsm/general/model_template/selected_output_OR_snapped_1p25eu_unique_obs_comas.xyn"
d.ini_external_forcing.ExtForceFile="/gpfs/home4/benitoli/papers/paper1/scripts/gtsm/general/model_template/gtsm_fine.ext"
d.ini_geometry.PartitionFile="/gpfs/home4/benitoli/papers/paper1/scripts/gtsm/general/model_template/step11_global_part.pol"
d.ini_numerics.Icgsolver=6

#Prescribe wind speed and pressure with omuse HollandHurricane model
tc_pres=HollandHurricane(d.flow_nodes,actf_file=tc_track)
tc_vel=HollandHurricane(d.flow_links,actf_file=tc_track)

pyplot.ion()
pyplot.show()
i=0

# prepare forcings
channel1=tc_vel.nodes.new_channel_to(d.flow_links_forcing)
channel2=tc_pres.nodes.new_channel_to(d.flow_nodes_forcing)
print('---preparing forcings---')
                                                                                              # set initial water level to 0. Would this also work if you have tides?
# --------------------------------------------------------------------------------------------------------------------------
#%% execute model
# --------------------------------------------------------------------------------------------------------------------------
# evolve the model in OMUSE
while d.model_time < tend-dt/2:
    i+=1
    # evolve model and apply forcings
    channel1.copy_attributes(["vx","vy"],target_names=["wind_vx","wind_vy"])
    channel2.copy_attributes(["pressure"],target_names=["atmospheric_pressure"])

    d.evolve_model(d.model_time+dt)

    tc_vel.evolve_model(d.model_time+dt) 
    tc_pres.evolve_model(d.model_time+dt) 
    
d.stop()
  
