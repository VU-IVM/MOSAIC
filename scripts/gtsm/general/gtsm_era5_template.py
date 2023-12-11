from functools import partial
from datetime import datetime, timedelta

from omuse.community.delft3d.interface import DFlowFM
from omuse.community.era5.interface import ERA5

from amuse.ext.grid_remappers import bilinear_2D_remapper, nearest_2D_remapper

from omuse.units import units
from matplotlib import pyplot
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import sys
import os
import numpy as np
import pandas as pd
from netCDF4 import num2date, Dataset, date2num

from time import perf_counter
bilinear_2D_remapper=partial(bilinear_2D_remapper, check_inside=False)

# --------------------------------------------------------------------------------------------------------------------------
#%% functions
# --------------------------------------------------------------------------------------------------------------------------
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

# --------------------------------------------------------------------------------------------------------------------------
#%% model inputs
# --------------------------------------------------------------------------------------------------------------------------

start_date=sys.argv[1]
start_datetime=datetime.strptime(start_date,"%Y%m%d")
tend=9. | units.day                                        # end period
dt=(1.| units.hour)                                        # timesteps
model_config= # TO DO
coordinates=sys.argv[2].split(',')      # bbox of netcdf file where water level's are being saved at grid level [lon min, lon max, lat min, lat max] 

# --------------------------------------------------------------------------------------------------------------------------
#%% defining the model
# --------------------------------------------------------------------------------------------------------------------------

# define the GTSM (Delft3d FM) model
d=DFlowFM(number_of_workers=32, ini_file="/gpfs/home4/benitoli/papers/paper1/scripts/gtsm/general/model_template/gtsm_fine.mdu", coordinates="spherical")

# set the model parameters
d.parameters.use_interface_wind=True
d.parameters.use_interface_patm=True
d.ini_time.RefDate=start_date
d.ini_physics.TidalForcing=1
d.ini_output.MapInterval=0.

if model_config=='TR':
    d.ini_output.HisInterval="600.  172800.  777600."  # this should be a system argument
    
if model_config=='OR':
    d.ini_output.ObsFile="/gpfs/home4/benitoli/papers/paper1/scripts/gtsm/output_resolution/observation_locations_snapped_1p25eu/selected_output_OR_haiyan_snapped_1p25eu_unique_obs.xyn"  # this should be a system argument

# point to the external forcings file and the partition file
d.ini_external_forcing.ExtForceFile="/gpfs/home4/benitoli/papers/paper1/scripts/gtsm/general/model_template/gtsm_fine.ext"
d.ini_geometry.PartitionFile="/gpfs/home4/benitoli/papers/paper1/scripts/gtsm/general/model_template/step11_global_part.pol"
d.ini_numerics.Icgsolver=6

# --------------------------------------------------------------------------------------------------------------------------
#%% defining meteorological forcing
# --------------------------------------------------------------------------------------------------------------------------

# setting up forcing - TO DO: make an if statement depending on the meteo forcing.. 
meteo=ERA5( variables=["10m_u_component_of_wind", 
                      "10m_v_component_of_wind",
                      "mean_sea_level_pressure"], 
            start_datetime=datetime.strptime(start_date,"%Y%m%d"),
            grid_resolution=0.25| units.deg,
            nwse_boundingbox=[90,-180,-90,180] | units.deg)

# prepare channels between the meteo forcing and the GTSM model
channel1=meteo.grid.new_remapping_channel_to(d.flow_links_forcing, bilinear_2D_remapper)
channel2=meteo.grid.new_remapping_channel_to(d.flow_nodes_forcing, bilinear_2D_remapper)
print('d.flow_nodes:', d.flow_nodes)
print('---preparing forcings---')

# --------------------------------------------------------------------------------------------------------------------------
#%% create map file only for the region of interest
# --------------------------------------------------------------------------------------------------------------------------

# create a file to save outputs in map format
try:
  os.mkdir('output')
except FileExistsError:
  print('Output directory already exists')
dest=Dataset(r'output/gtsm_map.nc', 'w', format='NETCDF4')  

# clip nodes and flow links to the bounding box - TO DO: Write this part of script more efficient
x_n= np.where((d.flow_nodes.lon.value_in(units.deg)>float(coordinates[0])) & (d.flow_nodes.lon.value_in(units.deg)<float(coordinates[1])) & (d.flow_nodes.lat.value_in(units.deg)<float(coordinates[3])) & (d.flow_nodes.lat.value_in(units.deg)>float(coordinates[2])), d.flow_nodes.lon.value_in(units.deg), np.nan)          # cropping lon,lat
x_n=x_n[~np.isnan(x_n)]
y_n= np.where((d.flow_nodes.lon.value_in(units.deg)>float(coordinates[0]))& (d.flow_nodes.lon.value_in(units.deg)<float(coordinates[1])) & (d.flow_nodes.lat.value_in(units.deg)<float(coordinates[3])) & (d.flow_nodes.lat.value_in(units.deg)>float(coordinates[2])), d.flow_nodes.lat.value_in(units.deg), np.nan)
y_n=y_n[~np.isnan(y_n)]

# define the netcdf file's coordinates
dest.createDimension('station',len(x_n))
lon=dest.createVariable('longitude',d.flow_nodes.lon.value_in(units.deg).dtype.str,('station',)) 
lon.units='degrees_east'
lon.standard_name='longitude'
lon[:]=x_n
lat=dest.createVariable('latitude',d.flow_nodes.lat.value_in(units.deg).dtype.str,('station',)) 
lat.units='degrees_north'
lat.standard_name='latitude'
lat[:]=y_n
# define the netcdf file's water levels
wl=dest.createVariable('waterlevel', np.float64,('station'))
wl.coordinates= 'longitude latitude'
wl.units='m'
wl.standard_name='sea_surface_height'

# save water level values at initial condition. This allows to remove dry cells later on.
wl_ini=dest.createVariable('waterlevel_ini', np.float64,('station'))
wl_ini.coordinates= 'longitude latitude'
wl_ini.units='m'
wl_ini.standard_name='initial_sea_surface_height'

# --------------------------------------------------------------------------------------------------------------------------
#%% execute model
# --------------------------------------------------------------------------------------------------------------------------
n=0
while d.model_time < tend-dt/2:
    n+=1

    # interpolate forcings
    time=-perf_counter()

    channel1.copy_attributes(["_10m_u_component_of_wind", 
                              "_10m_v_component_of_wind"], 
                              target_names=["wind_vx",
                                            "wind_vy"])
    channel2.copy_attributes(["_mean_sea_level_pressure"], 
                             target_names=["atmospheric_pressure"])

    time+=perf_counter()
    
    # evolve
    d.evolve_model(d.model_time+dt)
    meteo.evolve_model(d.model_time+dt)

    # -----------------------------------------------------------------------------------------------------------------------
    #%% save data to map file
    # -----------------------------------------------------------------------------------------------------------------------
    #Plotting waterlevels
    z=d.flow_nodes.water_level.value_in(units.m)
    z_bbox= np.where((d.flow_nodes.lon.value_in(units.deg)>float(coordinates[0])) & (d.flow_nodes.lon.value_in(units.deg)<float(coordinates[1])) & (d.flow_nodes.lat.value_in(units.deg)<float(coordinates[3])) & (d.flow_nodes.lat.value_in(units.deg)>float(coordinates[2])), z, np.nan)
    z_bbox=z_bbox[~np.isnan(z_bbox)]
    
    # save the initial water levels
    if d.model_time.value_in(units.s) == 3600.:
        wl_ini[:]=z_bbox
        wl[:]=wl_ini[:]
        #print('wl_ini[:]', wl_ini[:])
        
    # from second 172800 start saving water levels at grid level
    if 172800. < d.model_time.value_in(units.s) <= 518400. - dt.value_in(units.s):
        wl[:]=np.where(z_bbox>wl[:], z_bbox, wl[:])
        #print('wl[:]', wl[:])  

    '''
    #Plotting waterlevels
    z=d.flow_nodes.water_level.value_in(units.m)
    z_bbox= np.where((d.flow_nodes.lon.value_in(units.deg)>float(coordinates[0]))& (d.flow_nodes.lon.value_in(units.deg)<float(coordinates[1])) & (d.flow_nodes.lat.value_in(units.deg)<float(coordinates[3])) & (d.flow_nodes.lat.value_in(units.deg)> float(coordinates[2])), z, np.nan)
    z_bbox=z_bbox[~np.isnan(z_bbox)]
    print('waterlevel:', d.model_time, z_bbox.min(),z_bbox.max())
    name = "waterlevel" + str(i) + ".png"
    my_plot(x_n,y_n,z_bbox,-2.,2.,name)
    ''' 
    
d.stop()

