import sys
import geopandas as gpd
from hydromt_sfincs import SfincsModel
from hydromt.log import setuplog
import xarray as xr
import numpy as np
import geopandas as gpd
import pandas as pd


case=sys.argv[1]
meteo=sys.argv[2]
model_config=sys.argv[3]
event_bbox=[float(coord) for coord in sys.argv[4].split(',')]
tref=sys.argv[5]
tstart=sys.argv[5]
tstop=sys.argv[6]

print('event_bbox:', type(event_bbox[0]))

logger = setuplog('log_' + str(case) + '_' + str(meteo) + '_' + str(model_config))
path_to_yml=r'/projects/0/einf2224/paper1/scripts/hydromt-sfincs/general/data_catalog_v2.yml'

root = str(case) + '_' + str(meteo) + '_' + str(model_config)
mod = SfincsModel(root, mode='w+', data_libs=path_to_yml, logger=logger)

mod.setup_grid_from_region(region={'bbox': [124.5,10.5,125.5,11.5]}, res=150, rotated=False)

mod.setup_config(
    tref = tref,
    tstart = tstart,
    tstop = tstop,
    dtmaxout = 99999.0,
    dtout = 99999.0,
    dtwnd = 600.0,
    alpha = 0.5,
    zsini = 0.5,
    advection = 0, # before I had 0.0
    huthresh = 0.05,
)

datasets_dep = [{"elevtn":"fabdem_{case}", "zmin":0}, {"elevtn": "gebco"}]

mod.setup_dep(datasets_dep = datasets_dep)

mod.setup_mask_active(
    zmin = 0,                    # minimum elevation for valid cells
    exclude_mask = "osm_coastlines",
    drop_area = 1, # drops areas that are smaller than 1km2
)

mod.setup_mask_bounds(
    btype = "waterlevel",
    zmax = 10,
)

datasets_rgh = [{"lulc": "vito"}]
'''
mod.setup_subgrid(
    datasets_dep = datasets_dep,
    datasets_rgh=datasets_rgh,
    nr_subgrid_pixels=8, 
    write_dep_tif=True,
    write_man_tif=True,
)
'''

## Total water levels (GTSM)
# Import GTSM
gtsm_file=f'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/{case}/{meteo}_{model_config}/output/omuse_0000_his.nc'
gtsm = xr.open_dataset(gtsm_file).load()
# Optain the boundary point coordinates from GTSM
bnd = gpd.GeoDataFrame(
    index=np.atleast_1d(gtsm['stations'].values),
    geometry=gpd.points_from_xy(
        np.atleast_1d(gtsm['station_x_coordinate'].values), 
        np.atleast_1d(gtsm['station_y_coordinate'].values)
    ),
    crs=4326
).to_crs(mod.crs)

# Create a pandas dataframe to create the water level forcing   
df_timeseries=pd.DataFrame(index=gtsm.time, columns=bnd.index, data=gtsm.waterlevel)

mod.setup_waterlevel_forcing(
    timeseries=df_timeseries,
    locations=bnd,
    offset="dtu10mdt",
    merge=False,
)

mod.setup_observation_points(
    locations=f"/projects/0/einf2224/paper1/scripts/hydromt-sfincs/general/obs_points_{case}.geojson", merge=True
)

mod._write_gis = False
mod.write()