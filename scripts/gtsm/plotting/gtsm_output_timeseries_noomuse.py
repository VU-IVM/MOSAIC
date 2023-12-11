# -*- coding: utf-8 -*-
"""
Created on Wed Sept 28 11:30:28 2022

@author: -
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import matplotlib.dates as mdates

his_irma_BS = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_BS/output/omuse_0000_his.nc').load()
his_irma_TR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_TR/output/omuse_0000_his.nc').load()
his_irma_OR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_OR/output/omuse_0000_his.nc').load()
his_irma_era5 = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/era5_tides/output/omuse_0000_his.nc').load()
his_irma_holland_era5_BS = xr.open_dataset(r'/gpfs/home4/benitoli/gtsm_d3dfm/holland_era5_BS_irma/output/gtsm_fine_0000_his.nc').load()
his_irma_holland_era5_TR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/no_omuse/holland_era5_TR/output/gtsm_fine_0000_his.nc').load()

his_haiyan_BS = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_BS/output/omuse_0000_his.nc').load()
his_haiyan_TR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_TR/output/omuse_0000_his.nc').load()
his_haiyan_OR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_OR/output/omuse_0000_his.nc').load()

his_xynthia_BS = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/xynthia/era5_BS/output/omuse_0000_his.nc').load()
his_xynthia_TR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/xynthia/era5_TR/output/omuse_0000_his.nc').load()
his_xynthia_OR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/xynthia/era5_OR/output/omuse_0000_his.nc').load()



# Irma
irma_mask_lon = (his_irma_BS.station_x_coordinate >= -81.4) & (his_irma_BS.station_x_coordinate <= -81.3)
irma_mask_lat = (his_irma_BS.station_y_coordinate >= 31.0) & (his_irma_BS.station_y_coordinate <= 31.2)
his_irma_BS = his_irma_BS.where(irma_mask_lon & irma_mask_lat, drop = True)
his_irma_TR = his_irma_TR.where(irma_mask_lon & irma_mask_lat, drop = True)
his_irma_era5 = his_irma_era5.where(irma_mask_lon & irma_mask_lat, drop = True)
#print(his_irma_holland_era5)
irma_mask_lon2 = (his_irma_holland_era5_BS.station_x_coordinate >= -81.4) & (his_irma_holland_era5_BS.station_x_coordinate <= -81.3)
irma_mask_lat2 = (his_irma_holland_era5_BS.station_y_coordinate >= 31.0) & (his_irma_holland_era5_BS.station_y_coordinate <= 31.2)
his_irma_holland_era5_BS = his_irma_holland_era5_BS.where(irma_mask_lon2 & irma_mask_lat2, drop = True)
his_irma_holland_era5_TR = his_irma_holland_era5_TR.where(irma_mask_lon2 & irma_mask_lat2, drop = True)

start_date='2017-09-08T00:00:00.000000000'
end_date='2017-09-15T00:00:00.000000000'

his_irma_BS = his_irma_BS.sel(
    time=slice(start_date, end_date))
    
his_irma_TR = his_irma_TR.sel(
    time=slice(start_date, end_date))

his_irma_era5 = his_irma_era5.sel(
    time=slice(start_date, end_date))
    
his_irma_holland_era5_BS = his_irma_holland_era5_BS.sel(
    time=slice(start_date, end_date))

his_irma_holland_era5_TR = his_irma_holland_era5_TR.sel(
    time=slice(start_date, end_date))

# Haiyan    
haiyan_mask_lon = (his_haiyan_BS.station_x_coordinate >= 125) & (his_haiyan_BS.station_x_coordinate <= 125.2)
haiyan_mask_lat = (his_haiyan_BS.station_y_coordinate >= 11.2) & (his_haiyan_BS.station_y_coordinate <= 11.3)
his_haiyan_BS = his_haiyan_BS.where(haiyan_mask_lon & haiyan_mask_lat, drop = True)
his_haiyan_TR = his_haiyan_TR.where(haiyan_mask_lon & haiyan_mask_lat, drop = True)

start_date='2013-11-06T00:00:00.000000000'
end_date='2013-11-13T00:00:00.000000000'

his_haiyan_BS = his_haiyan_BS.sel(
    time=slice(start_date, end_date))
    
his_haiyan_TR = his_haiyan_TR.sel(
    time=slice(start_date, end_date))

#print('haiyan', his_haiyan_BS)

# Xynthia 
xynthia_mask_lon = (his_xynthia_BS.station_x_coordinate >= -1.5) & (his_xynthia_BS.station_x_coordinate <= -1.45)
xynthia_mask_lat = (his_xynthia_BS.station_y_coordinate >= 46.34) & (his_xynthia_BS.station_y_coordinate <= 46.4)   
#xynthia_mask_lon = (his_xynthia_BS.station_x_coordinate >= -1.27) & (his_xynthia_BS.station_x_coordinate <= -1.23)
#xynthia_mask_lat = (his_xynthia_BS.station_y_coordinate >= 46.25) & (his_xynthia_BS.station_y_coordinate <= 46.28)
his_xynthia_BS = his_xynthia_BS.where(xynthia_mask_lon & xynthia_mask_lat, drop = True)
his_xynthia_TR = his_xynthia_TR.where(xynthia_mask_lon & xynthia_mask_lat, drop = True)

start_date='2010-02-26T00:00:00.000000000'
end_date='2010-03-05T00:00:00.000000000'

his_xynthia_BS = his_xynthia_BS.sel(
    time=slice(start_date, end_date))
    
his_xynthia_TR = his_xynthia_TR.sel(
    time=slice(start_date, end_date))

#print('xynthia', his_xynthia_BS)

#fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, figsize=(20, 20))
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 4))
#plt.rcParams["figure.figsize"] = (15,5.5)

ax1.plot(his_irma_era5.time,his_irma_era5.waterlevel, 'g', linewidth=1.5, label='ERA5')
ax1.plot(his_irma_holland_era5_BS.time,his_irma_holland_era5_BS.waterlevel, 'r',linewidth=1.5, label='ERA5+Holland')

#ax1.plot(his_irma_TR.time,his_irma_TR.waterlevel,  'orange', linewidth=1.5, label='Holland Temporal resolution refined')
ax1.plot(his_irma_BS.time,his_irma_BS.waterlevel, 'b', linewidth=1.5, label='Holland')
#ax1.plot(his_irma_holland_era5_TR.time,his_irma_holland_era5_TR.waterlevel, 'm',linewidth=1.5, label='ERA5+Holland Temporal resolution refined')
ax1.set_title('Irma [lon: -81.313, lat: 31.187]')
ax1.tick_params(axis='x')
ax1.tick_params(axis='y')
ax1.xaxis.set_major_formatter(
    mdates.ConciseDateFormatter(ax1.xaxis.get_major_locator()))
ax1.legend()
#plt.legend(loc = ('lower center'), bbox_to_anchor=(-0.75, -0.3), ncol=3, frameon=False)
ax2.plot(his_haiyan_TR.time,his_haiyan_TR.waterlevel, 'orange', linewidth=1.5, label='Temporal resolution refined')
ax2.plot(his_haiyan_BS.time,his_haiyan_BS.waterlevel,'b', linewidth=1.5, label='Baseline configuration')
ax2.set_title('Haiyan [lon: 125.083, lat: 11.265]')
ax2.tick_params(axis='x')
ax2.tick_params(axis='y')
ax2.xaxis.set_major_formatter(
    mdates.ConciseDateFormatter(ax2.xaxis.get_major_locator()))
ax3.plot(his_xynthia_TR.time,his_xynthia_TR.waterlevel, 'orange', linewidth=1.5, label='Temporal resolution refined')
ax3.plot(his_xynthia_BS.time,his_xynthia_BS.waterlevel, 'b', linewidth=1.5,label='Baseline configuration')
ax3.set_title('Xynthia [lon: -1.472, lat: 46.384]')
ax3.tick_params(axis='x')
ax3.tick_params(axis='y')
ax3.xaxis.set_major_formatter(
    mdates.ConciseDateFormatter(ax3.xaxis.get_major_locator()))
ax1.set_ylabel('Water level [m]')
ax2.set_ylabel('Water level [m]')
ax3.set_ylabel('Water level [m]')
plt.legend(loc = ('lower center'), bbox_to_anchor=(-0.75, -0.3), ncol=3, frameon=False)
plt.gcf().autofmt_xdate()
plt.savefig(f"figures/gtsm_timeseries_irmaERA5_v2.png")


