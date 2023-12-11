# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 11:30:28 2022

@author: -
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 10:29:49 2022

@author: ibe202
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
cmap = plt.get_cmap('jet')
new_cmap = truncate_colormap(cmap, 0.2, 1.0)


#%% load track data
ibtracs = pd.read_csv(r'/gpfs/home4/benitoli/papers/paper1/data/gtsm_validation/ibtracs.csv', sep=';')
IRMAIB = ibtracs[ibtracs.NAME == 'IRMA']
data_to_extract = ['NAME','ISO_TIME','LAT','LON','WMO_WIND','WMO_PRES','USA_SSHS','USA_WIND','USA_PRES']
irma = IRMAIB[data_to_extract].copy()
HAIYANIB = ibtracs[ibtracs.NAME == 'HAIYAN']
data_to_extract = ['NAME','ISO_TIME','LAT','LON','WMO_WIND','WMO_PRES','USA_SSHS','USA_WIND','USA_PRES']
haiyan = HAIYANIB[data_to_extract].copy()
xynthia = pd.read_csv(r'/gpfs/home4/benitoli/papers/paper1/data/gtsm_validation/xynthia_track.csv', sep=';')


# load gtsm his data
his_irma_BS = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_BS/output/omuse_0000_his.nc').load() # I need to change it for IRMA!!
his_haiyan_BS = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_BS/output/omuse_0000_his.nc').load() 
his_xynthia_BS = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/xynthia/era5_BS/output/omuse_0000_his.nc').load() 

his_irma_TR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_TR/output/omuse_0000_his.nc').load() # I need to change it for IRMA!!
his_haiyan_TR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_TR/output/omuse_0000_his.nc').load() 
his_xynthia_TR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/xynthia/era5_TR/output/omuse_0000_his.nc').load() 

his_irma_OR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_OR/output/omuse_0000_his.nc').load() 
his_haiyan_OR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_OR_maxmap/output/omuse_0000_his.nc').load() 
his_xynthia_OR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/xynthia/era5_OR/output/omuse_0000_his.nc').load() 

# load gtsm map data
map_irma_BS = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_BS/output/gtsm_map.nc').load() # I need to change it for IRMA!!
#map_irma_BS = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/era5_tides/output/gtsm_map.nc').load() 
map_haiyan_BS = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_BS/output/gtsm_map.nc').load() 
map_xynthia_BS = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/xynthia/era5_BS/output/gtsm_map.nc').load() 

map_irma_TR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_TR/output/gtsm_map.nc').load()  # I need to change it for IRMA!!
map_haiyan_TR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_TR/output/gtsm_map.nc').load() 
map_xynthia_TR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/xynthia/era5_TR/output/gtsm_map.nc').load() 

map_irma_OR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_OR/output/gtsm_map.nc').load() 
map_haiyan_OR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_OR/output/gtsm_map.nc').load() 
map_xynthia_OR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/xynthia/era5_OR/output/gtsm_map.nc').load() 


#%%
# remove dry cells
#map_irma_BS_max = map_irma_BS.max('time', skipna=True)
#map_irma_BS_min = map_irma_BS.min('time', skipna=True)
#map_irma_BS_msk=map_irma_BS_max.where(map_irma_BS_max != map_irma_BS_min, drop=True)

#map_irma_BS_msk=map_irma_BS
map_irma_BS_msk=map_irma_BS.where(map_irma_BS.waterlevel != map_irma_BS.waterlevel_ini, drop=True)
his_irma_BS_max = his_irma_BS.max('time', skipna=True)

map_haiyan_BS_msk=map_haiyan_BS.where(map_haiyan_BS.waterlevel != map_haiyan_BS.waterlevel_ini, drop=True)
his_haiyan_BS_max = his_haiyan_BS.max('time', skipna=True)

map_xynthia_BS_msk=map_xynthia_BS.where(map_xynthia_BS.waterlevel != map_xynthia_BS.waterlevel_ini, drop=True)
his_xynthia_BS_max = his_xynthia_BS.max('time', skipna=True)

map_irma_TR_msk=map_irma_TR.where(map_irma_TR.waterlevel != map_irma_TR.waterlevel_ini, drop=True)
his_irma_TR_max = his_irma_TR.max('time', skipna=True)

map_haiyan_TR_msk=map_haiyan_TR.where(map_haiyan_TR.waterlevel != map_haiyan_TR.waterlevel_ini, drop=True)
his_haiyan_TR_max = his_haiyan_TR.max('time', skipna=True)

map_xynthia_TR_msk=map_xynthia_TR.where(map_xynthia_TR.waterlevel != map_xynthia_TR.waterlevel_ini, drop=True)
his_xynthia_TR_max = his_xynthia_TR.max('time', skipna=True)

map_irma_OR_msk=map_irma_OR.where(map_irma_OR.waterlevel != map_irma_OR.waterlevel_ini, drop=True)
his_irma_OR_max = his_irma_OR.max('time', skipna=True)

map_haiyan_OR_msk=map_haiyan_OR.where(map_haiyan_OR.waterlevel != map_haiyan_OR.waterlevel_ini, drop=True)
his_haiyan_OR_max = his_haiyan_OR.max('time', skipna=True)

map_xynthia_OR_msk=map_xynthia_OR.where(map_xynthia_OR.waterlevel != map_xynthia_OR.waterlevel_ini, drop=True)
his_xynthia_OR_max = his_xynthia_OR.max('time', skipna=True)


#%%
# cities for spatial reference in maps
jacksonville = np.asarray([-81.6,30.3])
tacloban = np.asarray([125.0,11.24])
rochelle = np.asarray([-1.2,46.2])
coast = cfeature.GSHHSFeature(scale='full')


fig, ((ax1, ax2, ax3),(ax4, ax5, ax6),(ax7, ax8, ax9), (ax10, ax11, ax12), (ax13, ax14, ax15)) = plt.subplots(5,3, figsize=(18, 25), subplot_kw={"projection": ccrs.PlateCarree()}, gridspec_kw = {'wspace':0.2, 'hspace':0.3})
ax1.set_title('Irma', fontsize=14, fontname='Arial')
ax1.plot(irma['LON'],irma['LAT'],'k',linestyle='--', linewidth=2.0,zorder=99)
ax1.scatter(his_irma_BS_max.station_x_coordinate,his_irma_BS_max.station_y_coordinate,c=his_irma_BS_max.waterlevel, marker='o', cmap=new_cmap,vmin=0.0,vmax=1.5,zorder=99,edgecolors='k')
#a1=ax1.scatter(map_irma_BS_msk.variables['longitude'][:],map_irma_BS_msk.variables['latitude'][:],c=map_irma_BS_msk.variables['waterlevel'][:], marker='o', cmap=new_cmap,vmin=0,vmax=4,zorder=99,edgecolors='k')
a1=ax1.tripcolor(map_irma_BS_msk.variables['longitude'][:],map_irma_BS_msk.variables['latitude'][:],map_irma_BS_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0.0,vmax=1.5, zorder=1)
#ax1.set_extent([-82.5, -79.0, 29.5, 33.0], crs=ccrs.PlateCarree())
ax1.set_extent([-82.2, -80.2, 30.0, 32.0], crs=ccrs.PlateCarree())
#ax1.set_extent([-85, -85+10*3/3.4, 24, 34], crs=ccrs.PlateCarree())
ax1.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
gl = ax1.gridlines(alpha=0, draw_labels=True)          
ax1.grid(False)
gl.top_labels=False   # suppress top labels
gl.right_labels=False # suppress right labels
gl.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
gl.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax1.scatter(jacksonville[0],jacksonville[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax1.text(-82.15,30.15,'Jacksonville',fontsize=12, family='Arial', zorder=99)
cbar1=fig.colorbar(a1,ax=ax1, fraction=0.046, pad=0.05)
cbar1.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)
ax1.text(-0.35, 0.55, 'Baseline configuration', fontsize=12, va='center', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax1.transAxes)

ax2.set_title('Haiyan', fontsize=14, fontname='Arial')
ax2.plot(haiyan['LON'],haiyan['LAT'],'k',linestyle='--', linewidth=2.0,zorder=4,label='Storm track')
ax2.scatter(his_haiyan_BS_max.station_x_coordinate,his_haiyan_BS_max.station_y_coordinate,c=his_haiyan_BS_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=7,zorder=99,edgecolors='k')
a2=ax2.tripcolor(map_haiyan_BS_msk.variables['longitude'][:],map_haiyan_BS_msk.variables['latitude'][:],map_haiyan_BS_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=7, zorder=1)
#ax2.set_extent([123, 126, 9.5, 12.9], crs=ccrs.PlateCarree())
ax2.set_extent([124.8, 125.3, 11.0, 11.5], crs=ccrs.PlateCarree())
ax2.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
gl2 = ax2.gridlines(alpha=0, draw_labels=True)
ax2.grid(False)
gl2.top_labels=False   # suppress top labels
gl2.right_labels=False # suppress right labels
gl2.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
gl2.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax2.scatter(tacloban[0],tacloban[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax2.text(124.85,11.22,'Tackloban',fontsize=12, family='Arial', zorder=99)
cbar2=fig.colorbar(a2,ax=ax2, fraction=0.046, pad=0.05)
cbar2.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)

ax3.set_title('Xynthia', fontsize=14, fontname='Arial')
ax3.plot(xynthia['LON'],xynthia['LAT'],'k',linestyle='--', linewidth=2.0,zorder=4)
ax3.scatter(his_xynthia_BS_max.station_x_coordinate,his_xynthia_BS_max.station_y_coordinate,c=his_xynthia_BS_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=5,zorder=99,edgecolors='k')
a3=ax3.tripcolor(map_xynthia_BS_msk.variables['longitude'][:],map_xynthia_BS_msk.variables['latitude'][:],map_xynthia_BS_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=5, zorder=1)
#ax3.set_extent([-3.5, -3.5+4*3/3.4, 43.5, 47.5], crs=ccrs.PlateCarree())
ax3.set_extent([-1.8, -0.8, 45.5, 46.5], crs=ccrs.PlateCarree())
ax3.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
gl3 = ax3.gridlines(alpha=0, draw_labels=True)
ax3.grid(False)
gl3.top_labels=False   # suppress top labels
gl3.right_labels=False # suppress right labels
gl3.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
gl3.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax3.scatter(rochelle[0],rochelle[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax3.text(-1.15,46.2,'La Rochelle',fontsize=12, family='Arial', zorder=99)
cbar3=fig.colorbar(a3,ax=ax3, fraction=0.046, pad=0.05)
cbar3.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)

#ax4.set_title('Irma', fontsize=14, fontname='Arial')
ax4.plot(irma['LON'],irma['LAT'],'k',linestyle='--', linewidth=2.0,zorder=99)
ax4.scatter(his_irma_TR_max.station_x_coordinate,his_irma_TR_max.station_y_coordinate,c=his_irma_TR_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=1.5,zorder=99,edgecolors='k')
#a4=ax4.scatter(map_irma_TR_msk.variables['longitude'][:],map_irma_TR_msk.variables['latitude'][:],c=map_irma_TR_msk.variables['waterlevel'][:], marker='o', cmap=new_cmap,vmin=0,vmax=4,zorder=99,edgecolors='k')
a4=ax4.tripcolor(map_irma_TR_msk.variables['longitude'][:],map_irma_TR_msk.variables['latitude'][:],map_irma_TR_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=1.5, zorder=1)
#ax1.set_extent([-82.5, -79.0, 29.5, 33.0], crs=ccrs.PlateCarree())
ax4.set_extent([-82.2, -80.2, 30.0, 32.0], crs=ccrs.PlateCarree())
#ax1.set_extent([-85, -85+10*3/3.4, 24, 34], crs=ccrs.PlateCarree())
ax4.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
g4 = ax4.gridlines(alpha=0, draw_labels=True)          
ax4.grid(False)
g4.top_labels=False   # suppress top labels
g4.right_labels=False # suppress right labels
g4.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
g4.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax4.scatter(jacksonville[0],jacksonville[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax4.text(-82.15,30.15,'Jacksonville',fontsize=12, family='Arial', zorder=99)
cbar4=fig.colorbar(a4,ax=ax4, fraction=0.046, pad=0.05)
cbar4.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)
ax4.text(-0.35, -0.75, 'Temporal resolution ref.', fontsize=12, va='center', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax1.transAxes)
        
#ax5.set_title('Haiyan', fontsize=14, fontname='Arial')
ax5.plot(haiyan['LON'],haiyan['LAT'],'k',linestyle='--', linewidth=2.0,zorder=4,label='Storm track')
ax5.scatter(his_haiyan_TR_max.station_x_coordinate,his_haiyan_TR_max.station_y_coordinate,c=his_haiyan_TR_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=7,zorder=99,edgecolors='k')
a5=ax5.tripcolor(map_haiyan_TR_msk.variables['longitude'][:],map_haiyan_TR_msk.variables['latitude'][:],map_haiyan_TR_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=7, zorder=1)
#ax2.set_extent([123, 126, 9.5, 12.9], crs=ccrs.PlateCarree())
ax5.set_extent([124.8, 125.3, 11.0, 11.5], crs=ccrs.PlateCarree())
ax5.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
gl5 = ax5.gridlines(alpha=0, draw_labels=True)
ax5.grid(False)
gl5.top_labels=False   # suppress top labels
gl5.right_labels=False # suppress right labels
gl5.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
gl5.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax5.scatter(tacloban[0],tacloban[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax5.text(124.85,11.22,'Tackloban',fontsize=12, family='Arial', zorder=99)
cbar5=fig.colorbar(a5,ax=ax5, fraction=0.046, pad=0.05)
cbar5.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)

#ax6.set_title('Xynthia', fontsize=14, fontname='Arial')
ax6.plot(xynthia['LON'],xynthia['LAT'],'k',linestyle='--', linewidth=2.0,zorder=4)
ax6.scatter(his_xynthia_TR_max.station_x_coordinate,his_xynthia_TR_max.station_y_coordinate,c=his_xynthia_TR_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=5,zorder=99,edgecolors='k')
a6=ax6.tripcolor(map_xynthia_TR_msk.variables['longitude'][:],map_xynthia_TR_msk.variables['latitude'][:],map_xynthia_TR_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=5, zorder=1)
#ax3.set_extent([-3.5, -3.5+4*3/3.4, 43.5, 47.5], crs=ccrs.PlateCarree())
ax6.set_extent([-1.8, -0.8, 45.5, 46.5], crs=ccrs.PlateCarree())
ax6.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
gl6 = ax6.gridlines(alpha=0, draw_labels=True)
ax6.grid(False)
gl6.top_labels=False   # suppress top labels
gl6.right_labels=False # suppress right labels
gl6.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
gl6.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax6.scatter(rochelle[0],rochelle[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax6.text(-1.15,46.2,'La Rochelle',fontsize=12, family='Arial', zorder=99)
cbar6=fig.colorbar(a6,ax=ax6, fraction=0.046, pad=0.05)
cbar6.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)

#ax4.set_title('Irma', fontsize=14, fontname='Arial')
ax7.plot(irma['LON'],irma['LAT'],'k',linestyle='--', linewidth=2.0,zorder=99)
#ax7.scatter(his_irma_BS_max.station_x_coordinate,his_irma_BS_max.station_y_coordinate,c=his_irma_BS_max.waterlevel, marker='v', cmap=new_cmap,vmin=0.5,vmax=1.5,zorder=99,edgecolors='k')
ax7.scatter(his_irma_OR_max.station_x_coordinate,his_irma_OR_max.station_y_coordinate,c=his_irma_OR_max.waterlevel, marker='o', cmap=new_cmap,vmin=0.0,vmax=1.5,zorder=99,edgecolors='k')
#ax7.scatter(his_irma_OR_max.station_x_coordinate,his_irma_OR_max.station_y_coordinate,c='none', marker='o', cmap=new_cmap,vmin=0.5,vmax=1.5,zorder=99,edgecolors='k')
#a7=ax7.scatter(map_irma_OR_msk.variables['longitude'][:],map_irma_OR_msk.variables['latitude'][:],c=map_irma_OR_msk.variables['waterlevel'][:], marker='o', cmap=new_cmap,vmin=0,vmax=4,zorder=99,edgecolors='k')
a7=ax7.tripcolor(map_irma_OR_msk.variables['longitude'][:],map_irma_OR_msk.variables['latitude'][:],map_irma_OR_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0.0,vmax=1.5, zorder=1)
#ax7.set_extent([-81.5, -80.5, 31, 32], crs=ccrs.PlateCarree())
ax7.set_extent([-82.2, -80.2, 30.0, 32.0], crs=ccrs.PlateCarree())
ax7.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
g7 = ax7.gridlines(alpha=0, draw_labels=True)          
ax7.grid(False)
g7.top_labels=False   # suppress top labels
g7.right_labels=False # suppress right labels
g7.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
g7.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax7.scatter(jacksonville[0],jacksonville[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax7.text(-82.15,30.15,'Jacksonville',fontsize=12, family='Arial', zorder=99)
cbar7=fig.colorbar(a7,ax=ax7, fraction=0.046, pad=0.05)
cbar7.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)
ax7.text(-0.35, -2.05, 'Output resolution ref.', fontsize=12, va='center', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax1.transAxes)
        
#ax5.set_title('Haiyan', fontsize=14, fontname='Arial')
ax8.plot(haiyan['LON'],haiyan['LAT'],'k',linestyle='--', linewidth=2.0,zorder=4,label='Storm track')
ax8.scatter(his_haiyan_OR_max.station_x_coordinate,his_haiyan_OR_max.station_y_coordinate,c=his_haiyan_OR_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=7,zorder=99,edgecolors='k')
a8=ax8.tripcolor(map_haiyan_OR_msk.variables['longitude'][:],map_haiyan_OR_msk.variables['latitude'][:],map_haiyan_OR_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=7, zorder=1)
#ax2.set_extent([123, 126, 9.5, 12.9], crs=ccrs.PlateCarree())
ax8.set_extent([124.8, 125.3, 11.0, 11.5], crs=ccrs.PlateCarree())
ax8.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
gl8 = ax8.gridlines(alpha=0, draw_labels=True)
ax8.grid(False)
gl8.top_labels=False   # suppress top labels
gl8.right_labels=False # suppress right labels
gl8.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
gl8.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax8.scatter(tacloban[0],tacloban[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax8.text(124.85,11.22,'Tackloban',fontsize=12, family='Arial', zorder=99)
cbar8=fig.colorbar(a8,ax=ax8, fraction=0.046, pad=0.05)
cbar8.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)

#ax6.set_title('Xynthia', fontsize=14, fontname='Arial')
ax9.plot(xynthia['LON'],xynthia['LAT'],'k',linestyle='--', linewidth=2.0,zorder=4)
ax9.scatter(his_xynthia_OR_max.station_x_coordinate,his_xynthia_OR_max.station_y_coordinate,c=his_xynthia_OR_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=5,zorder=99,edgecolors='k')
a9=ax9.tripcolor(map_xynthia_OR_msk.variables['longitude'][:],map_xynthia_OR_msk.variables['latitude'][:],map_xynthia_OR_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=5, zorder=1)
#ax3.set_extent([-3.5, -3.5+4*3/3.4, 43.5, 47.5], crs=ccrs.PlateCarree())
ax9.set_extent([-1.8, -0.8, 45.5, 46.5], crs=ccrs.PlateCarree())
ax9.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
gl9 = ax9.gridlines(alpha=0, draw_labels=True)
ax9.grid(False)
gl9.top_labels=False   # suppress top labels
gl9.right_labels=False # suppress right labels
gl9.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
gl9.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax9.scatter(rochelle[0],rochelle[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax9.text(-1.15,46.2,'La Rochelle',fontsize=12, family='Arial', zorder=99)
cbar9=fig.colorbar(a9,ax=ax9, fraction=0.046, pad=0.05)
cbar9.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)

#ax10.set_title('Irma', fontsize=12, fontname='Arial')
ax10.plot(irma['LON'],irma['LAT'],'k',linestyle='--', linewidth=2.0,zorder=99)
#ax10.scatter(his_irma_BS_max.station_x_coordinate,his_irma_BS_max.station_y_coordinate,c=his_irma_BS_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=2,zorder=99,edgecolors='k')
#a10=ax10.scatter(map_irma_BS_msk.variables['longitude'][:],map_irma_BS_msk.variables['latitude'][:],c=map_irma_BS_msk.variables['waterlevel'][:], marker='o', cmap=new_cmap,vmin=0,vmax=4,zorder=99,edgecolors='k')
#a10=ax10.tripcolor(map_irma_BS_msk.variables['longitude'][:],map_irma_BS_msk.variables['latitude'][:],map_irma_BS_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=2, zorder=1)
#ax10.set_extent([-82.5, -79.0, 29.5, 33.0], crs=ccrs.PlateCarree())
ax10.set_extent([-82.2, -80.2, 30.0, 32.0], crs=ccrs.PlateCarree())
#ax10.set_extent([-85, -85+10*3/3.4, 24, 34], crs=ccrs.PlateCarree())
ax10.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
gl10 = ax10.gridlines(alpha=0, draw_labels=True)          
ax10.grid(False)
gl10.top_labels=False   # suppress top labels
gl10.right_labels=False # suppress right labels
gl10.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
gl10.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax10.scatter(jacksonville[0],jacksonville[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax10.text(-82.15,30.15,'Jacksonville',fontsize=12, family='Arial', zorder=99)
cbar10=fig.colorbar(a1,ax=ax10, fraction=0.046, pad=0.05)
cbar10.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)
ax10.text(-0.35, -3.35, 'Grid resolution ref.', fontsize=12, va='center', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax1.transAxes)

#ax11.set_title('Haiyan', fontsize=12, fontname='Arial')
ax11.plot(haiyan['LON'],haiyan['LAT'],'k',linestyle='--', linewidth=2.0,zorder=4,label='Storm track')
#ax11.scatter(his_haiyan_BS_max.station_x_coordinate,his_haiyan_BS_max.station_y_coordinate,c=his_haiyan_BS_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=7,zorder=99,edgecolors='k')
#a11=ax11.tripcolor(map_haiyan_BS_msk.variables['longitude'][:],map_haiyan_BS_msk.variables['latitude'][:],map_haiyan_BS_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=7, zorder=1)
#ax11.set_extent([123, 126, 9.5, 12.9], crs=ccrs.PlateCarree())
ax11.set_extent([124.8, 125.3, 11.0, 11.5], crs=ccrs.PlateCarree())
ax11.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
gl11 = ax11.gridlines(alpha=0, draw_labels=True)
ax11.grid(False)
gl11.top_labels=False   # suppress top labels
gl11.right_labels=False # suppress right labels
gl11.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
gl11.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax11.scatter(tacloban[0],tacloban[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax11.text(124.85,11.22,'Tackloban',fontsize=12, family='Arial', zorder=99)
cbar11=fig.colorbar(a2,ax=ax11, fraction=0.046, pad=0.05)
cbar11.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)

#ax12.set_title('Xynthia', fontsize=12, fontname='Arial')
ax12.plot(xynthia['LON'],xynthia['LAT'],'k',linestyle='--', linewidth=2.0,zorder=4)
#ax12.scatter(his_xynthia_BS_max.station_x_coordinate,his_xynthia_BS_max.station_y_coordinate,c=his_xynthia_BS_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=5,zorder=99,edgecolors='k')
#a12=ax12.tripcolor(map_xynthia_BS_msk.variables['longitude'][:],map_xynthia_BS_msk.variables['latitude'][:],map_xynthia_BS_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=5, zorder=1)
#ax12.set_extent([-3.5, -3.5+4*3/3.4, 43.5, 47.5], crs=ccrs.PlateCarree())
ax12.set_extent([-1.8, -0.8, 45.5, 46.5], crs=ccrs.PlateCarree())
ax12.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
gl12 = ax12.gridlines(alpha=0, draw_labels=True)
ax12.grid(False)
gl12.top_labels=False   # suppress top labels
gl12.right_labels=False # suppress right labels
gl12.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
gl12.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax12.scatter(rochelle[0],rochelle[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax12.text(-1.15,46.2,'La Rochelle',fontsize=12, family='Arial', zorder=99)
cbar12=fig.colorbar(a3,ax=ax12, fraction=0.046, pad=0.05)
cbar12.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)

#ax13.set_title('Irma', fontsize=12, fontname='Arial')
ax13.plot(irma['LON'],irma['LAT'],'k',linestyle='--', linewidth=2.0,zorder=99)
#ax13.scatter(his_irma_BS_max.station_x_coordinate,his_irma_BS_max.station_y_coordinate,c=his_irma_BS_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=2,zorder=99,edgecolors='k')
#a13=ax13.scatter(map_irma_BS_msk.variables['longitude'][:],map_irma_BS_msk.variables['latitude'][:],c=map_irma_BS_msk.variables['waterlevel'][:], marker='o', cmap=new_cmap,vmin=0,vmax=4,zorder=99,edgecolors='k')
#a13=ax13.tripcolor(map_irma_BS_msk.variables['longitude'][:],map_irma_BS_msk.variables['latitude'][:],map_irma_BS_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=2, zorder=1)
#ax13.set_extent([-82.5, -79.0, 29.5, 33.0], crs=ccrs.PlateCarree())
ax13.set_extent([-82.2, -80.2, 30.0, 32.0], crs=ccrs.PlateCarree())
#ax13.set_extent([-85, -85+10*3/3.4, 24, 34], crs=ccrs.PlateCarree())
ax13.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
gl13 = ax13.gridlines(alpha=0, draw_labels=True)          
ax13.grid(False)
gl13.top_labels=False   # suppress top labels
gl13.right_labels=False # suppress right labels
gl13.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
gl13.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax13.scatter(jacksonville[0],jacksonville[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax13.text(-82.15,30.15,'Jacksonville',fontsize=12, family='Arial', zorder=99)
cbar13=fig.colorbar(a1,ax=ax13, fraction=0.046, pad=0.05)
cbar13.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)
ax13.text(-0.35, -4.65, 'Bed level ref.', fontsize=12, va='center', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax1.transAxes)

#ax14.set_title('Haiyan', fontsize=12, fontname='Arial')
ax14.plot(haiyan['LON'],haiyan['LAT'],'k',linestyle='--', linewidth=2.0,zorder=4,label='Storm track')
#ax14.scatter(his_haiyan_BS_max.station_x_coordinate,his_haiyan_BS_max.station_y_coordinate,c=his_haiyan_BS_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=7,zorder=99,edgecolors='k')
#a14=ax14.tripcolor(map_haiyan_BS_msk.variables['longitude'][:],map_haiyan_BS_msk.variables['latitude'][:],map_haiyan_BS_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=7, zorder=1)
#ax14.set_extent([123, 126, 9.5, 12.9], crs=ccrs.PlateCarree())
ax14.set_extent([124.8, 125.3, 11.0, 11.5], crs=ccrs.PlateCarree())
ax14.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
gl14 = ax14.gridlines(alpha=0, draw_labels=True)
ax14.grid(False)
gl14.top_labels=False   # suppress top labels
gl14.right_labels=False # suppress right labels
gl14.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
gl14.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax14.scatter(tacloban[0],tacloban[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax14.text(124.85,11.22,'Tackloban',fontsize=12, family='Arial', zorder=99)
cbar14=fig.colorbar(a2,ax=ax14, fraction=0.046, pad=0.05)
cbar14.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)

#ax15.set_title('Xynthia', fontsize=12, fontname='Arial')
ax15.plot(xynthia['LON'],xynthia['LAT'],'k',linestyle='--', linewidth=2.0,zorder=4)
#ax15.scatter(his_xynthia_BS_max.station_x_coordinate,his_xynthia_BS_max.station_y_coordinate,c=his_xynthia_BS_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=5,zorder=99,edgecolors='k')
#a15=ax15.tripcolor(map_xynthia_BS_msk.variables['longitude'][:],map_xynthia_BS_msk.variables['latitude'][:],map_xynthia_BS_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=5, zorder=1)
#ax15.set_extent([-3.5, -3.5+4*3/3.4, 43.5, 47.5], crs=ccrs.PlateCarree())
ax15.set_extent([-1.8, -0.8, 45.5, 46.5], crs=ccrs.PlateCarree())
ax15.add_feature(coast, facecolor='lightgray', linewidth=0)#, zorder 99)
gl15 = ax15.gridlines(alpha=0, draw_labels=True)
ax15.grid(False)
gl15.top_labels=False   # suppress top labels
gl15.right_labels=False # suppress right labels
gl15.xlabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
gl15.ylabel_style = {'color': 'k', 'size':12, 'fontname':'Arial'}
ax15.scatter(rochelle[0],rochelle[1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)
ax15.text(-1.15,46.2,'La Rochelle',fontsize=12, family='Arial', zorder=99)
cbar15=fig.colorbar(a3,ax=ax15, fraction=0.046, pad=0.05)
cbar15.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)

legend_elements= (Line2D([0],[0],color='k', linestyle='--', label='Storm track'),
                  Line2D([],[],color='k', marker='*', linestyle='None', markersize=8, label='Main cities'))#,
                  #Line2D([],[],color='k', marker='o', linestyle='None', markersize=8, markeredgewidth=1.5, markerfacecolor="none", label='GTSM Output stations'))
#lgd=plt.legend(handles=legend_elements,loc = ('best'), bbox_to_anchor=(1.3, 0.6), mode = "expand", prop={'size': 12, 'family':'Arial'}, frameon=False)
#lgd=plt.legend(handles=legend_elements, ncol=3, loc = ('lower center'), bbox_to_anchor=(-1.0, -0.2), prop={'size': 12, 'family':'Arial'}, frameon=False)
plt.savefig(f"figures/gtsm_output_v4.png")
#plt.savefig(f"figures/gtsm_BS.png", bbox_era_artists=lgd, bbox_inches='tight')
plt.clf()