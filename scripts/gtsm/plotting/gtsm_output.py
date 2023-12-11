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
his_irma = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_tides/output/omuse_0000_his.nc').load() 
his_haiyan = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_TR_maxmap/output/omuse_0000_his.nc').load() 
his_xynthia = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/xynthia/era5_tides/output/omuse_0000_his.nc').load() 


# load gtsm map data
map_irma = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_tides/output/gtsm_map.nc').load() 
#map_irma = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/era5_tides/output/gtsm_map.nc').load() 
map_haiyan = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_TR_maxmap/output/gtsm_map.nc').load() 
map_xynthia = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/xynthia/era5_tides/output/gtsm_map.nc').load() 


#%%
# remove dry cells

map_irma_max = map_irma.max('time', skipna=True)
map_irma_min = map_irma.min('time', skipna=True)
map_irma_msk=map_irma_max.where(map_irma_max != map_irma_min, drop=True)

#map_haiyan_max = map_haiyan.max('time', skipna=True)
#map_haiyan_min = map_haiyan.min('time', skipna=True)
#map_haiyan_msk=map_haiyan_max.where(map_haiyan_max != map_haiyan_min, drop=True)
map_haiyan_msk=map_haiyan.where(map_haiyan.waterlevel != map_haiyan.waterlevel_ini, drop=True)
his_haiyan_max = his_haiyan.max('time', skipna=True)

map_xynthia_max = map_xynthia.max('time', skipna=True)
map_xynthia_min = map_xynthia.min('time', skipna=True)
map_xynthia_msk=map_xynthia_max.where(map_xynthia_max != map_xynthia_min, drop=True)

#%%

# cities for spatial reference in maps
jacksonville = np.asarray([-81.6,30.3])
tacloban = np.asarray([125.0,11.24])
rochelle = np.asarray([-1.2,46.2])
coast = cfeature.GSHHSFeature(scale='full')


fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(20, 5), subplot_kw={"projection": ccrs.PlateCarree()}, gridspec_kw = {'wspace':0.5, 'hspace':0.3})
ax1.set_title('Irma', fontsize=14, fontname='Arial')
ax1.plot(irma['LON'],irma['LAT'],'k',linestyle='--', linewidth=2.0,zorder=99)
a1=ax1.tripcolor(map_irma_msk.variables['longitude'][:],map_irma_msk.variables['latitude'][:],map_irma_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=4, zorder=1)
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
ax1.text(-0.25, 0.55, 'Baseline configuration', fontsize=12, va='center', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax1.transAxes)

ax2.set_title('Haiyan', fontsize=14, fontname='Arial')
ax2.plot(haiyan['LON'],haiyan['LAT'],'k',linestyle='--', linewidth=2.0,zorder=4,label='Storm track')
ax2.scatter(his_haiyan_max.station_x_coordinate,his_haiyan_max.station_y_coordinate,c=his_haiyan_max.waterlevel, marker='o', cmap=new_cmap,vmin=0,vmax=7,zorder=99,edgecolors='k')
a2=ax2.tripcolor(map_haiyan_msk.variables['longitude'][:],map_haiyan_msk.variables['latitude'][:],map_haiyan_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=7, zorder=1)
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
a3=ax3.tripcolor(map_xynthia_msk.variables['longitude'][:],map_xynthia_msk.variables['latitude'][:],map_xynthia_msk.variables['waterlevel'][:],cmap=new_cmap,vmin=0,vmax=6, zorder=1)
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

legend_elements= (Line2D([0],[0],color='k', linestyle='--', label='Storm track'),
                  Line2D([],[],color='k', marker='*', linestyle='None', markersize=8, label='Main cities'))#,
                  #Line2D([],[],color='k', marker='o', linestyle='None', markersize=8, markeredgewidth=1.5, markerfacecolor="none", label='GTSM Output stations'))
#lgd=plt.legend(handles=legend_elements,loc = ('best'), bbox_to_anchor=(1.3, 0.6), mode = "expand", prop={'size': 12, 'family':'Arial'}, frameon=False)
lgd=plt.legend(handles=legend_elements, ncol=3, loc = ('lower center'), bbox_to_anchor=(-1.0, -0.2), prop={'size': 12, 'family':'Arial'}, frameon=False)
#plt.legend(handles=legend_elements,loc = ('center right'), prop={'size': 12, 'family':'Arial'}, frameon=False)
#plt.setp(axs[-1, :], xlabel='x axis label')
#ax4 = plt.subplot2grid((4,38),(1,35),rowspan=3)
#cbar=plt.colorbar(b, fraction=0.046, pad=0.04, anchor=(2.5, 0.5))
#cbar.ax.set_ylabel('Maximum water level [m]', rotation=90, fontsize=12)
plt.savefig(f"figures/gtsm_TR.png", bbox_era_artists=lgd, bbox_inches='tight')
plt.clf()