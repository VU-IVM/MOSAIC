import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import xarray as xr
'''
# GTSM data
# load gtsm his data
his_irma_BS = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_BS/output/omuse_0000_his.nc').load() # I need to change it for IRMA!!
his_haiyan_BS = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_BS/output/omuse_0000_his.nc').load() 
his_xynthia_BS = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/xynthia/era5_BS/output/omuse_0000_his.nc').load() 

his_irma_TR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_TR/output/omuse_0000_his.nc').load() # I need to change it for IRMA!!
his_haiyan_TR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_TR/output/omuse_0000_his.nc').load() 
his_xynthia_TR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/xynthia/era5_TR/output/omuse_0000_his.nc').load() 

his_irma_OR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/irma/holland_OR/output/omuse_0000_his.nc').load() 
his_haiyan_OR = xr.open_dataset(r'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/haiyan/holland_OR/output/omuse_0000_his.nc').load() 
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


# remove dry cells & calculate maximum water levels
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

map_xynthia_TR_msk=map_xynthia_TR.where(map_xynthia_BS.waterlevel != map_xynthia_BS.waterlevel_ini, drop=True)
#map_xynthia_TR_msk=map_xynthia_TR.where(map_xynthia_TR.waterlevel != map_xynthia_TR.waterlevel_ini, drop=True)
his_xynthia_TR_max = his_xynthia_TR.max('time', skipna=True)

map_irma_OR_msk=map_irma_OR.where(map_irma_OR.waterlevel != map_irma_OR.waterlevel_ini, drop=True)
his_irma_OR_max = his_irma_OR.max('time', skipna=True)

map_haiyan_OR_msk=map_haiyan_OR.where(map_haiyan_OR.waterlevel != map_haiyan_OR.waterlevel_ini, drop=True)
his_haiyan_OR_max = his_haiyan_OR.max('time', skipna=True)

map_xynthia_OR_msk=map_xynthia_OR.where(map_xynthia_OR.waterlevel != map_xynthia_OR.waterlevel_ini, drop=True)
his_xynthia_OR_max = his_xynthia_OR.max('time', skipna=True)

# Different extents for each row
extents = [
    [-82.2, -81.2, 30., 31.0],
    [124.8, 125.3, 11.0, 11.5],
    [-1.8, -0.8, 45.5, 46.5]
]

# Different colormaps for each column
cmaps = ['jet', 'bwr', 'bwr']

#def plot_waterlevel(gtsm, cmap):
def plot_waterlevel(ax, gtsm, cmap):
   # tripcolor = ax.tripcolor(gtsm.variables['longitude'][:], gtsm.variables['latitude'][:], gtsm.waterlevel, cmap=cmap, vmin=vmin, vmax=vmax, zorder=1)
    coast = cfeature.COASTLINE
    ax.add_feature(coast, facecolor='lightgray', linewidth=0)
    gl = ax.gridlines(alpha=0, draw_labels=True)
    ax.grid(False)
    gl.top_labels = False
    gl.right_labels = False
    gl.left_labels = False
    gl.xlabel_style = {'color': 'k', 'size': 12, 'fontname': 'Arial'}
    gl.ylabel_style = {'color': 'k', 'size': 12, 'fontname': 'Arial'}

    #ax.scatter(tacloban[0], tacloban[1], marker='*', s=50, c='k', edgecolor='k', linewidth=1.5, zorder=99)  # this should also be dependent on the event
    #ax.text(124.85, 11.22, 'Tackloban', fontsize=14, family='Arial', zorder=99) # this shoudl also be dependent on the event

    return #tripcolor

# Create a 3x3 grid of subplots
fig, axs = plt.subplots(3, 3, figsize=(12, 12), subplot_kw={"projection": ccrs.PlateCarree()}, constrained_layout=True)
coast = cfeature.GSHHSFeature(scale='full')

for i in range(3):
    for j in range(3):
        ax = axs[i, j]

        im = ax.tripcolor(map_irma_BS_msk.variables['longitude'][:], map_irma_BS_msk.variables['latitude'][:], map_irma_BS_msk.waterlevel, cmap=cmaps[j], vmin=0, vmax=1, zorder=1)
        ax.add_feature(coast, facecolor='lightgray', linewidth=1)
        ax.set_extent(extents[i], crs=ccrs.PlateCarree())
        gl = ax.gridlines(alpha=0, draw_labels=True)
        ax.grid(False)
        gl.top_labels = False
        gl.right_labels = False
        
        gl.xlabel_style = {'color': 'k', 'size': 12, 'fontname': 'Arial'}
        gl.ylabel_style = {'color': 'k', 'size': 12, 'fontname': 'Arial'}

        # Add colorbar for each column
        if j == 0:
            gl.left_labels = True
            cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02, fraction=0.05)
            cbar.set_label(f'Max. water level [m]', rotation=90, labelpad=15)
            
        else:
            gl.left_labels = False
            
        if j == 2:
            cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02, fraction=0.05)
            cbar.set_label(f'Max. water level difference [m]', rotation=90, labelpad=15)
            
           
def load_case(case, model_config):
    if case == 'xynthia':
        his = xr.open_dataset(f'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/{case}/era5_{model_config}/output/omuse_0000_his.nc').load()
        gtsm_map = xr.open_dataset(f'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/{case}/era5_{model_config}/output/gtsm_map.nc').load()
    else:
        his = xr.open_dataset(f'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/{case}/holland_{model_config}/output/omuse_0000_his.nc').load()
        gtsm_map = xr.open_dataset(f'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/{case}/holland_{model_config}/output/gtsm_map.nc').load()
    
    return his, gtsm_map
'''     
def load_case(case, model_config):
    if case == 'xynthia':
        his = xr.open_dataset(f'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/{case}/era5_{model_config}/output/omuse_0000_his.nc').load()
        gtsm_map_wl = xr.open_dataset(f'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/{case}/era5_{model_config}/output/gtsm_map.nc').load()
    else:
        his = xr.open_dataset(f'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/{case}/holland_{model_config}/output/omuse_0000_his.nc').load()
        gtsm_map_wl = xr.open_dataset(f'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/{case}/holland_{model_config}/output/gtsm_map.nc').load()
        
    gtsm_map_wl_msk = gtsm_map_wl.where(gtsm_map_wl.waterlevel != gtsm_map_wl.waterlevel_ini, drop=True)
    
    # If model_config is not 'BS', subtract gtsm_map for 'BS'
    if model_config != 'BS':
        if case == 'xynthia':
            gtsm_map_bs = xr.open_dataset(f'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/{case}/era5_BS/output/gtsm_map.nc').load()
            gtsm_map_bs_msk = gtsm_map_bs.where(gtsm_map_wl.waterlevel != gtsm_map_wl.waterlevel_ini, drop=True)
        else:
            gtsm_map_bs = xr.open_dataset(f'/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/{case}/holland_BS/output/gtsm_map.nc').load()
            gtsm_map_bs_msk = gtsm_map_bs.where(gtsm_map_wl.waterlevel != gtsm_map_wl.waterlevel_ini, drop=True)
        gtsm_map = gtsm_map_wl_msk.waterlevel - gtsm_map_bs_msk.waterlevel      
    else:
        gtsm_map = gtsm_map_wl_msk.waterlevel
        
    longitude = gtsm_map_wl_msk.longitude
    latitude = gtsm_map_wl_msk.latitude
    
    return his, gtsm_map, longitude, latitude
    
#def remove_dry_cells(map_data):
#    return map_data.where(map_data.waterlevel != map_data.waterlevel_ini, drop=True)

def calculate_max_water_level(his_data):
    return his_data.max('time', skipna=True)

# Set up subplots
fig, axs = plt.subplots(3, 3, figsize=(14, 12.5), subplot_kw={"projection": ccrs.PlateCarree()}, constrained_layout=True)

# Define coast feature
coast = cfeature.GSHHSFeature(scale='full')

# Different extents for each row
extents = [
    [-82.2, -80.2, 30., 32.0],
    [124.8, 125.3, 11.0, 11.5],
    [-1.8, -0.8, 45.5, 46.5]
]

# Define colormap
cmaps = ['jet', 'seismic', 'seismic']
cases = ['irma', 'haiyan', 'xynthia']
model_configs = ['BS', 'TR', 'OR']
model_names = ['Baseline configuration', 'Temporal Resolution', 'Improved bathymetry']
vmin_max_values = [
    [(0, 1.5), (-0.5, 0.5), (-0.5, 0.5)],  # vmin_max values for 'irma'
    [(0, 6), (-2, 2), (-2, 2)],  # vmin_max values for 'haiyan'
    [(0, 5), (-0.5, 0.5), (-0.5, 0.5)]  # vmin_max values for 'xynthia'
]

jacksonville = np.asarray([-81.6,30.3])
tacloban = np.asarray([125.0,11.24])
rochelle = np.asarray([-1.2,46.19])
city_locations = [jacksonville, tacloban, rochelle]
city_names = ['Jacksonville', 'Tacloban', 'La Rochelle']
city_name_locations = [np.asarray([-82.15,30.2]), np.asarray([124.89,11.22]), np.asarray([-1.16,46.17])]
case_names = ['Irma','Haiyan','Xynthia']
subplot_numbers = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)']

for i, case in enumerate(cases):
    for j, model_config in enumerate(model_configs):
        ax = axs[i, j]

        # Load datasets
        his, gtsm_map, longitude, latitude = load_case(case, model_config)
        
        # Remove dry cells
        #gtsm_map_msk = remove_dry_cells(gtsm_map)

        # Calculate max water level
        his_max = calculate_max_water_level(his)

        # Plotting
        vmin, vmax = vmin_max_values[i][j]
        
        ax.add_feature(coast, facecolor='lightgray', linewidth=1)
        ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='darkgray', zorder=99)
        ax.text(city_name_locations[i][0],city_name_locations[i][1],city_names[i],fontsize=12, family='Arial', zorder=99)

        ax.set_extent(extents[i], crs=ccrs.PlateCarree())
        gl = ax.gridlines(alpha=0, draw_labels=True)
        ax.grid(False)
        gl.top_labels = False
        gl.right_labels = False

        gl.xlabel_style = {'color': 'k', 'size': 12, 'fontname': 'Arial'}
        gl.ylabel_style = {'color': 'k', 'size': 12, 'fontname': 'Arial'}
        
        ax.scatter(city_locations[i][0],city_locations[i][1],marker='*',s=50,c='k',edgecolor='k',linewidth=1.5,zorder=99)

        square = patches.Rectangle((0.0, 0.9), 0.1, 0.1, transform=ax.transAxes, facecolor='white', edgecolor='black', linewidth=1, alpha=1, zorder=98)
        ax.add_patch(square)
        #square2 = patches.Rectangle((0.0, 0.9), 0.1, 0.1, transform=ax.transAxes, color='none', edgecolor='red', linewidth=10, linestyle='-', alpha=1, zorder=98)
        #ax.add_patch(square2)
        ax.text(0.03, 0.93, subplot_numbers[i * len(model_configs) + j],
                transform=ax.transAxes, fontsize=12, fontweight='bold', color='black', zorder=99)


        # Add colorbar for each column
        if j == 0:
            im = ax.tripcolor(longitude, latitude, gtsm_map, cmap=cmaps[j], vmin=vmin, vmax=vmax, zorder=1)
            gl.left_labels = True
            cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02, fraction=0.05)
            cbar.set_label(f'Max. water level [m]', rotation=90, fontsize=12, labelpad=8)
            cbar.ax.tick_params(labelsize=12)
            ax.text(-0.26, 0.55, case_names[i], fontsize=14, va='center', ha='center',
                rotation='vertical', rotation_mode='anchor',
                transform=ax.transAxes)

        else:
            im = ax.tripcolor(longitude, latitude, gtsm_map, cmap=cmaps[j], vmin=vmin, vmax=vmax, zorder=1)
            gl.left_labels = False

        if j == 2:
            cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.02, fraction=0.05)
            cbar.set_label(f'Max. water level difference [m]', rotation=90, fontsize=12, labelpad=8)
            cbar.ax.tick_params(labelsize=12)
            
        # Add titles to the first row
        if i == 0:
            ax.set_title(f'{model_names[j]}', fontsize=14)

fig.savefig(f"figures/test.png")