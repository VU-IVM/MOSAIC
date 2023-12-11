from hydromt_sfincs import SfincsModel
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import xarray as xr



def hmax(case, model_config):          # THIS NEEDS TO BE UPDATED!!! TO DO!!!
    #if case == 'xynthia':
    #    sfincs_root=f'/projects/0/einf2224/paper1/scripts/case_studies/hydromt-sfincs/{case}/era5_{model_config}'
    #else:
    #    sfincs_root=f'/projects/0/einf2224/paper1/scripts/case_studies/hydromt-sfincs/{case}/holland_{model_config}'
        
    sfincs_root=f'/projects/0/einf2224/paper1/scripts/case_studies/hydromt-sfincs/haiyan/test/haiyan_holland_BS'
    mod=SfincsModel(sfincs_root, mode="r")
    mod.write_raster("results.hmax", compress="LZW")
    mod.data_catalog.from_yml(r'/projects/0/einf2224/paper1/scripts/hydromt-sfincs/general/data_catalog.yml')
    gswo = mod.data_catalog.get_rasterdataset("gswo", geom=mod.region, buffer=10)
    gswo_mask = gswo.raster.reproject_like(mod.grid, method="max") <= 5                                        # permanent water where water occurence > 5%
    hmin = 0.1                                                                                                       # minimum flood depth [m] to plot
    da = mod.results["hmax"]    
    da_hmax = da.max(dim='timemax')
    da_hmax_fld = da_hmax.where(gswo_mask).where(da_hmax > hmin) 
    da_hmax_fld.rio.set_crs(mod.crs)
    flood_max = da_hmax_fld.raster.reproject_like(gswo)
    #flood_max_values = flood_max.compute()
    return flood_max, gswo_mask



# Set up subplots
fig, axs = plt.subplots(3, 4, figsize=(16, 11.5), subplot_kw={"projection": ccrs.PlateCarree()}, constrained_layout=True)

# Define coast feature
coast = cfeature.GSHHSFeature(scale='full')

# Different extents for each row
extents = [
    [-82.2, -80.2, 30., 32.0],
    [124.8, 125.3, 11.0, 11.5],
    [-1.8, -0.8, 45.5, 46.5]
]

# Define colormap & other features
cmaps = ['cool', 'seismic', 'seismic', 'seismic']
cases = ['irma', 'haiyan', 'xynthia']
model_configs = ['BS', 'TR', 'OR', 'IB']
model_names = ['Baseline configuration', 'Temporal resolution refined', 'Output resolution refined', 'Improved bathymetry']
vmin_max_values = [
    [(0, 3), (-0.5, 0.5), (-0.5, 0.5), (-0.5, 0.5)],  # vmin_max values for 'irma'
    [(0, 3), (-2, 2), (-2, 2), (-2, 2)],  # vmin_max values for 'haiyan'
    [(0, 5), (-0.5, 0.5), (-0.5, 0.5), (-0.5, 0.5)]  # vmin_max values for 'xynthia'
]

jacksonville = np.asarray([-81.6,30.3])
tacloban = np.asarray([125.0,11.24])
rochelle = np.asarray([-1.2,46.19])
city_locations = [jacksonville, tacloban, rochelle]
city_names = ['Jacksonville', 'Tacloban', 'La Rochelle']
city_name_locations = [np.asarray([-82.15,30.2]), np.asarray([124.89,11.22]), np.asarray([-1.16,46.17])]
case_names = ['Irma','Haiyan','Xynthia']
subplot_numbers = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)', 'i)', 'j)', 'k)', 'l)']

for i, case in enumerate(cases):
    print('CASE:', case)
    for j, model_config in enumerate(model_configs):
        print('MODEL CONFIG', model_config)
        ax = axs[i, j]

        # Flood max calculation
        flood_max, gswo_mask=hmax(case, model_config)
        print('flood_max:', flood_max)
        print('gswo_mask:', gswo_mask.compute())
        gswo_mask_true = gswo_mask.where(gswo_mask, drop=False)
        ax.imshow(gswo_mask_true, cmap='gray', vmin=0, vmax=1, zorder=97, alpha=0.5)


        # Plotting
        vmin, vmax = vmin_max_values[i][j]
        
        ax.set_facecolor('lightgray')
        ax.add_feature(coast, facecolor='white', linewidth=1)
        ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='darkgray', zorder=100)
        #ax.add_feature(cfeature.LAND, facecolor='lightgray')
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
        ax.text(0.02, 0.935, subplot_numbers[i * len(model_configs) + j],
                transform=ax.transAxes, fontsize=12, fontweight='bold', color='black', zorder=99)

        # Add colorbar for each column
        if j == 0:
            flood_max.plot(ax=ax, cmap=cmaps[j], vmin=vmin, vmax=vmax, zorder=97, cbar_kwargs={"label": "Max. Flood depth [m]", "shrink": 0.95, "pad": 0.02})
            gl.left_labels = True
            ax.text(-0.26, 0.55, case_names[i], fontsize=14, va='center', ha='center',
                rotation='vertical', rotation_mode='anchor',
                transform=ax.transAxes)
                
        elif j == 3: #maybe here i need to put elif:
            flood_max.plot(ax=ax, cmap=cmaps[j], vmin=vmin, vmax=vmax, zorder=97, cbar_kwargs={"label": "Max. Flood depth difference [m]", "shrink": 1, "pad": 0.02})
            gl.left_labels = False

        else:
            flood_max_plot=flood_max.plot(ax=ax, cmap=cmaps[j], vmin=vmin, vmax=vmax, zorder=97, add_colorbar=False)
            #colorbar = plt.gcf().colorbar(flood_max_plot, ax=ax, shrink=0.8, pad=0.02)
            #colorbar.remove()
            gl.left_labels = False

            
        # Add titles to the first row
        if i == 0:
            ax.set_title(f'{model_names[j]}', fontsize=14)   
        else:
            ax.set_title(' ')

fig.savefig(f"figures/test.png")
