import xarray as xr
import numpy as np
import math
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

mapfile=xr.open_dataset(r'/gpfs/home4/benitoli/gtsm_d3dfm/era5_BS_irma/output/gtsm_fine_0010_map.nc').load()

# set up a map
fig = plt.figure(figsize=(20, 8), tight_layout=False)
columns=4
rows=2
projex= ccrs.PlateCarree()

# Wind calculation
wind=np.sqrt((mapfile.windx)*(mapfile.windx)+(mapfile.windy)*(mapfile.windy))

def plotmymap(axs):
    # your plot specs of each map should replace this
    wind_t=wind[7 + i*6]
    #print(wind_t.time)
    #fig, (ax1) = plt.subplots(1, 1, subplot_kw={"projection": ccrs.PlateCarree()}, gridspec_kw = {'wspace':0.4, 'hspace':0.3})
    coast = cfeature.GSHHSFeature(scale='full')
    axs.add_feature(coast, facecolor='lightgray', linewidth=0, zorder=99)
    #axs.add_feature(coast, facecolor='none', linewidth=1.0)#, zorder 99)
    axs.set_extent([-90, -71, 18, 35], ccrs.PlateCarree())
    a=axs.tripcolor(wind_t.FlowElem_xcc,wind_t.FlowElem_ycc,wind_t, cmap='jet',vmin=0,vmax=60)
    axs.set_title(str(wind_t.time.values), fontsize=10, fontname='Arial')
    cbar=fig.colorbar(a,ax=axs, fraction=0.04, pad=0.05)
    cbar.ax.set_ylabel('Wind speed [m/s]', rotation=90, fontsize=10)
    return a  # for use by colorbar

for i in range(1, columns*rows +1):
    # add a subplot into the array of plots
    ax = fig.add_subplot(rows, columns, i, projection=projex)
    plims = plotmymap(ax)  # a simple maps is created on subplot


bottom, top = 0.1, 0.9
left, right = 0.1, 0.8
fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=0.15, wspace=0.25)
#cbar_ax = fig.add_axes([0.85, bottom, 0.05, top-bottom])
#fig.colorbar(plims, cax=cbar_ax)  # plot colorbar
plt.savefig(f"figures/gtsm_map_output_wind_era5.png")
plt.clf()

