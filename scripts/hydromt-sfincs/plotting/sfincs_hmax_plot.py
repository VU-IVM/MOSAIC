import xarray as xr
import numpy as np
from os.path import join
import matplotlib.pyplot as plt
import hydromt
from hydromt_sfincs import SfincsModel
from osgeo import gdal
import rioxarray
import cartopy.crs as ccrs
import cartopy.feature as cfeature
'''
#sfincs_root = "/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/hydromt-sfincs/haiyan/haiyan_track_data"
sfincs_root = "/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/hydromt-sfincs/xynthia/xynthia_ERA5_v2"  
#sfincs_root = "../../case_studies/hydromt-sfincs/xynthia/xynthia_ERA5"                                                        # (relative) path to sfincs root
mod = SfincsModel(sfincs_root, mode="r")
mod.read_results()
list(mod.results.keys())
mod.write_raster("results.hmax", compress="LZW")
# read global surface water occurance (GSWO) data to mask permanent water
mod.data_catalog.from_yml(r'/gpfs/home4/benitoli/papers/paper1/scripts/hydromt-sfincs/general/data_catalog.yml')
gswo = mod.data_catalog.get_rasterdataset("gswo", buffer=10)
gswo_mask = gswo.raster.reproject_like(mod.staticmaps, method="max") <= 5                                        # permanent water where water occurence > 5%
hmin = 0.1                                                                                                       # minimum flood depth [m] to plot
da_hmax = mod.results["hmax"]                                                                                    # hmax is computed from zsmax - zb
#da_hmax_fld=da_hmax
da_hmax_fld = da_hmax.where(gswo_mask).where(da_hmax > hmin)                                                     # get overland flood depth with GSWO and set minimum flood depth
# update attributes for colorbar label later
da_hmax.attrs.update(long_name="flood depth", unit="m")
#da_hmax_fld.rio.to_raster("irma_omuse_moreStations_larger.tif")                                          # save flood file
air = da_hmax_fld.rio.reproject("EPSG:4326")
'''
# calculate max water level & mask permanent water
def hmax(sfincs_root):
  mod=SfincsModel(sfincs_root, mode="r")
  mod.read_results()
  mod.write_raster("results.hmax", compress="LZW")
  mod.data_catalog.from_yml(r'/gpfs/home4/benitoli/papers/paper1/scripts/hydromt-sfincs/general/data_catalog.yml')
  gswo = mod.data_catalog.get_rasterdataset("gswo", buffer=10)
  gswo_mask = gswo.raster.reproject_like(mod.staticmaps, method="max") <= 5                                        # permanent water where water occurence > 5%
  hmin = 0.1                                                                                                       # minimum flood depth [m] to plot
  da_hmax = mod.results["hmax"]    
  da_hmax_fld = da_hmax.where(gswo_mask).where(da_hmax > hmin) 
  da_hmax.attrs.update(long_name="flood depth", unit="m") 
  flooding = da_hmax_fld.rio.reproject("EPSG:4326")
  return flooding

air = hmax("/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/hydromt-sfincs/xynthia/xynthia_ERA5")

# plotting
fig = plt.figure()
f, ((ax1, ax6, ax11),(ax2,ax7, ax12),(ax3,ax8,ax13),(ax4,ax9,ax14), (ax5,ax10,ax15)) = plt.subplots(5, 3, figsize=(12, 16), subplot_kw={"projection": ccrs.PlateCarree()}, gridspec_kw = {'wspace':0.5, 'hspace':0.3})
# Irma baseline
#a=hmax("/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/hydromt-sfincs/xynthia/xynthia_ERA5").plot(ax=ax1,  cmap=plt.cm.cool, vmin=0.1, vmax=2.0, cbar_kwargs={"label": "Maximum flood depth [m]"}, zorder=90)
ax1.set_title("Irma")
ax1.set_extent([-82.5, -79.5, 29.5, 32.5], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.LAND, facecolor='lightgray')
ax1.add_feature(cfeature.COASTLINE, linewidth=0.0)
ax1.text(-0.4, 0.55, 'Baseline configuration', fontsize=12, va='center', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax1.transAxes)
gl1 = ax1.gridlines(alpha=0, draw_labels=True)
ax1.grid(False)
gl1.top_labels=False   # suppress top labels
gl1.right_labels=False # suppress right labels
gl1.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl1.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Irma grid refinement
air.plot(ax=ax2, cmap=plt.cm.bwr, vmin=-0.5, vmax=0.5, cbar_kwargs={"label": "Flood depth difference [m]"})
ax2.set_title("")
ax2.set_extent([-82.5, -79.5, 29.5, 32.5], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.LAND, facecolor='lightgray')
ax2.add_feature(cfeature.COASTLINE, linewidth=0.0)
ax2.text(-0.4, -0.85, 'Grid refinement', fontsize=12, va='center', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax1.transAxes)
gl2 = ax2.gridlines(alpha=0, draw_labels=True)
ax2.grid(False)
gl2.top_labels=False   # suppress top labels
gl2.right_labels=False # suppress right labels
gl2.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl2.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Irma bed refinement
air.plot(ax=ax3, cmap=plt.cm.bwr, vmin=-0.5, vmax=0.5, cbar_kwargs={"label": "Flood depth difference [m]"})
ax3.set_title("")
ax3.set_extent([-82.5, -79.5, 29.5, 32.5], crs=ccrs.PlateCarree())
ax3.add_feature(cfeature.LAND, facecolor='lightgray')
ax3.add_feature(cfeature.COASTLINE, linewidth=0.0)
ax3.text(-0.4, -2.25, 'Bed refinement', fontsize=12, va='center', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax1.transAxes)
gl3 = ax3.gridlines(alpha=0, draw_labels=True)
ax3.grid(False)
gl3.top_labels=False   # suppress top labels
gl3.right_labels=False # suppress right labels
gl3.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl3.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Irma output refinement
air.plot(ax=ax4, cmap=plt.cm.bwr, vmin=-0.5, vmax=0.5, cbar_kwargs={"label": "Flood depth difference [m]"})
ax4.set_title("")
ax4.set_extent([-82.5, -79.5, 29.5, 32.5], crs=ccrs.PlateCarree())
ax4.add_feature(cfeature.LAND, facecolor='lightgray')
ax4.add_feature(cfeature.COASTLINE, linewidth=0.0)
ax4.set_xlabel("")
ax4.text(-0.4, -3.65, 'Output refinement', fontsize=12, va='center', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax1.transAxes)
gl4 = ax4.gridlines(alpha=0, draw_labels=True)
ax4.grid(False)
gl4.top_labels=False   # suppress top labels
gl4.right_labels=False # suppress right labels
gl4.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl4.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Irma time refinement
air.plot(ax=ax5, cmap=plt.cm.bwr, vmin=-0.5, vmax=0.5, cbar_kwargs={"label": "Flood depth difference [m]"})
ax5.set_title("")
ax5.set_extent([-82.5, -79.5, 29.5, 32.5], crs=ccrs.PlateCarree())
ax5.add_feature(cfeature.LAND, facecolor='lightgray')
ax5.add_feature(cfeature.COASTLINE, linewidth=0.0)
ax5.set_xlabel("")
ax5.text(-0.4, -5.05, 'Time refinement', fontsize=12, va='center', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax1.transAxes)
gl5 = ax5.gridlines(alpha=0, draw_labels=True)
ax5.grid(False)
gl5.top_labels=False   # suppress top labels
gl5.right_labels=False # suppress right labels
gl5.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl5.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Haiyan baseline
#b=hmax("/gpfs/home4/benitoli/papers/paper1/scripts/case_studies/hydromt-sfincs/haiyan/haiyan_track_data_smaller").plot(ax=ax6,  cmap=plt.cm.cool, vmin=0.1, vmax=2.0, cbar_kwargs={"label": "Maximum flood depth [m]"}, zorder=90)
ax6.set_title("Haiyan")
ax6.set_extent([123, 126, 9.75, 12.75], crs=ccrs.PlateCarree())
ax6.add_feature(cfeature.LAND, facecolor='lightgray')
ax6.add_feature(cfeature.COASTLINE, linewidth=0.0)
ax6.set_xlabel("")
#ax1.set_ylabel("test!!")
gl6 = ax6.gridlines(alpha=0, draw_labels=True)
ax6.grid(False)
gl6.top_labels=False   # suppress top labels
gl6.right_labels=False # suppress right labels
gl6.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl6.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Haiyan grid refinement
air.plot(ax=ax7, cmap=plt.cm.bwr, vmin=-0.5, vmax=0.5, cbar_kwargs={"label": "Flood depth difference [m]"})
ax7.set_title("")
ax7.set_extent([123, 126, 9.75, 12.75], crs=ccrs.PlateCarree())
ax7.add_feature(cfeature.LAND, facecolor='lightgray')
ax7.add_feature(cfeature.COASTLINE, linewidth=0.0)
ax7.set_xlabel("")
gl7 = ax7.gridlines(alpha=0, draw_labels=True)
ax7.grid(False)
gl7.top_labels=False   # suppress top labels
gl7.right_labels=False # suppress right labels
gl7.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl7.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Haiyan bed refinement
air.plot(ax=ax8, cmap=plt.cm.bwr, vmin=-0.5, vmax=0.5, cbar_kwargs={"label": "Flood depth difference [m]"})
ax8.set_title("")
ax8.set_extent([123, 126, 9.75, 12.75], crs=ccrs.PlateCarree())
ax8.add_feature(cfeature.LAND, facecolor='lightgray')
ax8.add_feature(cfeature.COASTLINE, linewidth=0.0)
ax8.set_xlabel("")
gl8 = ax8.gridlines(alpha=0, draw_labels=True)
ax8.grid(False)
gl8.top_labels=False   # suppress top labels
gl8.right_labels=False # suppress right labels
gl8.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl8.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Haiyan output refinement
air.plot(ax=ax9, cmap=plt.cm.bwr, vmin=-0.5, vmax=0.5, cbar_kwargs={"label": "Flood depth difference [m]"})
ax9.set_title("")
ax9.set_extent([123, 126, 9.75, 12.75], crs=ccrs.PlateCarree())
ax9.add_feature(cfeature.LAND, facecolor='lightgray')
ax9.add_feature(cfeature.COASTLINE, linewidth=0.0)
ax9.set_xlabel("")
gl9 = ax9.gridlines(alpha=0, draw_labels=True)
ax9.grid(False)
gl9.top_labels=False   # suppress top labels
gl9.right_labels=False # suppress right labels
gl9.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl9.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Haiyan time refinement
air.plot(ax=ax10, cmap=plt.cm.bwr, vmin=-0.5, vmax=0.5, cbar_kwargs={"label": "Flood depth difference [m]"})
ax10.set_title("")
ax10.set_extent([123, 126, 9.75, 12.75], crs=ccrs.PlateCarree())
ax10.add_feature(cfeature.LAND, facecolor='lightgray')
ax10.add_feature(cfeature.COASTLINE, linewidth=0.0)
ax10.set_xlabel("")
gl10 = ax10.gridlines(alpha=0, draw_labels=True)
ax10.grid(False)
gl10.top_labels=False   # suppress top labels
gl10.right_labels=False # suppress right labels
gl10.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl10.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Xynthia baseline
air.plot(ax=ax11, cmap=plt.cm.cool, vmin=0.1, vmax=2.0, cbar_kwargs={"label": "Maximum flood depth [m]"}, zorder=90)
ax11.set_title("Xynthia")
ax11.set_extent([-1.8, -0.8, 45.5, 46.5], crs=ccrs.PlateCarree())
coast = cfeature.GSHHSFeature(scale='full')
ax11.add_feature(coast, facecolor='lightgray', linewidth=0)
ax11.set_xlabel("")
#ax1.set_ylabel("test!!")
gl11 = ax11.gridlines(alpha=0, draw_labels=True)
ax11.grid(False)
gl11.top_labels=False   # suppress top labels
gl11.right_labels=False # suppress right labels
gl11.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl11.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Xynthia grid refinement
air.plot(ax=ax12, cmap=plt.cm.bwr, vmin=-0.5, vmax=0.5, cbar_kwargs={"label": "Flood depth difference [m]"}, zorder=99)
ax12.set_title("")
ax12.set_extent([-1.8, -0.8, 45.5, 46.5], crs=ccrs.PlateCarree())
#coast = cfeature.GSHHSFeature(scale='full')
ax12.add_feature(coast, facecolor='lightgray', linewidth=0.5)
ax12.set_xlabel("")
gl12 = ax12.gridlines(alpha=0, draw_labels=True)
ax12.grid(False)
gl12.top_labels=False   # suppress top labels
gl12.right_labels=False # suppress right labels
gl12.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl12.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Xynthia bed refinement
air.plot(ax=ax13, cmap=plt.cm.bwr, vmin=-0.5, vmax=0.5, cbar_kwargs={"label": "Flood depth difference [m]"}, zorder=99)
ax13.set_title("")
ax13.set_extent([-1.8, -0.8, 45.5, 46.5], crs=ccrs.PlateCarree())
ax13.add_feature(coast, facecolor='lightgray', linewidth=0.5)
ax13.set_xlabel("")
gl13 = ax13.gridlines(alpha=0, draw_labels=True)
ax13.grid(False)
gl13.top_labels=False   # suppress top labels
gl13.right_labels=False # suppress right labels
gl13.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl13.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Xynthia output refinement
air.plot(ax=ax14, cmap=plt.cm.bwr, vmin=-0.5, vmax=0.5, cbar_kwargs={"label": "Flood depth difference [m]"}, zorder=99)
ax14.set_title("")
ax14.set_extent([-1.8, -0.8, 45.5, 46.5], crs=ccrs.PlateCarree())
ax14.add_feature(coast, facecolor='lightgray', linewidth=0.5)
ax14.set_xlabel("")
gl14 = ax14.gridlines(alpha=0, draw_labels=True)
ax14.grid(False)
gl14.top_labels=False   # suppress top labels
gl14.right_labels=False # suppress right labels
gl14.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl14.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Xynthia time refinement
air.plot(ax=ax15, cmap=plt.cm.bwr, vmin=-0.5, vmax=0.5, cbar_kwargs={"label": "Flood depth difference [m]"}, zorder=99)
ax15.set_title("")
ax15.set_extent([-1.8, -0.8, 45.5, 46.5], crs=ccrs.PlateCarree())
ax15.add_feature(coast, facecolor='lightgray', linewidth=0.5)
ax15.set_xlabel("")
gl15 = ax15.gridlines(alpha=0, draw_labels=True)
ax15.grid(False)
gl15.top_labels=False   # suppress top labels
gl15.right_labels=False # suppress right labels
gl15.xlabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}
gl15.ylabel_style = {'color': 'k', 'size':10, 'fontname':'Arial'}

# Make it nice
plt.tight_layout()
plt.savefig(f"figures/sfincs_output_xynthia.png", bbox_inches='tight')
