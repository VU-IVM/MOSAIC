import xarray as xr
import numpy as np
from os.path import join
import matplotlib.pyplot as plt
import hydromt
from hydromt_sfincs import SfincsModel

'''
sfincs_root = "/projects/0/einf2224/paper1/scripts/case_studies/hydromt-sfincs/haiyan/test/haiyan_holland_BS"  # (relative) path to sfincs root
mod = SfincsModel(sfincs_root, mode="r")
mod.read_results()
h = mod.results["h"].where(mod.results["h"] > 0)
area_grid=h.raster.area_grid()
volume_grid=h*area_grid
print('volume_grid:', volume_grid)
total_volume = volume_grid.sum(dim=['x', 'y'])


print('total_volume:', total_volume)

plt.plot(total_volume.time, total_volume)



sfincs_root = "/projects/0/einf2224/paper1/scripts/case_studies/hydromt-sfincs/haiyan/test/haiyan_holland_BS"
mod = SfincsModel(sfincs_root, mode="r")
mod.read_results()

#------------------------------------------------------------------------------------------------------------------------------------
### Read the observation points to make the flood depth timeseries
#------------------------------------------------------------------------------------------------------------------------------------

h_point = mod.results["point_h"]#.rename({"stations": "station_id"})

#------------------------------------------------------------------------------------------------------------------------------------
### Read the map file to make the flood volume timeseries
#------------------------------------------------------------------------------------------------------------------------------------

h = mod.results["h"].where(mod.results["h"] > 0)
area_grid = h.raster.area_grid()
volume_grid = h * area_grid
total_volume = volume_grid.sum(dim=['x', 'y'])

#------------------------------------------------------------------------------------------------------------------------------------
### Plot
#------------------------------------------------------------------------------------------------------------------------------------

num_stations = len(h_point.stations)
num_columns_his = num_stations
num_columns_sfincs = 1  # Since there's only one subplot for sfincs volume timeseries

# Create a new figure with subplots for hydromt-sfincs and sfincs data
fig, axes = plt.subplots(1, num_columns_his + num_columns_sfincs, figsize=(15, 5 * (num_columns_his + num_columns_sfincs)))

# Loop through each station and plot in the corresponding subplot for hydromt-sfincs data
for station_idx in range(num_stations):
    # Extract data for the current station
    current_station_data = h_point.sel(stations=station_idx)

    # Plot your data for the current station
    # Replace the following lines with your specific plot commands
    axes[station_idx].plot(current_station_data.time, current_station_data, label='Flood depth 2 [m]')
    #axes[station_idx].set_title(f"Station: {current_station_data['station_name'].values.decode('utf-8')}")

    # Add labels, legends, etc. as needed
    axes[station_idx].set_xlabel('Time')
    axes[station_idx].set_ylabel('Flood depth [m]')
    axes[station_idx].legend()


# Plot sfincs data
axes[-1].plot(total_volume.time, total_volume)
axes[-1].set_xlabel('Time')
axes[-1].set_ylabel('Flood Volume [m3]')

# Adjust layout for better spacing
plt.tight_layout()
'''
cases = ['irma', 'haiyan', 'xynthia']
model_configs = ['BS', 'TR', 'OR', 'IB']

num_stations = 2
num_columns_his = num_stations
num_columns_sfincs = 1  # Since there's only one subplot for sfincs volume timeseries
num_rows=len(cases)

# Create a new figure with subplots for hydromt-sfincs and sfincs data
fig, axes = plt.subplots(num_rows, num_columns_his + num_columns_sfincs, figsize=(15, 3 * (num_columns_his + num_columns_sfincs)), constrained_layout=True)

for row_idx, case in enumerate(cases):
    for model_config in model_configs:
        # Construct sfincs_root for the current model_config
        #sfincs_root = f"/projects/0/einf2224/paper1/scripts/case_studies/hydromt-sfincs/{case}/test/{case}_holland_{model_config}"  # need to fix the holland, for the case of xynthia
        sfincs_root = f"/projects/0/einf2224/paper1/scripts/case_studies/hydromt-sfincs/haiyan/test/haiyan_holland_BS"
        
        # Read results for the current model_config & case
        mod = SfincsModel(sfincs_root, mode="r")
        mod.read_results()
    
        # Read the observation points to make the flood depth timeseries
        h_point = mod.results["point_h"]
    
        # Read the map file to make the flood volume timeseries
        h = mod.results["h"].where(mod.results["h"] > 0)
        area_grid = h.raster.area_grid()
        volume_grid = h * area_grid
        total_volume = volume_grid.sum(dim=['x', 'y'])
    
        # Loop through each station and plot in the corresponding subplot for hydromt-sfincs data
        for station_idx in range(len(h_point.stations)):
            # Extract data for a station
            current_station_data = h_point.sel(stations=station_idx) 
    
            # Plot flood depth timeseries for the station
            axes[row_idx, station_idx].plot(current_station_data.time, current_station_data)
            
            if row_idx == num_rows - 1:
                axes[row_idx, station_idx].set_xlabel('Time')
            axes[row_idx, station_idx].set_ylabel('Flood depth [m]')
            axes[row_idx, station_idx].set_title(f'{case} - Location {station_idx}')
            
            plt.setp(axes[row_idx, station_idx].xaxis.get_majorticklabels(), rotation=45, ha="right")
    
        # Plot sfincs data for the current model_config
        axes[row_idx, -1].set_title(f'{case} - Study area')
        axes[row_idx, -1].plot(total_volume.time, total_volume)
        if row_idx == num_rows - 1:
            axes[row_idx, -1].set_xlabel('Time')
        axes[row_idx, -1].set_ylabel('Flood Volume [m3]')
        
        plt.setp(axes[row_idx, -1].xaxis.get_majorticklabels(), rotation=45, ha="right")

fig.legend(model_configs, loc='lower center', bbox_to_anchor=(0.5, 0), ncol=len(model_configs), frameon=False)
plt.savefig('figures/sfincs_timeseries.png')
