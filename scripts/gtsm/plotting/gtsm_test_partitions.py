import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

base_file_path = r'/projects/0/einf2224/paper1/scripts/case_studies/haiyan/no_omuse/holland_GR/local_haiyan/output_twoweeks_charnock/gtsm_fine_00*_map.nc'
matching_files = glob.glob(base_file_path)

# Create a plot for each file
for file in matching_files:
    # Open the NetCDF file
    dataset = xr.open_dataset(file)

    # Calculate the maximum along the 'time' dimension
    test_max = dataset.max(dim='time')

    # Access the 's1' variable from the result
    s1_variable = test_max['s1']

    # Create a scatter plot for the current file
    plt.scatter(test_max['FlowElem_xcc'], test_max['FlowElem_ycc'], c=s1_variable, vmin=0,  cmap='jet', s=2, label=file)

# Add labels and legend
plt.xlabel('FlowElem_xcc')
plt.ylabel('FlowElem_ycc')
#plt.legend()
plt.title('Scatter Plots for Multiple Files')
plt.colorbar()
# Show the plot
plt.savefig('figures/test_charnock_v2.png')