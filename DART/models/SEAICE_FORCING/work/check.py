import xarray as xr
import numpy as np

data = xr.open_dataset('gridpt_cice_data.nc')
aice_org = data.si_data.values[:,1,:]
data.close()

data = xr.open_dataset('output_control.nc')
aice_compute = data.aice.values.T
data.close()

print(np.nanmean(aice_compute - aice_org))

