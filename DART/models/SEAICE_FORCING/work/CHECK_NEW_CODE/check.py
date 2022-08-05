import xarray as xr
import numpy as np
import matplotlib.pyplot as plt


data = xr.open_dataset('assim_output_EaKF_gaus.nc')
#data = xr.open_dataset('../output_control.nc')
output_aice = data.output_aice.values
output_hi = data.output_hi.values
output_hs = data.output_hs.values
data.close()

#data = xr.open_dataset('../CREATE_FORCING_FILE/gridpt_cice_data.nc')
data = xr.open_dataset('../output_control.nc')
check_aice = data.aice.values[:,:-1]
check_hi = data.hi.values[:,:-1]
check_hs = data.hs.values[:,:-1]

diff = (output_hi - check_hi).mean()
print('Diff HI:{0}'.format(diff))
diff = (output_hs - check_hs).mean()
print('Diff HS:{0}'.format(diff))

print('----------------------------')
data = xr.open_dataset('../DIFF_OBS_ERR/BOUNDEDRHF/CONTROL/assim_output_EaKF_gaus.nc')
old_aice = data.output_aice.values

print('AICE DIFF:{0}'.format((output_aice - old_aice).mean()))

