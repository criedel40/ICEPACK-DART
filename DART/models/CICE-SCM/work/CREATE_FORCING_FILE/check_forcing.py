import xarray as xr
import numpy as np

data = xr.open_dataset('gridpt_cice_data.nc')
check_hs = data.si_hs.values[87,1,:].T
data.close()

#data = xr.open_dataset('seaice_snow_forcings.nc')
#ics_hs = data.hs_ics.values
#hs_forcings = data.hs_forcings.values.T

hs_ics = check_hs[0]
forcings = check_hs[1:] - check_hs[:-1]
sys.exit()
#control = np.zeros((forcings.shape[0]+1))
#for t in range(0,control.shape[0]):
#  if t == 0:
#    control[t] = hs_ics
#  else:
    


#sys.exit()
#control = np.zeros((hs_forcings.shape[0]+1,hs_forcings.shape[1]))
#for t in range(0,hs_forcings.shape[0]+1):
#  if t == 0:
#    prior_snow = ics_hs
#  else:
#    #prior_snow = update_mem(post_snow,forcings_snow[t-1,:],1)
#    update_mem = post_snow + hs_forcings[t-1,:]*1
#    if np.any(update_mem[:] - check_hs[t,:])!=0:
#      print('STOP')
#      sys.exit()
#    update_mem[update_mem<0] = 0.0
#    prior_snow = update_mem
#  ##
#  post_snow = prior_snow
#  control[t,:] = post_snow
