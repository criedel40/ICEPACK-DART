import xarray as xr
import numpy as np

data = xr.open_dataset('gridpt_cice_data.nc')
aice_check = data.si_data.values[:,1,:]
data.close()

data = xr.open_dataset('ICs_Forcings.nc')
ics = data.ICs.values
forcings = data.forcings.values

test1 = np.zeros((100,3650))
test2 = np.zeros((100,3650))
test3 = np.zeros((100,3650))


prior = ics
for t in range(0,test1.shape[1]):
  if t == 0:
    prior = ics
  else:
    prior = post + forcings[:,t-1]*1.0
  ### DO ASSIM
  post = prior
  
  test1[:,t] = prior
  test2[:,t] = prior
  test3[:,t] = post
  #prior = post + forcings[:,t]*1.0
