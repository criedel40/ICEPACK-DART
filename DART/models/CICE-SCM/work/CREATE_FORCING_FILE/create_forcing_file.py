import xarray as xr
import netCDF4
import numpy as np
import sys
import os
import matplotlib.pyplot as plt

data = xr.open_dataset('gridpt_cice_data.nc')
aice = data.si_aice.values[:,1,:]
hi = data.si_hi.values[:,1,:]
hs = data.si_hs.values[:,1,:]


#ratio_hi = np.divide(hi,aice,out=np.zeros_like(aice).astype('float64'),where=aice!=0)
#ratio_hs = np.divide(hs,aice,out=np.zeros_like(aice).astype('float64'),where=aice!=0)

#new_hi = ratio_hi*aice
#new_hs = ratio_hs*aice

#tim = np.arange(0,hs[:,:720].shape[1])
#fig,ax = plt.subplots(3,1,sharex=True)
#ax[0].plot(tim,hs[:,:720].mean(axis=0),'-k')
#ax[1].plot(tim,new_hs[:,:720].mean(axis=0),'-k')
#ax[2].plot(tim,(new_hs-hs)[:,:720].mean(axis=0),'-k')
#plt.savefig('test.jpg',dpi=550)
#os.system('open test.jpg')

#sys.exit()
aice_ics = aice[:,0]
hi_ics = hi[:,0]
hs_ics = hs[:,0]

aice_forcings = aice[:,1:] - aice[:,:-1]
hi_forcings = hi[:,1:] - hi[:,:-1]
hs_forcings = hs[:,1:] - hs[:,:-1]


new_file = netCDF4.Dataset('seaice_snow_forcings.nc','w')
mem_dim = new_file.createDimension('Members',aice_ics.shape[0])
tim_dim = new_file.createDimension('Times',aice_forcings.shape[1])

var = new_file.createVariable('aice_ics','f8',('Members'))
var[:] = aice_ics
var2 = new_file.createVariable('aice_forcings','f8',('Members','Times'))
var2[:,:] = aice_forcings
##
var1 = new_file.createVariable('hi_ics','f8',('Members'))
var1[:] = hi_ics
var12 = new_file.createVariable('hi_forcings','f8',('Members','Times'))
var12[:,:] = hi_forcings
##
var2 = new_file.createVariable('hs_ics','f8',('Members'))
var2[:] = hs_ics
var22 = new_file.createVariable('hs_forcings','f8',('Members','Times'))
var22[:,:] = hs_forcings
new_file.close()
