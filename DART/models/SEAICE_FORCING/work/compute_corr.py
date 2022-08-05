import xarray as xr
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np
import os

data = xr.open_dataset('gridpt_cice_data.nc')
times = data.times.values
times = times[::365]

data = xr.open_dataset('output_control.nc')
aice = data.aice.values
hi = np.divide(data.hi.values,aice,where=aice!=0)


spearcorr = np.zeros((aice.shape[0]))
pearsoncorr = np.zeros((aice.shape[0]))

for t in range(0,aice.shape[0]):
  spearcorr[t] = scipy.stats.spearmanr(aice[t,:],hi[t,:],axis=None)[0]
  pearsoncorr[t] = scipy.stats.pearsonr(aice[t,:],hi[t,:])[0]



fig = plt.figure(figsize=(10,6))
plt.plot(spearcorr,label='Spear')
plt.plot(pearsoncorr,label='Pearson')
plt.xticks(np.arange(0,aice.shape[0])[::365],times,rotation=90)
plt.savefig('test.jpg',dpi=550)
os.system('open test.jpg')
