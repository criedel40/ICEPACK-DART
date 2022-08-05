import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
#import fortran_funcs
import os
import xarray as xr

data = xr.open_dataset('output_beta_stats.nc')
x = data.x.values
minef_results = data.cdfs.values
data.close()

#x = np.arange(0,1+0.01,0.01)
a_array = np.array([0.5,1,8,3,5])
b_array = np.array([5,3,8,1,0.5])

scipy_results = np.zeros((a_array.shape[0],x.shape[0]))
#minef_results = np.zeros((a_array.shape[0],x.shape[0]))

for i in range(0,a_array.shape[0]):
  for j in range(0,x.shape[0]):
    scipy_results[i,j] = scipy.stats.beta.cdf(x[j],a_array[i],b_array[i])
    #minef_results[i,j] = fortran_funcs.betacdf(x[j],a_array[i],b_array[i])

plt.rc('font', weight='bold')
fig,ax = plt.subplots(2,1,sharex=True)
colors = ['k','r','b','y','c']
for i in range(0,a_array.shape[0]):
  ax[0].plot(x,scipy_results[i,:],'-',color=colors[i],label='betacdf(x,{0},{1})'.format(a_array[i],b_array[i]))
  ax[0].plot(x,minef_results[i,:],'--',color=colors[i])
##
  ax[1].plot(x,minef_results[i,:]-scipy_results[i,:],'-',color=colors[i],label='betacdf(x,{0},{1})'.format(a_array[i],b_array[i]))

ax[0].legend(loc='upper left',ncol=2,fontsize=5)
ax[1].set_xlim(0,1)
ax[1].set_xlabel('X',weight='bold')
ax[0].set_ylabel('CDF',weight='bold')
ax[1].set_ylabel('CDF Differences',weight='bold')
ax[0].set_title('Beta CDFs - Solid Line: Scipy - Dashed line:My Fortran Coded Funcs.',weight='bold',fontsize=8)
ax[1].set_title('CDF differences between Scipy and My Fortran Functions',weight='bold',fontsize=8)
plt.savefig('compare_cdf.jpg',dpi=550,bbox_inches='tight')
os.system('open compare_cdf.jpg')


