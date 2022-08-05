import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import fortran_funcs
import os

x = np.arange(0,1+0.001,0.001)
a_array = np.array([0.5,1,8,3,5])
b_array = np.array([5,3,8,1,0.5])

scipy_results = np.zeros((a_array.shape[0],x.shape[0]))
minef_results = np.zeros((a_array.shape[0],x.shape[0]))

for i in range(0,a_array.shape[0]):
  for j in range(0,x.shape[0]):
    scipy_results[i,j] = scipy.stats.beta.ppf(x[j],a_array[i],b_array[i])
    minef_results[i,j] = fortran_funcs.betaincinv(x[j],a_array[i],b_array[i])

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
ax[1].set_xlabel('Quantiles',weight='bold')
ax[0].set_ylabel('Inverse CDF',weight='bold')
ax[1].set_ylabel('Inv. CDF Differences',weight='bold')
ax[0].set_title('Beta Inverse CDFs - Solid Line: Scipy - Dashed line:My Coded Funcs.',weight='bold',fontsize=8)
ax[1].set_title('Inverse CDF differences between Scipy and My Functions',weight='bold',fontsize=8)
plt.savefig('compare_invcdf.jpg',dpi=550,bbox_inches='tight')
os.system('open compare_invcdf.jpg')


