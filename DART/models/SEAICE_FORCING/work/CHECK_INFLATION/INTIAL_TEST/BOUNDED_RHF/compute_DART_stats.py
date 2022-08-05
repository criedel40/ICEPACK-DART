import xarray as xr
import netCDF4
import DART_functions
import numpy as np
import os
import matplotlib.pyplot as plt

#Get truth
data = xr.open_dataset('output_control.nc')
truth = data.aice.values[:,-1]
data.close()
#Get ensemble data
data = xr.open_dataset('assim_output_BoundedRHF_noinf.nc')
prior = data.input_aice.values
obs = data.obs.values
data.close()
obserr = np.zeros((obs.shape[0]))
for i in range(0,obs.shape[0]):
  if (truth[i]<0.01):
    obserr[i] = (0.01*0.15)**2
  else:
    obserr[i] = (truth[i]*0.15)**2
   

bias_truth = np.zeros((prior.shape[0]))
bias_obs = np.zeros((prior.shape[0]))
##
rmse_truth = np.zeros((prior.shape[0]))
rmse_obs = np.zeros((prior.shape[0]))
##
totspr_truth = np.zeros((prior.shape[0]))
totspr_obs = np.zeros((prior.shape[0]))
##
rank_truth = np.zeros((prior.shape[0],100))
rank_obs = np.zeros((prior.shape[0],100))


for i in range(0,prior.shape[0]):
    bias_truth[i],bias_obs[i],rmse_truth[i],rmse_obs[i],totspr_truth[i],totspr_obs[i],rank_truth[i,:],rank_obs[i,:] = DART_functions.compute_obsspace_diag(prior[i,:].mean(),prior[i,:].var(ddof=1),truth[i],obs[i],obserr[i],prior[i,:],99,i)


tim = np.arange(0,prior.shape[0])
fig,ax = plt.subplots(2,1,sharex=True,figsize=(8,6))
ax[0].plot(tim,bias_truth,'-k',label='Bias')
#ax[0].plot(tim,rmse_truth,'-r',label='RMSE')
#ax[0].plot(tim,totspr_truth,'-b',label='Total Spread')
ax[0].plot(tim,totspr_truth/rmse_truth,'-r',label='Total Spread')
##
ax[1].plot(tim,bias_obs,'-k',label='Bias')
ax[1].plot(tim,rmse_obs,'-r',label='RMSE')
ax[1].plot(tim,totspr_obs,'-b',label='Total Spread')
##
ax[1].set_xlim(tim[0],tim[-1])
plt.savefig('time_series.jpg',dpi=550,bbox_inches='tight')
os.system('open time_series.jpg')
##

tim = np.arange(0,prior.shape[0])
fig,ax = plt.subplots(2,1,sharex=True,figsize=(6,6))
ax[0].bar(np.arange(1,101),rank_truth.sum(axis=0))
##
ax[1].bar(np.arange(1,101),rank_obs.sum(axis=0))
##
ax[1].set_xlim(1-0.5,100+0.5)
plt.savefig('rank_hist.jpg',dpi=550,bbox_inches='tight')
os.system('open rank_hist.jpg')
