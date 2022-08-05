import xarray as xr
import netCDF4
import DART_functions
import numpy as np
import os
import matplotlib.pyplot as plt
import sys

#Get truth
data = xr.open_dataset('../../output_control.nc')
truth = data.aice.values[:,-1]
data.close()
#Get ensemble data
data = xr.open_dataset('assim_output_EaKF_beta.nc')
prior = data.input_aice.values
obs = data.obs.values
obserr = data.obs_var.values
data.close()

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
    bias_truth[i],bias_obs[i],rmse_truth[i],rmse_obs[i],totspr_truth[i],totspr_obs[i],rank_truth[i,:],rank_obs[i,:] = DART_functions.compute_obsspace_diag(prior[i,:].mean(),prior[i,:].var(ddof=1),truth[i],obs[i],obserr[i],prior[i,:],99,i,3)
    #sys.exit()


new_file = netCDF4.Dataset('DART_obsspace_diags_EaKF_beta.nc','w')
tim_dim = new_file.createDimension('Times',bias_obs.shape[0])
rank_dim = new_file.createDimension('ranks',rank_truth.shape[1])

var = new_file.createVariable('bias_truth','f8',('Times'))
var[:] = bias_truth
##
var1 = new_file.createVariable('bias_obs','f8',('Times'))
var1[:] = bias_obs
##
var2 = new_file.createVariable('rmse_truth','f8',('Times'))
var2[:] = rmse_truth
##
var3 = new_file.createVariable('rmse_obs','f8',('Times'))
var3[:] = rmse_truth
##
var4 = new_file.createVariable('totspread_truth','f8',('Times'))
var4[:] = totspr_truth
##
var5 = new_file.createVariable('totspread_obs','f8',('Times'))
var5[:] = totspr_obs
##
var6 = new_file.createVariable('rank_truth','f8',('Times','ranks'))
var6[:,:] = rank_truth
##
var7 = new_file.createVariable('rank_obs','f8',('Times','ranks'))
var7[:,:] = rank_obs
##
var6 = new_file.createVariable('truth','f8',('Times'))
var6[:] = truth
##
var6 = new_file.createVariable('obs','f8',('Times'))
var6[:] = obs
new_file.close()
sys.exit()
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
