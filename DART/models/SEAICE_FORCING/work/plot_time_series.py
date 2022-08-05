import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import sys


try:
  obs_err = sys.argv[1]
  obs_dist = sys.argv[2]
  filter_type = sys.argv[3]
except:
  print('WARNING!!!!')
  print('You need to input region and comparison type on command line')
  print('ex: python diag.py 1')
  sys.exit()

data = xr.open_dataset('gridpt_cice_data.nc')
dates = data.times.values
years = np.array([int(str(d)[:4]) for d in dates])
days = np.array([int(str(d)[6:8]) for d in dates])
mons = np.array([int(str(d)[4:6]) for d in dates])
data.close()
ind_years = np.where((days==1)&(mons==1))[0][:-1]
ind_years = np.insert(ind_years,0,0)
labels = years[ind_years]

data = xr.open_dataset('output_control.nc')
aice_control = data.aice.values
data.close()
##
data = xr.open_dataset('output_assim.nc')
aice_assim = data.aice.values
obs = data.obs.values
data.close()

tim = np.arange(0,aice_control.shape[0])
fig,ax = plt.subplots(2,1,sharex=True,figsize=(8,6))
for m in range(0,100):
  ax[0].plot(tim,aice_control[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.75)
ax[0].plot(tim,aice_control[:,:].mean(axis=1),'-k',linewidth=1.5)
ax[0].plot(tim,aice_control[:,-1],'-r',linewidth=1.5)
ax[0].set_xlim(tim[0],tim[-1])
###
for m in range(0,99):
  ax[1].plot(tim,aice_assim[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.75)
ax[1].plot(tim,aice_assim[:,:].mean(axis=1),'-k',linewidth=1.5)
ax[1].plot(tim,aice_control[:,-1],'-r',linewidth=1.5)
ax[1].set_xlim(tim[0],tim[-1])
ax[1].set_xticks(ind_years)
ax[1].set_xticklabels(labels,weight='bold')
plt.savefig('results_{0}_{1}_{2}.jpg'.format(filter_type,obs_err,obs_dist),dpi=550,bbox_inches='tight')
#####
tim = np.arange(0,aice_control.shape[0])
fig,ax = plt.subplots(2,1,sharex=True,figsize=(8,6))
ax[0].plot(tim,aice_control[:,:].mean(axis=1),'-k',linewidth=1.5)
ax[0].plot(tim,aice_control[:,-1],'-r',linewidth=1.5)
ax[0].set_xlim(tim[0],tim[-1])
###
ax[1].plot(tim,aice_assim[:,:].mean(axis=1),'-k',linewidth=1.5)
ax[1].plot(tim,aice_control[:,-1],'-r',linewidth=1.5)
ax[1].set_xlim(tim[0],tim[-1])
ax[1].set_xticks(ind_years)
ax[1].set_xticklabels(labels,weight='bold')
plt.savefig('results_mean_{0}_{1}_{2}.jpg'.format(filter_type,obs_err,obs_dist),dpi=550,bbox_inches='tight')
####
fig,ax = plt.subplots(1,1,sharex=True,figsize=(8,6))
ax.plot(tim,aice_control[:,:-1].mean(axis=1)-aice_control[:,-1],'-k',linewidth=1.5)
ax.plot(tim,aice_assim[:,:].mean(axis=1)-aice_control[:,-1],'-r',linewidth=1.5)
#ax[0].plot(tim,aice_control[:,-1],'-r',linewidth=1.5)
ax.set_xlim(tim[0],tim[-1])
###
#ax[1].plot(tim,aice_assim[:,:].mean(axis=1)-aice_control[:,-1],'-k',linewidth=1.5)
##ax[1].plot(tim,aice_control[:,-1],'-r',linewidth=1.5)
ax.set_xticks(ind_years)
ax.set_xticklabels(labels,weight='bold')
plt.savefig('results_bias_{0}_{1}_{2}.jpg'.format(filter_type,obs_err,obs_dist),dpi=550,bbox_inches='tight')


