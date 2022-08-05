import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import xarray as xr
from datetime import datetime

data = xr.open_dataset('../gridpt_cice_data.nc')
dates = data.times.values
years = np.array([int(str(d)[:4]) for d in dates])
days = np.array([int(str(d)[6:8]) for d in dates])
mons = np.array([int(str(d)[4:6]) for d in dates])
data.close()
ind_years = np.where((days==1)&(mons==1))[0][:-1]
ind_years = np.insert(ind_years,0,0)
labels = years[ind_years]

ind = np.where((years>=2006)&(years<=2007))[0]
dates = dates[ind]
days = days[ind]
mons = mons[ind]
ind_mons = np.where(days==1)[0]
labels = []
for x in ind_mons:
  name = datetime(years[x],mons[x],days[x]).strftime('%b')
  labels.append(name)

data = xr.open_dataset('../output_control.nc')
truth_aice = data.aice.values[ind,-1]
control_aice = data.aice.values[ind,:-1]
truth_hi = data.hi.values[ind,-1]
control_hi = data.hi.values[ind,:-1]
data.close()

#files1 = sorted(glob.glob('*/assim_output_*.nc'))
#files2 = sorted(glob.glob('RHF/assim_output_RHF_*.nc'))
#files3 = sorted(glob.glob('BOUNDED_RHF/assim_output_BoundedRHF_*.nc'))

data = xr.open_dataset('assim_output_EaKF_gaus.nc')
output_aice = data.output_aice.values[ind,:]
output_hi = data.output_hi.values[ind,:]
input_aice = data.input_aice.values[ind,:]
input_hi = data.input_hi.values[ind,:]
######################
plt.rc('font', weight='bold')
tim = np.arange(0,dates.shape[0])
fig,ax = plt.subplots(2,1,sharex=True,figsize=(8,6))
for m in range(0,99):
  if m == 0:
    ax[0].plot(tim,control_aice[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5,label='Individual Members')
  else:
    ax[0].plot(tim,control_aice[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
ax[0].plot(tim,control_aice.mean(axis=1),'-k',linewidth=1.5,label='Ensemble Mean')
ax[0].plot(tim,truth_aice,'-r',linewidth=1.5,label='Truth')
ax[0].set_title('Free Forecast',loc='left',weight='bold',fontsize=12,pad=1.1)
ax[0].legend(loc='upper center',ncol=3,fontsize=8,bbox_to_anchor=(0.6,1.12))
####
#
for m in range(0,99):
  ax[1].plot(tim,output_aice[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
ax[1].plot(tim,output_aice[:,:].mean(axis=1),'-k',linewidth=1.5)
ax[1].plot(tim,truth_aice,'-r',linewidth=1.5)
ax[1].set_title('SI Concentration',loc='left',weight='bold',fontsize=12,pad=1.9)
##
ax[1].set_xlim(tim[0],tim[-1])
ax[1].set_xticks(ind_mons)
ax[1].set_xticklabels(labels,fontsize=10,rotation=90)
for x in range(0,2):
  ax[x].set_yticks(np.arange(0.2,1.0+0.2,0.2))
  ax[x].tick_params(axis='y', labelsize=9)
  ax[x].set_ylim(0,1.05)
  ax[x].grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3,zorder=1)
  ax[x].set_ylabel('Concentration',weight='bold',fontsize=10)
tit = plt.suptitle('Idealized Sea Ice OSSE\nTrunc. Gauss. Obs Error Distribution',weight='bold',fontsize=12,y=0.97)
#fig.subplots_adjust(hspace=0.35)
plt.savefig('idealized_si_BoundedRHF_aice.jpg',dpi=550,bbox_inches='tight')
os.system('open idealized_si_BoundedRHF_aice.jpg')
####################
plt.rc('font', weight='bold')
tim = np.arange(0,dates.shape[0])
fig,ax = plt.subplots(2,1,sharex=True,figsize=(8,6))
for m in range(0,99):
  if m == 0:
    ax[0].plot(tim,control_hi[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5,label='Individual Members')
  else:
    ax[0].plot(tim,control_hi[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
ax[0].plot(tim,control_hi.mean(axis=1),'-k',linewidth=1.5,label='Ensemble Mean')
ax[0].plot(tim,truth_hi,'-r',linewidth=1.5,label='Truth')
ax[0].set_title('Free Forecast',loc='left',weight='bold',fontsize=12,pad=1.1)
ax[0].legend(loc='upper center',ncol=3,fontsize=8,bbox_to_anchor=(0.6,1.12))
####
#
for m in range(0,99):
  ax[1].plot(tim,output_hi[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
ax[1].plot(tim,output_hi[:,:].mean(axis=1),'-k',linewidth=1.5)
ax[1].plot(tim,truth_hi,'-r',linewidth=1.5)
ax[1].set_title('SI Volume',loc='left',weight='bold',fontsize=12,pad=1.9)
##
ax[1].set_xlim(tim[0],tim[-1])
ax[1].set_xticks(ind_mons)
ax[1].set_xticklabels(labels,fontsize=10,rotation=90)
for x in range(0,2):
  #ax[x].set_yticks(np.arange(0.2,1.0+0.2,0.2))
  ax[x].tick_params(axis='y', labelsize=9)
  #ax[x].set_ylim(0,6)
  ax[x].grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3,zorder=1)
  ax[x].set_ylabel('Volume',weight='bold',fontsize=10)
tit = plt.suptitle('Idealized Sea Ice OSSE\nTrunc. Gauss. Obs Error Distribution',weight='bold',fontsize=12,y=0.97)
#fig.subplots_adjust(hspace=0.35)
plt.savefig('idealized_si_BoundedRHF_hi.jpg',dpi=550,bbox_inches='tight')
os.system('open idealized_si_BoundedRHF_hi.jpg')
