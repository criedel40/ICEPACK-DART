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
truth = data.aice.values[ind,-1]
control = data.aice.values[ind,:-1]
data.close()

data = xr.open_mfdataset('assim_output_EaKF_noinf.nc')
aice_eakf_noinf = data.output_aice.values[ind,:]
data.close()
##
data = xr.open_mfdataset('assim_output_EaKF_inf.nc')
aice_eakf_inf = data.output_aice.values[ind,:]
data.close()
##
data = xr.open_mfdataset('assim_output_BoundedRHF_noinf.nc')
aice_bdrhf_noinf = data.output_aice.values[ind,:]
data.close()
##
data = xr.open_mfdataset('assim_output_BoundedRHF_inf.nc')
aice_bdrhf_inf = data.output_aice.values[ind,:]
data.close()

plt.rc('font', weight='bold')
tim = np.arange(0,dates.shape[0])
fig,ax = plt.subplots(5,1,sharex=True,figsize=(8,8))
for m in range(0,99):
  if m == 0:
    ax[0].plot(tim,control[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5,label='Individual Members')
  else:
    ax[0].plot(tim,control[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
ax[0].plot(tim,control.mean(axis=1),'-k',linewidth=1.5,label='Ensemble Mean')
ax[0].plot(tim,truth,'-r',linewidth=1.5,label='Truth')
ax[0].set_title('Control',loc='left',weight='bold',fontsize=8,pad=1.1)
ax[0].legend(loc='upper center',ncol=3,fontsize=7,bbox_to_anchor=(0.5,1.225))
####
for m in range(0,99):
  ax[1].plot(tim,aice_eakf_noinf[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
ax[1].plot(tim,aice_eakf_noinf.mean(axis=1),'-k',linewidth=1.5)
ax[1].plot(tim,truth,'-r',linewidth=1.5)
ax[1].set_title('EaKF - Gauss. Prior Inflation Off',loc='left',weight='bold',fontsize=8,pad=1.1)
##
for m in range(0,99):
  ax[2].plot(tim,aice_eakf_inf[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
ax[2].plot(tim,aice_eakf_inf.mean(axis=1),'-k',linewidth=1.5)
ax[2].plot(tim,truth,'-r',linewidth=1.5)
ax[2].set_title('EaKF - Gauss. Prior Inflation On',loc='left',weight='bold',fontsize=8,pad=1.9)
##
for m in range(0,99):
  ax[3].plot(tim,aice_bdrhf_noinf[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
ax[3].plot(tim,aice_bdrhf_noinf.mean(axis=1),'-k',linewidth=1.5)
ax[3].plot(tim,truth,'-r',linewidth=1.5)
ax[3].set_title('Bounded Rank Histogram Filter - Beta Prior Inflation Off',loc='left',weight='bold',fontsize=8,pad=1.1)
##
for m in range(0,99):
  ax[4].plot(tim,aice_bdrhf_inf[:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
ax[4].plot(tim,aice_bdrhf_inf.mean(axis=1),'-k',linewidth=1.5)
ax[4].plot(tim,truth,'-r',linewidth=1.5)
ax[4].set_title('Bounded Rank Histogram Filter - Beta Prior Inflation On',loc='left',weight='bold',fontsize=8,pad=1.4)
####
ax[4].set_xlim(tim[0],tim[-1])
ax[4].set_xticks(ind_mons)
ax[4].set_xticklabels(labels,fontsize=9,rotation=90)
for x in range(0,5):
  ax[x].set_yticks(np.arange(0.2,1.0+0.2,0.2))
  ax[x].set_ylim(0,1.05)
  ax[x].grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3,zorder=1)
  ax[x].set_ylabel('Concentration',weight='bold',fontsize=8)
tit = plt.suptitle('Idealized Sea Ice OSSE - Obs Generated using Truncated Normal Dist.\nObs Err = truth*0.15',weight='bold',fontsize=10,y=0.95)
plt.savefig('idealized_si_varyinflation.jpg',dpi=550,bbox_inches='tight')
os.system('open idealized_si_varyinflation.jpg')
############################
