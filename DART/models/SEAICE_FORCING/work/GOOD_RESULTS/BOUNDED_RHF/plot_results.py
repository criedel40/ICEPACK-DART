import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import xarray as xr
from datetime import datetime

data = xr.open_dataset('../../gridpt_cice_data.nc')
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

data = xr.open_dataset('../../output_control.nc')
truth_aice = data.aice.values[ind,-1]
control_aice = data.aice.values[ind,:-1]
truth_hi = data.hi.values[ind,-1]
control_hi = data.hi.values[ind,:-1]
data.close()

files = sorted(glob.glob('*/*.nc'))
output_aice = np.zeros((len(files),ind.shape[0],99))
output_hi = np.zeros((len(files),ind.shape[0],99))
for f in range(0,len(files)):
  print(files[f])
  data = xr.open_dataset(files[f])
  output_aice[f,:] = data.output_aice.values[ind,:]
  output_hi[f,:] = data.output_hi.values[ind,:]

######################
plt.rc('font', weight='bold')
tim = np.arange(0,dates.shape[0])
fig,ax = plt.subplots(5,1,sharex=True,figsize=(8,7))
for m in range(0,99):
  if m == 0:
    ax[0].plot(tim,output_aice[0,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5,label='Ind. Members')
    ax[1].plot(tim,output_aice[1,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5,label='Individual Members')
    ax[2].plot(tim,output_aice[2,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5,label='Individual Members')
    ax[3].plot(tim,output_aice[3,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5,label='Individual Members')
    ax[4].plot(tim,output_aice[4,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5,label='Individual Members')
  else:
    ax[0].plot(tim,output_aice[0,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
    ax[1].plot(tim,output_aice[1,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
    ax[2].plot(tim,output_aice[2,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
    ax[3].plot(tim,output_aice[3,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
    ax[4].plot(tim,output_aice[4,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
ax[0].plot(tim,output_aice[0,:,:].mean(axis=1),'-k',linewidth=1.5,label='Ensemble Mean')
ax[1].plot(tim,output_aice[1,:,:].mean(axis=1),'-k',linewidth=1.5,label='Ensemble Mean')
ax[2].plot(tim,output_aice[2,:,:].mean(axis=1),'-k',linewidth=1.5,label='Ensemble Mean')
ax[3].plot(tim,output_aice[3,:,:].mean(axis=1),'-k',linewidth=1.5,label='Ensemble Mean')
ax[4].plot(tim,output_aice[4,:,:].mean(axis=1),'-k',linewidth=1.5,label='Ensemble Mean')
###
ax[0].plot(tim,truth_aice,'-r',linewidth=1.5,label='Truth')
ax[1].plot(tim,truth_aice,'-r',linewidth=1.5,label='Truth')
ax[2].plot(tim,truth_aice,'-r',linewidth=1.5,label='Truth')
ax[3].plot(tim,truth_aice,'-r',linewidth=1.5,label='Truth')
ax[4].plot(tim,truth_aice,'-r',linewidth=1.5,label='Truth')

ax[0].set_title('State Vars:AICE - Obs: AICE',loc='left',weight='bold',fontsize=12,pad=1.1)
ax[1].set_title('State Vars:AICE,HI - Obs: AICE',loc='left',weight='bold',fontsize=12,pad=1.1)
ax[2].set_title('State Vars:AICE,HI - Obs: AICE,THICKNESS',loc='left',weight='bold',fontsize=12,pad=1.1)
ax[3].set_title('State Vars:AICE,HI - Obs: THICKNESS',loc='left',weight='bold',fontsize=12,pad=1.1)
ax[4].set_title('State Vars:HI - Obs: THICKNESS',loc='left',weight='bold',fontsize=12,pad=1.1)
ax[0].legend(loc='upper center',ncol=3,fontsize=8,bbox_to_anchor=(0.71,1.285))
#ax[1].legend(loc='upper center',ncol=3,fontsize=8)#,bbox_to_anchor=(0.6,1.12))
#ax[2].legend(loc='upper center',ncol=3,fontsize=8)#,bbox_to_anchor=(0.6,1.12))
#ax[3].legend(loc='upper center',ncol=3,fontsize=8)#,bbox_to_anchor=(0.6,1.12))
#ax[4].legend(loc='upper center',ncol=3,fontsize=8)#,bbox_to_anchor=(0.6,1.12))
##
ax[4].set_xlim(tim[0],tim[-1])
ax[4].set_xticks(ind_mons)
ax[4].set_xticklabels(labels,fontsize=10,rotation=90)
for x in range(0,5):
  #ax[x].set_yticks(np.arange(0.2,1.0+0.2,0.2))
  ax[x].tick_params(axis='y', labelsize=9)
  #ax[x].set_ylim(0,6)
  ax[x].grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3,zorder=1)
  ax[x].set_ylabel('Concen.',weight='bold',fontsize=10)
tit = plt.suptitle('Idealized Sea Ice OSSE\nTrunc. Gauss. Obs Error Distribution',weight='bold',fontsize=12,y=0.97)
#fig.subplots_adjust(hspace=0.35)
plt.savefig('idealized_si_BoundedRHF_aice.jpg',dpi=550,bbox_inches='tight')
os.system('open idealized_si_BoundedRHF_aice.jpg')
##############################
##############################
tim = np.arange(0,dates.shape[0])
fig,ax = plt.subplots(5,1,sharex=True,figsize=(8,7))
for m in range(0,99):
  if m == 0:
    ax[0].plot(tim,output_hi[0,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5,label='Ind. Members')
    ax[1].plot(tim,output_hi[1,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5,label='Individual Members')
    ax[2].plot(tim,output_hi[2,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5,label='Individual Members')
    ax[3].plot(tim,output_hi[3,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5,label='Individual Members')
    ax[4].plot(tim,output_hi[4,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5,label='Individual Members')
  else:
    ax[0].plot(tim,output_hi[0,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
    ax[1].plot(tim,output_hi[1,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
    ax[2].plot(tim,output_hi[2,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
    ax[3].plot(tim,output_hi[3,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
    ax[4].plot(tim,output_hi[4,:,m],'-',color='darkgrey',linewidth=0.5,alpha=0.5)
ax[0].plot(tim,output_hi[0,:,:].mean(axis=1),'-k',linewidth=1.5,label='Ensemble Mean')
ax[1].plot(tim,output_hi[1,:,:].mean(axis=1),'-k',linewidth=1.5,label='Ensemble Mean')
ax[2].plot(tim,output_hi[2,:,:].mean(axis=1),'-k',linewidth=1.5,label='Ensemble Mean')
ax[3].plot(tim,output_hi[3,:,:].mean(axis=1),'-k',linewidth=1.5,label='Ensemble Mean')
ax[4].plot(tim,output_hi[4,:,:].mean(axis=1),'-k',linewidth=1.5,label='Ensemble Mean')
###
ax[0].plot(tim,truth_hi,'-r',linewidth=1.5,label='Truth')
ax[1].plot(tim,truth_hi,'-r',linewidth=1.5,label='Truth')
ax[2].plot(tim,truth_hi,'-r',linewidth=1.5,label='Truth')
ax[3].plot(tim,truth_hi,'-r',linewidth=1.5,label='Truth')
ax[4].plot(tim,truth_hi,'-r',linewidth=1.5,label='Truth')

ax[0].set_title('State Vars:AICE - Obs: AICE',loc='left',weight='bold',fontsize=12,pad=1.1)
ax[1].set_title('State Vars:AICE,HI - Obs: AICE',loc='left',weight='bold',fontsize=12,pad=1.1)
ax[2].set_title('State Vars:AICE,HI - Obs: AICE,THICKNESS',loc='left',weight='bold',fontsize=12,pad=1.1)
ax[3].set_title('State Vars:AICE,HI - Obs: THICKNESS',loc='left',weight='bold',fontsize=12,pad=1.1)
ax[4].set_title('State Vars:HI - Obs: THICKNESS',loc='left',weight='bold',fontsize=12,pad=1.1)
ax[0].legend(loc='upper center',ncol=3,fontsize=8,bbox_to_anchor=(0.71,1.285))
#ax[1].legend(loc='upper center',ncol=3,fontsize=8)#,bbox_to_anchor=(0.6,1.12))
#ax[2].legend(loc='upper center',ncol=3,fontsize=8)#,bbox_to_anchor=(0.6,1.12))
#ax[3].legend(loc='upper center',ncol=3,fontsize=8)#,bbox_to_anchor=(0.6,1.12))
#ax[4].legend(loc='upper center',ncol=3,fontsize=8)#,bbox_to_anchor=(0.6,1.12))
##
ax[4].set_xlim(tim[0],tim[-1])
ax[4].set_xticks(ind_mons)
ax[4].set_xticklabels(labels,fontsize=10,rotation=90)
for x in range(0,5):
  #ax[x].set_yticks(np.arange(0.2,1.0+0.2,0.2))
  ax[x].tick_params(axis='y', labelsize=9)
  #ax[x].set_ylim(0,6)
  ax[x].grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3,zorder=1)
  ax[x].set_ylabel('Volume',weight='bold',fontsize=10)
tit = plt.suptitle('Idealized Sea Ice OSSE\nTrunc. Gauss. Obs Error Distribution',weight='bold',fontsize=12,y=0.97)
#fig.subplots_adjust(hspace=0.35)
plt.savefig('idealized_si_BoundedRHF_hi.jpg',dpi=550,bbox_inches='tight')
os.system('open idealized_si_BoundedRHF_hi.jpg')
