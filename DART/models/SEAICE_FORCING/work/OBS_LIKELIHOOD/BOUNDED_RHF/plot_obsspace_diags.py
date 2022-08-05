import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.offsetbox import AnchoredText

data = xr.open_dataset('DART_obsspace_diags_BoundedRHF_beta.nc')
truth = data.truth.values
obs = data.obs.values
rank_truth = data.rank_truth.values
rank_obs = data.rank_obs.values

textstr ='Filter Type: Bounded RHF\nObs Error Dist.: Beta\nObs Err: Truth*0.15\nPrior Inflation: Off'
save_name = 'rank_hist_beta_boundedRHF'

bin_centers = np.arange(0.05,0.95+0.1,0.1)
bin_edges = np.zeros((bin_centers.shape[0],2))
for b in range(0,bin_centers.shape[0]):
  bin_edges[b,0] = bin_centers[b] - 0.05
  bin_edges[b,1] = bin_centers[b] + 0.05


rank_truth_sort = np.zeros((bin_centers.shape[0],100))
rank_obs_sort = np.zeros((bin_centers.shape[0],100))
ave_obs = np.zeros((bin_centers.shape[0]))
for b in range(0,bin_centers.shape[0]):
  ind = np.where((truth>=bin_edges[b,0])&(truth<bin_edges[b,1]))[0]
  rank_truth_sort[b,:] = np.nansum(rank_truth[ind,:],axis=0)
  rank_obs_sort[b,:] = np.nansum(rank_obs[ind,:],axis=0)
  ave_obs[b] = np.nanmean(obs[ind])

bar_width = 0.75
plt.rc('font', weight='bold')
fig,ax = plt.subplots(3,4,sharex=False,sharey=False)
count = 0
for a in ax.flatten():
  if count == 10:
    a.bar(np.arange(1,101),rank_truth.sum(axis=0),bar_width,color='blue',align='center')
    a.set_xlim(1-2,100+2)
    a.set_xticks(np.insert(np.arange(20,100+20,20),0,1))
    a.tick_params(axis='y', labelsize=6)
    a.tick_params(axis='x', labelsize=5)
    a.set_title('All Bins',fontsize=7,weight='bold',pad=0.95)
    break
  elif count == 11:
    break
  a.bar(np.arange(1,101),rank_truth_sort[count,:],bar_width,color='b',align='center')
  a.set_xlim(1-2,100+2)
  a.set_xticks(np.insert(np.arange(20,100+20,20),0,1))
  a.set_title('{0}<Truth<{1}'.format(round(bin_edges[count,0],2),round(bin_edges[count,1],2)),fontsize=7,weight='bold',pad=0.95)
  a.tick_params(axis='y', labelsize=6)
  a.tick_params(axis='x', labelsize=5)
  count += 1

ax[-1,-1].axis('off')
ax[-1,-1].text(0.0,0.92,'Experiment Details',va='top',fontsize=8)
ax[-1,-1].text(0.0,0.8,textstr,va='top',fontsize=8)
plt.suptitle('Truth Rank Histogram',weight='bold',y=1.025)
plt.tight_layout()
plt.savefig('{0}_truth.jpg'.format(save_name),dpi=550,bbox_inches='tight')
os.system('open {0}_truth.jpg'.format(save_name))
############
plt.rc('font', weight='bold')
fig,ax = plt.subplots(3,4,sharex=False,sharey=False)
count = 0
for a in ax.flatten():
  if count == 10:
    a.bar(np.arange(1,101),rank_obs.sum(axis=0),bar_width,color='b',align='center')
    a.set_xlim(1-2,100+2)
    a.set_xticks(np.insert(np.arange(20,100+20,20),0,1))
    a.tick_params(axis='y', labelsize=6)
    a.tick_params(axis='x', labelsize=5)
    a.set_title('All Bins',fontsize=7,weight='bold',pad=0.95)
    break
  elif count == 11:
    break
  a.bar(np.arange(1,101),rank_obs_sort[count,:],bar_width,color='b',align='center')
  anchored_text = AnchoredText('Ave. Obs. Value:{0}'.format(round(ave_obs[count],2)), loc=2,prop=dict(fontsize=5),frameon=True)
  a.add_artist(anchored_text)
  a.set_xlim(1-2,100+2)
  a.set_xticks(np.insert(np.arange(20,100+20,20),0,1))
  a.set_title('{0}<Truth<{1}'.format(round(bin_edges[count,0],2),round(bin_edges[count,1],2)),fontsize=7,weight='bold',pad=0.95)
  a.tick_params(axis='y', labelsize=6)
  a.tick_params(axis='x', labelsize=5)
  count += 1

ax[-1,-1].axis('off')
ax[-1,-1].text(0.0,0.92,'Experiment Details',va='top',fontsize=8)
ax[-1,-1].text(0.0,0.8,textstr,va='top',fontsize=8)
#ax[-1,-2].axis('off')
plt.suptitle('Observation Rank Histogram',weight='bold',y=1.025)
plt.tight_layout()
plt.savefig('{0}_obs.jpg'.format(save_name),dpi=550,bbox_inches='tight')
os.system('open {0}_obs.jpg'.format(save_name))
###
comd = 'convert +append {0}_truth.jpg {0}_obs.jpg complete_{0}.jpg'.format(save_name)
os.system(comd)

