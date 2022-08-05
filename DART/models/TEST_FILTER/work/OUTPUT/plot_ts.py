import numpy as np
import matplotlib.pyplot as plt


mems = np.arange(1,81)
###########################
mem_data = np.array(open('ORG_RHF/fort.42').readlines()[:-1])
total_ens = []
for m in range(0,len(mem_data)):
  if mem_data[m] == ' ##############\n':
    continue
  else:
    line = list(map(float,mem_data[m].split()))
    total_ens.append(line)
total_ens = np.array(total_ens)

sep_data_org = np.zeros((5000,80))
for m in range(0,mems.shape[0]):
  ind = np.where(total_ens[:,0]==mems[m])[0]
  sep_data_org[:,m] = total_ens[ind,1]

obs_data = np.array(open('ORG_RHF/fort.41').readlines())
obs_array_org = np.array([float(d.split()[-1]) for d in obs_data])

mean_org = np.nanmean(sep_data_org,axis=1)
#############################
mem_data = np.array(open('BOUNDED_RHF/fort.42').readlines()[:-1])
total_ens = []
for m in range(0,len(mem_data)):
  if mem_data[m] == ' ##############\n':
    continue
  else:
    line = list(map(float,mem_data[m].split()))
    total_ens.append(line)
total_ens = np.array(total_ens)

sep_data_bounded = np.zeros((5000,80))
for m in range(0,mems.shape[0]):
  ind = np.where(total_ens[:,0]==mems[m])[0]
  sep_data_bounded[:,m] = total_ens[ind,1]

obs_data = np.array(open('BOUNDED_RHF/fort.41').readlines())
obs_array_bounded = np.array([float(d.split()[-1]) for d in obs_data])

mean_bounded = np.nanmean(sep_data_bounded,axis=1)
############################


truth = np.zeros((5000))
hold = np.zeros((5000))
hold[:] = obs_array_org.mean()
truth[:] = 0.9899
tim = np.arange(0,5000)
plt.rc('font', weight='bold')
fig,ax = plt.subplots(2,1,sharex=True,figsize=(8,10))
#################
ax[0].plot(tim,hold,'-r',label='Average Observation')
ax[0].plot(tim,truth,'-k',label='Truth',linewidth=3.0)
#for m in range(0,80):
#  plt.plot(tim,sep_data[:,m],'-',color='grey',alpha=0.5,linewidth=0.25)
ax[0].plot(tim,mean_org,'-b',label='Prior Ens. Mean',linewidth=3.0)
ax[0].legend(loc='best')
ax[0].set_title('Sea Ice Concentration - Orginal Rank Histogram Filter',weight='bold')
##################
ax[1].plot(tim,hold,'-r',label='Average Observation')
ax[1].plot(tim,truth,'-k',label='Truth',linewidth=3.0)
#for m in range(0,80):
#  plt.plot(tim,sep_data[:,m],'-',color='grey',alpha=0.5,linewidth=0.25)
ax[1].plot(tim,mean_bounded,'-b',label='Prior Ens. Mean',linewidth=3.0)
ax[1].legend(loc='best')
ax[1].set_title('Sea Ice Concentration - Bounded Rank Histogram Filter',weight='bold')
##################
ax[1].set_xlim(tim[0]-0.3,tim[-1])
plt.legend(loc='best')
ax[0].set_ylim(0.7,1.001)
ax[1].set_ylim(0.7,1.001)
plt.xlabel('Cycling Times',weight='bold')
fig.tight_layout()
plt.savefig('rhf_simple.jpg',dpi=550,bbox_inches='tight')



