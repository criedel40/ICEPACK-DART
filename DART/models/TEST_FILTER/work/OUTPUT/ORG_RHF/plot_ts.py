import numpy as np
import matplotlib.pyplot as plt


mems = np.arange(1,81)

mem_data = np.array(open('fort.42').readlines()[:-1])
total_ens = []
for m in range(0,len(mem_data)):
  if mem_data[m] == ' ##############\n':
    continue
  else:
    line = list(map(float,mem_data[m].split()))
    total_ens.append(line)
total_ens = np.array(total_ens)

sep_data = np.zeros((5000,80))
for m in range(0,mems.shape[0]):
  ind = np.where(total_ens[:,0]==mems[m])[0]
  sep_data[:,m] = total_ens[ind,1]

obs_data = np.array(open('fort.41').readlines())
obs_array = np.array([float(d.split()[-1]) for d in obs_data])

mean = np.nanmean(sep_data,axis=1)


truth = np.zeros((5000))
hold = np.zeros((5000))
hold[:] = obs_array.mean()
truth[:] = 0.9899
tim = np.arange(0,5000)
plt.rc('font', weight='bold')
fig = plt.figure(figsize=(10,6))
plt.plot(tim,hold,'-r',label='Average Observation')
plt.plot(tim,truth,'-k',label='Truth',linewidth=3.0)
#for m in range(0,80):
#  plt.plot(tim,sep_data[:,m],'-',color='grey',alpha=0.5,linewidth=0.25)
plt.plot(tim,mean,'-b',label='Prior Ens. Mean',linewidth=3.0)
plt.xlim(tim[0]-0.2,tim[-1])
plt.legend(loc='best')
plt.ylim(0.7,1.001)
plt.xlabel('Cycling Times',weight='bold')
plt.title('Sea Ice Concentration - Orginal Rank Histogram Filter',weight='bold')
plt.savefig('org_rhf_simple.jpg',dpi=550,bbox_inches='tight')



