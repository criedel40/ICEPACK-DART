import matplotlib.pyplot as plt
import numpy as np

data = open('fort.41').readlines()
data2 = open('fort.42').readlines()[:-1]
data3 = open('fort.42').readlines()[-1]


members = []
mem_data = []
for l in range(0,len(data2)):
  line = data2[l].split()
  members.append(int(line[0]))
  mem_data.append(float(line[1]))

members = np.array(members)
mem_data = np.array(mem_data)
ens_data = np.zeros((80,10000))
mems = np.arange(1,81)
for i in range(0,mems.shape[0]):
  ind = np.where(mems[i]==members)[0]
  ens_data[i,:] = mem_data[ind]

mean = []
var = []

for l in data:
  line = l.split()
  #print(line)
  mean.append(float(line[0]))
  

truth = np.zeros((len(mean)))
truth[:] = float(data3.split()[1])

plt.plot(mean,'-k')
plt.plot(truth,'-r')
max_truth = truth[0]+truth[0]*0.2
if max_truth > 1.0:
  max_truth = 1.001
min_truth = truth[0]-truth[0]*0.2
if min_truth < 0.0:
  min_truth = -0.001
plt.ylim(min_truth,max_truth)
plt.show()

tim = np.arange(0,ens_data.shape[1])
fig = plt.figure()
ax = plt.gca()
for m in range(0,ens_data.shape[0]):
  if m == 0:
    ax.plot(tim,ens_data[m,:],'-',color='grey',alpha=0.9,label='Individual Members',linewidth=0.5)
  else:
    ax.plot(tim,ens_data[m,:],'-',color='grey',alpha=0.9,linewidth=0.5)
ax.plot(tim,np.nanmean(ens_data,axis=0),'-',color='black',linewidth=2.0,label='Prior Ens. Mean')
ax.plot(tim,truth,'-',color='red',linewidth=3.0,label='Truth')
ax.legend(loc='best')
max_truth = truth[0]+truth[0]*0.2
if max_truth > 1.0:
  max_truth = 1.01
min_truth = truth[0]-truth[0]*0.2
if min_truth < 0.0:
  min_truth = -0.001
ax.set_ylim(min_truth,max_truth)
ax.set_xlim(tim[0],tim[-1])
plt.show()


caf = plt.contourf(np.arange(1,81),np.arange(0,10000),ens_data.T,np.linspace(0,1,50))
plt.colorbar(caf)
plt.show()







