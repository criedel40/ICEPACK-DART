import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import CubicSpline


data = open('sic_err_list.txt').readlines()
sic_err_array = np.zeros((len(data),2))
for l in range(0,len(data)):
  line = data[l].strip().split(',')
  sic_err_array[l,:] = np.array(list(map(float,line))).astype('float64')

cs = CubicSpline(sic_err_array[:,0],sic_err_array[:,1])


sic = np.linspace(0.0,1.0,100)
error_vals = np.zeros((4,sic.shape[0]))
for i in range(0,sic.shape[0]):
  if sic[i] < 0.01:
    error_vals[0,i] = 0.01*0.15
  else:
    error_vals[0,i] = sic[i]*0.15
##
for i in range(0,sic.shape[0]):
  error_vals[1,i] = cs(sic[i])
##
error_vals[2,:] = 0.5*0.15
##
ratio = 0.24997550249974496/0.0225
for i in range(0,sic.shape[0]):
  if sic[i] < 0.01:
    error_vals[3,i] = 0.01*0.15
  elif sic[i] > 0.99:
    error_vals[3,i] = (0.15 - (0.99*0.15))
  else:
    error_vals[3,i] = np.sqrt((sic[i]*(1.0-sic[i]))/ratio)

fig = plt.figure()
plt.plot(sic,error_vals[0,:],'-k',label='Control')
plt.plot(sic,error_vals[1,:],'-r',label='Control Inverse')
plt.plot(sic,error_vals[2,:],'-b',label='constant value')
plt.plot(sic,error_vals[3,:],'-c',label='Beta Bell Shape')
plt.xlim(0,1)
plt.xlabel('Sea Ice Concentration',weight='bold')
plt.ylabel('Obs. Error (STDDEV)',weight='bold')
plt.ylim(0,0.26)
plt.legend(loc='best')
plt.grid(True,axis='y',linestyle='dashed',alpha=0.5)
plt.savefig('sic_err_dist.jpg',dpi=550,bbox_inches='tight')
os.system('open sic_err_dist.jpg')

