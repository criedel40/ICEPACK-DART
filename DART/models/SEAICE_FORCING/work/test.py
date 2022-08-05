import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os


data = xr.open_dataset('assim_output_EaKF_gaus.nc')
output_hi = data.output_hi.values[2,:]
input_hi = data.input_hi.values[2,:]
obs_hi = data.obs_hi.values[2]


fig,ax = plt.subplots(2,1,sharex=True)
ax[0].bar(np.arange(1,100),input_hi)
ax[0].plot(np.array([0,100]),np.array([obs_hi,obs_hi]),'--k',label='Observation')
ax[1].bar(np.arange(1,100),output_hi)
ax[1].plot(np.array([0,100]),np.array([obs_hi,obs_hi]),'--k',label='Observation')
ax[0].legend(loc='lower left')
ax[1].legend(loc='lower left')
ax[1].set_xlim(0,100)
ax[1].set_xlabel('Ensemble Members')
ax[0].set_title('SI Thickness Prior Ensemble',weight='bold')
ax[1].set_title('SI Thickness Posterior Ensemble',weight='bold')
ax[0].set_ylabel('Thickness (m)',weight='bold')
ax[1].set_ylabel('Thickness (m)',weight='bold')
plt.savefig('test.jpg',dpi=550,bbox_inches='tight')
os.system('open test.jpg')


new_file = open('prior_si_thick_case.txt','w')
new_file2 = open('posterior_si_thick_case.txt','w')

for l in range(0,input_hi.shape[0]):
  new_file.writelines('{0}\n'.format(input_hi[l]))
  new_file2.writelines('{0}\n'.format(output_hi[l]))
new_file.writelines('--Observation--\n')
new_file2.writelines('--Observation--\n')
new_file.writelines('{0}'.format(obs_hi))
new_file2.writelines('{0}'.format(obs_hi))
new_file.close()
new_file2.close()


