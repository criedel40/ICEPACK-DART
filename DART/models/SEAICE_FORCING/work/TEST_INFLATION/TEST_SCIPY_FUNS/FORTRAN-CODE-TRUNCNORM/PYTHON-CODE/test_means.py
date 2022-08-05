import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import os

#top = scipy.stats.norm.pdf(0.0) - scipy.stats.norm.pdf(1.0)
#bot = scipy.stats.norm.cdf(1.0) - scipy.stats.norm.cdf(0.0)

STDDEV = np.linspace(0.0,0.5,50)
MEANS = np.linspace(0,1,100)

adjusted_mean = np.zeros((MEANS.shape[0],STDDEV.shape[0]))
adjusted_mean_diff = np.zeros((MEANS.shape[0],STDDEV.shape[0]))
for s in range(0,STDDEV.shape[0]):
  for m in range(0,MEANS.shape[0]):
    if STDDEV[s] == 0.0:
      adjusted_mean[m,s] = MEANS[m] + (1.0)*STDDEV[s]
    else:
      alpha = (0.0-MEANS[m])/STDDEV[s]
      beta = (1.0-MEANS[m])/STDDEV[s]
      top = scipy.stats.norm.pdf(alpha) - scipy.stats.norm.pdf(beta)
      bot = scipy.stats.norm.cdf(beta) - scipy.stats.norm.cdf(alpha)
      adjusted_mean[m,s] = MEANS[m] + (top/bot)*STDDEV[s]
  adjusted_mean_diff[:,s] = adjusted_mean[:,s] - MEANS

plt.rc('font', weight='bold')

fig,ax = plt.subplots(2,1,figsize=(8,8),sharex=True)
caf1 = ax[0].contourf(MEANS,STDDEV,adjusted_mean.T,np.linspace(0,1,50),cmap=plt.get_cmap('seismic'),extend='both')
plt.colorbar(caf1,ax=ax[0],ticks=np.arange(0,1+0.1,0.1),pad=0.008)
##
caf1 = ax[1].contourf(MEANS,STDDEV,adjusted_mean_diff.T,np.linspace(-0.5,0.5,50),cmap=plt.get_cmap('seismic'))
plt.colorbar(caf1,ax=ax[1],pad=0.008,ticks=np.arange(-0.5,0.5+0.1,0.1))
ax[1].set_xlim(0,1)
ax[1].set_xlabel('Input Mean',weight='bold')
ax[0].set_ylabel('STDDEV',weight='bold')
ax[1].set_ylabel('STDDEV',weight='bold')
ax[0].set_title('Inflation Mean',weight='bold')
ax[1].set_title('Inflation Mean - Uninflated Mean',weight='bold')
plt.savefig('test_doublebound.jpg',dpi=550,bbox_inches='tight')
os.system('open test_doublebound.jpg')
##################################
##################################
STDDEV = np.linspace(0.0,1.0,50)
MEANS = np.linspace(0,10,100)

adjusted_mean = np.zeros((MEANS.shape[0],STDDEV.shape[0]))
adjusted_mean_diff = np.zeros((MEANS.shape[0],STDDEV.shape[0]))
for s in range(0,STDDEV.shape[0]):
  for m in range(0,MEANS.shape[0]):
    if STDDEV[s] == 0.0:
      adjusted_mean[m,s] = MEANS[m] + (1.0)*STDDEV[s]
    else:
      alpha = (0.0-MEANS[m])/STDDEV[s]
      #beta = (1.0-MEANS[m])/STDDEV[s]
      top = scipy.stats.norm.pdf(alpha) #- scipy.stats.norm.pdf(beta)
      bot = 1.0 - scipy.stats.norm.cdf(alpha)
      adjusted_mean[m,s] = MEANS[m] + (top/bot)*STDDEV[s]
  adjusted_mean_diff[:,s] = adjusted_mean[:,s] - MEANS


fig,ax = plt.subplots(2,1,figsize=(8,8),sharex=True)
caf1 = ax[0].contourf(MEANS,STDDEV,adjusted_mean.T,np.linspace(0,10,50),cmap=plt.get_cmap('seismic'),extend='both')
plt.colorbar(caf1,ax=ax[0],ticks=np.arange(0,10+1,1),pad=0.008)
##
caf1 = ax[1].contourf(MEANS,STDDEV,adjusted_mean_diff.T,np.linspace(-1,1,50),cmap=plt.get_cmap('seismic'),extend='both')
plt.colorbar(caf1,ax=ax[1],pad=0.008,ticks=np.arange(-1,1+0.5,0.5))
ax[1].set_xlim(0,10)
ax[1].set_xlabel('Input Mean',weight='bold')
ax[0].set_ylabel('STDDEV',weight='bold')
ax[1].set_ylabel('STDDEV',weight='bold')
ax[0].set_title('Inflation Mean',weight='bold')
ax[1].set_title('Inflation Mean - Uninflated Mean',weight='bold')
plt.savefig('test_singlebound.jpg',dpi=550,bbox_inches='tight')
os.system('open test_singlebound.jpg')
###########################################
###########################################
inflation_factor = np.arange(1,3+0.25,0.25)
MEANS = np.linspace(0.0000001,0.9999999,100)
STDDEV = np.sqrt(MEANS*(1.0-MEANS))

adjusted_mean = np.zeros((MEANS.shape[0],inflation_factor.shape[0]))
adjusted_mean_diff = np.zeros((MEANS.shape[0],inflation_factor.shape[0]))
for s in range(0,inflation_factor.shape[0]):
  for m in range(0,MEANS.shape[0]):
    if STDDEV[s] == 0.0:
      adjusted_mean[m,s] = MEANS[m] + (1.0)*STDDEV[s]
    else:
      hold_std = STDDEV[m]*inflation_factor[s]
      alpha = (0.0-MEANS[m])/hold_std
      beta = (1.0-MEANS[m])/hold_std
      top = scipy.stats.norm.pdf(alpha) - scipy.stats.norm.pdf(beta)
      bot = scipy.stats.norm.cdf(beta) - scipy.stats.norm.cdf(alpha)
      adjusted_mean[m,s] = MEANS[m] + (top/bot)*hold_std
  adjusted_mean_diff[:,s] = adjusted_mean[:,s] - MEANS

fig = plt.figure()
plt.plot(MEANS,STDDEV,'-k',label='STDDEV')
plt.xlim(0,1)
plt.xlabel('Mean',weight='bold')
plt.ylabel('STDDEV',weight='bold')
plt.title('STDDEV vs. Input Mean',weight='bold')
plt.savefig('stddev_values.jpg',dpi=550,bbox_inches='tight')
os.system('open stddev_values.jpg')

fig,ax = plt.subplots(2,1,figsize=(8,8),sharex=True)
caf1 = ax[0].contourf(MEANS,inflation_factor,adjusted_mean.T,np.linspace(0,1,50),cmap=plt.get_cmap('seismic'),extend='both')
plt.colorbar(caf1,ax=ax[0],ticks=np.arange(0,1+0.1,0.1),pad=0.008)
##
caf1 = ax[1].contourf(MEANS,inflation_factor,adjusted_mean_diff.T,np.linspace(-0.5,0.5,50),cmap=plt.get_cmap('seismic'),extend='both')
plt.colorbar(caf1,ax=ax[1],pad=0.008,ticks=np.arange(-0.5,0.5+0.1,0.1))
ax[1].set_xlim(0,1)
ax[1].set_xlabel('Input Mean',weight='bold')
ax[0].set_ylabel('Inf. Factor',weight='bold')
ax[1].set_ylabel('Inf. Factor',weight='bold')
ax[0].set_title('Inflation Mean',weight='bold')
ax[1].set_title('Inflation Mean - Uninflated Mean',weight='bold')
plt.savefig('test_inflationfact.jpg',dpi=550,bbox_inches='tight')
os.system('open test_inflationfact.jpg')
