import truncnorm_dist_funcs
import scipy.stats
import numpy as np
import os
import matplotlib.pyplot as plt

dist = scipy.stats.norm.rvs(size=10000)
dist = np.sort(dist)

scipy_cdf = np.zeros((dist.shape[0]))
scipy_cdfinv = np.zeros((dist.shape[0]))

mine_cdf = np.zeros((dist.shape[0]))
mine_cdfinv = np.zeros((dist.shape[0]))

for z in range(0,dist.shape[0]):
  scipy_cdf[z] = scipy.stats.norm.cdf(dist[z])
  scipy_cdfinv[z] = scipy.stats.norm.ppf(scipy_cdf[z])
  ####
  mine_cdf[z] = truncnorm_dist_funcs.norm_cdf(dist[z])
  mine_cdfinv[z] = truncnorm_dist_funcs.norm_cdf_inv(mine_cdf[z])

fig,ax = plt.subplots(2,1,figsize=(6,4),sharex=True)
ax[0].plot(dist,scipy_cdf,'-k',label='Scipy')
ax[0].plot(dist,mine_cdf,'-r',label='Mine')
##
ax[1].plot(dist,mine_cdf-scipy_cdf,'k')
plt.savefig('cdf_test.jpg',dpi=550,bbox_inches='tight')
os.system('open cdf_test.jpg')
#########
fig,ax = plt.subplots(2,1,figsize=(6,4),sharex=True)
ax[0].plot(dist,dist,'-k',label='Input')
ax[0].plot(dist,scipy_cdfinv,'-r',label='Scipy')
ax[0].plot(dist,mine_cdfinv,'-b',label='Mine')
ax[0].legend()
##
ax[1].plot(dist,scipy_cdfinv-dist,'k',label='Scipy-Truth')
ax[1].plot(dist,mine_cdfinv-dist,'r',label='Mine-Truth')
ax[1].plot(dist,mine_cdfinv-scipy_cdfinv,'b',label='Mine-Scipy')
ax[1].legend()
plt.savefig('cdfinv_test.jpg',dpi=550,bbox_inches='tight')
os.system('open cdfinv_test.jpg')
###########################################################
###########################################################
dist = scipy.stats.truncnorm.rvs(loc=0.0,scale=1.0,a=-1,b=1,size=10000)
dist = np.sort(dist)

scipy_cdf = np.zeros((dist.shape[0]))
scipy_cdfinv = np.zeros((dist.shape[0]))

mine_cdf = np.zeros((dist.shape[0]))
mine_cdfinv = np.zeros((dist.shape[0]))

for z in range(0,dist.shape[0]):
  scipy_cdf[z] = scipy.stats.truncnorm.cdf(dist[z],loc=0.0,scale=1.0,a=-1,b=1)
  scipy_cdfinv[z] = scipy.stats.truncnorm.ppf(scipy_cdf[z],loc=0.0,scale=1.0,a=-1,b=1)
  ####
  mine_cdf[z] = truncnorm_dist_funcs.truncnorm_cdf(dist[z],0.0,1.0,a=-1,b=1)
  mine_cdfinv[z] = truncnorm_dist_funcs.truncnorm_cdf_inv(mine_cdf[z],0.0,1.0,a=-1,b=1)

fig,ax = plt.subplots(2,1,figsize=(6,4),sharex=True)
ax[0].plot(dist,scipy_cdf,'-k',label='Scipy')
ax[0].plot(dist,mine_cdf,'-r',label='Mine')
##
ax[1].plot(dist,mine_cdf-scipy_cdf,'k')
plt.savefig('trunccdf_test.jpg',dpi=550,bbox_inches='tight')
os.system('open trunccdf_test.jpg')
#########
fig,ax = plt.subplots(2,1,figsize=(6,4),sharex=True)
ax[0].plot(dist,dist,'-k',label='Input')
ax[0].plot(dist,scipy_cdfinv,'-r',label='Scipy')
ax[0].plot(dist,mine_cdfinv,'-b',label='Mine')
ax[0].legend()
##
ax[1].plot(dist,scipy_cdfinv-dist,'k',label='Scipy-Truth')
ax[1].plot(dist,mine_cdfinv-dist,'r',label='Mine-Truth')
ax[1].plot(dist,mine_cdfinv-scipy_cdfinv,'b',label='Mine-Scipy')
ax[1].legend()
plt.savefig('trunccdfinv_test.jpg',dpi=550,bbox_inches='tight')
os.system('open trunccdfinv_test.jpg')


