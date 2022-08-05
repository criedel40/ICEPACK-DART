import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import os
import xarray as xr
import sys

test = False
if test:
  mean = np.linspace(1.0e-8,0.999,50)
  variance = np.linspace(1.0e-8,0.26,50)
  a = np.zeros((50,50))
  b = np.zeros((50,50))
  for i in range(0,50):
    for j in range(0,50):
      a[i,j] = mean[i]*((mean[i]*(1.0-mean[i]) - variance[j])/variance[j])
      b[i,j] = (1.0-mean[i])*((mean[i]*(1.0-mean[i]) - variance[j])/variance[j])
  fig,ax = plt.subplots(2,1,sharex=True)
  caf = ax[0].contourf(mean,variance,a.T,np.linspace(0.001,100))
  caf2 = ax[1].contourf(mean,variance,b.T,np.linspace(0.001,100))
  plt.colorbar(caf,ax=ax[0])
  plt.colorbar(caf2,ax=ax[1])
  plt.savefig('test_params.jpg',bbox_inches='tight',dpi=550)
  os.system('open test_params.jpg')
  sys.exit()

##############################################
inflated_val = 2.0

data = xr.open_dataset('gridpt_cice_data.nc')
sic = data.si_data.values[:,1,:]

inflate_sic = np.zeros((sic.shape[0],sic.shape[1]))

for t in range(0,sic.shape[1]):
  mean = sic[:,t].mean()
  variance = sic[:,t].var(ddof=1)
  if variance == 0 or variance*inflated_val > mean*(1.0-mean):
    #print('Skipping')
    inflate_sic[:,t] = sic[:,t]
    continue
#if variance == 0 or variance > mean*(1.0-mean):
#  print('Need to adjust input mean and variance')
  

  a = mean*(((mean*(1.0-mean))/variance) - 1.0)
  b = (1.0-mean)*(((mean*(1.0-mean))/variance) - 1.0)

  #org_dist = scipy.stats.beta.rvs(a,b,size=100)

  ## Inflated Variance and new beta params
  variance2 = variance*inflated_val
  a2 = mean*(((mean*(1.0-mean))/variance2) - 1.0)
  b2 = (1.0-mean)*(((mean*(1.0-mean))/variance2) - 1.0)
  sys.exit()
  new_dist = scipy.stats.beta.ppf(scipy.stats.beta.cdf(sic[:,t],a,b),a2,b2)
  inflate_sic[:,t] = new_dist

tim = np.arange(0,sic.shape[1])
fig,ax = plt.subplots(2,1,sharex=True)
ax[0].plot(tim,sic.mean(axis=0),'-k',label='Orginal')
ax[0].plot(tim,inflate_sic.mean(axis=0),'-r',label='Inflated')
##
ax[1].plot(tim,sic.var(ddof=1,axis=0),'-k',label='Orginal')
ax[1].plot(tim,inflate_sic.var(ddof=1,axis=0),'-r',label='Inflated')
##
ax[1].set_xlim(tim[0],tim[-1])
plt.savefig('timseries_mean_var.jpg',dpi=550,bbox_inches='tight')
os.system('open timseries_mean_var.jpg')
####
fig,ax = plt.subplots(2,1,sharex=True)
ax[0].plot(tim,inflate_sic.mean(axis=0)-sic.mean(axis=0),'-k',label='Orginal')
#ax[0].plot(tim,inflate_sic.mean(axis=0),'-r',label='Inflated')
##
ax[1].plot(tim,inflate_sic.var(ddof=1,axis=0)-sic.var(ddof=1,axis=0),'-k',label='Orginal')
#ax[1].plot(tim,inflate_sic.var(ddof=1,axis=0),'-r',label='Inflated')
##
ax[1].set_xlim(tim[0],tim[-1])
ax[1].set_xticks(tim[::365])
plt.savefig('timseries_bias_mean_var.jpg',dpi=550,bbox_inches='tight')
os.system('open timseries_bias_mean_var.jpg')





sys.exit()

fig,ax = plt.subplots(2,1,sharex=True)
ax[0].hist(sic[:,566],bins=75)
textstr = '\n'.join((
    r'$\mu=%.3f$' % (sic[:,566].mean(), ),
    r'$\sigma^2=%.3f$' % (sic[:,566].var(ddof=1), )))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax[0].text(0.05, 0.95, textstr, transform=ax[0].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
ax[1].hist(inflate_sic[:,566],bins=75)
textstr = '\n'.join((
    r'$\mu=%.3f$' % (inflate_sic[:,566].mean(), ),
    r'$\sigma^2=%.3f$' % (inflate_sic[:,566].var(ddof=1), )))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax[1].text(0.05, 0.95, textstr, transform=ax[1].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

ax[1].set_xlim(0,1)
plt.savefig('test2.jpg',dpi=550,bbox_inches='tight')
os.system('open test2.jpg')




plt.rc('font', weight='bold')

fig,ax = plt.subplots(2,1,sharex=True)
ax[0].hist(org_dist,bins=50)
textstr = '\n'.join((
    r'$\mu=%.3f$' % (org_dist.mean(), ),
    r'$\sigma^2=%.3f$' % (org_dist.var(ddof=1), )))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax[0].text(0.05, 0.95, textstr, transform=ax[0].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
ax[1].hist(new_dist,bins=50)
textstr = '\n'.join((
    r'$\mu=%.3f$' % (new_dist.mean(), ),
    r'$\sigma^2=%.3f$' % (new_dist.var(ddof=1), )))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax[1].text(0.05, 0.95, textstr, transform=ax[1].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
ax[1].set_xlim(0,1)

plt.savefig('test.jpg',dpi=550,bbox_inches='tight')
os.system('open test.jpg')
