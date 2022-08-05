import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import os
import xarray as xr
import sys

data = xr.open_dataset('gridpt_cice_data.nc')
sic = data.si_data.values[:,1,566]

#sic = (sic - np.min(sic))/(np.max(sic) - np.min(sic))

inflation_vals = np.arange(1,2+0.1,0.1)

mean = sic.mean()
variance = sic.var(ddof=1)

a = mean*(((mean*(1.0-mean))/variance) - 1.0)
b = (1.0-mean)*(((mean*(1.0-mean))/variance) - 1.0)
#a1,b1,loc1,scale1 = scipy.stats.beta.fit(sic)
rand_HOLD = scipy.stats.beta.rvs(a,b,size=10000)
sys.exit()
#############################################
#fig,ax = plt.subplots(2,1,sharex=True)
#ax[0].hist(rand_dists,bins=50)
#textstr = '\n'.join((
#    r'$\mu=%.3f$' % (rand_dists.mean(), ),
#    r'$\sigma^2=%.3f$' % (rand_dists.var(ddof=1), )))
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#ax[0].text(0.05, 0.95, textstr, transform=ax[0].transAxes, fontsize=10,
#        verticalalignment='top', bbox=props)
#ax[0].set_xlim(0,1)
#ax[1].hist(inflated_dist[0,:],bins=50)
#textstr = '\n'.join((
#    r'$\mu=%.3f$' % (inflated_dist[0,:].mean(), ),
#    r'$\sigma^2=%.3f$' % (inflated_dist[0,:].var(ddof=1), )))
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#ax[1].text(0.05, 0.95, textstr, transform=ax[1].transAxes, fontsize=10,
#        verticalalignment='top', bbox=props)
#ax[1].set_xlim(0,1)
#ax[0].set_title('Randomly Generated Distribution',weight='bold')
#ax[1].set_title('Data from Sea Ice Model',weight='bold')
#plt.savefig('no_inflation_rand_real.jpg',dpi=550,bbox_inches='tight')
#os.system('open no_inflation_rand_real.jpg')
#sys.exit()
#############################################
#sys.exit()
rand_dists = np.zeros((inflation_vals.shape[0],10000))
for i in range(0,inflation_vals.shape[0]):
  if variance*inflation_vals[i] > mean*(1.0-mean):
    print('Variance is too large...STOP!')
    sys.exit()
  variance2 = variance*inflation_vals[i]
  a2 = mean*(((mean*(1.0-mean))/variance2) - 1.0)
  b2 = (1.0-mean)*(((mean*(1.0-mean))/variance2) - 1.0)
  rand_dists[i,:] = scipy.stats.beta.rvs(a2,b2,size=10000)
  #rand_dists[i,:] = rand_dists[i,:]
plt.rc('font', weight='bold')

fig,ax = plt.subplots(5,2,sharex=True,sharey=False)
count = 0
for aa in ax.flatten():
  hist = aa.hist(rand_dists[count,:],bins=50)
  aa.set_xlim(0,1)
  aa.plot(np.array([mean,mean]),np.array([0,hist[0].max()]),'--k',label='Org. Mean')
  aa.plot(np.array([rand_dists[count,:].mean(),rand_dists[count,:].mean()]),np.array([0,hist[0].max()]),'-k',label='Inflated Mean')
  aa.legend(loc='best',fontsize=4)
  aa.set_title('Inflation Val:{0:.1f} - Inflated Var:{1:.3f} - Desired Var:{1:.3f}'.format(inflation_vals[count],rand_dists[count,:].var(ddof=1),variance*inflation_vals[count]),fontsize=5,weight='bold',pad=0.8)
  count+=1
  aa.set_xticks(np.arange(0,1+0.1,0.1))
  aa.tick_params(axis='both', which='major', labelsize=5)
  #aa.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
plt.savefig('total_hist_rand.jpg',dpi=550,bbox_inches='tight')
os.system('open total_hist_rand.jpg')
#sys.exit()
#############################################################
inflated_dist = np.zeros((inflation_vals.shape[0],sic.shape[0]))
for i in range(0,inflation_vals.shape[0]):
  if variance*inflation_vals[i] > mean*(1.0-mean):
    print('Variance is too large...STOP!')
    sys.exit()
  variance2 = variance*inflation_vals[i]
  org_theo_max = mean*(1.0-mean)
  upper_bound = 0.5*(np.sqrt(1.0 - 4.0*variance) + 1.0)
  diff_mean = upper_bound - mean
  inflate_upper_bound = 0.5*(np.sqrt(1.0 - 4.0*variance2) + 1.0)
  new_mean = inflate_upper_bound - diff_mean
  
  a2 = new_mean*(((new_mean*(1.0-new_mean))/variance2) - 1.0)
  b2 = (1.0-new_mean)*(((new_mean*(1.0-new_mean))/variance2) - 1.0)
  inflated_dist[i,:] = scipy.stats.beta.ppf(scipy.stats.beta.cdf(sic[:],a,b),a2,b2)
  #inflated_dist[i,:] = inflated_dist[i,:] - (mean-0.5)
  inflated_mean = inflated_dist[i,:].mean()
  #diff = inflated_mean - mean
  #inflated_dist[i,:] = inflated_dist[i,:] - diff
fig,ax = plt.subplots(5,2,sharex=True,sharey=False)
count = 1
for aa in ax.flatten():
  hist = aa.hist(inflated_dist[count,:],bins=50)
  aa.plot(np.array([mean,mean]),np.array([0,hist[0].max()]),'--k',label='Org. Mean')
  aa.plot(np.array([inflated_dist[count,:].mean(),inflated_dist[count,:].mean()]),np.array([0,hist[0].max()]),'-k',label='Inflated Mean')
  aa.legend(loc='best',fontsize=4)
  aa.set_title('Inflation Val:{0:.1f} - Inflated Var:{1:.3f} - Desired Var:{1:.3f}'.format(inflation_vals[count],inflated_dist[count,:].var(ddof=1),variance*inflation_vals[count]),fontsize=5,weight='bold',pad=0.8)
  aa.set_xticks(np.arange(0,1+0.1,0.1))
  aa.set_xlim(0,1)
  aa.tick_params(axis='both', which='major', labelsize=5)
  #aa.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
  count+=1
plt.savefig('total_hist_real.jpg',dpi=550,bbox_inches='tight')
os.system('open total_hist_real.jpg')
#sys.exit()
#############################################################
fig,ax = plt.subplots(2,1,sharex=True)
ax[0].hist(rand_HOLD,bins=50)
textstr = '\n'.join((
    r'$\mu=%.3f$' % (rand_HOLD.mean(), ),
    r'$\sigma^2=%.3f$' % (rand_HOLD.var(ddof=1), )))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax[0].text(0.05, 0.95, textstr, transform=ax[0].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
ax[0].set_xlim(0,1)
ax[1].hist(inflated_dist[0,:],bins=50)
textstr = '\n'.join((
    r'$\mu=%.3f$' % (inflated_dist[0,:].mean(), ),
    r'$\sigma^2=%.3f$' % (inflated_dist[0,:].var(ddof=1), )))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax[1].text(0.05, 0.95, textstr, transform=ax[1].transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
ax[1].set_xlim(0,1)
ax[0].set_title('Randomly Generated Distribution',weight='bold')
ax[1].set_title('Data from Sea Ice Model',weight='bold')
plt.savefig('no_inflation_rand_real.jpg',dpi=550,bbox_inches='tight')
os.system('open no_inflation_rand_real.jpg')
#sys.exit()
#############################################################
plt.rc('font', weight='bold')
tim = np.arange(0,inflation_vals.shape[0])
fig,ax = plt.subplots(2,1,sharex=True)
ax[0].bar(tim,rand_dists.mean(axis=1))
hold = np.zeros((inflation_vals.shape[0]))
hold[:] = mean
ax[0].plot(tim,hold,'ok')
ax[1].bar(tim,rand_dists.var(ddof=1,axis=1))
ax[1].plot(tim,variance*inflation_vals,'ok')
ax[0].set_ylabel('Mean')
ax[1].set_ylabel('Variance')
ax[1].set_xlabel('Inflation Values')
ax[1].set_xticks(tim)
ax[1].set_xticklabels(np.array([1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]),rotation=90)
plt.savefig('test_inflation_vals_rand.jpg',dpi=550,bbox_inches='tight')
os.system('open test_inflation_vals_rand.jpg')
###
fig,ax = plt.subplots(2,1,sharex=True)
ax[0].bar(tim,rand_dists.mean(axis=1)-mean)
ax[1].bar(tim,rand_dists.var(ddof=1,axis=1)-(variance*inflation_vals))
#ax[1].plot(tim,variance*inflation_vals,'ok')
ax[0].set_ylabel('Mean')
ax[1].set_ylabel('Variance')
ax[1].set_xlabel('Inflation Values')
ax[1].set_xticks(tim)
ax[1].set_xticklabels(np.array([1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]),rotation=90)
plt.savefig('test_inflation_vals_rand_bias.jpg',dpi=550,bbox_inches='tight')
os.system('open test_inflation_vals_rand_bias.jpg')
##
###################
fig,ax = plt.subplots(2,1,sharex=True)
ax[0].bar(tim,inflated_dist.mean(axis=1))
hold = np.zeros((inflation_vals.shape[0]))
hold[:] = mean
ax[0].plot(tim,hold,'ok')
ax[1].bar(tim,inflated_dist.var(ddof=1,axis=1))
ax[1].plot(tim,variance*inflation_vals,'ok')
ax[0].set_ylabel('Mean')
ax[1].set_ylabel('Variance')
ax[1].set_xlabel('Inflation Values')
ax[1].set_xticks(tim)
ax[1].set_xticklabels(np.array([1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]),rotation=90)
plt.savefig('test_inflation_vals_actual.jpg',dpi=550,bbox_inches='tight')
os.system('open test_inflation_vals_actual.jpg')
###
fig,ax = plt.subplots(2,1,sharex=True)
ax[0].bar(tim,inflated_dist.mean(axis=1)-mean)
ax[1].bar(tim,inflated_dist.var(ddof=1,axis=1)-(variance*inflation_vals))
#ax[1].plot(tim,variance*inflation_vals,'ok')
ax[0].set_ylabel('Mean')
ax[1].set_ylabel('Variance')
ax[1].set_xlabel('Inflation Values')
ax[1].set_xticks(tim)
ax[1].set_xticklabels(np.array([1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]),rotation=90)
plt.savefig('test_inflation_vals_actual_bias.jpg',dpi=550,bbox_inches='tight')
os.system('open test_inflation_vals_actual_bias.jpg')



