import truncnorm_dist_funcs
import fortran_funcs
import scipy.stats
import numpy as np
import os
import matplotlib.pyplot as plt
import xarray as xr
import os

data = xr.open_dataset('assim_output_EaKF_gaus.nc')
input_aice = data.input_aice.values[1:,:]

inflated_values_trunc = np.zeros((input_aice.shape[0],input_aice.shape[1]))
inflated_values_beta = np.zeros((input_aice.shape[0],input_aice.shape[1]))

inflate_val = 1.5
for t in range(0,input_aice.shape[0]):
  print(t)
  prior = input_aice[t,:]
  prior[prior==1] = 0.99998888
  #### Start with trunc
  mean = prior.mean()
  std = prior.std(ddof=1)
  inflated_std = std*inflate_val
  cdfs = np.array([truncnorm_dist_funcs.truncnorm_cdf(d,mean,std,0.0,1.0) for d in prior])
  inv_cdfs = np.array([truncnorm_dist_funcs.truncnorm_cdf_inv(d,mean,inflated_std,0.0,1.0) for d in cdfs])
  inflated_values_trunc[t,:] = inv_cdfs
  ##### Start Beta
  mean = prior.mean()
  std = prior.std(ddof=1)
  var = std**2
  inflated_std = std*inflate_val
  inflated_var = inflated_std**2
  if inflated_var > mean*(1.0-mean):
    theo_max_var = mean*(1.0 - mean)
    if (theo_max_var > var):
      print('Use Theo max')
      inflate_prior_var = theo_max_var - 1.0e-5
    else:
      print('Can not inflate prior ensemble....continue')
      inflate_prior_var = var
  else:
    inflate_prior_var = inflated_var
  a = mean*(((mean*(1.0-mean))/var) - 1.0)
  b = (1.0-mean)*(((mean*(1.0-mean))/var) - 1.0)
  inflate_a = mean*(((mean*(1.0-mean))/inflate_prior_var) - 1.0)
  inflate_b = (1.0-mean)*(((mean*(1.0-mean))/inflate_prior_var) - 1.0)
  cdfs = np.array([fortran_funcs.betacdf(d,a,b) for d in prior])
  inv_cdfs = np.array([fortran_funcs.betaincinv(d,inflate_a,inflate_b) for d in cdfs])
  inflated_values_beta[t,:] = inv_cdfs


mean_bias_trunc = inflated_values_trunc.mean(axis=1) - input_aice.mean(axis=1)
mean_bias_beta = inflated_values_beta.mean(axis=1) - input_aice.mean(axis=1)
###
std_bias_trunc = inflated_values_trunc.std(ddof=1,axis=1) - input_aice.std(ddof=1,axis=1)
std_bias_beta = inflated_values_beta.std(ddof=1,axis=1) - input_aice.std(ddof=1,axis=1)
mean_control = input_aice.mean(axis=1)
mean_trunc = inflated_values_trunc.mean(axis=1)
mean_beta =inflated_values_beta.mean(axis=1)

plt.rc('font', weight='bold')
tim = np.arange(0,input_aice.shape[0])
fig,ax = plt.subplots(4,1,sharex=True,figsize=(10,8))
ax[0].plot(tim,mean_control,'-k',label='Uninflated')
ax[0].plot(tim,mean_trunc,'-r',label='Inflated-Truncnorm')
ax[0].plot(tim,mean_beta,'-b',label='Inflated-Beta')
ax[0].set_title('Ensemble Mean',weight='bold')
ax[0].legend(ncol=3,loc='lower right')
#####
ax[1].plot(tim,input_aice.std(ddof=1,axis=1),'-c',label='Uninflated Prior')
ax[1].plot(tim,input_aice.std(ddof=1,axis=1)*inflate_val ,'-k',label='Desired Inflated Spread')
ax[1].plot(tim,inflated_values_trunc.std(ddof=1,axis=1),'-r',label='Inflated-Truncnorm')
ax[1].plot(tim,inflated_values_beta.std(ddof=1,axis=1),'-b',label='Inflated-Beta')
ax[1].set_title('Ensemble Inflated STDDEV',weight='bold')
ax[1].legend(ncol=3)
#####
ax[2].plot(tim,mean_bias_trunc,'-r',label='Inflated-Truncnorm')
ax[2].plot(tim,mean_bias_beta,'-b',label='Inflated-Beta')
ax[2].set_title('Ensemble Mean Bias',weight='bold')
ax[2].legend(ncol=3)
#####
ax[3].plot(tim,input_aice.std(ddof=1,axis=1)*inflate_val-input_aice.std(ddof=1,axis=1),'-k',label='Desired Inflated Spread')
ax[3].plot(tim,inflated_values_trunc.std(ddof=1,axis=1)-input_aice.std(ddof=1,axis=1),'-r',label='Inflated-Truncnorm')
ax[3].plot(tim,inflated_values_beta.std(ddof=1,axis=1)-input_aice.std(ddof=1,axis=1),'-b',label='Inflated-Beta')
ax[3].set_title('Ensemble Inflated STDDEV Difference to uninflated',weight='bold')
ax[3].legend(ncol=3)
#####
ax[3].set_xlim(tim[0],tim[-1])
tit = plt.suptitle('Sea Ice Concentration - Inflation Factor: {0}'.format(inflate_val),weight='bold',y=0.95)
plt.savefig('inflated_ts_{0}.jpg'.format(inflate_val),dpi=550,bbox_inches='tight')
os.system('open inflated_ts_{0}.jpg'.format(inflate_val))
###############################################################
###############################################################
plt.rc('font', weight='bold')
tim = np.arange(0,input_aice.shape[0])
fig,ax = plt.subplots(3,1,sharex=True,figsize=(10,8))
for m in range(0,input_aice.shape[1]):
  if m == 0:
    ax[0].plot(tim,input_aice[:,m],'-',color='grey',linewidth=1,label='Individual Members',alpha=0.5)
    ax[1].plot(tim,inflated_values_trunc[:,m],'-',color='grey',linewidth=1,label='Individual Members',alpha=0.5)
    ax[2].plot(tim,inflated_values_beta[:,m],'-',color='grey',linewidth=1,label='Individual Members',alpha=0.5)
  else:
    ax[0].plot(tim,input_aice[:,m],'-',color='grey',linewidth=1,alpha=0.5)
    ax[1].plot(tim,inflated_values_trunc[:,m],'-',color='grey',linewidth=1,alpha=0.5)
    ax[2].plot(tim,inflated_values_beta[:,m],'-',color='grey',linewidth=1,alpha=0.5)
ax[0].plot(tim,np.nanmean(input_aice[:,:],axis=1),'-',color='k',linewidth=2,label='Mean')
ax[1].plot(tim,np.nanmean(inflated_values_trunc[:,:],axis=1),'-',color='k',linewidth=2,label='Mean')
ax[2].plot(tim,np.nanmean(inflated_values_beta[:,:],axis=1),'-',color='k',linewidth=2,label='Mean')
####
ax[2].set_xlim(tim[0],tim[-1])
plt.savefig('inflated_ts_mems_{0}.jpg'.format(inflate_val),dpi=550,bbox_inches='tight')
os.system('open inflated_ts_mems_{0}.jpg'.format(inflate_val))
