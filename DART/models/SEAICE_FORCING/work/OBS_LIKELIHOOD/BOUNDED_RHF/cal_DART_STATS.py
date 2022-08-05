import xarray as xr
import numpy as np
import glob
import netCDF4
import os
import sys

print(remove_negative_sic)

#files = sorted(glob.glob('2013*/obs_epoch_001.nc'))
# LOAD IN DATA
data = xr.open_dataset('assim_output_BoundedRHF_gaus.nc')



marginal_stats = np.zeros((2,3,len(files)))
marginal_qcs = np.zeros((len(files),9))
marginal_rankvalues = np.zeros((len(files),2,80))

####
interm_stats = np.zeros((2,3,len(files)))
interm_qcs = np.zeros((len(files),9))

interm_rankvalues = np.zeros((len(files),2,80))
####
com_stats = np.zeros((2,3,len(files)))
com_qcs = np.zeros((len(files),9))

com_rankvalues = np.zeros((len(files),2,80))
##########
all_stats = np.zeros((2,3,len(files)))
all_qcs = np.zeros((len(files),9))

all_stats_thic = np.zeros((2,3,len(files)))
all_qcs_thic = np.zeros((len(files),9))

all_stats_snow = np.zeros((2,3,len(files)))
all_qcs_snow = np.zeros((len(files),9))

all_rankvalues = np.zeros((len(files),2,80))
all_rankvalues_thic = np.zeros((len(files),2,80))
all_rankvalues_snow = np.zeros((len(files),2,80))
###########
times = []
negative_count = 0
for d in range(0,len(files)):
  print(files[d])
  data = xr.open_dataset(files[d])
  times.append(int(data.time.values[0].strftime('%Y%m%d%H')))
  types = data.obs_type.values
  ind_conc = np.where(types==12)[0]
  ind_thic = np.where(types==15)[0]
  ind_snow = np.where(types==16)[0]
  conc_obs = data.observations.values[ind_conc,0]
  thic_obs = data.observations.values[ind_thic,0]
  snow_obs = data.observations.values[ind_snow,0]
  #
  conc_truth = data.observations.values[ind_conc,1]
  thic_truth = data.observations.values[ind_thic,1]
  snow_truth = data.observations.values[ind_snow,1]
  #
  conc_prmean = data.observations.values[ind_conc,2]
  thic_prmean = data.observations.values[ind_thic,2]
  snow_prmean = data.observations.values[ind_snow,2]
  #
  conc_prspr = data.observations.values[ind_conc,4]
  thic_prspr = data.observations.values[ind_thic,4]
  snow_prspr = data.observations.values[ind_snow,4]
  #
  conc_oberr = data.observations.values[ind_conc,-1]
  thic_oberr = data.observations.values[ind_thic,-1]
  snow_oberr = data.observations.values[ind_snow,-1]
  #
  locs_conc = data.location.values[ind_conc,:2]
  locs_thic = data.location.values[ind_thic,:2]
  locs_snow = data.location.values[ind_snow,:2]
  #
  conc_qc = data.qc.values[ind_conc,1]
  thic_qc = data.qc.values[ind_thic,1]
  snow_qc = data.qc.values[ind_snow,1]
  #
  conc_prmems = data.observations.values[ind_conc,6:-1][:,::2]
  thic_prmems = data.observations.values[ind_thic,6:-1][:,::2]
  snow_prmems = data.observations.values[ind_snow,6:-1][:,::2]
  ####
  conc_bias = []
  conc_rmse = []
  conc_totspr = []
  tracker_conc = []
  conc_qcs = []
  thic_qcs = []
  thic_bias = []
  thic_rmse = []
  thic_totspr = []
  tracker_thic = []
  #
  snow_qcs = []
  snow_bias = []
  snow_rmse = []
  snow_totspr = []
  tracker_snow = []
  #
  qc_zeros = np.zeros((ind_conc.shape[0],3,9))
  qc_zeros_thic = np.zeros((ind_conc.shape[0],3,9))
  qc_zeros_snow = np.zeros((ind_conc.shape[0],3,9))
  rank_vals_zeros = []
  rank_vals_zeros_thic = []
  rank_vals_zeros_snow = []
  for p in range(0,ind_conc.shape[0]):
    if remove_negative_sic:
      if np.any(conc_prmems[p,:]<0):
        print('Negative concentrations...skip!')
        print('DART QC: '+str(conc_qc[p]))
        negative_count += 1
        continue
    if conc_qc[p] == 0.0 or conc_qc[p] == 1.0:
      if conc_truth[p] > 0.0 and conc_truth[p] < 0.15:
        tracker = 1
        if conc_qc[p] == 0.0:
          qc_zeros[p,0,0] += 1
        else:
          qc_zeros[p,0,1] += 1
      elif conc_truth[p] > 0.15 and conc_truth[p] < 0.9:
        tracker = 2
        if conc_qc[p] == 0.0:
          qc_zeros[p,1,0] += 1
        else:
          qc_zeros[p,1,1] += 1
      else:
        tracker = 3
        if conc_qc[p] == 0.0:
          qc_zeros[p,2,0] += 1
        else:
          qc_zeros[p,2,1] += 1
      tracker_conc.append(tracker)
      #if conc_qc[p] == 0.0:
      hold_b = []
      hold_b.append(conc_prmean[p] - conc_truth[p])
      hold_b.append(conc_prmean[p] - conc_obs[p])
      conc_bias.append(hold_b)
      hold_rmse = []
      hold_rmse.append((conc_prmean[p] - conc_truth[p])**2)
      hold_rmse.append((conc_prmean[p] - conc_obs[p])**2)
      conc_rmse.append(hold_rmse)
      hold_totspr = []
      hold_totspr.append(conc_prspr[p]**2)
      hold_totspr.append(conc_prspr[p]**2 + conc_oberr[p])
      conc_totspr.append(hold_totspr)
      conc_qcs.append(conc_qc[p])
      hold_rank = []
      np.random.seed(12*p*d)
      rank_vals = np.sort(conc_prmems[p,:])
      count = np.zeros((80))
      rank = -1
      for i in range(0,rank_vals.shape[0]):
        if conc_truth[p] <= rank_vals[i]:
          rank = i
          break
      if rank == -1:
        rank = 79
      count[rank] += 1
      noise = np.random.normal(0.0,np.sqrt(conc_oberr[p]),size=79)
      rank_vals = np.sort(conc_prmems[p,:] + noise)
      hold_rank.append(count)
      count = np.zeros((80))
      rank = -1
      for i in range(0,rank_vals.shape[0]):
        if conc_obs[p] <= rank_vals[i]:
          rank = i
          break
      if rank == -1:
        rank = 79
      count[rank] += 1
      hold_rank.append(count)
      rank_vals_zeros.append(hold_rank)
    else:
      if conc_truth[p] > 0.0 and conc_truth[p] < 0.15:
        val = conc_qc[p]
        qc_zeros[p,0,val] += 1
      elif conc_truth[p] > 0.15 and conc_truth[p] < 0.9:
        val = conc_qc[p]
        qc_zeros[p,1,val] += 1
      else:
        val = conc_qc[p]
        qc_zeros[p,2,val] += 1
    ###########################
    ind = np.where((locs_conc[p,0]==locs_thic[:,0])&(locs_conc[p,1]==locs_thic[:,1]))[0]
    if ind.shape[0] == 0:
      print('no thickness ob at this concentration pt...')
      continue
    else:
      ind = ind[0]
    if thic_qc[ind] == 0.0:
      if conc_truth[ind] > 0.0 and conc_truth[ind] < 0.15:
        tracker = 1
        qc_zeros_thic[p,0,0] += 1
      elif conc_truth[ind] > 0.15 and conc_truth[ind] < 0.9:
        tracker = 2
        qc_zeros_thic[p,1,0] += 1
      else:
        tracker = 3
        qc_zeros_thic[p,2,0] += 1
      tracker_thic.append(tracker)
      #if conc_qc[p] == 0.0:
      hold_b = []
      hold_b.append(thic_prmean[p] - thic_truth[p])
      hold_b.append(thic_prmean[p] - thic_obs[p])
      thic_bias.append(hold_b)
      hold_rmse = []
      hold_rmse.append((thic_prmean[p] - thic_truth[p])**2)
      hold_rmse.append((thic_prmean[p] - thic_obs[p])**2)
      thic_rmse.append(hold_rmse)
      hold_totspr = []
      hold_totspr.append(thic_prspr[p]**2)
      hold_totspr.append(thic_prspr[p]**2 + thic_oberr[p])
      thic_totspr.append(hold_totspr)
      thic_qcs.append(thic_qc[p])
      hold_rank = []
      np.random.seed(15*p*d)
      rank_vals = np.sort(thic_prmems[p,:])
      count = np.zeros((80))
      rank = -1
      for i in range(0,rank_vals.shape[0]):
        if thic_truth[p] <= rank_vals[i]:
          rank = i
          break
      if rank == -1:
        rank = 79
      count[rank] += 1
      hold_rank.append(count)
      noise = np.random.normal(0.0,np.sqrt(thic_oberr[p]),size=79)
      rank_vals = np.sort(thic_prmems[p,:] + noise)
      count = np.zeros((80))
      rank = -1
      for i in range(0,rank_vals.shape[0]):
        if thic_obs[p] <= rank_vals[i]:
          rank = i
          break
      if rank == -1:
        rank = 79
      count[rank] += 1
      hold_rank.append(count)
      rank_vals_zeros_thic.append(hold_rank)
    else:
      if conc_truth[ind] > 0.0 and conc_truth[ind] < 0.15:
        val = thic_qc[ind]
        qc_zeros_thic[p,0,val] += 1
      elif conc_truth[ind] > 0.15 and conc_truth[ind] < 0.9:
        val = thic_qc[ind]
        qc_zeros_thic[p,1,val] += 1
      else:
        val = thic_qc[ind]
        qc_zeros_thic[p,2,val] += 1
    ###################
    ind = np.where((locs_conc[p,0]==locs_snow[:,0])&(locs_conc[p,1]==locs_snow[:,1]))[0]
    if ind.shape[0] == 0:
      print('no snow thickness ob at this concentration pt...')
      continue
    else:
      ind = ind[0]
    if thic_qc[ind] == 0.0:
      if conc_truth[ind] > 0.0 and conc_truth[ind] < 0.15:
        tracker = 1
        qc_zeros_snow[p,0,0] += 1
      elif conc_truth[ind] > 0.15 and conc_truth[ind] < 0.9:
        tracker = 2
        qc_zeros_snow[p,1,0] += 1
      else:
        tracker = 3
        qc_zeros_snow[p,2,0] += 1
      tracker_snow.append(tracker)
      #if conc_qc[p] == 0.0:
      hold_b = []
      hold_b.append(snow_prmean[p] - snow_truth[p])
      hold_b.append(snow_prmean[p] - snow_obs[p])
      snow_bias.append(hold_b)
      hold_rmse = []
      hold_rmse.append((snow_prmean[p] - snow_truth[p])**2)
      hold_rmse.append((snow_prmean[p] - snow_obs[p])**2)
      snow_rmse.append(hold_rmse)
      hold_totspr = []
      hold_totspr.append(snow_prspr[p]**2)
      hold_totspr.append(snow_prspr[p]**2 + snow_oberr[p])
      snow_totspr.append(hold_totspr)
      snow_qcs.append(snow_qc[p])
      hold_rank = []
      np.random.seed(15*p*d)
      rank_vals = np.sort(snow_prmems[p,:])
      count = np.zeros((80))
      rank = -1
      for i in range(0,rank_vals.shape[0]):
        if snow_truth[p] <= rank_vals[i]:
          rank = i
          break
      if rank == -1:
        rank = 79
      count[rank] += 1
      hold_rank.append(count)
      noise = np.random.normal(0.0,np.sqrt(snow_oberr[p]),size=79)
      rank_vals = np.sort(snow_prmems[p,:] + noise)
      count = np.zeros((80))
      rank = -1
      for i in range(0,rank_vals.shape[0]):
        if snow_obs[p] <= rank_vals[i]:
          rank = i
          break
      if rank == -1:
        rank = 79
      count[rank] += 1
      hold_rank.append(count)
      rank_vals_zeros_snow.append(hold_rank)
    else:
      if conc_truth[ind] > 0.0 and conc_truth[ind] < 0.15:
        val = snow_qc[ind]
        qc_zeros_snow[p,0,val] += 1
      elif conc_truth[ind] > 0.15 and conc_truth[ind] < 0.9:
        val = snow_qc[ind]
        qc_zeros_snow[p,1,val] += 1
      else:
        val = snow_qc[ind]
        qc_zeros_snow[p,2,val] += 1
    ####################
  conc_bias = np.array(conc_bias)
  conc_rmse = np.array(conc_rmse)
  conc_totspr = np.array(conc_totspr)
  tracker_conc = np.array(tracker_conc)
  rank_vals_zeros = np.array(rank_vals_zeros)
  all_stats[:,0,d] = np.nanmean(conc_bias[:,:],axis=0)
  all_stats[:,1,d] = np.sqrt(np.nanmean(conc_rmse[:,:],axis=0))
  all_stats[:,2,d] = np.sqrt(np.nanmean(conc_totspr[:,:],axis=0))
  all_qcs[d,:] = np.nansum(qc_zeros[:,:,:],axis=(0,1))
  all_rankvalues[d,:] = np.nansum(rank_vals_zeros[:,:],axis=0)
  ind = np.where(tracker_conc==1)[0]
  marginal_stats[:,0,d] = np.nanmean(conc_bias[ind,:],axis=0)
  marginal_stats[:,1,d] = np.sqrt(np.nanmean(conc_rmse[ind,:],axis=0))
  marginal_stats[:,2,d] = np.sqrt(np.nanmean(conc_totspr[ind,:],axis=0))
  marginal_qcs[d,:] = np.nansum(qc_zeros[:,0,:],axis=0)
  marginal_rankvalues[d,:] = np.nansum(rank_vals_zeros[ind,:],axis=0)
  ind = np.where(tracker_conc==2)[0] 
  interm_stats[:,0,d] = np.nanmean(conc_bias[ind,:],axis=0)
  interm_stats[:,1,d] = np.sqrt(np.nanmean(conc_rmse[ind,:],axis=0))
  interm_stats[:,2,d] = np.sqrt(np.nanmean(conc_totspr[ind,:],axis=0))
  interm_qcs[d,:] = np.nansum(qc_zeros[:,1,:],axis=0)    
  interm_rankvalues[d,:] = np.nansum(rank_vals_zeros[ind,:],axis=0)
  ind = np.where(tracker_conc==3)[0]
  com_stats[:,0,d] = np.nanmean(conc_bias[ind,:],axis=0)
  com_stats[:,1,d] = np.sqrt(np.nanmean(conc_rmse[ind,:],axis=0))
  com_stats[:,2,d] = np.sqrt(np.nanmean(conc_totspr[ind,:],axis=0))
  com_qcs[d,:] = np.nansum(qc_zeros[:,2,:],axis=0)
  com_rankvalues[d,:] = np.nansum(rank_vals_zeros[ind,:],axis=0)
  #####################################
  thic_bias = np.array(thic_bias)
  thic_rmse = np.array(thic_rmse)
  thic_totspr = np.array(thic_totspr)
  tracker_thic = np.array(tracker_thic)
  rank_vals_zeros_thic = np.array(rank_vals_zeros_thic)
  all_stats_thic[:,0,d] = np.nanmean(thic_bias[:,:],axis=0)
  all_stats_thic[:,1,d] = np.sqrt(np.nanmean(thic_rmse[:,:],axis=0))
  all_stats_thic[:,2,d] = np.sqrt(np.nanmean(thic_totspr[:,:],axis=0))
  all_qcs_thic[d,:] = np.nansum(qc_zeros_thic[:,:,:],axis=(0,1))
  all_rankvalues_thic[d,:] = np.nansum(rank_vals_zeros_thic[:,:],axis=0)
  ind = np.where(tracker_thic==1)[0]
  marginal_stats_thic[:,0,d] = np.nanmean(thic_bias[ind,:],axis=0)
  marginal_stats_thic[:,1,d] = np.sqrt(np.nanmean(thic_rmse[ind,:],axis=0))
  marginal_stats_thic[:,2,d] = np.sqrt(np.nanmean(thic_totspr[ind,:],axis=0))
  marginal_qcs_thic[d,:] = np.nansum(qc_zeros_thic[:,0,:],axis=0)
  marginal_rankvalues_thic[d,:] = np.nansum(rank_vals_zeros_thic[ind,:],axis=0)
  ind = np.where(tracker_thic==2)[0]
  interm_stats_thic[:,0,d] = np.nanmean(thic_bias[ind,:],axis=0)
  interm_stats_thic[:,1,d] = np.sqrt(np.nanmean(thic_rmse[ind,:],axis=0))
  interm_stats_thic[:,2,d] = np.sqrt(np.nanmean(thic_totspr[ind,:],axis=0))
  interm_qcs_thic[d,:] = np.nansum(qc_zeros_thic[:,1,:],axis=0)
  interm_rankvalues_thic[d,:] = np.nansum(rank_vals_zeros_thic[ind,:],axis=0)
  ind = np.where(tracker_thic==3)[0]
  com_stats_thic[:,0,d] = np.nanmean(thic_bias[ind,:],axis=0)
  com_stats_thic[:,1,d] = np.sqrt(np.nanmean(thic_rmse[ind,:],axis=0))
  com_stats_thic[:,2,d] = np.sqrt(np.nanmean(thic_totspr[ind,:],axis=0))
  com_qcs_thic[d,:] = np.nansum(qc_zeros_thic[:,2,:],axis=0) 
  com_rankvalues_thic[d,:] = np.nansum(rank_vals_zeros_thic[ind,:],axis=0)
  #################################
  snow_bias = np.array(snow_bias)
  snow_rmse = np.array(snow_rmse)
  snow_totspr = np.array(snow_totspr)
  tracker_snow = np.array(tracker_snow)
  rank_vals_zeros_snow = np.array(rank_vals_zeros_snow)
  all_stats_snow[:,0,d] = np.nanmean(snow_bias[:,:],axis=0)
  all_stats_snow[:,1,d] = np.sqrt(np.nanmean(snow_rmse[:,:],axis=0))
  all_stats_snow[:,2,d] = np.sqrt(np.nanmean(snow_totspr[:,:],axis=0))
  all_qcs_snow[d,:] = np.nansum(qc_zeros_snow[:,:,:],axis=(0,1))
  all_rankvalues_snow[d,:] = np.nansum(rank_vals_zeros_snow[:,:],axis=0)
  ind = np.where(tracker_snow==1)[0]
  marginal_stats_snow[:,0,d] = np.nanmean(snow_bias[ind,:],axis=0)
  marginal_stats_snow[:,1,d] = np.sqrt(np.nanmean(snow_rmse[ind,:],axis=0))
  marginal_stats_snow[:,2,d] = np.sqrt(np.nanmean(snow_totspr[ind,:],axis=0))
  marginal_qcs_snow[d,:] = np.nansum(qc_zeros_snow[:,0,:],axis=0)
  marginal_rankvalues_snow[d,:] = np.nansum(rank_vals_zeros_snow[ind,:],axis=0)
  ind = np.where(tracker_snow==2)[0]
  interm_stats_snow[:,0,d] = np.nanmean(snow_bias[ind,:],axis=0)
  interm_stats_snow[:,1,d] = np.sqrt(np.nanmean(snow_rmse[ind,:],axis=0))
  interm_stats_snow[:,2,d] = np.sqrt(np.nanmean(snow_totspr[ind,:],axis=0))
  interm_qcs_snow[d,:] = np.nansum(qc_zeros_snow[:,1,:],axis=0)
  interm_rankvalues_snow[d,:] = np.nansum(rank_vals_zeros_snow[ind,:],axis=0)
  ind = np.where(tracker_snow==3)[0]
  com_stats_snow[:,0,d] = np.nanmean(snow_bias[ind,:],axis=0)
  com_stats_snow[:,1,d] = np.sqrt(np.nanmean(snow_rmse[ind,:],axis=0))
  com_stats_snow[:,2,d] = np.sqrt(np.nanmean(snow_totspr[ind,:],axis=0))
  com_qcs_snow[d,:] = np.nansum(qc_zeros_snow[:,2,:],axis=0)
  com_rankvalues_snow[d,:] = np.nansum(rank_vals_zeros_snow[ind,:],axis=0)

print(negative_count)

if remove_negative_sic:
  save_name = 'seaice_dartstats_zones_removeneg.nc'
else:
  save_name = 'seaice_dartstats_zones.nc'

new_file = netCDF4.Dataset(save_name,'w')
tim_dim = new_file.createDimension('Times',None)
ver_dim = new_file.createDimension('verif',2)
stats_dim = new_file.createDimension('stats',3)
qc_dim = new_file.createDimension('qc_dim',9)
rhmems_dim = new_file.createDimension('rankmems',80)

tim_var = new_file.createVariable('times','i8',('Times'))
tim_var[:] = times

all_stats_var = new_file.createVariable('all_stats','f8',('verif','stats','Times'))
all_stats_var[:,:,:] = all_stats

all_qc_var = new_file.createVariable('all_qcs','f8',('Times','qc_dim'))
all_qc_var[:,:] = all_qcs

all_rh_var = new_file.createVariable('all_rankhist','f8',('Times','verif','rankmems'))
all_rh_var[:,:,:] = all_rankvalues

marg_stats_var = new_file.createVariable('marginal_stats','f8',('verif','stats','Times'))
marg_stats_var[:,:,:] = marginal_stats

marg_qc_var = new_file.createVariable('marginal_qcs','f8',('Times','qc_dim'))
marg_qc_var[:,:] = marginal_qcs

marg_rh_var = new_file.createVariable('marginal_rankhist','f8',('Times','verif','rankmems'))
marg_rh_var[:,:,:] = marginal_rankvalues

in_stats_var = new_file.createVariable('interm_stats','f8',('verif','stats','Times'))
in_stats_var[:,:,:] = interm_stats

in_qc_var = new_file.createVariable('interm_qcs','f8',('Times','qc_dim'))
in_qc_var[:,:] = interm_qcs

in_rh_var = new_file.createVariable('interm_rankhist','f8',('Times','verif','rankmems'))
in_rh_var[:,:,:] = interm_rankvalues

com_stats_var = new_file.createVariable('com_stats','f8',('verif','stats','Times'))
com_stats_var[:,:,:] = com_stats

com_qc_var = new_file.createVariable('com_qcs','f8',('Times','qc_dim'))
com_qc_var[:,:] = com_qcs

com_rh_var = new_file.createVariable('com_rankhist','f8',('Times','verif','rankmems'))
com_rh_var[:,:,:] = com_rankvalues
#############
all_stats_var2 = new_file.createVariable('all_stats_thic','f8',('verif','stats','Times'))
all_stats_var2[:,:,:] = all_stats_thic

all_qc_var2 = new_file.createVariable('all_qcs_thic','f8',('Times','qc_dim'))
all_qc_var2[:,:] = all_qcs_thic

all_rh_var2 = new_file.createVariable('all_rankhist_thic','f8',('Times','verif','rankmems'))
all_rh_var2[:,:,:] = all_rankvalues_thic

marg_stats_var2 = new_file.createVariable('marginal_stats_thic','f8',('verif','stats','Times'))
marg_stats_var2[:,:,:] = marginal_stats_thic

marg_qc_var2 = new_file.createVariable('marginal_qcs_thic','f8',('Times','qc_dim'))
marg_qc_var2[:,:] = marginal_qcs_thic

marg_rh_var2 = new_file.createVariable('marginal_rankhist_thic','f8',('Times','verif','rankmems'))
marg_rh_var2[:,:,:] = marginal_rankvalues_thic

in_stats_var2 = new_file.createVariable('interm_stats_thic','f8',('verif','stats','Times'))
in_stats_var2[:,:,:] = interm_stats_thic

in_qc_var2 = new_file.createVariable('interm_qcs_thic','f8',('Times','qc_dim'))
in_qc_var2[:,:] = interm_qcs_thic

in_rh_var2 = new_file.createVariable('interm_rankhist_thic','f8',('Times','verif','rankmems'))
in_rh_var2[:,:,:] = interm_rankvalues_thic

com_stats_var2 = new_file.createVariable('com_stats_thic','f8',('verif','stats','Times'))
com_stats_var2[:,:,:] = com_stats_thic

com_qc_var2 = new_file.createVariable('com_qcs_thic','f8',('Times','qc_dim'))
com_qc_var2[:,:] = com_qcs_thic

com_rh_var2 = new_file.createVariable('com_rankhist_thic','f8',('Times','verif','rankmems'))
com_rh_var2[:,:,:] = com_rankvalues_thic
#############
all_stats_var2 = new_file.createVariable('all_stats_snow','f8',('verif','stats','Times'))
all_stats_var2[:,:,:] = all_stats_snow

all_qc_var2 = new_file.createVariable('all_qcs_snow','f8',('Times','qc_dim'))
all_qc_var2[:,:] = all_qcs_snow

all_rh_var2 = new_file.createVariable('all_rankhist_snow','f8',('Times','verif','rankmems'))
all_rh_var2[:,:,:] = all_rankvalues_snow

marg_stats_var2 = new_file.createVariable('marginal_stats_snow','f8',('verif','stats','Times'))
marg_stats_var2[:,:,:] = marginal_stats_snow

marg_qc_var2 = new_file.createVariable('marginal_qcs_snow','f8',('Times','qc_dim'))
marg_qc_var2[:,:] = marginal_qcs_snow

marg_rh_var2 = new_file.createVariable('marginal_rankhist_snow','f8',('Times','verif','rankmems'))
marg_rh_var2[:,:,:] = marginal_rankvalues_snow

in_stats_var2 = new_file.createVariable('interm_stats_snow','f8',('verif','stats','Times'))
in_stats_var2[:,:,:] = interm_stats_snow

in_qc_var2 = new_file.createVariable('interm_qcs_snow','f8',('Times','qc_dim'))
in_qc_var2[:,:] = interm_qcs_snow

in_rh_var2 = new_file.createVariable('interm_rankhist_snow','f8',('Times','verif','rankmems'))
in_rh_var2[:,:,:] = interm_rankvalues_snow

com_stats_var2 = new_file.createVariable('com_stats_snow','f8',('verif','stats','Times'))
com_stats_var2[:,:,:] = com_stats_snow

com_qc_var2 = new_file.createVariable('com_qcs_snow','f8',('Times','qc_dim'))
com_qc_var2[:,:] = com_qcs_snow

com_rh_var2 = new_file.createVariable('com_rankhist_snow','f8',('Times','verif','rankmems'))
com_rh_var2[:,:,:] = com_rankvalues_snow
###########
new_file.close()

