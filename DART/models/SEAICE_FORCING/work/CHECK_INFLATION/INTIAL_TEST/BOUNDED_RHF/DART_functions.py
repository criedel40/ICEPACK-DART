import numpy as np



def compute_obsspace_diag(ensmean,ensspr,truth,obs,obserr,mems,ens_size,seed):
  bias_truth = ensmean - truth
  bias_obs = ensmean - obs
  ##
  rmse_truth = np.sqrt(np.nanmean((ensmean-truth)**2))
  rmse_obs = np.sqrt(np.nanmean((ensmean - obs)**2))
  ##
  totspr_truth = np.sqrt(np.nanmean(ensspr**2))
  totspr_obs = np.sqrt(np.nanmean(ensspr**2 + obserr))
  ##
  np.random.seed(12*seed)
  rank_vals = np.sort(mems)
  count_truth = np.zeros((ens_size+1))
  rank = -1
  for i in range(0,rank_vals.shape[0]):
    if truth <= rank_vals[i]:
      rank = i
      break
  if rank == -1:
     rank = 99
  count_truth[rank] += 1
  noise = np.random.normal(0.0,np.sqrt(obserr),size=ens_size)
  rank_vals = np.sort(mems + noise)
  count_obs = np.zeros((ens_size+1))
  rank = -1
  for i in range(0,rank_vals.shape[0]):
    if obs <= rank_vals[i]:
      rank = i
      break
  if rank == -1:
    rank = 99
  count_obs[rank] += 1
  ###
  return bias_truth,bias_obs,rmse_truth,rmse_obs,totspr_truth,totspr_obs,count_truth,count_obs
  
