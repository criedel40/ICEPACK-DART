import numpy as np
import scipy.stats

def randg(a):
  if a < 1:
    u = scipy.stats.uniform.rvs(size=1)
    g = scipy.stats.norm.rvs(loc=1.0+a,scale=1.0,size=1) * u**(1.0/a)
    return g
  d = a - 1.0/3.0
  c = (1.0/3.0)/np.sqrt(d)
  v = -1.0
  
  STOP=False
  while STOP != True:
    v = -1.0
    while v <= 0.0:
      x = scipy.stats.norm.rvs(size=1)
      v = 1.0 + c*x
      
    v = v*v*v
    u = scipy.stats.uniform.rvs(size=1)
    if u < 1.0 - 0.0331 * x * x * x * x: STOP=True
    if np.log(u) < 0.5 * x * x + d * (1 - v + np.log(v)): STOP=True
  g = 1.0*d*v
  return g
##########
def betarnd(mean,variance):
  mu = mean
  if (variance > mu*(1.0-mu)):
    mu = 0.5*(np.sqrt(1.0 - 4.0*(variance+1.0e-5)) + 1.0)


  a = mu*(((mu*(1.0-mu))/variance) - 1.0)
  b = (1.0-mu)*(((mu*(1.0-mu))/variance) - 1.0)

  if a<0 or b<0:
    print('Bad beta input parameters..STOP')
    sys.exit()
    
  #g1 = scipy.stats.gamma.rvs(a,size=1)
  #g2 = scipy.stats.gamma.rvs(b,size=1)
  g1 = randg(a)
  g2 = randg(b)
  
  
  if g1 == 0 or g2 == 0:
    mu = 0.5*(np.sqrt(1.0 - 4.0*(variance*2.0)) + 1.0)
    a = mu*(((mu*(1.0-mu))/variance) - 1.0)
    b = (1.0-mu)*(((mu*(1.0-mu))/variance) - 1.0)
    g1 = randg(a)
    g2 = randg(b)
    if g1 == 0 or g2 == 0:
      p = a/(a+b)
      #r = scipy.stats.binom.rvs(n=1,p=p,size=1)
      #r = np.sum(scipy.stats.uniform.rvs(size=1) < p)
      if scipy.stats.uniform.rvs(size=1) < p:
        r = 1.0
      else:
        r = 0.0
    else:
      r = g1/(g1+g2)
  else:
    r = g1/(g1+g2)
  return r
###########

def compute_obsspace_diag(ensmean,ensvar,truth,obs,obsvar,mems,ens_size,seed,dist_type):
  bias_truth = ensmean - truth
  bias_obs = ensmean - obs
  ##
  rmse_truth = np.sqrt(np.nanmean((ensmean-truth)**2))
  rmse_obs = np.sqrt(np.nanmean((ensmean - obs)**2))
  ##
  totspr_truth = np.sqrt(np.nanmean(ensvar))
  totspr_obs = np.sqrt(np.nanmean(ensvar + obsvar))
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
  if dist_type == 1:
    noise = np.random.normal(0.0,np.sqrt(obsvar),size=ens_size)
    rank_vals = np.sort(mems + noise)
  if dist_type == 2:
    rank_vals = np.zeros((ens_size))
    for m in range(0,ens_size):
      a = (0.0-mems[m])/np.sqrt(obsvar)
      b = (1.0-mems[m])/np.sqrt(obsvar)
      rank_vals[m] = scipy.stats.truncnorm.rvs(a,b,loc=mems[m],scale=np.sqrt(obsvar),size=1)
    rank_vals = np.sort(rank_vals)
  if dist_type == 3:
    rank_vals = np.zeros((ens_size))
    for m in range(0,ens_size):
      rank_vals[m] = betarnd(mems[m],obsvar)
    rank_vals = np.sort(rank_vals)
  #print(mems)
  #print(rank_vals)
  #sys.exit()
  #rank_vals = np.sort(mems + noise)
  #rank_vals = np.sort(rank_vals)
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
  
