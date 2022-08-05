import numpy as np
import scipy.stats
from tabulate import tabulate
import matplotlib.pyplot as plt

def bounded_rhf(ens,ens_size,prior_var,obs,obs_var,is_bounded,bounds):
  obs_inc = np.ones((ens_size))*np.nan
  sort_ind = np.argsort(ens)
  x = ens[sort_ind]
  #computed likelihoods
  like = np.zeros((ens_size))
  for i in range(0,ens_size):
    if (is_bounded[0] == True or is_bounded[1]==True):
      if i == 0: print('BOUND')
      a,b = (bounds-obs)/np.sqrt(obs_var)
      like[i] = scipy.stats.truncnorm.pdf(x[i],loc=obs,scale=np.sqrt(obs_var),a=a,b=b)
    else:
      if i == 0:print('Not Bounded')
      like[i] = scipy.stats.norm.pdf(x[i],loc=obs,scale=np.sqrt(obs_var))
  
  like_dense = (like[1:]+like[:-1])/2.0
  prior_sd = np.sqrt(prior_var)
  
  #compute range
  if (is_bounded[0] == True or is_bounded[1]==True):
    hold = 1.0/(ens_size+1.0)
    a,b = (bounds - x[0])/prior_sd
    left_mean = scipy.stats.truncnorm.ppf(1.0-hold,loc=x[0],scale=prior_sd,a=a,b=b)
    a,b = (bounds-x[-1])/prior_sd
    right_mean = scipy.stats.truncnorm.ppf(hold,loc=x[-1],scale=prior_sd,a=a,b=b)
  else:
    dist_for_unit_sd = scipy.stats.norm.ppf(1.0/(ens_size+1.0),loc=0.0,scale=1.0)
    dist_for_unit_sd = dist_for_unit_sd*-1.0
    print(dist_for_unit_sd)
    #print(dist_for_unit_sd)
    left_mean = x[0] + dist_for_unit_sd*prior_sd
    right_mean = x[-1] - dist_for_unit_sd*prior_sd
  #print(left_mean,right_mean)
  left_var = prior_var
  left_sd = prior_sd
  right_var = prior_var
  right_sd = prior_sd
  
  mass = np.zeros((ens_size+1))
  #Assume flat tails for now
  new_var_left = left_var
  new_sd_left = left_sd
  new_mean_left = left_mean
  prod_weight_left = like[0]
  mass[0] = like[0]/(ens_size + 1.0)
  
  new_var_right = right_var
  new_sd_right = right_sd
  new_mean_right = right_mean
  prod_weight_right = like[-1]
  mass[ens_size] = like[-1]/(ens_size + 1.0)
  print(new_mean_right,prod_weight_right)
  
  for i in range(1,ens_size):
    mass[i] = like_dense[i-1]/(ens_size + 1.0)
  mass_sum = np.sum(mass)
  nmass = mass/mass_sum
  print(mass_sum)
  #print(nmass)
  #print(prod_weight_right)
  right_amp = prod_weight_right/mass_sum
  left_amp = prod_weight_left/mass_sum
  #print(new_mean_left,new_mean_right)
  #print(left_amp,right_amp)
  
  cumul_mass = np.zeros((ens_size+2))
  for i in range(1,ens_size+2):
    cumul_mass[i] = cumul_mass[i-1] + nmass[i-1]
  #print(nmass)
  
  lowest_box = 0
  #FIRST MEMBER 0.90006956606374400
  new_ens = np.zeros((ens_size))
  for i in range(0,ens_size):
    umass = (1.0*(i+1))/(ens_size+1.0)
    if umass < cumul_mass[1]:
      print('Left Tail')
      if (is_bounded[0]==True or is_bounded[1]==True):
         a,b = (bounds-new_mean_left)/new_sd_left
         new_ens[i] = scipy.stats.truncnorm.ppf(umass/left_amp,loc=new_mean_left,scale=new_sd_left,a=a,b=b)
      else:
         new_ens[i] = scipy.stats.norm.ppf(umass/left_amp,loc=new_mean_left,scale=new_sd_left)
         print('HOLD UP:{0} {1}'.format(new_ens[i],new_mean_left))
    elif umass > cumul_mass[ens_size]:
      print('Right Tail')
      if (is_bounded[0] == True or is_bounded[1]==True):
         a,b = (bounds-new_mean_right)/new_sd_right
         new_ens[i] = scipy.stats.truncnorm.ppf((1.0-umass)/right_amp,loc=new_mean_right,scale=new_sd_right,a=a,b=b)
         #print(right_amp,new_mean_right,new_sd_right,umass,a,b)
         alpha_cdf = scipy.stats.norm.cdf(a)
         beta_cdf = scipy.stats.norm.cdf(b)
         xi = (beta_cdf-alpha_cdf)*(1.0-umass) + alpha_cdf
         #new_ens[i] = scipy.stats.norm.ppf(xi/right_amp,loc=new_mean_right,scale=new_sd_right)
      else:
         new_ens[i] = scipy.stats.norm.ppf((1.0 - umass)/right_amp,loc=new_mean_right,scale=new_sd_right)
         print(right_amp,new_mean_right,new_sd_right,1.0 - umass)
      #new_ens[i] = new_mean_right + (new_mean_right - new_ens[i])
    else:
      for j in range(lowest_box,ens_size):
        if umass >= cumul_mass[j] and umass <= cumul_mass[j+1]:
          new_ens[i] = x[j-1] + ((umass - cumul_mass[j]) /
                  (cumul_mass[j+1] - cumul_mass[j])) * (x[j] - x[j-1])
          lowest_box = j
          break
  for i in range(0,ens_size):
    obs_inc[sort_ind[i]] = new_ens[i] - x[i]
  return obs_inc,new_ens
  
#ens = np.array([0.89999997615814209, 0.91200000047683716, 0.93000000715255737, 0.98400002717971802, 0.99199998378753662, 0.94999998807907104, 0.94700002670288086, 0.95999997854232788, 0.97000002861022949, 0.97299998998641968],dtype=np.float64)
obs = 0.93634204556225842
#obs = 0.80000000000
obs_sd = np.sqrt(0.15)
obs_var = 0.15
bounds = np.array([0.0,1.0])
is_bounds = np.array([True,True])

mems = np.arange(1,81)
ens = (mems*1.0)/(mems.shape[0]*1.0)


a,b = (np.array([0.0,1.0])-0.99)/0.15

#obs = scipy.stats.truncnorm.rvs(size=5000,a=a,b=b,loc=0.99,scale=0.15)

#mass_array = np.zeros((obs.shape[0],ens.shape[0]+1))
#for i in range(0,obs.shape[0]):
hold,new_ens = bounded_rhf(ens,ens.shape[0],ens.var(ddof=1),obs,obs_var,is_bounds,bounds)

sys.exit()
obs_array = np.zeros((hold.shape[0]))
obs_array[:] = obs

total = np.zeros((4,hold.shape[0]),dtype=np.float64)
total[0,:] = obs_array
total[1,:] = ens
total[2,:] = ens + hold
total[3,:] = hold

#print(tabulate(total.T,["Observation", "Prior", "Post","Obs Inc"], tablefmt="grid",floatfmt=".16f"))






sys.exit()
sort_post = np.sort(total[2,:])
sort_cdf = scipy.stats.norm.cdf(sort_post,loc=sort_post.mean(),scale=sort_post.std(ddof=1))
for l in range(0,sort_post.shape[0]):
  if l == 0:
    print(sort_cdf[l] - 0.0)
  elif l == sort_post.shape[0]-1:
    print(1.0 - sort_cdf[l])
  else:
    print(sort_cdf[l] - sort_cdf[l-1])


#############TEST#############
#check = np.zeros((total.shape[1])+2)
#check[1:-1] = total[2,:]
#check[0] = 0.9404724869567163
#check[-1] = 0.9515274729889625
umass = 1.0/(10.0+1.0)
ens_size = ens.shape[0]
x = np.sort(ens)
prior_sd = ens.std(ddof=1)
prior_var = ens.var(ddof=1)
dist_for_unit_sd = scipy.stats.norm.ppf(1.0/(ens_size+1.0),loc=0.0,scale=1.0)
dist_for_unit_sd = dist_for_unit_sd*-1.0
  #print(dist_for_unit_sd)
left_mean = x[0] + dist_for_unit_sd*prior_sd
right_mean = x[-1] - dist_for_unit_sd*prior_sd
left_var = prior_var
left_sd = prior_sd
right_var = prior_var
right_sd = prior_sd
#################################
mass_sum = 0.6293898816697646
new_mean_right = 0.9515274729889625
prod_weight_right = 0.6306696000663637
right_amp = prod_weight_right/mass_sum

first_norm = scipy.stats.norm.rvs(size=10000,loc=x[-1],scale=prior_sd)
first_mean_right = scipy.stats.norm.ppf(1.0/(10.+1.0),loc=x[-1],scale=prior_sd)
first_norm_sort = np.sort(first_norm)
first_cdf = scipy.stats.norm.cdf(first_norm_sort,loc=x[-1],scale=prior_sd)
point_cdf = scipy.stats.norm.cdf(first_mean_right,loc=x[-1],scale=prior_sd)
ind_cdf1 = abs(point_cdf - first_cdf).argmin()
#####
first_norm = scipy.stats.norm.rvs(size=10000,loc=first_mean_right,scale=prior_sd)
first_mean_right = scipy.stats.norm.ppf(1.0/(10.+1.0),loc=first_mean_right,scale=prior_sd)
first_norm_sort = np.sort(first_norm)
first_cdf = scipy.stats.norm.cdf(first_norm_sort,loc=first_mean_right,scale=prior_sd)
point_cdf = scipy.stats.norm.cdf(first_mean_right,loc=first_mean_right,scale=prior_sd)
ind_cdf1 = abs(point_cdf - first_cdf).argmin()



fig = plt.figure()
plt.plot(first_norm_sort,first_cdf,'-k',label='Orginal')
plt.plot(first_norm_sort[ind_cdf1],first_cdf[ind_cdf1],'ok')
#plt.plot(sort_second_mean,second_cdf,'-r',label='Need Adjust')
#plt.plot(sort_second_mean[ind_cdf2],second_cdf[ind_cdf2],'or')
#plt.plot(sort_second_mean[ind_cdf22],second_cdf[ind_cdf22],'xr')
plt.legend()
plt.show()

sys.exit()
fig = plt.figure()
plt.plot(sort_first_mean,first_pdf/max(first_pdf),'-k',label='Orginal')
plt.plot(sort_first_mean[index_first],first_pdf[index_first]/max(first_pdf),'ok')
plt.plot(sort_second_mean,second_pdf/max(second_pdf),'-r',label='Need Adjust')
plt.plot(sort_second_mean[index_second],second_pdf[index_second]/max(second_pdf),'or')
plt.legend()
plt.show()

