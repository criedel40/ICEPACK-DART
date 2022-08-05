import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import sys

def trunc_normal_b_cdf_inv(cdf,mu,sigma,b):
  beta = (b-mu)/sigma
  beta_cdf = scipy.stats.norm.cdf(beta)
  xi_cdf = beta_cdf*cdf
  xi = scipy.stats.norm.ppf(xi_cdf)
  x = mu + sigma*xi
  return x
def trunc_normal_a_cdf_inv(cdf,mu,sigma,a):
  alpha = (a - mu)/sigma
  alpha_cdf = scipy.stats.norm.cdf(alpha)
  xi_cdf = (1.0-alpha_cdf)*cdf + alpha_cdf
  #print(xi_cdf)
  xi = scipy.stats.norm.ppf(xi_cdf)
  x = mu + sigma*xi
  return x


ens_size = 80
ens = (np.arange(1,ens_size+1)*1.0)/(ens_size+1.0)
#print(ens)
prior_sd = ens.std(ddof=1)
prior_var = ens.var(ddof=1)
truth = 0.99
obs_sd = truth*0.15
obs_var = obs_sd**2
obs = 0.92697862510663576
#obs = 1.0
bounds = np.array([0.0,1.0])
is_bounded = np.array([True,True])

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
    #like[i] = scipy.stats.norm.pdf(x[i],loc=obs,scale=np.sqrt(obs_var))
    a,b = (bounds-obs)/np.sqrt(obs_var)
    like[i] = scipy.stats.norm.pdf(x[i],loc=obs,scale=np.sqrt(obs_var))
like_dense = (like[1:]+like[:-1])/2.0
prior_sd = np.sqrt(prior_var)
  
  #compute range
if (is_bounded[0] == True or is_bounded[1]==True):
  hold = (1.0/(ens_size+1.0))
#  dist = trunc_normal_a_cdf_inv(1.0-hold,x[0],prior_sd,0.0)
#  #dist = dist*-1.0
#  #left_mean = x[0] + prior_sd*dist
#  left_mean  = dist
#  dist = trunc_normal_b_cdf_inv(hold,x[-1],prior_sd,1.0)
#  #dist = dist*-1.0
#  #right_mean = x[-1] - prior_sd*dist
#  right_mean = dist
#  sys.exit()
  hold = (1.0/(ens_size+1.0))
  alpha_cdf = scipy.stats.norm.cdf(0.0)
  xi_cdf = (1.0-alpha_cdf)*hold + alpha_cdf
  amp = 1.0*(hold/xi_cdf)
  dist_for_unit_sd = scipy.stats.norm.ppf((1.0/(ens_size+1.0))/amp,loc=0.0,scale=1.0)
  dist_for_unit_sd = dist_for_unit_sd*-1.0
  left_mean = x[0] + dist_for_unit_sd*prior_sd

  beta_cdf = scipy.stats.norm.cdf(1.0)
  xi_cdf = beta_cdf*hold
  amp = 1.0*(hold/xi_cdf)
  dist_for_unit_sd = scipy.stats.norm.ppf((1.0/(ens_size+1.0))/amp,loc=0.0,scale=1.0)
  dist_for_unit_sd = dist_for_unit_sd*-1.0
  right_mean = x[-1] - dist_for_unit_sd*prior_sd
  
else:
  dist_for_unit_sd = scipy.stats.norm.ppf((1.0/(ens_size+1.0)),loc=0.0,scale=1.0)
  dist_for_unit_sd = dist_for_unit_sd*-1.0
  left_mean = x[0] + dist_for_unit_sd*prior_sd
  right_mean = x[-1] - dist_for_unit_sd*prior_sd
print(left_mean,right_mean)
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
#print(new_mean_right,prod_weight_right)

for i in range(1,ens_size):
  mass[i] = like_dense[i-1]/(ens_size + 1.0)
mass_sum = np.sum(mass)
nmass = mass/mass_sum
#print(mass_sum)
#print(nmass)
#print(prod_weight_right)
if (is_bounded[0] == True or is_bounded[1]==True):
  a,b = (bounds-new_mean_right)/prior_sd
  right_amp = (prod_weight_right/mass_sum)# - scipy.stats.truncnorm.pdf(1.0,loc=obs,scale=np.sqrt(obs_var),a=a,b=b)/mass_sum
  a,b = (bounds-new_mean_left)/prior_sd
  left_amp = (prod_weight_left/mass_sum)# - scipy.stats.truncnorm.pdf(0.0,loc=obs,scale=np.sqrt(obs_var),a=a,b=b)/mass_sum
else:
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
  #print(umass)
  if umass < cumul_mass[1]:
    print('Left Tail')
    #if (is_bounded[0]==True or is_bounded[1]==True):
    #  #a,b = (bounds-new_mean_left)/new_sd_left
    #  #new_ens[i] = scipy.stats.truncnorm.ppf(umass/left_amp,loc=new_mean_left,scale=new_sd_left,a=a,b=b)
    #  new_ens[i] = scipy.stats.norm.ppf(umass/left_amp,loc=new_mean_left,scale=new_sd_left)
    #else:
#    alpha_cdf = scipy.stats.norm.cdf(0.0,loc=new_mean_left,scale=new_sd_left)
#    xi_cdf = (1.0-alpha_cdf)*umass + alpha_cdf
#    amp = left_amp*(hold/xi_cdf)
    new_ens[i] = scipy.stats.norm.ppf(umass/amp,loc=new_mean_left,scale=new_sd_left)
    #new_ens[i] = trunc_normal_a_cdf_inv(umass/left_amp,new_mean_left,new_sd_left,0.0)
      #print('HOLD UP:{0} {1}'.format(new_ens[i],new_mean_left))
  elif umass > cumul_mass[ens_size]:
    print('Right Tail')
    #if (is_bounded[0] == True or is_bounded[1]==True):
#    beta_cdf = scipy.stats.norm.cdf(1.0,loc=new_mean_right,scale=new_sd_right)
#    xi_cdf = beta_cdf*hold
#    amp = right_amp*(hold/xi_cdf)
    new_ens[i] = scipy.stats.norm.ppf((1.0-umass)/right_amp,loc=new_mean_right,scale=new_sd_right)
    #  print((1.0-umass)/right_amp,new_ens[i])
    new_ens[i] = new_mean_right + (new_mean_right - new_ens[i])
    #else:
    #new_ens[i] = trunc_normal_b_cdf_inv((umass)/right_amp,new_mean_right,new_sd_right,1.0)
    #sys.exit()
    #new_ens[i] = new_mean_right + (new_mean_right - new_ens[i])
    #print(new_ens[i])
  else:
    for j in range(lowest_box,ens_size):
      if umass >= cumul_mass[j] and umass <= cumul_mass[j+1]:
        new_ens[i] = x[j-1] + ((umass - cumul_mass[j]) /
                  (cumul_mass[j+1] - cumul_mass[j])) * (x[j] - x[j-1])
        lowest_box = j
        break
for i in range(0,ens_size):
  obs_inc[sort_ind[i]] = new_ens[i] - x[i]

print(ens)
print(ens+obs_inc)

sys.exit()
umass = (1.0*(80))/(ens_size+1.0)
#plt.plot(scipy.stats.norm.ppf(np.sort(np.linspace(0.0001,0.9999,100)),loc=x[-1],scale=prior_sd),np.sort(np.linspace(0.0001,0.9999,100)),'-k')
#plt.plot(scipy.stats.norm.ppf(np.sort(np.linspace(0.0001,0.9999,100)),loc=new_mean_right,scale=prior_sd),np.sort(np.linspace(0.0001,0.9999,100)),'-r')
#plt.plot(scipy.stats.norm.ppf(np.array([1.0 - (umass)]),loc=x[-1],scale=prior_sd),np.array([1.0 - (umass)]),'ok')
#plt.plot(scipy.stats.norm.ppf(np.array([1.0 - (umass)]),loc=new_mean_right,scale=prior_sd),np.array([1.0 - (umass)]),'or')
#plt.plot(scipy.stats.norm.ppf(np.array([1.0 - (umass)])/right_amp,loc=new_mean_right,scale=prior_sd),np.array([1.0 - (umass)])/right_amp,'xr')

plt.plot(x,scipy.stats.norm.pdf(1.0 - x,loc=x[-1],scale=prior_sd),'-k')
plt.plot(x,scipy.stats.norm.pdf(1.0- x,loc=new_mean_right,scale=prior_sd),'-r')


org = scipy.stats.norm.ppf(np.array([1.0 - (umass)]),loc=x[-1],scale=prior_sd)
new = scipy.stats.norm.ppf(np.array([1.0 - (umass)]),loc=new_mean_right,scale=prior_sd)
new_adj = scipy.stats.norm.ppf(np.array([1.0 - (umass)])/right_amp,loc=new_mean_right,scale=prior_sd)

print(scipy.stats.norm.ppf((1.0 - umass),loc=x[-1],scale=new_sd_right))
print(scipy.stats.norm.ppf((1.0 - umass),loc=new_mean_right,scale=new_sd_right))
print(scipy.stats.norm.ppf((1.0 - umass)/right_amp,loc=new_mean_right,scale=new_sd_right))

print(scipy.stats.norm.cdf(0.3475385764128236,loc=x[-1],scale=new_sd_right))
print(scipy.stats.norm.cdf(-0.3049228471743535,loc=new_mean_right,scale=new_sd_right))
print(scipy.stats.norm.cdf(-0.43099667181318435,loc=new_mean_right,scale=new_sd_right))


org =  scipy.stats.norm.ppf(1.0/(ens_size+1.0),loc=0.0,scale=1.0)
org = org*-1.0
right_test = x[-1] - org*prior_sd
test = scipy.stats.norm.rvs(size=10000,loc=x[-1],scale=prior_sd)
hold = scipy.stats.norm.cdf(np.sort(test),loc=x[-1],scale=prior_sd)
test2 = scipy.stats.norm.rvs(size=10000,loc=right_test,scale=prior_sd)
hold2 = scipy.stats.norm.cdf(np.sort(test2),loc=right_test,scale=prior_sd)



plt.plot(np.sort(test),hold/max(hold),'-k')
plt.plot(np.sort(test2),hold2/max(hold2),'-r')
plt.show()


