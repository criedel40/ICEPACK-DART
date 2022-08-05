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
  print(xi_cdf)
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


sort_ind = np.argsort(ens)
x = ens[sort_ind]
#computed likelihoods
like_normal = np.zeros((ens_size))
like_trunc = np.zeros((ens_size))
for i in range(0,ens_size):
    a,b = (bounds-obs)/np.sqrt(obs_var)
    like_trunc[i] = scipy.stats.truncnorm.pdf(x[i],loc=obs,scale=np.sqrt(obs_var),a=a,b=b)
    like_normal[i]= scipy.stats.truncnorm.pdf(x[i],loc=obs,scale=np.sqrt(obs_var),a=a,b=b) #= scipy.stats.norm.pdf(x[i],loc=obs,scale=np.sqrt(obs_var))


#################################
hold = 1.0/(ens_size+1.0)
a,b = (bounds - x[0])/prior_sd
left_mean_trunc = scipy.stats.truncnorm.ppf((1.0-hold),loc=x[0],scale=prior_sd,a=a,b=np.inf)
a,b = (bounds-x[-1])/prior_sd
right_mean_tunc = scipy.stats.truncnorm.ppf(hold,loc=x[-1],scale=prior_sd,b=b,a=np.inf)

alpha = (0.0-x[0])/prior_sd
xi_cdf = (1.0-scipy.stats.norm.cdf(alpha))*(1.0-hold) + scipy.stats.norm.cdf(alpha)
xi = scipy.stats.norm.ppf(xi_cdf)
xi = xi*-1.
left_mean_normal = x[0] + xi*prior_sd

sys.exit()
alpha = (0.0 - x[-1])/prior_sd
beta = (1.0 - x[-1])/prior_sd
test = scipy.stats.norm.cdf(beta)*hold
#dist_for_unit_sd = scipy.stats.norm.ppf(test,loc=0.0,scale=1.0)
#dist_for_unit_sd = dist_for_unit_sd*-1.0
#right_mean_normal = x[-1] - dist_for_unit_sd*prior_sd
right_mean_normal = scipy.stats.norm.ppf(test,loc=x[-1],scale=prior_sd)
#################################
print('Normal Prior Tail Values')
print(left_mean_normal,right_mean_normal)
print('Trunc Prior Tail Values')
print(left_mean_trunc,right_mean_tunc)
sys.exit()
###################
like_dense_normal = (like_normal[1:]+like_normal[:-1])/2.0
like_dense_trunc = (like_trunc[1:]+like_trunc[:-1])/2.0
###########################
prior_sd = np.sqrt(prior_var)

left_var = prior_var
left_sd = prior_sd
right_var = prior_var
right_sd = prior_sd
###############################
mass_normal = np.zeros((ens_size+1))
mass_trunc = np.zeros((ens_size+1))

prod_weight_left_normal = like_normal[0]
mass_normal[0] = like_normal[0]/(ens_size + 1.0)

prod_weight_left_trunc = like_trunc[0]
mass_trunc[0] = like_trunc[0]/(ens_size + 1.0)

prod_weight_right_normal = like_normal[-1]
mass_normal[ens_size] = like_normal[-1]/(ens_size + 1.0)

prod_weight_right_trunc = like_trunc[-1]
mass_trunc[ens_size] = like_trunc[-1]/(ens_size + 1.0)
#################################

for i in range(1,ens_size):
  mass_normal[i] = like_dense_normal[i-1]/(ens_size + 1.0)
  mass_trunc[i] = like_dense_trunc[i-1]/(ens_size + 1.0)

mass_sum_normal = np.sum(mass_normal)
nmass_normal = mass_normal/mass_sum_normal

mass_sum_trunc = np.sum(mass_trunc)
nmass_trunc = mass_trunc/mass_sum_trunc
###################################
right_amp_normal = prod_weight_right_normal/mass_sum_normal
left_amp_normal = prod_weight_left_normal/mass_sum_normal
  
right_amp_trunc = prod_weight_right_trunc/mass_sum_trunc
left_amp_trunc = prod_weight_left_trunc/mass_sum_trunc
#####################################
print('Normal Amps')
print(left_amp_normal,right_amp_trunc)
print('Trunc Amps')
print(left_amp_trunc,right_amp_trunc)


umass = (80*1.0)/(80+1.0)
print('-----------------------------------------')
print('NORMAL-Adjustment')
alpha = (0.0 - right_mean_normal)/prior_sd
beta = (1.0 - right_mean_normal)/prior_sd
#test2 = (scipy.stats.norm.cdf(beta) - scipy.stats.norm.cdf(alpha)) + scipy.stats.norm.cdf(alpha)/(1.0 - umass)
hold_right2 = right_amp_normal
hold_right = (1/right_amp_trunc)*((scipy.stats.norm.cdf(beta) - scipy.stats.norm.cdf(alpha))*(1.0-umass) + scipy.stats.norm.cdf(alpha)*right_amp_trunc)
new_ens_normal = scipy.stats.norm.ppf((1.0-umass)/hold_right2,loc=right_mean_normal,scale=prior_sd)
print(new_ens_normal)
new_ens_normal = right_mean_normal + (right_mean_normal - new_ens_normal)
print(new_ens_normal)

print('TRUNC')
a,b = (bounds-right_mean_tunc)/prior_sd
right_HOLD = right_amp_trunc
new_ens_trunc = scipy.stats.truncnorm.ppf((1.0 - umass)/right_HOLD,loc=right_mean_tunc,scale=prior_sd,a=a,b=b)
print(new_ens_trunc)
invers_cdf = 1.0 - scipy.stats.truncnorm.cdf(new_ens_trunc,loc=right_mean_tunc,scale=prior_sd,a=a,b=b)
new_ens_trunc = scipy.stats.truncnorm.ppf(invers_cdf,loc=right_mean_tunc,scale=prior_sd,a=a,b=b)
print(new_ens_trunc)
print('-----------------------------------------')





sys.exit()
#########################
mean = 0.9
sd = mean*0.15
print(scipy.stats.norm.ppf(0.25,loc=mean,scale=sd),scipy.stats.norm.ppf(0.5,loc=mean,scale=sd),scipy.stats.norm.ppf(0.75,loc=mean,scale=sd))
print('--------------------------------------')
a,b = (np.array([0.0,1.0]) - mean)/sd
print(scipy.stats.truncnorm.ppf(0.25,loc=mean,scale=sd,a=a,b=b),scipy.stats.truncnorm.ppf(0.5,loc=mean,scale=sd,a=a,b=b),scipy.stats.truncnorm.ppf(0.75,loc=mean,scale=sd,a=a,b=b))
