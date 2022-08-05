import matplotlib.pyplot as plt
import scipy.stats
import numpy as np

ens_size = 80.
ens = (np.arange(1,ens_size+1)*1.0)/(ens_size+1.0)
x = np.sort(ens)
prior_sd = ens.std(ddof=1)

alpha = (0.0 - 0.0)/1.0
beta = (1.0 - 0.0)/1.0

cdf = 1.0/(80+1.0)

input_xi = (1.0 - scipy.stats.norm.cdf(alpha))*cdf + scipy.stats.norm.cdf(alpha)
dist = scipy.stats.norm.ppf(input_xi)
print(dist)
dist = dist*-1.0
left_mean = x[0] + dist*prior_sd

input_xi = scipy.stats.norm.cdf(beta)*cdf
dist = scipy.stats.norm.ppf(input_xi)
print(dist)
dist = dist*-1.0
right_mean = x[-1] - dist*prior_sd

print(left_mean,right_mean)
print('-------------------')
left_mean2 = scipy.stats.truncnorm.ppf(1.0-cdf,loc=x[0],scale=prior_sd,a=(0.0-x[0])/prior_sd,b=np.inf)
right_mean2 = scipy.stats.truncnorm.ppf(cdf,loc=x[-1],scale=prior_sd,b=(1.0-x[-1])/prior_sd,a=-np.inf)
print(left_mean2,right_mean2)

sys.exit()
ens_size = 80
ens = (np.arange(1,ens_size+1)*1.0)/(ens_size+1.0)
x = np.sort(ens)
prior_sd = ens.std(ddof=1)
prior_var = ens.var(ddof=1)
bounded_above = False
if bounded_above:
  mean_loc = x[-1]
  a = -np.inf
  b = (1.0-mean_loc)/prior_sd
  cdf = 1.0/(ens_size+1.0)
else:
  mean_loc = x[0]
  a = (0.0-mean_loc)/prior_sd
  b = np.inf
  
prior_sd = ens.std(ddof=1)
prior_var = ens.var(ddof=1)

cdf = 1.0/(ens_size+1.0)

sorted_rand_norm = np.sort(scipy.stats.norm.rvs(size=10000,loc=mean_loc,scale=prior_sd))
#sorted_rand_norm = mean_loc + prior_sd*sorted_rand_norm
pdf = scipy.stats.norm.pdf(sorted_rand_norm,loc=mean_loc,scale=prior_sd)


sorted_rand_trunc = np.sort(scipy.stats.truncnorm.rvs(size=10000,loc=mean_loc,scale=prior_sd,a=a,b=b))
#sorted_rand_trunc = mean_loc + prior_sd*sorted_rand_trunc
pdf_trunc = scipy.stats.truncnorm.pdf(sorted_rand_trunc,loc=mean_loc,scale=prior_sd,a=a,b=b)

print('------Normal Dist------')
lown =scipy.stats.norm.ppf(cdf,loc=mean_loc,scale=prior_sd)
midn =scipy.stats.norm.ppf(0.5,loc=mean_loc,scale=prior_sd)
highn =scipy.stats.norm.ppf(1.0-cdf,loc=mean_loc,scale=prior_sd)
print(lown,midn,highn)
print(scipy.stats.norm.cdf(lown,loc=mean_loc,scale=prior_sd),scipy.stats.norm.cdf(mean_loc,loc=mean_loc,scale=prior_sd),scipy.stats.norm.cdf(highn,loc=mean_loc,scale=prior_sd))
print('------Trunc Dist------')
lowt =scipy.stats.truncnorm.ppf(cdf,loc=mean_loc,scale=prior_sd,a=a,b=b)
midt =scipy.stats.truncnorm.ppf(0.5,loc=mean_loc,scale=prior_sd,a=a,b=b)
hight =scipy.stats.truncnorm.ppf(1.0-cdf,loc=mean_loc,scale=prior_sd,a=a,b=b)
print(lowt,midt,hight)
print(scipy.stats.truncnorm.cdf(lowt,loc=mean_loc,scale=prior_sd,a=a,b=b),scipy.stats.truncnorm.cdf(mean_loc,loc=mean_loc,scale=prior_sd,a=a,b=b),scipy.stats.truncnorm.cdf(hight,loc=mean_loc,scale=prior_sd,a=a,b=b))

if bounded_above:
  beta = scipy.stats.norm.cdf(b)
  adjusted_cdf = beta*cdf
else:
  alpha = scipy.stats.norm.cdf(a)
  adjusted_cdf = (1.0-alpha)*(1.0-cdf) + alpha
adjusted_val =scipy.stats.norm.ppf(adjusted_cdf,loc=mean_loc,scale=prior_sd)

fig = plt.figure()
plt.plot(sorted_rand_norm,pdf,'-k')
plt.plot(sorted_rand_trunc,pdf_trunc,'-r')
plt.plot(highn,scipy.stats.norm.pdf(highn,loc=mean_loc,scale=prior_sd),'ok')
plt.plot(hight,scipy.stats.truncnorm.pdf(hight,loc=mean_loc,scale=prior_sd,a=a,b=b),'or')
plt.plot(adjusted_val,scipy.stats.norm.pdf(adjusted_val,loc=mean_loc,scale=prior_sd),'xk')
#plt.xlim(min(sorted_rand_norm)-0.1,max(sorted_rand_norm)+0.1)
plt.show()



