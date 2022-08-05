import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import sys

def truncated_normal_ab_cdf_inv( cdf, mu, sigma, a, b):
  alpha = (a-mu)/sigma
  beta = (b-mu)/sigma
  
  alpha_cdf = scipy.stats.norm.cdf(alpha)
  beta_cdf = scipy.stats.norm.cdf(beta)
  
  xi_cdf = (beta_cdf - alpha_cdf)*cdf + alpha_cdf
  xi = scipy.stats.norm.ppf(xi_cdf)
  x = mu + sigma*xi
  return x
  
def truncated_normal_a_cdf_inv( cdf, mu, sigma, a):
  alpha = (a - mu)/sigma
  
  alpha_cdf = scipy.stats.norm.cdf(alpha)
  xi_cdf = (1.0-alpha_cdf)*cdf + alpha_cdf
  xi = scipy.stats.norm.ppf(xi_cdf)
  x = mu + sigma*xi
  return x

def truncated_normal_b_cdf_inv( cdf, mu, sigma, b):
  beta = ( b - mu ) / sigma
  beta_cdf = scipy.stats.norm.cdf(beta)
  xi_cdf = beta_cdf*cdf
  xi = scipy.stats.norm.ppf(xi_cdf)
  x = mu + sigma*xi
  return x
def pdf_ab(x,mu,sigma,a,b):
  alpha = (a - mu)/sigma
  beta = (b-mu)/sigma
  xi = (x - mu)/sigma

  alpha_cdf = scipy.stats.norm.cdf(alpha)
  beta_cdf = scipy.stats.norm.cdf(beta)
  xi_pdf = np.exp ( -0.5 * xi * xi ) / np.sqrt ( 2.0 * np.pi )
  pdf = xi_pdf/(beta_cdf - alpha_cdf)/sigma
  return pdf
def pdf_b(x,mu,sigma,a,b):
  beta = (b-mu)/sigma
  xi = (x - mu)/sigma
  beta_cdf = scipy.stats.norm.cdf(beta)
  xi_pdf = np.exp ( -0.5 * xi * xi ) / np.sqrt ( 2.0 * np.pi )
  pdf = xi_pdf/beta_cdf/sigma
  return pdf
def pdf_a(x,mu,sigma,a,b):
  alpha = (a - mu)/sigma
  xi = (x - mu)/sigma
  alpha_cdf = scipy.stats.norm.cdf(alpha)
  xi_pdf = np.exp ( -0.5 * xi * xi ) / np.sqrt ( 2.0 * np.pi )
  pdf = xi_pdf / ( 1.0 - alpha_cdf ) / sigma
  return pdf

hold = np.array([0.9, 0.912, 0.93, 0.984, 0.992, 0.95, 0.947, 0.96, 0.97, 0.973])
post = np.array([0.90006956606374400,0.91218833013610845,0.93021552742651448,0.98403313070642928,0.99203418698744195,0.95011722612156169,0.94703811566138196,0.96010318786349280,0.97002529462339304,0.97306935559645813])
x = np.sort(hold)
mean = hold.mean()
std = hold.std(ddof=1)
obs = 0.98000001907348633
obs_var = 0.4
obs_std = np.sqrt(obs_var)
new_mean_right = 0.95152747301110041
new_sd_right = 3.0312451820997753E-002
right_amp = 1.0020332700558889

new_mean_left = 0.94047248693457830
new_sd_left = 3.0312451820997753E-002
left_amp = 0.99422791907545582


low_mass = (1.0/(10.+1.0))
high_mass = (1.0*10.)/(10.+1.0)
########## NORMAL NO BOUNDS###############
#low = scipy.stats.norm.ppf(low_mass,loc=new_mean_left,scale=new_sd_left)
#high = scipy.stats.norm.ppf((1.0-high_mass),loc=new_mean_right,scale=new_sd_right)
#test_high = new_mean_right + (new_mean_right - high)
#high = scipy.stats.norm.ppf(high_mass,loc=new_mean_right,scale=new_sd_right)
#
#print(low,high,test_high)
#low_cdf =scipy.stats.norm.cdf(low,loc=new_mean_left,scale=new_sd_left)
#high_cdf =scipy.stats.norm.cdf(high,loc=new_mean_right,scale=new_sd_right)
#high_cdf_test =scipy.stats.norm.cdf(test_high,loc=new_mean_right,scale=new_sd_right)
#print(low_cdf,1.0-high_cdf,1.0-high_cdf_test)

########## NORMAL DOUBLE BOUNDS############
a,b = (np.array([0.0,1.0]) - new_mean_right)/new_sd_right
low = scipy.stats.truncnorm.ppf(low_mass,loc=new_mean_left,scale=new_sd_left,a=a,b=b)
high = scipy.stats.truncnorm.ppf((1.0-high_mass),loc=new_mean_right,scale=new_sd_right,a=a,b=b)
test_high = new_mean_right + (new_mean_right - high)
right_amp = high_mass/(1.0 - ((1.0 - high_mass)/left_amp))
high = scipy.stats.truncnorm.ppf(high_mass,loc=new_mean_right,scale=new_sd_right,a=a,b=b)
print(low,high,test_high)
low_cdf =scipy.stats.truncnorm.cdf(low,loc=new_mean_left,scale=new_sd_left,a=a,b=b)
high_cdf =scipy.stats.truncnorm.cdf(high,loc=new_mean_right,scale=new_sd_right,a=a,b=b)
high_cdf_test =scipy.stats.truncnorm.cdf(test_high,loc=new_mean_right,scale=new_sd_right,a=a,b=b)
print(low_cdf,1.0-high_cdf,1.0-high_cdf_test)



sys.exit()
#cdf = 1.0/(79+1.0)
cdf = 9.0909090909090912E-002
print(pdf_a(x[0],obs,np.sqrt(obs_var),0,1))
a,b = (np.array([0.0,np.inf]) - obs)/np.sqrt(obs_var)
print(scipy.stats.truncnorm.pdf(x[0],loc=obs,scale=np.sqrt(obs_var),a=a,b=b))
sys.exit()
print(scipy.stats.norm.ppf(cdf,loc=mean,scale=std),scipy.stats.norm.ppf(1.0-cdf,loc=mean,scale=std))
#print(scipy.stats.truncnorm.pdf()
bounds = np.array([0.0,1.0])
a,b = (bounds-mean)/std

print('-------------------------------')
print(truncated_normal_ab_cdf_inv(cdf,mean,std,bounds[0],bounds[1]),scipy.stats.truncnorm.ppf(cdf,a,b,loc=mean,scale=std))
print(truncated_normal_ab_cdf_inv(1.0-cdf,mean,std,bounds[0],bounds[1]),scipy.stats.truncnorm.ppf(1.0-cdf,a,b,loc=mean,scale=std))
print('-------------------------------')
bounds = np.array([0.0,np.inf])
a,b = (bounds-mean)/std
print(truncated_normal_a_cdf_inv(cdf,mean,std,bounds[0]),scipy.stats.truncnorm.ppf(cdf,a,b,loc=mean,scale=std))
print(truncated_normal_a_cdf_inv(1.0-cdf,mean,std,bounds[0]),scipy.stats.truncnorm.ppf(1.0-cdf,a,b,loc=mean,scale=std))
print('-------------------------------')

bounds = np.array([-np.inf,1.0])
a,b = (bounds-mean)/std
print(truncated_normal_b_cdf_inv(cdf,mean,std,bounds[1]),scipy.stats.truncnorm.ppf(cdf,a,b,loc=mean,scale=std))
print(truncated_normal_b_cdf_inv(1.0-cdf,mean,std,bounds[1]),scipy.stats.truncnorm.ppf(1.0-cdf,a,b,loc=mean,scale=std))
print('-------------------------------')


#####################################
normal_rand = scipy.stats.norm.rvs(size=10000)
cdf_val = scipy.stats.norm.ppf(cdf)
####
hold = np.array([0.9, 0.912, 0.93, 0.984, 0.992, 0.95, 0.947, 0.96, 0.97, 0.973])
mean = hold.mean()
std = hold.std(ddof=1)
x = np.sort(hold)
bounds = np.array([0.0,1.0])
cdf = 1.0/(10.+1.0)
obs = 0.98
obs_var = 0.4
#xi = scipy.stats.norm.ppf(cdf,loc=0.0,scale=1.0)
#xi = xi*-1.0
#
#left = 0.0+1.0*xi
#right = 0.0-1.0*xi
#print(left,right)
#print(scipy.stats.norm.cdf(right)-scipy.stats.norm.cdf(left))
#al,bl = (bounds-x[0])/std
left = scipy.stats.truncnorm.ppf(cdf,a=0,b=1,loc=x[0],scale=std)
right = scipy.stats.truncnorm.ppf(cdf,a=0,b=1,loc=x[-1],scale=std)
#ar,br = (bounds-x[-1])/std
#xi2 = scipy.stats.truncnorm.ppf(1.0 - (1.0/(hold.shape[0]-1.0)),a=0,b=1,loc=0.0,scale=1.0)
#xi = scipy.stats.norm.ppf()

#xi = xi
##xi2 = xi2*-1.0
#left = x[0] + xi*std
#right = x[-1] - xi*std

#print(scipy.stats.truncnorm.cdf(left,a=0,b=1,loc=0.0,scale=1.0) - scipy.stats.truncnorm.cdf(right,a=0,b=1,loc=0.0,scale=1.0))
#print(scipy.stats.truncnorm.cdf(left,a=0,b=1,loc=0.0,scale=1.0) - scipy.stats.truncnorm.cdf(right,a=0,b=1,loc=0.0,scale=1.0))

#truncnorm_rand = scipy.stats.truncnorm.rvs(size=10000,a=0.0,b=1.0,loc=hold.mean(),scale=hold.std())
#trunccdf_val = scipy.stats.truncnorm.ppf(cdf,a=0.0,b=1.0,loc=hold.mean(),scale=hold.std())


#plt.hist(normal_rand,bins=50)
#plt.plot(cdf_val)
#plt.show()
#
#plt.hist(truncnorm_rand,bins=50)
#plt.plot(trunccdf_val)
#plt.show()
