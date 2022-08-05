import scipy.stats
import scipy.special
import os
import numpy as np
import matplotlib.pyplot as plt
a = 75.0
b = 66334470.0

mean,variance = scipy.stats.beta.stats(a,b)
print('Distribution mean and variance: {0} - {1}'.format(mean,variance))

rand_dist = scipy.stats.beta.rvs(a,b,size=10000)
x = np.linspace(scipy.stats.beta.ppf(0.01, a, b),
                scipy.stats.beta.ppf(0.99, a, b), 100)

fig = plt.figure()
rv = scipy.stats.beta(a,b)
plt.plot(x, rv.pdf(x), 'k-', lw=2, label='pdf')
plt.hist(rand_dist,density=True,bins=50)
plt.xlim(0,2.0e-6)
plt.savefig('dist_example.jpg',dpi=550,bbox_inches='tight')
os.system('open dist_example.jpg')

reg_beta = scipy.stats.beta.ppf(0.9999999999997369,a,b)
cephens_beta = scipy.special.btdtri(a, b,0.9999999999997369)
boost_beta = scipy.stats.beta.isf(1.0-0.9999999999997369, a, b)
matlab_beta = scipy.special.betaincinv(a,b,0.9999999999997369)

print('-----------------------------')
print('Different Methods for Computing Beta Inverse CDF')
print('Default Scipy Method: {0}'.format(reg_beta))
print('Cephes Mathematical Function: {0}'.format(cephens_beta))
print('Boost C++ Function: {0}'.format(boost_beta))
print('Matlab Function: {0}'.format(matlab_beta))
