import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import CubicSpline

x = np.arange(0,1+0.01,0.01)

max_var = x*(1.0-x)


new_need = False

inflate_factor = 3.0
org_mean = 0.6
org_var = 0.05
inflate_var = org_var*inflate_factor
org_theo_max = org_mean*(1.0-org_mean)
upper_bound = 0.5*(np.sqrt(1.0 - 4.0*org_var) + 1.0)
diff_mean = upper_bound - org_mean
diff_var = org_theo_max - org_var
hypot = np.sqrt(diff_mean**2 + diff_var**2)
inflate_upper_bound = 0.5*(np.sqrt(1.0 - 4.0*inflate_var) + 1.0)
new_mean = inflate_upper_bound - diff_mean

diff_mean2 = inflate_upper_bound - new_mean
new_theo_max = new_mean*(1.0-new_mean)
diff_var2 = new_theo_max - inflate_var
hypot2 = np.sqrt(diff_mean2**2+diff_var2**2)
print(hypot,hypot2)
#ratio = hypot/hypot2
#new_mean = new_mean/ratio

if new_mean < 0.5:
  new_mean = 0.5

#sys.exit()
if org_var > org_theo_max:
  org_var = org_theo_max
  org_mean = 0.5*(np.sqrt(1.0 - 4.0*org_theo_max) + 1.0)
  new_var = inflate_var
  new_mean = 0.5*(np.sqrt(1.0 - 4.0*inflate_var) + 1.0)
  new_need = True
#
upper_bound = 0.5*(np.sqrt(1.0 - 4.0*org_var) + 1.0)
lower_bound = 0.5*(1.0 - np.sqrt(1.0 - 4.0*org_var))
#####
inflate_upper_bound = 0.5*(np.sqrt(1.0 - 4.0*inflate_var) + 1.0)
inflate_lower_bound = 0.5*(1.0 - np.sqrt(1.0 - 4.0*inflate_var))
####################################################
#if not new_need:
#  if org_mean > 0.65:
#    diff_var = org_theo_max - org_var
#    new_theo = (diff_var + inflate_var)
#    new_var = org_theo_max*(1.0 - (org_mean-0.7)/(1.0-0.7)) + (org_mean-0.7)/(1.0-0.7)*new_theo
#    if new_var > 0.25:
#      new_var = 0.25
#      new_mean = 0.25
#    else:
#      new_mean = 0.5*(np.sqrt(1.0 - 4.0*new_var) + 1.0)
#  else:
#    new_mean = org_mean
#    new_var = inflate_var
####################################################

plt.rc('font', weight='bold')
fig = plt.figure()
plt.plot(x,max_var,'-k',label='Theo. Max',zorder=3.0)
plt.plot(np.array([lower_bound,upper_bound]),np.array([org_var,org_var]),'-r',label='Mean Range',zorder=2.0)
plt.plot(np.array([inflate_lower_bound,inflate_upper_bound]),np.array([inflate_var,inflate_var]),'--r',label='Inflated Mean Range')
plt.plot(org_mean,org_var,'k*',label='Org. Mean')
plt.plot(new_mean,inflate_var,'ko',label='Shifted Mean',markersize=4)
plt.xticks(np.arange(0,1+0.1,0.1))
plt.xlim(0,1)
plt.grid(True,axis='y',alpha=0.5,zorder=1)
plt.legend(loc='upper left',ncol=1,fontsize=7)
plt.xlabel('Mean',weight='bold')
plt.ylabel('Variance',weight='bold')
plt.savefig('test.jpg',dpi=550,bbox_inches='tight')
os.system('open test.jpg')




