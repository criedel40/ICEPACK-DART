import xarray as xr
import numpy as np


data = xr.open_dataset('assim_output_EaKF_gaus.nc')
prior = data.input_aice.values[0,:]
#post = data.output_aice.values[0,:]

#obs_incs = post - prior

state_mean = prior.mean()
obs_state_cov = np.sum((prior - state_mean)*(prior - state_mean))/(99-1)
obs_prior_var = prior.var(ddof=1)

reg_coef = obs_state_cov/obs_prior_var

state_incs = reg_coef*obs_incs

new_state_post = prior + state_incs


