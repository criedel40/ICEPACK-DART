from scipy.stats import lognorm
import numpy as np
import matplotlib.pyplot as plt

def lognorm_params(mode, stddev):
    """
    Given the mode and std. dev. of the log-normal distribution, this function
    returns the shape and scale parameters for scipy's parameterization of the
    distribution.
    """
    p = np.poly1d([1, -1, 0, 0, -(stddev/mode)**2])
    r = p.roots
    sol = r[(r.imag == 0) & (r.real > 0)].real
    shape = np.sqrt(np.log(sol))
    scale = mode * sol
    return shape, scale


mode = 0.2
sd = mode*0.15

sigma, scale = lognorm_params(mode, sd)

sample = lognorm.rvs(sigma, 0, scale, size=1000000)
