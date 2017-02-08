"""
Various tool functions used throughout the package.

Author: Panos Tsilifis
Date: 2/7/2017

"""


__all__ = ['dirac1D', 'diracND', 'inverse_transform_sampling']


import numpy as np
import scipy.stats as st
import orthpol
from scipy import interpolate

def dirac1D(a,b):
	"""
	Univariate Dirac function.
	"""
    val = 0 
    if a == b:
        val = 1
    return val

def diracND(a,b):
	"""
	Multivariate Dirac function.
	"""
    assert len(a) == len(b)
    n = len(a)
    val = 1
    for i in range(n):
        val = val * dirac1D(a[i], b[i])
    return val

def inverse_transform_sampling(data, eval_points = 'None', n_bins=40, n_samples=1000):
    hist, bin_edges = np.histogram(data, bins=n_bins, density=True)
    cum_values = np.zeros(bin_edges.shape)
    cum_values[1:] = np.cumsum(hist*np.diff(bin_edges))
    cum_values[-1] += 1.e-15
    inv_cdf = interpolate.interp1d(cum_values, bin_edges)
    if eval_points == 'None':
        eval_points = np.random.rand(n_samples)
    return inv_cdf(eval_points)