"""
Various tool functions used throughout the package.

Author: Panos Tsilifis
Date: 2/7/2017

"""


__all__ = ['dirac1D', 'diracND', 'grad_inner_prod_Legendre', 'gradgrad_inner_prod_Legendre', 'inverse_transform_sampling']


import numpy as np
import scipy.stats as st
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

def grad_inner_prod_Legendre(a,b):
    """
    Computes the inner product (expectation) E[psi'_a(xi) psi_b(xi)],
    where psi'_a(xi) the derivative of a Legendre polynomial of order a
    and psi_b(xi) a Legendre polynomial of order b. Both polynomials as assumed
    normalized and xi follows Uniform(-1,1).
    """
    a = int(a)
    b = int(b)
    val = 0
    if a == 1:
        val = np.sqrt(3) * dirac1D(0,b)
    elif a > 1:
        c =  a - 1
        while c > 0:
            val = val + dirac1D(c, b)
            c = c - 2
        val = val * np.sqrt((2*a+1)*(2*b+1))
    return val


def gradgrad_inner_prod_Legendre(a,b):
    """
    Computes the inner product (expectation) E[psi'_a(xi) psi'_b(xi)],
    where psi'_a(xi) the derivative of a Legendre polynomial of order a
    and psi'_b(xi) the derivative of a Legendre polynomial of order b. Both 
    polynomials as assumed normalized and xi follows Uniform(-1,1).
    """
    a = int(a)
    b = int(b)
    val = 0.
    if a == 1 and b > 0:
        if b == 1:
            val = 3.
        elif b > 1 and b%2 == 0:
            val = 0.
        elif b > 1 and b%2 != 0:
            val = val + np.sqrt(3*(2*b+1))
    elif a > 1 and b > 0:
        c = range(a-1, 0, -2)
        d = range(b-1, 0, -2)
        for i in range(len(c)):
            for j in range(len(d)):
                val = val + (2*c[i]+1) * dirac1D(c[i], d[j])
        val = val * np.sqrt((2*a+1)*(2*b+1))
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
