"""
This script computes the projection vector of the quadratic function 
by making use of the active subspace method. The polynomial chaos 
coefficients are computed analytically and are passed to an 
ActiveSubspaceAdaptation object that computes the 1D active subspace 
and its corresponding eigenvalue.

Author : Panagiotis Tsilifis
Date : 2/8/2017
"""

import numpy as np 
import scipy.stats as st 
import orthpol
import chaos_basispy as cb
import matplotlib.pyplot as plt

def f(xi, a, b, c, W):
    assert xi.shape[0] == 10
    assert W.shape[0] == 10
    return a + b * np.dot(W.T, xi) + c * np.dot(xi.reshape(1,xi.shape[0]), np.dot(np.dot(W, W.T) , xi))

dim = 10
np.random.seed(1234)
W = np.random.normal(size = (dim,1))
W = W / np.linalg.norm(W)

a = np.random.normal()
b = np.random.normal()
c = np.random.normal()
print 'POLYNOMIAL COEFFICIENTS : '
print a, b, c

rvs = [st.uniform(loc = -1., scale = 2.)] * dim
pol = orthpol.ProductBasis(rvs, degree = 2)

xi = st.uniform.rvs(loc = -1. ,scale = 2., size = (100,dim))


coeffs = np.zeros(len(pol._terms))
coeffs[0] = a + c*np.sum(W[:,0]**2) / 3.
for i in range(1, 11):
    coeffs[i] = W[i-1,0] * b / np.sqrt(3.)
for i in range(11,66):
    alpha = pol._terms[i]
    coeffs[i] = (2/3.) * c * np.prod( W[:,0] ** alpha)
    if alpha.max() == 2:
        coeffs[i] = coeffs[i] / np.sqrt(5)


inps = st.uniform.rvs(loc = -1., scale = 2., size = (100, dim))
out1 = np.zeros(100)
for i in range(100):
    out1[i] = f(inps[i,:], a, b, c, W)
out2 = np.dot(pol(inps), coeffs)


basis = cb.ActiveSubspaceAdaptation(dim)
basis._poly_type = 'Legendre'


# Test the stiffness matrix is symmetric 
print "Test that the stiffness matrix is symmetric"
print "Stiff_ij - Stiff_ji"
for i in range(10):
    for j in range(i):
        stiff1 = basis._stiffness_K(2, i, j)
        stiff2 = basis._stiffness_K(2, j, i)
        S = stiff1-stiff2.T
        print np.sum(S.flatten()), np.sum(stiff1.flatten()), np.sum(stiff2.flatten())


for i in range(10):
    stiff = basis._stiffness_K(2, i, i)
    print np.sum(stiff.flatten())

G = basis._grad_covar(2, coeffs)
[l, v] = np.linalg.eigh(G)

print 'EIGENVALUES'
print 'Numerical'+ ' '*5 + '|' + ' '*5 + 'Analytical'
print '-'*40
print str(l[-1]) + ' | ' +str((b**2 + 4 * c**2 / 3.))

plt.plot(l[::-1], 'x')

A = v[:,::-1].T
a = A[0,:]
print 'PROJECTION VECTOR : '
print '-'*30
print 'True'
print '-'*30
print W
print 'Estimated'
print '-'*30
print a.reshape(A.shape[0],1)

plt.plot(a, 'bx', ms = 10)
plt.plot(W[:,0], 'ro', ms = 5)
plt.show()
